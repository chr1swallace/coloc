## NB - susieR is not (yet?) on cran. Until then it must live in the Enhances
## catergory, and we need workarounds in these functions and their docs.

## '@importFrom susieR susie_rss susie_get_cs
## NULL
globalVariables(c("pp4", "i", "j"))



##' convert alpha matrix to log BF matrix
##'
##' @title alpha_to_logbf
##' @param alpha an L by p+1 matrix of posterior inclusion probabilites
##' @param pi per snp prior probability of causality
##' @return L by p+1 matrix of log BF
##' @author Chris Wallace
alpha_to_logbf=function(alpha,pi) {
  n=ncol(alpha)-1 # number of snps
  if(length(pi)==1 && pi > 1/n) # scalar pi, but too big
    pi=1/n
  if(length(pi)==1) { # scalar pi
    pi=c(rep(pi,n),1-n*pi)
  }
  if(any(pi == 0)) { # avoid division by zero
    pi[ pi==0 ] = 1e-16
    pi=pi/sum(pi)
  }
  if(any(alpha == 0)) { # avoid division by zero
    alpha[ alpha==0 ] = min(1e-64,min(alpha[alpha>0]))
  }
  ## cope with approx from runsusie
  approx=length(pi) < ncol(alpha)
  if(approx) { # vector pi, add on zeroes to cope with legacy runsusie version
    nullcols=1:(ncol(alpha)-length(pi))
    alpha0=alpha[,nullcols,drop=FALSE]
    alpha=alpha[,-nullcols,drop=FALSE]
    ## pi=c(rep(0, ncol(alpha)-length(pi)), pi)
  }
  bf_unscaled=log(alpha) - matrix(log(pi),nrow(alpha),length(pi),byrow=TRUE)
  bf=bf_unscaled - bf_unscaled[,"null"]
  colnames(bf)=colnames(alpha)
  if(approx) {
    bf0=matrix(apply(bf,1,min),nrow(bf),length(nullcols),dimnames=list(NULL,colnames(alpha0)))
    bf=cbind(bf0,bf)
  }
  bf
}

##' convert logbf matrix to PP matrix
##'
##' @title logbf 2 pp
##' @param bf an L by p or p+1 matrix of log Bayes factors
##' @param pi *either* a scalar representing the prior probability for any snp
##'   to be causal, *or* a full vector of per snp / null prior probabilities
##' @param last_is_null TRUE if last value of the bf vector or last column of a
##'   bf matrix relates to the null hypothesis of no association. This is
##'   standard for SuSiE results, but may not be for BF constructed in other
##'   ways.
##' @return matrix of posterior probabilies, same dimensions as bf
##' @author Chris Wallace
logbf_to_pp=function(bf,pi, last_is_null) {
  n=if(last_is_null) {
      ncol(bf)-1 # number of snps, null is at n+1
    } else {
      n=ncol(bf) # number of snps
    }
  if(length(pi)==1) { # scalar pi
    if(pi > 1/n)
      pi=1/n
    pi=if(last_is_null) {
         c(rep(pi,n),1-n*pi)
       } else {
         rep(pi,n)
       }
  }
  if(any(pi == 0)) { # avoid division by zero
    pi[ pi==0 ] = 1e-16
    pi=pi/sum(pi)
  }
  if(last_is_null) {
    bf=bf - bf[,ncol(bf)] ## scale bf in case needed, so log BF for null = 0
  }
  priors=matrix(log(pi),nrow(bf),ncol(bf),byrow=TRUE)
  denom=matrix(apply(bf + priors,1,logsum), nrow(bf), ncol(bf))
  exp(bf + priors - denom)
}

##' colocalisation with multiple causal variants via SuSiE
##'
##' @title run coloc using susie to detect separate signals
##' @return a list, containing elements
##' * summary a data.table of posterior
##'   probabilities of each global hypothesis, one row per pairwise comparison
##'   of signals from the two traits
##' * results a data.table of detailed results giving the posterior probability
##'   for each snp to be jointly causal for both traits *assuming H4 is true*.
##'   Please ignore this column if the corresponding posterior support for H4
##'   is not high.
##' * priors a vector of the priors used for the analysis
##' @export
##' @author Chris Wallace
##' @param dataset1 *either* a coloc-style input dataset (see
##'   \link{check_dataset}), or the result of running \link{runsusie} on such a
##'   dataset
##' @param dataset2 *either* a coloc-style input dataset (see
##'   \link{check_dataset}), or the result of running \link{runsusie} on such a
##'   dataset
##' @param nref number of individuals from whom the LD matrix was estimated
##' @param susie.args a named list of additional arguments to be passed to
##'   \link{runsusie}
##' @param ... other arguments passed to \link{coloc.bf_bf}, in particular prior
##'   values for causal association with one trait (p1, p2) or both (p12)
coloc.susie=function(dataset1,dataset2, nref, back_calculate_lbf=FALSE, susie.args=list(),  ...) {
  if(!requireNamespace("susieR", quietly = TRUE)) {
    message("please install susieR https://github.com/stephenslab/susieR")
    return(NULL)
  }
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1,nref=nref),susie.args))
  if("susie" %in% class(dataset2))
    s2=dataset2
  else
    s2=do.call("runsusie", c(list(d=dataset2,suffix=1,nref=nref),susie.args))
  cs1=s1$sets
  cs2=s2$sets
  ## cs1=susie_get_cs(s1)
  ## cs2=susie_get_cs(s2)
  if(is.null(cs1$cs) || is.null(cs2$cs) || length(cs1$cs)==0 || length(cs2$cs)==0 )
    return(data.table(nsnps=NA))
  ## names(cs) = paste0("L", which(include_idx))
  ## get_idx=function(x)
  ##   sub("L","",x) %>% as.numeric()
  isnps=setdiff(intersect(colnames(s1$alpha),colnames(s2$alpha)),
                "null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  message("Using ",length(isnps),"/ ",ncol(s1$alpha)-1," and ",ncol(s2$alpha)-1," available")

  idx1=cs1$cs_index #sapply(names(cs1$cs), get_idx) ## use sapply here to keep set names
  idx2=cs2$cs_index #sapply(names(cs2$cs), get_idx)
  ## calculate bayes factors, using susie's priors
  ## if(back_calculate_lbf) {
  ##   bf1=alpha_to_logbf(s1$alpha[idx1,,drop=FALSE], s1$pi)
  ##   bf2=alpha_to_logbf(s2$alpha[idx2,,drop=FALSE], s2$pi)
  ## } else {
    bf1=s1$lbf_variable[idx1,,drop=FALSE][,setdiff(colnames(s1$lbf_variable),"")]
    bf2=s2$lbf_variable[idx2,,drop=FALSE][,setdiff(colnames(s2$lbf_variable),"")]
  ## }

  ret=coloc.bf_bf(bf1,bf2,...)
  ## renumber index to match
  ret$summary[,idx1:=cs1$cs_index[idx1]]
  ret$summary[,idx2:=cs2$cs_index[idx2]]
  ret
}

##' coloc for susie output + a separate BF matrix
##'
##' @title run coloc using susie to detect separate signals
##' @inheritParams coloc.signals
##' @param bf2 named vector of BF, names are snp ids and will be matched to column names of susie object's alpha
##' @param ... other arguments passed to \link{coloc.bf_bf}, in particular prior
##'   values for causal association with one trait (p1, p2) or both (p12)
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
coloc.susie_bf=function(dataset1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, ...) {
  if(!requireNamespace("susieR", quietly = TRUE)) {
    message("please install susieR https://github.com/stephenslab/susieR")
    return(NULL)
  }
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=runsusie(dataset1,suffix=1,...)
  ## cs1=susie_get_cs(s1)
  cs1=s1$sets
  if(is.null(cs1$cs) || length(cs1$cs)==0 )
    return(data.table(nsnps=NA))
  idx1=cs1$cs_index
  ## alpha: an L by p matrix of posterior inclusion probabilites
  isnps=setdiff(intersect(colnames(s1$alpha),names(bf2)),
                "null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  bf1=alpha_to_logbf(alpha=s1$alpha[idx1,,drop=FALSE], pi=s1$pi)

  ret=coloc.bf_bf(bf1,bf2, ...)
  ## renumber index to match
  ret$summary[,idx1:=cs1$cs_index[idx1]]
  ret
}

susie_get_cs_with_names=function(s) {
  if(!requireNamespace("susieR", quietly = TRUE)) {
    message("please install susieR https://github.com/stephenslab/susieR")
    return(NULL)
  }
  sets=susieR::susie_get_cs(s)
  sets$cs=lapply(sets$cs, function(x) structure(x, names=colnames(s$alpha)[x]))
  sets
}

##' Colocalise two datasets represented by Bayes factors
##'
##' This is the workhorse behind many coloc functions
##'
##' @title Coloc data through Bayes factors
##' @inheritParams coloc.signals
##' @param bf1 named vector of BF, or matrix of BF with colnames (cols=snps, rows=signals)
##' @param bf2 named vector of BF, or matrix of BF with colnames (cols=snps, rows=signals)
##' @param trim_by_posterior it is important that the signals to be colocalised
##'   are covered by adequate numbers of snps in both datasets. If TRUE, signals
##'   for which snps in common do not capture least overlap.min proportion of
##'   their posteriors support are dropped and colocalisation not attempted.
##' @param overlap.min see trim_by_posterior
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
coloc.bf_bf=function(bf1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5,trim_by_posterior=TRUE) {
  if(is.vector(bf1))
    bf1=matrix(bf1,nrow=1,dimnames=list(NULL,names(bf1)))
  if(is.vector(bf2))
    bf2=matrix(bf2,nrow=1,dimnames=list(NULL,names(bf2)))
  todo <- as.data.table(expand.grid(i=1:nrow(bf1),j=1:nrow(bf2)))
  todo[,pp4:=0]
  isnps=setdiff(intersect(colnames(bf1),colnames(bf2)),
                "null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  ## scale bf in case needed
  if("null" %in% colnames(bf1))
    bf1=bf1 - matrix(bf1[,"null"],nrow(bf1),ncol(bf1))
  if("null" %in% colnames(bf2))
    bf2=bf2 - matrix(bf2[,"null"],nrow(bf2),ncol(bf2))

  ## check whether isnps covers the signal for each trait
  pp1=logbf_to_pp(bf1,p1, last_is_null=TRUE)
  pp2=logbf_to_pp(bf2,p2, last_is_null=TRUE)
  ph0.1=if("null" %in% colnames(pp1)) { pp1[,"null"] } else { 1 - rowSums(pp1) }
  ph0.2=if("null" %in% colnames(pp2)) { pp2[,"null"] } else { 1 - rowSums(pp2) }
  prop1=rowSums(pp1[,c(isnps),drop=FALSE]) / rowSums(pp1[,setdiff(colnames(pp1),"null"),drop=FALSE])
  prop2=rowSums(pp2[,c(isnps),drop=FALSE]) / rowSums(pp2[,setdiff(colnames(pp2),"null"),drop=FALSE])
  if(trim_by_posterior==TRUE) {
    drop=sapply(1:nrow(todo), function(k) {
      prop1[todo$i[k]] < overlap.min | prop2[todo$j[k]] < overlap.min
    })
    if(all(drop)) {
      warning("snp overlap too small between datasets: too few snps with high posterior in one trait represented in other")
      return(
        list(summary=cbind(
               data.table(nsnps=length(isnps),
                          hit1=colnames(pp1)[apply(pp1,1,which.max)][todo$i],
                          hit2=colnames(pp2)[apply(pp2,1,which.max)][todo$j],
                          PP.H0.abf=pmin(ph0.1[todo$i],ph0.2[todo$j]),
                          PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA),
               todo[,.(idx1=i,idx2=j)]))
      )
    }
    if(any(drop))
      todo=todo[!drop,,drop=FALSE]
  }

  bf1=bf1[,isnps,drop=FALSE]
  bf2=bf2[,isnps,drop=FALSE]
  results <- PP <- vector("list",nrow(todo))
  ## results=lapply(1:nrow(todo), function(k) {
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf1=bf1[todo$i[k], ], bf2=bf2[todo$j[k], ])
    df$internal.sum.lABF <- with(df, bf1 + bf2)
    my.denom.log.abf <- logsum(df$internal.sum.lABF)
    df$SNP.PP.H4 <- exp(df$internal.sum.lABF - my.denom.log.abf)
    pp.abf <- combine.abf(df$bf1, df$bf2, p1, p2, p12)
    PP[[k]]=df$SNP.PP.H4
    if(all(is.na(df$SNP.PP.H4))) {
      df$SNP.PP.H4=0
      pp.abf[1:5]=c(1,0,0,0,0)
    }
    common.snps <- nrow(df)
    hit1=names(which.max(bf1[todo$i[k],])) #df$snp[ which.max(abs(df$bf1)) ]
    if(is.null(hit1)) {
      hit1="-"
      pp.abf[c(1,3)]=c(0,1)
    }
    hit2=names(which.max(bf2[todo$j[k],])) #df$snp[ which.max(abs(df$bf2)) ]
    if(is.null(hit2)) {
      hit2="-"
      pp.abf[c(1,2)]=c(0,1)
    }
    results[[k]]=do.call("data.frame",c(list(nsnps=common.snps, hit1=hit1, hit2=hit2),
                                        as.list(pp.abf)))
  }
  results <- as.data.table(do.call("rbind",results))
  PP <- as.data.table(do.call("cbind",PP))
  if(nrow(todo)>1)
    colnames(PP)=paste0("SNP.PP.H4.row",1:nrow(todo))
  else
    colnames(PP)="SNP.PP.H4.abf"

  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  ## rarely, susie puts the same cred set twice. check, and prune if found
  hits=paste(results$hit1,results$hit2,sep=".")
  if(any(duplicated(hits))) {
    results=results[!duplicated(hits)]
    PP=PP[,!duplicated(hits)]
  }
  PP=cbind(data.table(snp=isnps),PP)
  list(summary=results,
       results=PP,
       priors=c(p1=p1,p2=p2,p12=p12))
}

##' run susie_rss storing some additional information for coloc
##'
##' @title Run susie on a single coloc-structured dataset
##' @param d coloc dataset, must include LD (signed correlation matrix)
##' @param suffix suffix label that will be printed with any error messages
##' @param p prior probability a snp is causal (equivalent to p1 or p2 in
##'   coloc.abf). By default, this is set to NULL, upon which we will set a
##'   small null_weight to pass to susie_rss() (see vignette a06 for details
##'   why). You can override this by setting p as you would p1 or p2 in a coloc
##'   function, but note that you may miss some true signals that way. Also note
##'   that neither of these options correspond to the susie_rss() defaults,
##'   because our goal here is not fine mapping alone.
##' @param trimz used to trim datasets for development purposes
##' @param r2.prune sometimes SuSiE can return multiple signals in high LD. if
##'   you set r2.prune to a value between 0 and 1, sets with index SNPs with LD
##'   greater than r2.prune
##' @param maxit maximum number of iterations for the first run of susie_rss().
##'   If susie_rss() does not report convergence, runs will be extended assuming
##'   repeat_until_convergence=TRUE. Most users will not need to change this
##'   default.
##' @param repeat_until_convergence keep running until susie_rss() indicates
##'   convergence. Default TRUE. If FALSE, susie_rss() will run with maxit
##'   iterations, and if not converged, runsusie() will error. Most users will
##'   not need to change this default.
##' @param s_init used internally to extend runs that haven't converged. don't
##'   use.
##' @param nref number of samples in the dataset used to estimate the LD matrix
##'   in d. If you supply this argument, it will be used to set the z_ld_weight
##'   parameter in susie_rss. This is no longer recommended by the susieR
##'   authors, so this parameter is deprecated.
##' @param ... arguments passed to susie_rss. In particular, if you want to
##'   match some coloc defaults, set
##'
##'   * prior_variance=0.2^2 (if a case-control trait) or (0.15/sd(Y))^2 if a
##'   quantitative trait
##'   * estimate_prior_variance=FALSE
##'
##'   otherwise susie_rss will estimate the prior variance itself
##' @return results of a susie_rss run, with some added dimnames
##' @export
##' @examples
##' library(coloc)
##' data(coloc_test_data)
##' if(requireNamespace("susieR", quietly = TRUE)) {
##'   result=runsusie(coloc_test_data$D1,nref=500)
##'   summary(result)
##' }
##' @author Chris Wallace
runsusie=function(d,suffix=1,p=NULL,
                  trimz=NULL,r2.prune=NULL,
                  maxit=100,repeat_until_convergence=TRUE,
                  s_init=NULL,nref=NULL, ...) {
  if(!requireNamespace("susieR", quietly = TRUE)) {
    message("please install susieR https://github.com/stephenslab/susieR")
    return(NULL)
  }
  ## if(!is.null(ld.prune) && !is.null(ld.merge))
  ##   stop("please specicify at most one of ld.prune and ld.merge")
  check_dataset(d,suffix,req=c("beta","varbeta","LD","snp"))
  check_ld(d,d$LD)
  ##make copies of z and LD so we can subset if needed
  if(!"z" %in% names(d))
    z=d$beta/sqrt(d$varbeta)
  else
    z=d$z
  LD=d$LD[d$snp,d$snp] # just in case
  names(z)=d$snp
  snp=d$snp
  if(!is.null(trimz) && trimz!=0) {
    dbak=d
    keep=abs(z) > abs(trimz)
    message("trimming to subset of SNPs with |z| > ",trimz," : ",sum(keep)," / ",length(keep))
    z=z[keep]
    LD=LD[which(keep),which(keep)]
    snp=snp[keep]
  }
  converged=FALSE;
  ## set some defaults for susie arguments
  susie_args=list(...)
  if(!("z_ld_weight" %in% names(susie_args))) {
    if(!is.null(nref))
      ## stop("Please give nref, the number of samples used to estimate the LD matrix")
      susie_args$z_ld_weight=1/nref
  }
  if(!("null_weight" %in% names(susie_args))) { # set it ourselves
    susie_args$null_weight=if(!is.null(p)) {
                             max(1 - length(d$snp)*p, p)
                           } else {
                             susie_args$null_weight=1/(length(d$snp)+1) #max(1 - length(d$snp)*p, p)
                           }
  }
  while(!converged) {
    message("running max iterations: ",maxit)
  ## if(!is.null(s_init)) {
  ##   res=do.call(susieR::susie_rss,
  ##               c(list(z=z, R=LD, max_iter=maxit, s_init=s_init), susie_args))
  ##   ## res0=susie_rss(z=z,R=LD,z_ld_weight=1/nref,max_iter=maxit,s_init=s_init)
  ##   ## res1=susie_rss(z=z,R=LD,z_ld_weight=1/nref,max_iter=maxit,s_init=s_init,
  ##   ##                null_weight=susie_args$null_weight)
  ##   ## res1$sets
  ## } else {
    res=do.call(susieR::susie_rss,
                c(list(z=z, R=LD, max_iter=maxit), susie_args))
    converged=res$converged; s_init=res; maxit=maxit*2
    message("\tconverged: ",converged)
    if(!converged && repeat_until_convergence==FALSE)
      stop("susie_rss() did not converge in ",maxit," iterations. Try running with run_until_convergence=TRUE")
    if(!converged)
      maxit=maxit * 100 # no point in half measures!
  }
  colnames(res$lbf_variable) = colnames(res$alpha) =
    c(snp,"null")[1:ncol(res$alpha)]
  names(res$pip)=snp
  if(length(res$sets$cs))
    res$sets$cs = lapply(res$sets$cs, function(x) { names(x) = snp[x]; x })
  res$sld=.susie_setld(res$sets$cs,LD)
  res$pruned=FALSE

  ## prune sets in high LD
  if(!is.null(r2.prune))
    res=.susie_prune(res,r2.prune)

  ## ## merge sets in high LD
  ## if(!is.null(r2.merge)) {
  ##   s=res$sets$cs
  ##   if(length(s)>1) {
  ##     ## trim null if present
  ##     sld=.susie_setld(s,LD)
  ##     wh=which(sld>r2.prune,arr.ind=TRUE)
  ##     while(length(wh)) {
  ##       message("merging sets in high LD")
  ##       wh=which(sld==max(sld),arr.ind=TRUE) # make sure we merge one pair at a time
  ##       res = .susie_mergesets(res, wh)
  ##       wh=which(sld>r2.merge,arr.ind=TRUE)
  ##     }
  ##   }
  ## }

  if(!is.null(trimz) && trimz!=0) {
    res$trimz=trimz
    snps_to_add=setdiff(dbak$snp,snp)
    res$pip=c(structure(rep(0,length(snps_to_add)), names=snps_to_add),
              res$pip)
    res$alpha=cbind(matrix(apply(res$alpha,1,min,na.rm=TRUE),
                           nrow(res$alpha),length(snps_to_add),
                           dimnames=list(rownames(res$alpha),snps_to_add)),
                    res$alpha)
    res$alpha=res$alpha/matrix(rowSums(res$alpha),nrow(res$alpha),ncol(res$alph))
  }
  res
}

.susie_prune=function(res,r2.prune) {
  s=res$sets$cs
  if(length(s)>1) {
    wh=which(res$sld>r2.prune,arr.ind=TRUE)
    if(length(wh)) {
      message("pruning sets in high LD")
      res = .susie_dropsets(res, wh)
      res$pruned=TRUE
    }
  }
  res
}

.susie_setld=function(s,ld) {
  ## stmp=lapply(s, setdiff, ncol(ld)+1)
  stmp=lapply(s, names)
  if(!length(stmp))
    return(0)
  if(length(stmp)==1)
    return(matrix(0,1,1,dimnames=list(names(stmp),names(stmp))))
  sld=matrix(0,length(stmp),length(stmp),dimnames=list(names(stmp),names(stmp)))
  for(i in 2:length(stmp)) {
    for(j in 1:(i-1)) {
      if(length(stmp[[i]]) && length(stmp[[j]])) {
        sld[i,j]=max(ld[stmp[[i]],stmp[[j]]]^2, na.rm=TRUE)
      } else {
        sld[i,j]=0
      }
    }
  }
  sld
}

.susie_dropsets=function(res, wh) {
  ## FIX BUG : cs_index and names of cs don't match alpha after drop - just drop from the sets
  ## drop=res$sets$cs_index[as.vector(wh)] # convert to index
  ## for(v in c("alpha","mu","mu2"))
  ##   res[[v]]=res[[v]][-drop,,drop=FALSE]
  ## for(v in c("KL","lbf","V"))
  ##   res[[v]]=res[[v]][-drop]
  drop=as.vector(wh) # raw form
  ## if("sld" %in% names(res))
  ##   res$sld=res$sld[-drop,-drop,drop=FALSE]
  res$sets$cs=res$sets$cs[-drop]
  res$sets$cs_index=res$sets$cs_index[-drop]
  res$sets$purity=res$sets$purity[-drop,,drop=FALSE]
  return(res)
}

.funrows=function(x,m,FUN) {
  x[m[1],] = apply(x[m,,drop=FALSE],2,FUN)
  x[-m[2],,drop=FALSE]
}

.sum1=function(x) {
  pmin(1, sum(x))
}

.susie_mergesets=function(res, wh) {
  m=res$sets$cs_index[as.vector(wh)] # convert to index
  res[["alpha"]] = .funrows(res[["alpha"]], m, .sum1)
  for(v in c("mu","mu2")) {
    res[[v]] = .funrows(res[[v]], m, sum)
  }
  for(v in c("KL","lbf","V"))
    res[[v]]=res[[v]][-m[2]]
  drop=as.vector(wh) # raw form
  res$sets$cs=unlist(res$sets$cs[m])
  res$sets$cs_index=res$sets$cs_index[-m[2]]
  res$sets$purity=res$sets$purity[-drop,,drop=FALSE]
  return(res)
}

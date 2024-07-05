globalVariables(c("pp4", "i", "j","idx"))


##' generic convenience function to convert logbf matrix to PP matrix
##'
##' @title logbf 2 pp
##' @param bf an L by p or p+1 matrix of log Bayes factors
##' @param pi *either* a scalar representing the prior probability for any snp
##'   to be causal, *or* a full vector of per snp / null prior probabilities
##' @param last_is_null TRUE if last value of the bf vector or last column of a
##'   bf matrix relates to the null hypothesis of no association. This is
##'   standard for SuSiE results, but may not be for BF constructed in other
##'   ways.
##' @return matrix of posterior probabilities, same dimensions as bf
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
##' @return a list, containing elements * summary a data.table of posterior
##'   probabilities of each global hypothesis, one row per pairwise comparison
##'   of signals from the two traits * results a data.table of detailed results
##'   giving the posterior probability for each snp to be jointly causal for
##'   both traits *assuming H4 is true*. Please ignore this column if the
##'   corresponding posterior support for H4 is not high. * priors a vector of
##'   the priors used for the analysis
##' @export
##' @author Chris Wallace
##' @param dataset1 *either* a coloc-style input dataset (see
##'   \link{check_dataset}), or the result of running \link{runsusie} on such a
##'   dataset
##' @param dataset2 *either* a coloc-style input dataset (see
##'   \link{check_dataset}), or the result of running \link{runsusie} on such a
##'   dataset
##' @param back_calculate_lbf by default, use the log Bayes factors returned by
##'   susie_rss. It is also possible to back-calculate these from the posterior
##'   probabilities. It is not advised to set this to TRUE, the option exists
##'   really for testing purposes only.
##' @param susie.args a named list of additional arguments to be passed to
##'   \link{runsusie}
##' @param ... other arguments passed to \link{coloc.bf_bf}, in particular prior
##'   values for causal association with one trait (p1, p2) or both (p12) and
##'   and prior weights.
coloc.susie=function(dataset1, dataset2, back_calculate_lbf=FALSE, susie.args=list(), ...) {
  ## if(!requireNamespace("susieR", quietly = TRUE)) {
  ##   message("please install susieR https://github.com/stephenslab/susieR")
  ##   return(NULL)
  ## }
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1),susie.args))
  if("susie" %in% class(dataset2))
    s2=dataset2
  else
    s2=do.call("runsusie", c(list(d=dataset2,suffix=2),susie.args))
  cs1=s1$sets
  cs2=s2$sets
  if(is.null(cs1$cs) || is.null(cs2$cs) || length(cs1$cs)==0 || length(cs2$cs)==0 ) {
	  warning("at least one dataset has no credible sets, nothing to colocalise")
    return(data.table(nsnps=NA))
	  }

  idx1=cs1$cs_index
  idx2=cs2$cs_index
  bf1=s1$lbf_variable[idx1,,drop=FALSE]
  bf2=s2$lbf_variable[idx2,,drop=FALSE]

  ret=coloc.bf_bf(bf1,bf2, ...)
  ## renumber index to match
  ret$summary[,idx1:=cs1$cs_index[idx1]]
  ret$summary[,idx2:=cs2$cs_index[idx2]]
  ret
}

finemap.susie=function(dataset1, susie.args=list(),  ...) {
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1),susie.args))
  cs1=s1$sets
  if(is.null(cs1$cs) || length(cs1$cs)==0 )
    return(data.table(nsnps=NA))

  idx1=cs1$cs_index
  bf1=s1$lbf_variable[idx1,,drop=FALSE]

  ret=finemap.bf(bf1)
  ret$summary[,idx:=cs1$cs_index]
  ret
}


##' coloc for susie output + a separate BF matrix
##'
##' @title run coloc using susie to detect separate signals
##' @inheritParams coloc.signals
##' @param bf2 named vector of log BF, names are snp ids and will be matched to column names of susie object's alpha
##' @param susie.args named list of arguments to be passed to susieR::susie_rss()
##' @param ... other arguments passed to \link{coloc.bf_bf}, in particular prior
##'   values for causal association with one trait (p1, p2) or both (p12)
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
coloc.susie_bf=function(dataset1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, susie.args=list(), ...) {
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=do.call("runsusie", c(list(d=dataset1,suffix=1),susie.args))
  cs1=s1$sets
  if(is.null(cs1$cs) || length(cs1$cs)==0 )
    return(data.table(nsnps=NA))
  idx1=cs1$cs_index
  ## alpha: an L by p matrix of posterior inclusion probabilites
  bf1=s1$lbf_variable[idx1,,drop=FALSE][,setdiff(colnames(s1$lbf_variable),"")]

  ret=coloc.bf_bf(bf1,bf2, ...)
  ## renumber index to match
  ret$summary[,idx1:=cs1$cs_index[idx1]]
  ret
}

##' Finemap one dataset represented by Bayes factors
##'
##' This is the workhorse behind many finemap functions
##'
##' @title Finemap data through Bayes factors
##' @inheritParams finemap.abf
##' @param bf1 named vector of log BF, or matrix of log BF with colnames (cols=snps, rows=signals)
##' @return finemap.signals style result
##' @export
##' @author Chris Wallace
finemap.bf=function(bf1, p1=1e-4) {
  if(is.vector(bf1))
    bf1=matrix(bf1,nrow=1,dimnames=list(NULL,names(bf1)))
  todo <- data.table(i=1:nrow(bf1))
  todo[,pp4:=0]
  isnps=setdiff(colnames(bf1), "null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  ## scale bf in case needed
  if("null" %in% colnames(bf1))
    bf1=bf1 - matrix(bf1[,"null"],nrow(bf1),ncol(bf1))

  ## check whether isnps covers the signal for each trait
  pp1=logbf_to_pp(bf1,p1, last_is_null=TRUE)
  ph0.1=if("null" %in% colnames(pp1)) { pp1[,"null"] } else { 1 - rowSums(pp1) }
  prop1=rowSums(pp1[,c(isnps),drop=FALSE]) / rowSums(pp1[,setdiff(colnames(pp1),"null"),drop=FALSE])
  bf1=bf1[,isnps,drop=FALSE]
  results <- PP <- vector("list",nrow(todo))
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf1=bf1[todo$i[k], ])
    my.denom.log.abf <- logsum(df$bf1)
    df$SNP.PP <- exp(df$bf1 - my.denom.log.abf)
    PP[[k]]=df$SNP.PP
    if(all(is.na(df$SNP.PP))) {
      df$SNP.PP=0
      pp.abf[1:5]=c(1,0,0,0,0)
    }
    common.snps <- nrow(df)
    hit1=names(which.max(bf1[todo$i[k],])) #df$snp[ which.max(abs(df$bf1)) ]
    if(is.null(hit1)) {
      hit1="-"
    }
    results[[k]]=do.call("data.frame",c(list(nsnps=common.snps, hit=hit1)))
  }
  results <- as.data.table(do.call("rbind",results))
  PP <- as.data.table(do.call("cbind",PP))
  if(nrow(todo)>1)
    colnames(PP)=paste0("SNP.PP.row",1:nrow(todo))
  else
    colnames(PP)="SNP.PP.abf"

  ## rarely, susie puts the same cred set twice. check, and prune if found
  hits=results$hit1
  if(any(duplicated(hits))) {
    PP=PP[,!duplicated(hits)]
  }
  PP=cbind(data.table(snp=isnps),PP)
  ## print(results)
  list(summary=results,
       results=PP,
       priors=c(p1=p1))
}

##' Colocalise two datasets represented by Bayes factors
##'
##' This is the workhorse behind many coloc functions
##'
##' @title Coloc data through Bayes factors
##' @inheritParams coloc.signals
##' @param bf1 named vector of log BF, or matrix of BF with colnames (cols=snps, rows=signals)
##' @param bf2 named vector of log BF, or matrix of BF with colnames (cols=snps, rows=signals)
##' @param trim_by_posterior it is important that the signals to be colocalised
##'   are covered by adequate numbers of snps in both datasets. If TRUE, signals
##'   for which snps in common do not capture least overlap.min proportion of
##'   their posteriors support are dropped and colocalisation not attempted.
##' @param overlap.min see trim_by_posterior
##' @param prior_weights1 Non-negative weights for the prior probability a SNP is associated with trait 1 
##' @param prior_weights2 Non-negative weights for the prior probability a SNP is asscoiated with trait 2
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
coloc.bf_bf=function(bf1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5,
                     trim_by_posterior=TRUE, prior_weights1 = NULL, prior_weights2 = NULL) {
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
  ## check p12
  if(length(p12)>1) { # only proceed if didn't need to trim snps
    if(length(isnps)!=ncol(bf1) || length(isnps)!=ncol(bf2))
      stop("p12 must have length 1 or equal length to p1 and p2, but different numbers of snps between datasets")
  }

  ## restrict bf1/2, p1, p2, for incomplete snp overlap
  if(length(isnps)!=ncol(bf1)) {
    keep=match(isnps,colnames(bf1))
    bf1=bf1[,keep,drop=FALSE]
    if(length(p1)>1)
      p1=p1[keep]
  }
  if(length(isnps)!=ncol(bf2)) {
    keep=match(isnps,colnames(bf2))
    bf2=bf2[,keep,drop=FALSE]
    if(length(p2)>2)
      p2=p2[keep]
  }
  ## sort p12 if length(p1)>1 || length(p2)>1
  if(length(p12)==1) {
    if(length(p1)>1 || length(p2)>1) {
      p1_at_dist=min(p1)
      p2_at_dist=min(p2)
      c12=p12 / p1_at_dist / p2_at_dist
      p12=p1*p2*c12
    }
  }

  results <- PP <- vector("list",nrow(todo))
  ## results=lapply(1:nrow(todo), function(k) {
  for(k in 1:nrow(todo)) {
    df <- data.frame(snp=isnps, bf1=bf1[todo$i[k], ], bf2=bf2[todo$j[k], ])
    df$internal.sum.lABF <- with(df, bf1 + bf2)
    my.denom.log.abf <- logsum(df$internal.sum.lABF)
    df$SNP.PP.H4 <- exp(df$internal.sum.lABF - my.denom.log.abf)

    if (!is.null(prior_weights1) || !is.null(prior_weights2)) {
      prior_weights1 <- prior_weights1[match(isnps, colnames(bf1))]
      prior_weights2 <- prior_weights2[match(isnps, colnames(bf2))]
      pp.abf <- combine_abf_weighted(df$bf1, df$bf2,
                                     p1, p2, p12,
                                     prior_weights1,
                                     prior_weights2,
                                     quiet = TRUE)
    } else {
      pp.abf <- combine.abf(df$bf1, df$bf2, p1, p2, p12, quiet=TRUE)
    }

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
  PP=cbind(data.table(snp=isnps),PP)
  list(summary=results,
       results=PP,
       priors=if(length(p1)==1 && length(p2)==1 && length(p12)==1) {
                c(p1=p1,p2=p2,p12=p12)
              } else {
                list(p1=p1,p2=p2,p12=p12)
              })
}




##' run susie_rss storing some additional information for coloc
##'
##' @title Run susie on a single coloc-structured dataset
##' @param d coloc dataset, must include LD (signed correlation matrix) and N
##'   (sample size)
##' @param suffix suffix label that will be printed with any error messages
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
##' result=runsusie(coloc_test_data$D1)
##' summary(result)
##' @author Chris Wallace
runsusie=function(d,suffix=1,
                  maxit=100,repeat_until_convergence=TRUE,
                  s_init=NULL, ...) {
  check_dataset(d,suffix,req=c("beta","varbeta","LD","snp","N"))
  check_ld(d,d$LD)
  ##make copies of z and LD so we can subset if needed
  if(!("z" %in% names(d)))
    z=d$beta/sqrt(d$varbeta)
  else
    z=d$z
  LD=d$LD[d$snp,d$snp] # just in case
  names(z)=d$snp
  snp=d$snp

  converged=FALSE;
  ## set some defaults for susie arguments
  susie_args=list(...)
  if("max_iter" %in% names(susie_args)) {
    maxit=susie_args$max_iter
    susie_args = susie_args[ setdiff(names(susie_args), "max_iter") ]
  }
  ## at 0.12.6 susieR introduced need for n = sample size
  if(!("n" %in% names(susie_args)))
	  susie_args=c(list(n=d$N), susie_args)

  while(!converged) {
    message("running max iterations: ",maxit)
    res=do.call(susie_rss,
                c(list(z=z, R=LD, max_iter=maxit), susie_args))
    converged=res$converged; #s_init=res; maxit=maxit*2
    message("\tconverged: ",converged)
    if(!converged && repeat_until_convergence==FALSE)
      stop("susie_rss() did not converge in ",maxit," iterations. Try running with run_until_convergence=TRUE")
    if(!converged)
      maxit=maxit * 100 # no point in half measures!
  }
 res=annotate_susie(res, snp, LD)
  ## ## prune sets in high LD
  ## if(!is.null(r2.prune))
  ##   res=.susie_prune(res,r2.prune)

  ## if(!all(keep)) {
  ##   res$trimz=trimz
  ##   res$trimpeaks=trimpeaks
  ##   snps_to_add=setdiff(dbak$snp,snp)
  ##   res$pip=c(structure(rep(0,length(snps_to_add)), names=snps_to_add),
  ##             res$pip)
  ##   res$alpha=cbind(matrix(apply(res$alpha,1,min,na.rm=TRUE),
  ##                          nrow(res$alpha),length(snps_to_add),
  ##                          dimnames=list(rownames(res$alpha),snps_to_add)),
  ##                   res$alpha)
  ##   res$alpha=res$alpha/matrix(rowSums(res$alpha),nrow(res$alpha),ncol(res$alph))
  ##   res$lbf_variable=cbind(matrix(apply(res$lbf_variable,1,min,na.rm=TRUE),
  ##                                 nrow(res$lbf_variable),length(snps_to_add),
  ##                                 dimnames=list(rownames(res$lbf_variable),snps_to_add)),
  ##                          res$lbf_variable)
  ## }
  res
}
##' annotate susie_rss output for use with coloc_susie
##'
##' coloc functions need to be able to link summary stats from two
##' different datasets and they do this through snp identifiers.  This
##' function takes the output of susie_rss() and adds snp
##' identifiers. It is entirely the user's responsibility to ensure
##' snp identifiers are in the correct order, coloc cannot make any
##' sanity checks.
##'
##' Note: this annotation step is not needed if you use runsusie() -
##' this is only required if you use the susieR functions directly
##' @param res output of susie_rss()
##' @param snp vector of snp identifiers
##' @param LD matrix of LD (r) between snps in snp
##'     identifiers. Columns, rows should be named by a string that
##'     exists in the vector snp
##' @return res with column names added to some components
##' @export
##' @author Chris Wallace
annotate_susie=function(res,snp, LD) {
    ## if(ncol(res$lbf_variable) != length(snp)+1)
    ##     stop("length of snp vector should be 1 less than ncol(res$lbf_variable)")
    colnames(res$lbf_variable) = c(snp,"null")[1:ncol(res$lbf_variable)]
    colnames(res$alpha) = c(snp,"null")[1:ncol(res$alpha)]
    names(res$pip)=snp
    if(length(res$sets$cs))
        res$sets$cs = lapply(res$sets$cs, function(x) { names(x) = snp[x]; x })
    res$sld=.susie_setld(res$sets$cs,LD)
    res$pruned=FALSE
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

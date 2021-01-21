#'@importFrom susieR susie_rss susie_get_cs
NULL
  ## alpha:
  ## pi:
  ## bf:
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param alpha an L by p+1 matrix of posterior inclusion probabilites
##' @param pi per snp prior probability of causality
##' @param last_is_null TRUE if last column corresponds to NULL (normally TRUE)
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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param bf an L by p or p+1 matrix of log Bayes factors
##' @param pi *either* a scalar representing the prior probability for any snp
##'   to be causal, *or* a full vector of per snp / null prior probabilities
##' @param last_is_null
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

##' colocalisation with multiple causal variants via SUSIE
##'
##' .. content for \details{} ..
##' @title run coloc using susie to detect separate signals
##' @inheritParams coloc.signals
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
##' @param dataset1 *either* a coloc-style input dataset, or the result of
##'   running runsusie on such a dataset
##' @param dataset2 *either* a coloc-style input dataset, or the result of
##'   running runusie on such a dataset
##' @param overlap.min when SNPs differ between datasets, only run coloc for
##'   signals for which at least overlap.min of the posterior support is
##'   captured by overlapping snps
coloc.susie=function(dataset1,dataset2, ...) {
  if("susie" %in% class(dataset1))
    s1=dataset1
  else
    s1=runsusie(dataset1,suffix=1)
  if("susie" %in% class(dataset2))
    s2=dataset2
  else
    s2=runsusie(dataset2,suffix=2)
  cs1=s1$sets
  cs2=s2$sets
  ## cs1=susie_get_cs(s1)
  ## cs2=susie_get_cs(s2)
  if(is.null(cs1$cs) || is.null(cs2$cs) || length(cs1$cs)==0 || length(cs2$cs)==0 )
    return(data.table(nsnps=NA))
  ## names(cs) = paste0("L", which(include_idx))
  ## get_idx=function(x)
  ##   sub("L","",x) %>% as.numeric()
  isnps=intersect(colnames(s1$alpha),colnames(s2$alpha)) %>%setdiff(.,"null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  message("Using ",length(isnps),"/ ",ncol(s1$alpha)-1," and ",ncol(s2$alpha)-1," available")

  idx1=cs1$cs_index #sapply(names(cs1$cs), get_idx) ## use sapply here to keep set names
  idx2=cs2$cs_index #sapply(names(cs2$cs), get_idx)
  ## calculate bayes factors, using susie's priors
  bf1=alpha_to_logbf(s1$alpha[idx1,,drop=FALSE], s1$pi)
  bf2=alpha_to_logbf(s2$alpha[idx2,,drop=FALSE], s2$pi)

  ret=coloc.bf_bf(bf1,bf2,...)
  ## renumber index to match
  ret[,idx1:=cs1$cs_index[idx1]]
  ret[,idx2:=cs2$cs_index[idx2]]
  ret
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title run coloc using susie to detect separate signals
##' @inheritParams coloc.signals
##' @param bf2 named vector of BF, names are snp ids and will be matched to column names of susie object's alpha
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
coloc.susie_bf=function(dataset1,bf2, ...) {
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
  isnps=intersect(colnames(s1$alpha),names(bf2)) %>%setdiff(.,"null")
  if(!length(isnps))
    return(data.table(nsnps=NA))
  bf1=alpha_to_logbf(alpha=s1$alpha[idx1,,drop=FALSE], pi=s1$pi)

  ret=coloc.bf_bf(bf1,bf2,...)
  ## renumber index to match
  ret[,idx1:=cs1$cs_index[idx1]]
  ret
}

susie_get_cs_with_names=function(s) {
  sets=susie_get_cs(s)
  sets$cs=lapply(sets$cs, function(x) structure(x, names=colnames(s$alpha)[x]))
  sets
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title run coloc using susie to detect separate signals
##' @inheritParams coloc.signals
##' @param bf1 named vector of BF, or matrix of BF (cols=snps, rows=signals)
##' @param bf2 named vector of BF, or matrix of BF (cols=snps, rows=signals)
##' @return coloc.signals style result
##' @export
##' @author Chris Wallace
coloc.bf_bf=function(bf1,bf2, p1=1e-4, p2=1e-4, p12=5e-6, overlap.min=0.5,trim_by_posterior=TRUE) {
  if(is.vector(bf1))
    bf1=matrix(bf1,nrow=1,dimnames=list(NULL,names(bf1)))
  if(is.vector(bf2))
    bf2=matrix(bf2,nrow=1,dimnames=list(NULL,names(bf2)))
  todo <- expand.grid(i=1:nrow(bf1),j=1:nrow(bf2)) %>%as.data.table()
  todo[,susie.pp4:=0]
  isnps=intersect(colnames(bf1),colnames(bf2)) %>%setdiff(.,"null")
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
        cbind(
          data.table(nsnps=length(isnps),
                     hit1=colnames(pp1)[apply(pp1,1,which.max)][todo$i],
                     hit2=colnames(pp2)[apply(pp2,1,which.max)][todo$j],
                     PP.H0.abf=pmin(ph0.1[todo$i],ph0.2[todo$j]),
                     PP.H1.abf=NA, PP.H2.abf=NA, PP.H3.abf=NA, PP.H4.abf=NA),
          todo[,.(idx1=i,idx2=j)])
      )
    }
    if(any(drop))
      todo=todo[!drop,,drop=FALSE]
  }

  bf1=bf1[,isnps,drop=FALSE]
  bf2=bf2[,isnps,drop=FALSE]
  results=lapply(1:nrow(todo), function(k) {
    df <- data.frame(snp=isnps, bf1=bf1[todo$i[k], ], bf2=bf2[todo$j[k], ])
    df$internal.sum.lABF <- with(df, bf1 + bf2)
    my.denom.log.abf <- logsum(df$internal.sum.lABF)
    df$SNP.PP.H4 <- exp(df$internal.sum.lABF - my.denom.log.abf)
    pp.abf <- combine.abf(df$bf1, df$bf2, p1, p2, p12)
    if(all(is.na(df$SNP.PP.H4))) {
      df$SNP.PP.H4=0
      pp.abf[1:5]=c(1,0,0,0,0)
    }
    common.snps <- nrow(df)
      hit1=which.max(bf1[todo$i[k],]) %>% names() #df$snp[ which.max(abs(df$bf1)) ]
      if(is.null(hit1)) {
        hit1="-"
        pp.abf[c(1,3)]=c(0,1)
      }
    hit2=which.max(bf2[todo$j[k],]) %>% names() #df$snp[ which.max(abs(df$bf2)) ]
      if(is.null(hit2)) {
        hit2="-"
        pp.abf[c(1,2)]=c(0,1)
      }
    do.call("data.frame",c(list(nsnps=common.snps, hit1=hit1, hit2=hit2), as.list(pp.abf)))
  }) %>% do.call("rbind",.) %>% as.data.table()

  results=cbind(results,todo[,.(idx1=i,idx2=j)])
  ## rarely, susie puts the same cred set twice. check, and prune if found
  hits=paste(results$hit1,results$hit2,sep=".")
  if(any(duplicated(hits)))
    results=results[!duplicated(hits)]

  results
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Run susie on a single coloc-structured dataset
##' @param d dataset
##' @param suffix suffix label that will be printed with any error messages
##' @return results of a susie_rss run, with some added dimnames
##' @author Chris Wallace
runsusie=function(d,suffix=1,nref=503,r2.prune=NULL,r2.merge=NULL,p=1e-4,trimz=NULL,L=10,
                  s_init=NULL,estimate_prior_variance=FALSE) {
  ## if(!is.null(ld.prune) && !is.null(ld.merge))
  ##   stop("please specicify at most one of ld.prune and ld.merge")
  check.dataset(d,suffix,req=c("beta","varbeta","LD","snp"))
  check.ld(d,d$LD)
  if(!"z" %in% names(d))
    d$z=d$beta/sqrt(d$varbeta)
  LD=d$LD[d$snp,d$snp] # just in case
  z=d$z
  names(z)=d$snp
  snp=d$snp
  if(!is.null(trimz)) {
    dbak=d
    keep=abs(z) > abs(trimz)
    message("trimming to subset of SNPs with |z| > ",trimz," : ",sum(keep)," / ",length(keep))
    z=z[keep]
    LD=LD[which(keep),which(keep)]
    snp=snp[keep]
  }
  converged=FALSE; maxit=100;
  while(!converged) {
    message("running iterations: ",maxit)
    if(!is.null(s_init)) {
      res=susie_rss(z, LD, z_ld_weight = 1/nref,
                 ,null_weight=max(1 - length(d$snp)*p, p),
                  ,estimate_prior_method="EM"
                  ## ,estimate_prior_variance=estimate_prior_variance
                  ##  ,prior_variance=0.2^2
                    ## ,estimate_residual_variance=FALSE
                  ## ,verbose=TRUE
                  ,max_iter=maxit
                  ,s_init=s_init )
    } else {
      res=susie_rss(z, LD, z_ld_weight = 1/nref,
                   ,null_weight=max(1 - length(d$snp)*p, p),
                                        ,estimate_prior_method="EM",
                    L=L
                  ## ,estimate_prior_variance=estimate_prior_variance
                  ##  ,prior_variance=0.2^2
                    ## ,estimate_residual_variance=FALSE
                    ## ,verbose=TRUE
                   ,max_iter=maxit)
    }
    converged=res$converged; s_init=res; maxit=maxit*2
    message("\tconverged: ",converged)
  }
  colnames(res$alpha)=c(snp,"null")
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

  if(!is.null(trimz)) {
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
  stmp=lapply(s, setdiff, ncol(ld)+1)
  if(!length(stmp))
    return(0)
  if(length(stmp)==1)
    return(matrix(0,1,1))
  sld=matrix(0,length(stmp),length(stmp))
  for(i in 2:length(stmp))
    for(j in 1:(i-1))
      sld[i,j]=max(ld[stmp[[i]],stmp[[j]]]^2,na.rm=TRUE)
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

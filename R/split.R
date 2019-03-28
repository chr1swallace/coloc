##' Internal function, approx.bf.p
##'
##' Calculate approximate Bayes Factors
##' @title Internal function, approx.bf.p
##' @param p p value
##' @param f MAF
##' @param type "quant" or "cc"
##' @param N sample size
##' @param s proportion of samples that are cases, ignored if type=="quant"
##' @param suffix suffix to append to column names of returned data.frame
##' @return data.frame containing lABF and intermediate calculations
##' @author Claudia Giambartolomei, Chris Wallace
approx.bf.p <- function(p,f,type, N, s, suffix=NULL) {
  if(type=="quant") {
    sd.prior <- 0.15
    V <- Var.data(f, N)
  } else {
    sd.prior <- 0.2
    V <- Var.data.cc(f, N, s)
  }
  z <- qnorm(0.5 * p, lower.tail = FALSE)
  ## Shrinkage factor: ratio of the prior variance to the total variance
  r <- sd.prior^2 / (sd.prior^2 + V)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  lABF = 0.5 * (log(1-r) + (r * z^2))
  ret <- data.frame(V,z,r,lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)  
}

##' Internal function, approx.bf.estimates
##'
##' Calculate approximate Bayes Factors using supplied variance of the regression coefficients
##' @title Internal function, approx.bf.estimates
##' @param z normal deviate associated with regression coefficient and its variance
##' @param V its variance
##' @param sdY standard deviation of the trait. If not supplied, will be estimated.
##' @inheritParams approx.bf.p
##' @return data.frame containing lABF and intermediate calculations
##' @author Vincent Plagnol, Chris Wallace
approx.bf.estimates <- function (z, V, type, suffix=NULL, sdY=1) {
  sd.prior <- if (type == "quant") { 0.15*sdY } else { 0.2 }
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  if(!is.null(suffix))
    colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}


##' Internal function, calculate posterior probabilities for configurations, given logABFs for each SNP and prior probs
##'
##' @title combine.abf
##' @param l1 merged.df$lABF.df1
##' @param l2 merged.df$lABF.df2
##' @inheritParams coloc.abf
##' @return named numeric vector of posterior probabilities
##' @author Claudia Giambartolomei, Chris Wallace
combine.abf <- function(l1, l2, p1, p2, p12) {
  lsum <- l1 + l2
  lH0.abf <- 0
  lH1.abf <- log(p1) + logsum(l1)
  lH2.abf <- log(p2) + logsum(l2)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(l1) + logsum(l2), logsum(lsum))
  lH4.abf <- log(p12) + logsum(lsum)

  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf)
  names(pp.abf) <- paste("PP.H", (1:length(pp.abf)) - 1, ".abf", sep = "")
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  return(pp.abf)
}
##' Internal function, process each dataset list for coloc.abf
##'
##' @title process.dataset
##' @param d list
##' @param suffix "df1" or "df2"
##' @return data.frame with log(abf) or log(bf)
##' @author Chris Wallace
process.dataset <- function(d, suffix) {
  #message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

  if(!(d$type %in% c("quant","cc")))
      stop("dataset ",suffix,": ","type must be quant or cc")
  
  if(d$type=="cc") {
      if(! "s" %in% nd)
          stop("dataset ",suffix,": ","please give s, proportion of samples who are cases")
      if("pvalues" %in% nd && !( "MAF" %in% nd))
          stop("dataset ",suffix,": ","please give MAF if using p values")
      if(d$s<=0 || d$s>=1)
          stop("dataset ",suffix,": ","s must be between 0 and 1")
  }
  
  if(d$type=="quant") {
      if(!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd )))
          stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
  }
  
  if("beta" %in% nd && "varbeta" %in% nd) {  ## use beta/varbeta.  sdY should be estimated by now for quant
    if(length(d$beta) != length(d$varbeta))
      stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")
 
    if(d$type=="quant" && !('sdY' %in% nd)) 
          d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
   df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    return(df)
  }

  if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) { ## no beta/varbeta: use p value / MAF approximation
    if (length(d$pvalues) != length(d$MAF))
      stop('Length of the P-value vectors and MAF vector must match')
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
    df <- data.frame(pvalues = d$pvalues,
                     MAF = d$MAF,
                     snp=as.character(d$snp))    
    colnames(df)[-3] <- paste(colnames(df)[-3], suffix, sep=".")
    df <- subset(df, df$MAF>0 & df$pvalues>0) # all p values and MAF > 0
    abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=d$N, s=d$s, suffix=suffix)
    df <- cbind(df, abf)
    return(df)  
  }

  stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}


##' Bayesian colocalisation analysis, detailed output
##'
##' This function replicates coloc.abf, but outputs more detail for
##' further processing using coloc.process
##'
##' @title Bayesian colocalisation analysis with detailed output
##' @inheritParams coloc.abf
##' @return a list of three \code{data.tables}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item df is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.H4 of the SNP being causal for the shared signal
##' \item df3 is the same for all 2 SNP H3 models
##' }
##' @author Chris Wallace
##' @seealso \code{\link{coloc.process}}, \code{\link{coloc.abf}} 
##' @export
coloc.detail <- function(dataset1, dataset2, MAF=NULL, 
                      p1=1e-4, p2=1e-4, p12=1e-5) {

  if(!is.list(dataset1) || !is.list(dataset2))
    stop("dataset1 and dataset2 must be lists.")
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
    
  df1 <- as.data.table(process.dataset(d=dataset1, suffix="df1"))
  df2 <- as.data.table(process.dataset(d=dataset2, suffix="df2"))
  df <- merge(df1,df2,by="snp")

   if(!nrow(df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  df$lABF.h4 <- with(df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(df$lABF.h4)
  df$H4 <- exp(df$lABF.h4 - my.denom.log.abf)
    df3 <- as.data.table(expand.grid(snp1=df$snp,snp2=df$snp))
    df3 <- merge(df3,df1[,.(snp,lABF.df1)],by.x="snp1",by.y="snp")
    df3 <- merge(df3,df2[,.(snp,lABF.df2)],by.x="snp2",by.y="snp")
    df3[,lABF.h3:=lABF.df1+lABF.df2]
    setnames(df,c("lABF.df1","lABF.df2"),c("lABF.h1","lABF.h2"),skip_absent=TRUE)
    
############################## 

  pp.abf <- combine.abf(df$lABF.h1, df$lABF.h2, p1, p2, p12)  
  common.snps <- nrow(df)
  results <- c(nsnps=common.snps, pp.abf)
  
    list(summary=results,
         df=df,
         df3=df3
         )
}

finemap.indep.signals <- function(D,LD,
                                  r2thr=0.1,
                                  zthr=qnorm(1e-6/2,lower.tail = FALSE)) {
    ## make nicer data table
    x <- as.data.table(D[c("beta","varbeta","snp")])
    x[,z:=beta/sqrt(varbeta)]
    ## check LD in matching order
    LD <- LD[x$snp,x$snp]
    hits <- NULL
    expectedz <- rep(0,nrow(x))
    .f <- function(x,mask=NULL) {
        use <- rep(TRUE,nrow(x))
        if(!is.null(mask)) {
            friends <- apply(abs(LD[,mask,drop=FALSE])>sqrt(r2thr),1,any,na.rm=TRUE)
            use <- use & !friends
            imask <- match(mask,x$snp)
            expectedz <- LD[,mask,drop=FALSE] %*% x$z[imask]
        }
        zdiff <- abs(x$z[use] - expectedz[use])
        wh <- which.max(zdiff)
        if(zdiff[wh] > zthr)
            x$snp[use][wh]
        else
            NULL
    }
    while(length(hits)<10) {
       newhit=.f(x,mask=hits)
        if(is.null(newhit))
            break
       hits <- c(hits,newhit)
    }
    hits
}


group.indep.signals <- function(...,LD,r2thr=0.5) {
    hits <- list(...)
    nhits <- sapply(hits,length)
    hs <- unique(unlist(hits[nhits>1]))
    if(!length(hs)) {
        message("no multiple hits found for any trait")
        return(NULL)
    }
    hld <- LD[hs,hs] > r2thr
    diag(hld) <- FALSE
    g <- graph_from_adjacency_matrix(hld)
    ## plot(g)
    cl <- clusters(g)$membership
    split(names(cl),cl)
}


coloc.process.signals <- function(result,signals,LD=NULL,r2thr=0.1,p1=1e-4,p2=1e-4,p12=1e-5) {
    ## which format?
    if("results" %in% names(result)) {
        result$df <- result$results
        setnames(result$df3,c("snpA","snpB"),c("snp1","snp2"))
    }
    setnames(result$df,paste0("lABF.h",c(1,2,4)),paste0("lbf",c(1,2,4)),skip_absent=TRUE)
    setnames(result$df3,paste0("lABF.h",3),paste0("lbf",3),skip_absent=TRUE)
    ## overall coloc on a given set of snp results
    .f <- function(df,df3) {
        lH0.abf <- 0
        lH1.abf <- log(p1) + logsum(df$lbf1)
        lH2.abf <- log(p2) + logsum(df$lbf2)
        lH3.abf <- log(p1) + log(p2) + logsum(df3$lbf3)
        lH4.abf <- log(p12) + logsum(df$lbf4)
        all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
        best1 <- df[which.max(lbf1),]$snp
        best2 <- df[which.max(lbf2),]$snp
        best4 <- df[which.max(lbf4),]$snp
        my.denom.log.abf <- logsum(all.abf)
        pp.abf <- c(as.list(exp(all.abf - my.denom.log.abf)),list(best1=best1,best2=best2,best4=best4))
        names(pp.abf) <- c(paste("PP.H", (1:(length(pp.abf)-3)) - 1, ".abf", sep = ""),
                           "best1","best2","best4")
        cbind(nsnps=nrow(df),as.data.frame(pp.abf,stringsAsFactors=FALSE))
    }

    if(is.list(signals))
        signals <- unlist(signals)

    ## 
    ldfriends <-  abs(LD[signals,,drop=FALSE])
    masks <- lapply(seq_along(signals),function(i) signals[-i])
    todo <- expand.grid(i=seq_along(signals),j=seq_along(signals))
    todo <- t(todo) #[i<=j,])
    ret <- lapply(1:ncol(todo), function(r) {
        ij <- todo[,r]
        ldin <- apply(ldfriends[ signals[ ij ], ,drop=FALSE ],2,max)
        ldout <- apply(ldfriends[ signals[ -ij ], , drop=FALSE ],2,max)
        dropsnps <- colnames(ldfriends)[ldin < ldout]
        if(length(setdiff(result$df$snp, dropsnps)) <= 1)
            return(NULL)
        ## mask <- masks[[ todo[r,1] ]], mask2[[ todo[r,2] ]]),"")
        ## friends <-  apply(abs(LD[masks[,,drop=FALSE]) > sqrt(r2thr),2,any)
        ## dropsnps <- rownames(LD)[which(friends)]
        cbind(.f(df=result$df[ !(snp %in% dropsnps), ],
                 df3=result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ]),
              hit1=signals[todo[1,r] ],
              hit2=signals[ todo[2,r] ])
    })  %>% rbindlist()

    ret
}
## 599-701
## 928-976

## ret[,e:=PP.H4.abf/(PP.H4.abf+PP.H3.abf)]
## a <- dcast(ret,hit1 ~ hit2, value.var="PP.H4.abf")
## am <- as.matrix(a[,-1])
## rownames(am) <- a$hit1
## pheatmap(am)

## 0
## ## ret <- myouter(mods[[1]], mods[[2]], do.coloc)
## ## ret$g1 <- args$g1
## ## ret$g2 <- args$g2
## ## ret

    
##     if(is.null(hitgroups))
##         return(cbind(.f(result$df,result$df3),unmask=""))
##     sresult <- lapply(seq_along(hitgroups), function(i) {
##         mask <- unlist(hitgroups[-i])
##         dropsnps <- rownames(LD)[ apply(LD[mask,]>r2thr,2,any) ]
##         .f(df=result$df[ !(snp %in% dropsnps), ],
##                  df3=result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ])
##     })
##     sresult <- do.call("rbind",sresult)
##     sresult$unmask <- sapply(hitgroups,paste,collapse="%")
##     sresult 
## }
coloc.process <- function(result,hits1=NULL,hits2=NULL,LD=NULL,r2thr=0.1,p1=1e-4,p2=1e-4,p12=1e-5) {
    ## which format?
    if("results" %in% names(result)) {
        result$df <- result$results
        setnames(result$df3,c("snpA","snpB"),c("snp1","snp2"))
    }
        setnames(result$df,paste0("lABF.h",c(1,2,4)),paste0("lbf",c(1,2,4)),skip_absent=TRUE)
        setnames(result$df3,paste0("lABF.h",3),paste0("lbf",3),skip_absent=TRUE)
    ## overall coloc on a given set of snp results
    .f <- function(df,df3) {
        lH0.abf <- 0
        lH1.abf <- log(p1) + logsum(df$lbf1)
        lH2.abf <- log(p2) + logsum(df$lbf2)
        lH3.abf <- log(p1) + log(p2) + logsum(df3$lbf3)
        lH4.abf <- log(p12) + logsum(df$lbf4)
        all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
        best1 <- df[which.max(lbf1),]$snp
        best2 <- df[which.max(lbf2),]$snp
        best4 <- df[which.max(lbf4),]$snp
        my.denom.log.abf <- logsum(all.abf)
        pp.abf <- c(as.list(exp(all.abf - my.denom.log.abf)),list(best1=best1,best2=best2,best4=best4))
        names(pp.abf) <- c(paste("PP.H", (1:(length(pp.abf)-3)) - 1, ".abf", sep = ""),
                           "best1","best2","best4")
        cbind(nsnps=nrow(df),as.data.frame(pp.abf,stringsAsFactors=FALSE))
    }

    if(is.null(hits1))
        hits1 <- ""
    if(is.null(hits2))
        hits2 <- ""
    
    ## one signal per trait
    if(length(hits1)<=1 & length(hits2)<=1)
        return(cbind(.f(result$df,result$df3),
                     hit1=hits1,hit2=hits2))

    ## ## multi signals per trait, but same masks
    ## if(identical(sort(hits1), sort(hits2))) {
    ## mask1 <- lapply(seq_along(hits1),function(i) hits1[-i])
    ## ret <- lapply(todo), function(r) {
    ##     mask <- setdiff(c(mask1[[ todo[r,1] ]], mask2[[ todo[r,2] ]]),"")
    ##     friends <-  apply(abs(LD[mask,,drop=FALSE]) > sqrt(r2thr),2,any)
    ##     dropsnps <- rownames(LD)[which(friends)]
    ##     cbind(.f(df=result$df[ !(snp %in% dropsnps), ],
    ##              df3=result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ]),
    ##           mask=paste(mask,collapse="/"),
    ##           hit1=hits1[ todo[r,1] ],
    ##           hit2=hits2[ todo[r,2] ])
    ## })  %>% rbindlist()

    
    ## otherwise mask all but one signal from each trait with  >1 signal
      ## 
    ldfriends <-  LD[setdiff(unique(c(hits1,hits2)),""),,drop=FALSE]^2
    todo <- expand.grid(i=seq_along(hits1),j=seq_along(hits2))
    todo <- t(todo) #[i<=j,])
    ret <- lapply(1:ncol(todo), function(r) {
        i <- todo[1,r]
        j <- todo[2,r]
        ldin <- apply(ldfriends[ setdiff(c(hits1[ i ],hits2[j]),""), ,drop=FALSE ],2,max)
        ldout <- apply(ldfriends[ setdiff(c(hits1[-i],hits2[-j]),""), , drop=FALSE ],2,max)
        dropsnps <- colnames(ldfriends)[ldout > r2thr & ldin < ldout]
        if(length(setdiff(result$df$snp, dropsnps)) <= 1)
            return(NULL)
        ## mask <- masks[[ todo[r,1] ]], mask2[[ todo[r,2] ]]),"")
        ## friends <-  apply(abs(LD[masks[,,drop=FALSE]) > sqrt(r2thr),2,any)
        ## dropsnps <- rownames(LD)[which(friends)]
        cbind(.f(df=result$df[ !(snp %in% dropsnps), ],
                 df3=result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ]),
              hit1=hits1[todo[1,r] ],
              hit2=hits2[ todo[2,r] ])
    })  %>% rbindlist()
    
    ret
}

## ret <- myouter(mods[[1]], mods[[2]], do.coloc)
## ret$g1 <- args$g1
## ret$g2 <- args$g2
## ret

    
##     if(is.null(hitgroups))
##         return(cbind(.f(result$df,result$df3),unmask=""))
##     sresult <- lapply(seq_along(hitgroups), function(i) {
##         mask <- unlist(hitgroups[-i])
##         dropsnps <- rownames(LD)[ apply(LD[mask,]>r2thr,2,any) ]
##         .f(df=result$df[ !(snp %in% dropsnps), ],
##                  df3=result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ])
##     })
##     sresult <- do.call("rbind",sresult)
##     sresult$unmask <- sapply(hitgroups,paste,collapse="%")
##     sresult 
## }

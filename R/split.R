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
    
  df1 <- as.data.table(process.dataset(d=dataset1, suffix="df1"))[!is.na(lABF.df1),]
  df2 <- as.data.table(process.dataset(d=dataset2, suffix="df2"))[!is.na(lABF.df2),]
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
         df3=df3,
         priors=c(p1=p1,p2=p2,p12=p12)
         )
}

map_mask <- function(D,LD,r2thr=0.01,sigsnps=NULL) {
    ## make nicer data table
    x <- as.data.table(D[c("beta","varbeta","snp","MAF")])
    x[,z:=beta/sqrt(varbeta)]
    use <- rep(TRUE,nrow(x))
    if(!is.null(sigsnps)) {
         expectedz <- rep(0,nrow(x))
         friends <- apply(abs(LD[,sigsnps,drop=FALSE])>sqrt(r2thr),1,any,na.rm=TRUE)
         use <- use & !friends
         if(!any(use))
             return(NULL)
         imask <- match(sigsnps,x$snp)
         expectedz <- LD[,sigsnps,drop=FALSE] %*% x$z[imask]
         zdiff <- abs(x$z[use]) - abs(expectedz[use])
     } else {
         zdiff <- abs(x$z)
     }
     wh <- which.max(zdiff)
     ## if(zdiff[wh] > zthr)
         structure(x$z[use][wh],names=x$snp[use][wh])
     ## else
     ##     NULL
}

est_cond <- function(x,LD,N,YY,sigsnps,cc=FALSE) {
    LD <- LD[x$snp,x$snp]
    use <- !(rownames(LD) %in% sigsnps) #setdiff(rownames(LD),sigsnps)
    nuse <- match(sigsnps,x$snp)
    VX <- 2 * x$MAF * (1-x$MAF) # expected variance of X, h_buf in GCTA, Dj/N in Yang et al
    D <- VX * N
    ## if case control, LD is from controls, need to adjust
    XX1 <-  N * LD[nuse,nuse,drop=FALSE] * sqrt(matrix(VX[nuse],ncol=1) %*% matrix(VX[nuse],nrow=1))
    if("MAF0" %in% names(x)) {
        VW <- 2 * x$MAF0 * (1-x$MAF0)
        XX1 <-  XX1 / N * sqrt(matrix(VW[nuse],ncol=1) %*% matrix(VW[nuse],nrow=1))
    }
    ## message("YY input: ",YY)
    ## YY <- median( D * x$varbeta * (N-1) + D * x$beta^2 )
    ## ## print(summary( (D * x$varbeta + D * x$beta^2) / (N - 1) ))
    ## message("YY calc = ",YY)
    ## cat("->\t",YY,"\n")

    ## joint effects at sigsnps 
    D1 <-  if(length(nuse)==1) {
               matrix(D[nuse],nrow=1,ncol=1)
           } else {
               diag(D[nuse])
           }
    b1 <- solve(XX1) %*% D1 %*% matrix(x$beta[nuse],ncol=1)

    ## conditional effects 2 | 1
    XX2 <- N * LD[use,use,drop=FALSE] * sqrt(matrix(VX[use],ncol=1) %*% matrix(VX[use],nrow=1)) # cov(X[,use],X[,use])
    XX21 <- N * LD[use,nuse,drop=FALSE] * sqrt(matrix(VX[use],ncol=1) %*% matrix(VX[nuse],nrow=1)) # cov(X[,use],X[,nuse])
    if("MAF0" %in% names(x)) {
        XX2 <-  XX2 / N * sqrt(matrix(VW[use],ncol=1) %*% matrix(VW[use],nrow=1))
        XX21 <-  XX21 / N * sqrt(matrix(VW[use],ncol=1) %*% matrix(VW[nuse],nrow=1))
    }
    D2 <-  if(sum(use)==1) {
               diag(D[use],nrow=1,ncol=1)
           } else {
               diag(D[use])
           }
    b2 <- x$beta[use] - XX21 %*% solve(XX1) %*% D1 %*% matrix(x$beta[nuse],ncol=1) / (D[use])

    ## residual and conditional b2 variance
    Sc <- c(YY -
            matrix(b1,nrow=1) %*% D1 %*% matrix(x$beta[nuse],ncol=1))
    Sc <- ( Sc -  b2 * D[use] * x$beta[use] ) / (N - length(sigsnps) - 1)
    vb2 <- c(Sc) * ( D[use] - diag(XX21 %*% solve(XX1) %*% t(XX21))) / (D[use]^2)
    vb2 <- c(vb2)
    ## na <- apply(abs(LD[sigsnps,use,drop=FALSE])==1,2,any)
    data.table(snp=x$snp[use],beta=c(b2),varbeta=vb2)
    ## if(any(na))
    ##     ret[na,c("beta","varbeta"):=list(NA,NA)]
}

map_cond <- function(D,LD,N,YY,sigsnps=NULL) {
    ## make nicer data table
    x <- as.data.table(D[c("beta","varbeta","snp","MAF")])
    x[,z:=beta/sqrt(varbeta)]
    LD <- LD[x$snp,x$snp]
    if(is.null(sigsnps)) { # take max unconditional abs(z)
        wh <- which.max(x$z)
        ## if(x$z[wh] > zthr)
            return(structure(x$z[wh],names=x$snp[wh]))
        ## else
        ##     return(NULL)
    }
    ## ignore any SNPs in complete LD with sigsnps - results in non invertible matrices
    wh <- setdiff(which(apply(abs(LD[sigsnps,,drop=FALSE])>0.99,2,any)),
                  which(rownames(LD) %in% sigsnps))
    if(length(wh)) {
        drop <- rownames(LD)[wh]
        x <- x[!(snp %in% drop),]
        LD <- LD[-wh,-wh]
    }
    est <- est_cond(x,LD,N,YY,sigsnps)
    Z <- est$beta/sqrt(est$varbeta)
    wh <- which.max(Z)
    ## if(Z[wh] > zthr)
    structure(Z[wh],names=est$snp[wh])
    ## else
    ##     NULL
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param D list of summary stats for a single disease, as defined
##'     for coloc.abf
##' @param LD matrix of signed r values (not rsq!) giving correlation
##'     between SNPs
##' @param aligned if TRUE, the LD matrix is considered aligned to the
##'     effect estimates in D, so that conditional analysis can be
##'     used
##' @param mask use masking if TRUE, otherwise conditioning. defaults
##'     to !aligned
##' @param r2thr if mask==TRUE, all snps will be masked with r2 >
##'     r2thr with any sigsnps. Otherwise ignored
##' @param sigsnps SNPs already deemed significant, to condition on
##' @param pthr when p > pthr, stop successive conditioning
##' @param zthr alternative to specifying pthr, stop when abs(z) <
##'     zthr
##' @param maxhits maximum depth of conditioning. procedure will stop
##'     if p > pthr OR abs(z)<zthr OR maxhits hits have been found.
##' @return list of successive significance fine mapped SNPs, named by
##'     the SNPs
##' @author Chris Wallace
finemap.indep.signals <- function(D,LD=D$LD,aligned=FALSE,
                                  method=if(aligned) { "cond" } else { "mask" },
                                  r2thr=0.01,
                                  sigsnps=NULL,
                                  pthr=1e-6,
                                  zthr=qnorm(pthr/2,lower.tail = FALSE),
                                  maxhits=10) {
    ## make nicer data table
    ## x <- as.data.table(D[c("beta","varbeta","snp","MAF")])
    ## x[,z:=beta/sqrt(varbeta)]
    if(D$type=="cc" & method=="cond") {
        message("approximating linear analysis of binary trait")
        D <- bin2lin(D)
    }
    YY=if(D$type=="quant") {
           D$N * D$sdY^2
       } else {
           D$N * D$s * (1 - D$s)
       }
    ## check LD in matching order
    LD <- LD[D$snp,D$snp]
    hits <- NULL
    while(length(hits)<maxhits) {
        newhit=if(method=="mask") {
                   map_mask(D,LD,r2thr,sigsnps=names(hits))
               } else {
                   map_cond(D,LD,D$N,
                            YY, # sum(Y^2)
                            sigsnps=names(hits))
               }
        if(is.null(newhit) || abs(newhit) < zthr )
            break
        hits <- c(hits,newhit)
    }
    hits
}


## group.indep.signals <- function(...,LD,r2thr=0.5) {
##     hits <- list(...)
##     nhits <- sapply(hits,length)
##     hs <- unique(unlist(hits[nhits>1]))
##     if(!length(hs)) {
##         message("no multiple hits found for any trait")
##         return(NULL)
##     }
##     hld <- LD[hs,hs] > r2thr
##     diag(hld) <- FALSE
##     g <- graph_from_adjacency_matrix(hld)
##     ## plot(g)
##     cl <- clusters(g)$membership
##     split(names(cl),cl)
## }


## coloc.process.signals <- function(result,signals,LD=NULL,r2thr=0.1,p1=1e-4,p2=1e-4,p12=1e-5) {
##     ## which format?
##     if("results" %in% names(result)) {
##         result$df <- result$results
##         setnames(result$df3,c("snpA","snpB"),c("snp1","snp2"))
##     }
##     setnames(result$df,paste0("lABF.h",c(1,2,4)),paste0("lbf",c(1,2,4)),skip_absent=TRUE)
##     setnames(result$df3,paste0("lABF.h",3),paste0("lbf",3),skip_absent=TRUE)
##     ## overall coloc on a given set of snp results
##     .f <- function(df,df3) {
##         lH0.abf <- 0
##         lH1.abf <- log(p1) + logsum(df$lbf1)
##         lH2.abf <- log(p2) + logsum(df$lbf2)
##         lH3.abf <- log(p1) + log(p2) + logsum(df3$lbf3)
##         lH4.abf <- log(p12) + logsum(df$lbf4)
##         all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
##         best1 <- df[which.max(lbf1),]$snp
##         best2 <- df[which.max(lbf2),]$snp
##         best4 <- df[which.max(lbf4),]$snp
##         my.denom.log.abf <- logsum(all.abf)
##         pp.abf <- c(as.list(exp(all.abf - my.denom.log.abf)),list(best1=best1,best2=best2,best4=best4))
##         names(pp.abf) <- c(paste("PP.H", (1:(length(pp.abf)-3)) - 1, ".abf", sep = ""),
##                            "best1","best2","best4")
##         cbind(nsnps=nrow(df),as.data.frame(pp.abf,stringsAsFactors=FALSE))
##     }

##     if(is.list(signals))
##         signals <- unlist(signals)

##     ## 
##     ldfriends <-  abs(LD[signals,,drop=FALSE])
##     masks <- lapply(seq_along(signals),function(i) signals[-i])
##     todo <- expand.grid(i=seq_along(signals),j=seq_along(signals))
##     todo <- t(todo) #[i<=j,])
##     ret <- lapply(1:ncol(todo), function(r) {
##         ij <- todo[,r]
##         ldin <- apply(ldfriends[ signals[ ij ], ,drop=FALSE ],2,max)
##         ldout <- apply(ldfriends[ signals[ -ij ], , drop=FALSE ],2,max)
##         ## dropsnps <- colnames(ldfriends)[ldin < ldout]
##         dropsnps <- colnames(ldfriends)[ldout > r2thr]
##         if(length(setdiff(result$df$snp, dropsnps)) <= 1)
##             return(NULL)
##         ## mask <- masks[[ todo[r,1] ]], mask2[[ todo[r,2] ]]),"")
##         ## friends <-  apply(abs(LD[masks[,,drop=FALSE]) > sqrt(r2thr),2,any)
##         ## dropsnps <- rownames(LD)[which(friends)]
##         cbind(.f(df=result$df[ !(snp %in% dropsnps), ],
##                  df3=result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ]),
##               hit1=signals[todo[1,r] ],
##               hit2=signals[ todo[2,r] ])
##     })  %>% rbindlist()

##     ret
## }
## ## 599-701
## ## 928-976

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

##' Post process a coloc.details result using masking
##'
##' .. content for \details{} ..
##' @title
##' @param result
##' @param hits1
##' @param hits2
##' @param LD
##' @param r2thr
##' @param p1
##' @param p2
##' @param p12
##' @return 
##' @author Chris Wallace
coloc.process <- function(result,hits1=NULL,hits2=NULL,LD=NULL,r2thr=0.01,p1=1e-4,p2=1e-4,p12=1e-6,LD1=LD,LD2=LD) {
    ## which format?
    ## if("results" %in% names(result)) { # for legacy objects
    ##     result$df <- result$results
    ##     setnames(result$df3,c("snpA","snpB"),c("snp1","snp2"))
    ## }
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
        as.data.table(cbind(nsnps=nrow(df),as.data.frame(pp.abf,stringsAsFactors=FALSE)))
    }

    if(is.null(hits1))
        hits1 <- ""
    if(is.null(hits2))
        hits2 <- ""
    
    ## one signal per trait
    if(length(hits1)<=1 & length(hits2)<=1) {
        tmp <- .f(result$df,result$df3)
        tmp$hit1=hits1
        tmp$hit2=hits2
        return(tmp)
    }

    ## otherwise mask all but one signal from each trait with  >1 signal
    ldfriends1 <-  LD1[setdiff(unique(c(hits1)),""),,drop=FALSE]^2
    ldfriends2 <-  LD2[setdiff(unique(c(hits2)),""),,drop=FALSE]^2
    todo <- expand.grid(i=seq_along(hits1),j=seq_along(hits2))
    todo <- t(todo) #[i<=j,])
    ret <- lapply(1:ncol(todo), function(r) {
        i <- todo[1,r]
        j <- todo[2,r]
        message(r)
        ## ldin <- apply(ldfriends[ setdiff(c(hits1[ i ],hits2[j]),""), ,drop=FALSE ],2,max)
        drop1 <- drop2 <- NULL
        if(length(hits1)>1) {
            ldout1 <- apply(ldfriends1[ setdiff(c(hits1[-i]),""), , drop=FALSE ],2,max)
            drop1 <- colnames(ldfriends1)[ldout1 > r2thr]
        }
        if(length(hits2)>1) {
            ldout2 <- apply(ldfriends2[ setdiff(c(hits2[-j]),""), , drop=FALSE ],2,max)
            drop2 <- colnames(ldfriends1)[ldout2 > r2thr]
        }
        ## ld1 <- apply(ldfriends[ setdiff(hits1[-i],""), , drop=FALSE ],2,max)
        ## ld2 <- apply(ldfriends[ setdiff(hits2[-i],""), , drop=FALSE ],2,max)
        ## drop1 <- colnames(ldfriends)[ld1 > r2thr]
        ## drop2 <- colnames(ldfriends)[ld2 > r2thr]
        dropsnps <- unique(c(drop1,drop2))

        message("dropping ",length(dropsnps),"/",nrow(result$df)," : ",length(drop1)," + ",length(drop2))
        if(length(setdiff(result$df$snp, dropsnps)) <= 1)
            return(NULL)
        df <- copy(result$df)
        df3 <- copy(result$df3)
        df[snp %in% drop1, c("lbf1","lbf4"):=list(0,0)]
        df[snp %in% drop2, c("lbf2","lbf4"):=list(0,0)]
        ## df[snp %in% drop2, c("lbf2","lbf4"):=list(log(1/1000),log(1/1000))]
        ## df3[snp1 %in% drop1, lABF.df1:=0]
        ## df3[snp2 %in% drop2, lABF.df2:=0]
        df3[snp1 %in% drop1 | snp2 %in% drop2,lbf3:=0]
        ## mask <- masks[[ todo[r,1] ]], mask2[[ todo[r,2] ]]),"")
        ## friends <-  apply(abs(LD[masks[,,drop=FALSE]) > sqrt(r2thr),2,any)
        ## dropsnps <- rownames(LD)[which(friends)]
        ## cbind(.f(result$df[ !(snp %in% dropsnps), ],
        ##          result$df3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ]),
        tmp <- cbind(.f(df, df3))
        tmp$hit1=hits1[todo[1,r] ]
        tmp$hit2=hits2[ todo[2,r] ]
        ## zz <- signif(coloc.res[2 : 6], 2)
	## 			ii <- zz < 0.01
	## 			zz[ii] <- "< 0.01"
	## 			zz[!ii] <- paste0("= ", zz[!ii])		
	## 			tt <- paste(paste0("H", 0 : 4, " ", zz), collapse = ", ")
						
        ## DF <- melt(df, measure.vars = c("trait", "eqtl"), value.name = "pval")
        ## DF[, dropped := snp %in% dropsnps]
        ## dummy <- rbind(DF[snp %in% hits1 & variable == 'trait', ], DF[snp %in% hits2 & variable == 'eqtl', ]) ​
        ## gg <- ggplot(data = DF, aes(x = pos, y = -log10(pval), colour = dropped)) + geom_point(data = DF, aes(alpha = 0.5)) + geom_vline(data = dummy, aes(xintercept = pos, colour = dropped), linetype = "dashed") + facet_wrap(variable ~ ., nrow = 2, scales = "free_y") + scale_colour_manual(values = c("firebrick1", "grey50")) + theme(legend.position="none") + labs(title = tt)
        ## gg
        tmp
    })  %>% rbindlist()
   ret 
    ## c(result,ret)
}

##' Convert binomial to linear regression
##' 
##' Estimate beta and varbeta if a linear regression had been run on a
##' binary outcome, given log OR and their variance + MAF in controls
##'
##' sets beta = cov(x,y)/var(x)
##' varbeta = (var(y)/var(x) - cov(x,y)^2/var(x)^2)/N
##' @title binomial to linear regression conversion
##' @param D
##' @return D, with original beta and varbeta in beta.bin, varbeta.bin, and beta and varbeta updated to linear estimates
##' @author Chris Wallace
bin2lin <- function(D) {
    if(D$type!="cc")
        stop("type != cc")
    ## est maf in cases
    g0 <- vestgeno.1.ctl(D$MAF)
    g1 <- vestgeno.1.cse(g0,D$beta)
    ex0 <- tcrossprod(c(0,1,2),g0)
    ex1 <- tcrossprod(c(0,1,2),g1)
    ex <- (1 - D$s) * ex0 + D$s * ex1
    vx0 <- tcrossprod(c(0,1,4),g0)
    vx1 <- tcrossprod(c(0,1,4),g1)
    vx <- (1 - D$s) * vx0 + D$s * vx1 - ex^2
    vy <- D$s * (1-D$s)
    vxy <- D$s*ex1 - ex * D$s

    D1 <- D
    ## D1$beta.bin <- D1$beta
    ## D1$varbeta.bin <- D1$varbeta
    D1$beta <- c(vxy/vx)
    D1$varbeta <- c(vy/vx - vxy^2/vx^2)/D$N
    D1
}
   
    
est_all_cond <- function(D,FM) {
    if(D$type=="cc")
        D <- bin2lin(D)
    x <- as.data.table(D[c("beta","varbeta","snp","MAF")])
    x[,z:=beta/sqrt(varbeta)]
    if(length(FM)==1)
        return(structure(list(x),names=names(FM)))
    YY=if(D$type=="quant") {
                        D$N * D$sdY^2
                    } else {
                        D$N * D$s * ( 1 - D$s )
                    }
    cond <- lapply(seq_along(FM), function(i) {
        tmp <- est_cond(x,D$LD,D$N,YY=YY,
                        sigsnps=names(FM)[-i])
        tmp$MAF <- x$MAF[ match(tmp$snp, x$snp) ] 
        tmp
    })
    names(cond) <- names(FM)
    cond
}


coloc.signals <- function(dataset1, dataset2,
                          MAF=NULL, LD=NULL, method="",
                          p1=1e-4, p2=1e-4, p12=NULL, maxhits=10, r2thr=0.01,
                          pthr = 1e-06, zthr = qnorm(pthr/2,lower.tail=FALSE),
                          ...) {
    if(is.null(p12))
        stop("default value for p12 has been removed. please read ... and choose a value appropriate for your study")
    ## error checking
    if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
        dataset1$MAF <- MAF
    if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
        dataset2$MAF <- MAF
    if(!("LD" %in% names(dataset1)) & !is.null(LD)) {
        dataset1$LD <- LD
    }
    if(!("LD" %in% names(dataset2)) & !is.null(LD))
        dataset2$LD <- LD
    if(!("method" %in% names(dataset1)) & !is.null(method))
        dataset1$method <- method
    if(!("method" %in% names(dataset2)) & !is.null(method))
        dataset2$method <- method
    check.dataset(dataset1,1)
    check.dataset(dataset2,2)
    
    ## fine map
    fm1 <- finemap.indep.signals(dataset1,
                                 method=dataset1$method,
                                 maxhits=maxhits,r2thr=r2thr,zthr=zthr,...)
    fm2 <- finemap.indep.signals(dataset2,
                                 method=dataset2$method,
                                 maxhits=maxhits,r2thr=r2thr,zthr=zthr,...)
    print(fm1)
    print(fm2)
    if(is.null(fm1) || is.null(fm2)) {
        warning("signals with p <",pthr,"not present in both datasets")
    }
    if(dataset1$method=="")
        fm1 <- fm1[1]
    if(dataset2$method=="")
        fm2 <- fm2[1]

    ## conditionals if needed
    if(dataset1$method=="cond") {
        cond1 <- est_all_cond(dataset1,fm1)
        X1 <- dataset1[ setdiff(names(dataset1),names(cond1[[1]])) ]
    }
    if(dataset2$method=="cond")  {
        cond2 <- est_all_cond(dataset2,fm2)
        X2 <- dataset2[ setdiff(names(dataset2),names(cond2[[1]])) ]
    }

    
    ## double mask
    if(dataset1$method=="mask" & dataset2$method=="mask") {
        col <- coloc.detail(dataset1,dataset2, p1=p1,p2=p2,p12=p12)
        res <- coloc.process(col, hits1=names(fm1), hits2=names(fm2),
                             LD1=dataset1$LD, LD2=dataset2$LD, r2thr=r2thr, ...)
    } 

    ## double cond
    if(dataset1$method=="cond" & dataset2$method=="cond") {
        todo <- expand.grid(i=seq_along(cond1),j=seq_along(cond2))
        res <- vector("list",nrow(todo))
        for(k in 1:nrow(todo)) {
            col <- coloc.detail(c(cond1[[ todo[k,"i"] ]],X1),
                                c(cond2[[ todo[k,"j"] ]],X2),
                                p1=p1,p2=p2,p12=p12)
            res[[k]] <- coloc.process(col,
                                 hits1=names(fm1)[ todo[k,"i"] ],
                                 hits2=names(fm2)[ todo[k,"j"] ],
                                 LD1=dataset1$LD, LD2=dataset2$LD) #, ...)
        }
        res <- rbindlist(res)
    }

    ## cond-
    if(dataset1$method=="cond" & dataset2$method=="") {
        res <- vector("list",length(cond1))
        for(k in 1:length(res)) {
            col <- coloc.detail(c(cond1[[ k ]],X1),
                                dataset2,
                                p1=p1,p2=p2,p12=p12)
            res[[k]] <- coloc.process(col,
                                 hits1=names(fm1)[ k ],
                                 hits2=names(fm2)[1],
                                 LD1=dataset1$LD, LD2=dataset2$LD,r2thr=r2thr) #, ...)
        }
        res <- rbindlist(res)
    }

    ## cond-mask
    if(dataset1$method=="cond" & dataset2$method=="mask") {
        res <- vector("list",length(cond1))
        for(k in 1:length(res)) {
            col <- coloc.detail(c(cond1[[ k ]],X1),
                                dataset2,
                                p1=p1,p2=p2,p12=p12)
            res[[k]] <- coloc.process(col,
                                 hits1=names(fm1)[ k ],
                                 hits2=names(fm2),
                                 LD1=dataset1$LD, LD2=dataset2$LD,r2thr=r2thr) #, ...)
        }
        res <- rbindlist(res)
    }

    ## -cond
    if(dataset1$method=="" & dataset2$method=="cond") {
        res <- vector("list",length(cond2))
        for(k in 1:length(res)) {
            col <- coloc.detail(dataset1,
                                c(cond2[[ k ]],X2),
                                p1=p1,p2=p2,p12=p12)
            res[[k]] <- coloc.process(col,
                                 hits1=names(fm1)[1],
                                 hits2=names(fm2)[ k ],
                                 LD1=dataset1$LD, LD2=dataset2$LD,r2thr=r2thr) #, ...)
        }
        res <- rbindlist(res)
    }

    ## mask-cond
    if(dataset1$method=="mask" & dataset2$method=="cond") {
        res <- vector("list",length(cond2))
        for(k in 1:length(res)) {
            col <- coloc.detail(dataset1,
                                c(cond2[[ k ]],X2),
                                p1=p1,p2=p2,p12=p12)
            res[[k]] <- coloc.process(col,
                                 hits1=names(fm1),
                                 hits2=names(fm2)[ k ],
                                 LD1=dataset1$LD, LD2=dataset2$LD,r2thr=r2thr) #, ...)
        }
        res <- rbindlist(res)
    }

     p1 <- data.table(hit1=names(fm1),hit1.margz=c(fm1))
    res <- merge(res,p1,by="hit1")
    p2 <- data.table(hit2=names(fm2),hit2.margz=c(fm2))
    res <- merge(res,p2,by="hit2")
    
   res 
}


    


## NB - LD in separate datasets now
        
    

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

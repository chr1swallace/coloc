


##' Bayesian colocalisation analysis, detailed output
##'
##' This function replicates coloc.abf, but outputs more detail for
##' further processing using coloc.process
##'
##' Intended to be called internally by coloc.signals
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
    if("position" %in% names(df1) && "position" %in% names(df2))
        df <- merge(df1,df2,by=c("snp","position"))
    else
        df <- merge(df1,df2,by="snp")

   if(!nrow(df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  df$lABF.h4 <- with(df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(df$lABF.h4)
  df$SNP.PP.H4 <- exp(df$lABF.h4 - my.denom.log.abf)
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
         results=df,
         results.H3=df3,
         priors=c(p1=p1,p2=p2,p12=p12)
         )
}
##' Internal helper function for finemap.signals
##'
##' @title find the next most significant SNP, masking a list of sigsnps
##' @param D dataset in standard coloc format
##' @param LD named matrix of r
##' @param r2thr mask all snps with r2 > r2thr with any in sigsnps
##' @param sigsnps names of snps to mask
##' @return named numeric - Z score named by snp
##' @author Chris Wallace
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

##' Internal helper function for map_cond
##' 
##' @title generate conditional summary stats
##' @param x data.table created from standard coloc dataset
##' @inheritParams map_cond
##' @param cc TRUE if case control
##' @return data.table giving snp, beta and varbeta on remaining snps after conditioning
##' @author Chris Wallace
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
    rbind(data.table(snp=x$snp[use],beta=c(b2),varbeta=vb2),
          data.table(snp=x$snp[!use],beta=rep(0,sum(!use)),varbeta=x$varbeta[!use]))
    ## if(any(na))
    ##     ret[na,c("beta","varbeta"):=list(NA,NA)]
}

##' Internal helper function for finemap.signals
##'
##' @title find the next most significant SNP, conditioning on a list
##'     of sigsnps
##' @param D dataset in standard coloc format
##' @param LD named matrix of r
##' @param N sample size
##' @param YY sum(y^2)
##' @param sigsnps names of snps to mask
##' @return named numeric - Z score named by snp
##' @author Chris Wallace
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


##' Internal helper function
##'
##' @title Pick out snp with most extreme Z score
##' @param D standard format coloc dataset
##' @return z at most significant snp, named by that snp id
##' @author Chris Wallace
find.best.signal <- function(D) {
    ## make nicer data table
    ## x <- as.data.table(D[c("beta","varbeta","snp","MAF")])
    ## x[,z:=beta/sqrt(varbeta)]
    if(D$type=="cc" & D$method=="cond") {
        message("approximating linear analysis of binary trait")
        D <- bin2lin(D)
    }
    YY=if(D$type=="quant") {
           D$N * D$sdY^2
       } else {
           D$N * D$s * (1 - D$s)
       }
    z <- D$beta/sqrt(D$varbeta)
    zdiff <- abs(z)
    wh <- which.max(zdiff)
    structure(z[wh],names=D$snp[wh])
}


##' This is an analogue to finemap.abf, adapted to find multiple
##' signals where they exist, via conditioning or masking - ie a
##' stepwise procedure
##'
##' @title Finemap multiple signals in a single dataset
##' @param D list of summary stats for a single disease, as defined
##'     for coloc.abf
##' @param LD matrix of signed r values (not rsq!) giving correlation
##'     between SNPs
##' @param mask use masking if TRUE, otherwise conditioning. defaults
##'     to TRUE
##' @param r2thr if mask==TRUE, all snps will be masked with r2 >
##'     r2thr with any sigsnps. Otherwise ignored
##' @param sigsnps SNPs already deemed significant, to condition on or
##'     mask
##' @param pthr when p > pthr, stop successive searching
##' @param maxhits maximum depth of conditioning. procedure will stop
##'     if p > pthr OR abs(z)<zthr OR maxhits hits have been found.
##' @export
##' @return list of successively significant fine mapped SNPs, named by
##'     the SNPs
##' @author Chris Wallace
finemap.signals <- function(D,LD=D$LD,aligned=FALSE,
                                  method=if(aligned) { "cond" } else { "mask" },
                                  r2thr=0.01,
                                  sigsnps=NULL,
                                  pthr=1e-6,
                                  maxhits=3) {
    check.dataset(D)
    zthr=qnorm(pthr/2,lower.tail = FALSE)
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
    if(ncol(LD)<2)
        stop("<2 SNPs in D found in LD matrix - do the SNP names match?")
    hits <- NULL
    while(length(hits)<maxhits) {
        newhit=if(method=="mask") {
                   map_mask(D,LD,r2thr,sigsnps=names(hits))
               } else {
                   map_cond(D,LD,D$N,
                            YY, # sum(Y^2)
                            sigsnps=names(hits))
               }
        if(is.null(newhit) || !length(newhit) || abs(newhit) < zthr )
            break
        hits <- c(hits,newhit)
    }
    hits
}


##'  Internal helper function
##'
##' @title Post process a coloc.details result using masking
##' @param result object returned by coloc.detail()
##' @param hits1 lead snps for trait 1. If length > 1, will use masking
##' @param hits2 lead snps for trait 2. If length > 1, will use masking
##' @param LD named LD matrix (for masking)
##' @param LD1 named LD matrix (for masking) for trait 1 only
##' @param LD2 named LD matrix (for masking) for trait 2 only
##' @param r2thr r2 threshold at which to mask
##' @inheritParams coloc.abf
##' @return data.table of coloc results
##' @author Chris Wallace
coloc.process <- function(obj,hits1=NULL,hits2=NULL,LD=NULL,r2thr=0.01,p1=1e-4,p2=1e-4,p12=1e-6,LD1=LD,LD2=LD,detail=FALSE) {
    ## which format?
    ## if("results" %in% names(result)) { # for legacy objects
    ##     obj$results <- obj$results
    ##     setnames(obj$results3,c("snpA","snpB"),c("snp1","snp2"))
    ## }
    setnames(obj$results,paste0("lABF.h",c(1,2,4)),paste0("lbf",c(1,2,4)),skip_absent=TRUE)
    setnames(obj$results.H3,paste0("lABF.h",3),paste0("lbf",3),skip_absent=TRUE)
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
        tmp <- .f(obj$results,obj$results.H3)
        tmp$hit1=hits1
        tmp$hit2=hits2
        ## setnames(obj$results,"H4","SNP.PP.H4")
        return(list(summary=tmp,
                    results=obj$results,
                    priors=c(p1,p2,p12)))
    }

    ## otherwise mask all but one signal from each trait with  >1 signal
    ldfriends1 <-  LD1[setdiff(unique(c(hits1)),""),,drop=FALSE]^2
    ldfriends2 <-  LD2[setdiff(unique(c(hits2)),""),,drop=FALSE]^2
    todo <- expand.grid(i=seq_along(hits1),j=seq_along(hits2))
    todo <- t(todo) #[i<=j,])
    newresult <- obj$results[,intersect(c("snp","position"),
                                      colnames(obj$results)),
                           drop=FALSE,with=FALSE]
    ret <- vector("list",ncol(todo))
    for(r in 1:ncol(todo)) {
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

        message("dropping ",length(dropsnps),"/",nrow(obj$results)," : ",length(drop1)," + ",length(drop2))
        if(length(setdiff(obj$results$snp, dropsnps)) <= 1)
            return(NULL)
        df <- copy(obj$results)
        df3 <- copy(obj$results.H3)
        df[snp %in% drop1, c("lbf1","lbf4"):=list(0,0)]
        df[snp %in% drop2, c("lbf2","lbf4"):=list(0,0)]
        ## PP | H4
        my.denom.log.abf <- logsum(df$lbf4)
        newresult[[ paste0("SNP.PP.H4.row",r) ]] <- exp(df$lbf4 - my.denom.log.abf)
        newresult[[ paste0("z.df1.row",r) ]] <- ifelse(df$snp %in% drop1, 0, df$z.df1)
        newresult[[ paste0("z.df2.row",r) ]] <- ifelse(df$snp %in% drop2, 0, df$z.df2)

        ## df[snp %in% drop2, c("lbf2","lbf4"):=list(log(1/1000),log(1/1000))]
        ## df3[snp1 %in% drop1, lABF.df1:=0]
        ## df3[snp2 %in% drop2, lABF.df2:=0]
        df3[snp1 %in% drop1 | snp2 %in% drop2,lbf3:=0]
        ## mask <- masks[[ todo[r,1] ]], mask2[[ todo[r,2] ]]),"")
        ## friends <-  apply(abs(LD[masks[,,drop=FALSE]) > sqrt(r2thr),2,any)
        ## dropsnps <- rownames(LD)[which(friends)]
        ## cbind(.f(obj$results[ !(snp %in% dropsnps), ],
        ##          obj$results3[ !(snp1 %in% dropsnps | snp2 %in% dropsnps), ]),
        tmp <- cbind(.f(df, df3))
        tmp$hit1=hits1[todo[1,r] ]
        tmp$hit2=hits2[ todo[2,r] ]
        ret[[i]] <- copy(tmp)
    }
    return(list(summary=rbindlist(ret),
                results=newresult,
                priors=c(p1=p1,p2=p2,p12=p12)))
}

##' Convert binomial to linear regression
##' 
##' Estimate beta and varbeta if a linear regression had been run on a
##' binary outcome, given log OR and their variance + MAF in controls
##'
##' sets beta = cov(x,y)/var(x)
##' varbeta = (var(y)/var(x) - cov(x,y)^2/var(x)^2)/N
##' @title binomial to linear regression conversion
##' @param D standard format coloc dataset
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
    x <- as.data.table(D[intersect(c("beta","varbeta","snp","MAF","position"),names(D))])
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
        ## tmp$MAF <- x$MAF[ match(tmp$snp, x$snp) ]
        ## if("position" %in% names(D))
        ##     tmp$position <- x$position[ match(tmp$snp, x$snp) ]
        tmp
    })
    names(cond) <- names(FM)
    cond
}

## merge_cond <- function(cond,D) {
##     m <- match(cond$snp,D$snp)
##     for(nm in intersect(names(D), c("MAF","position")))
        

##' New coloc function, builds on coloc.abf() by allowing for multiple
##' independent causal variants per trait through conditioning or
##' masking.
##'
##' @title Coloc with multiple signals per trait
##' @inheritParams coloc.abf
##' @param LD required if method="cond". matrix of genotype
##'     *correlation* (ie r, not r^2) between SNPs. If dataset1 and
##'     dataset2 may have different LD, you can instead add LD=LD1 to
##'     the list of dataset1 and a different LD matrix for dataset2
##' @param method default "" means do no conditioning, should return
##'     similar to coloc.abf.  if method="cond", then use conditioning
##'     to coloc multiple signals.  if method="mask", use masking to
##'     coloc multiple signals. if different datasets need different
##'     methods (eg LD is only available for one of them) you can set
##'     method on a per-dataset basis by adding method="..." to the
##'     list for that dataset.
##' @param maxhits maximum number of levels to condition/mask
##' @param r2thr if masking, the threshold on r2 should be used to
##'     call two signals independent.  our experience is that this
##'     needs to be set low to avoid double calling the same strong
##'     signal.
##' @param pthr if masking or conditioning, what p value threshold to
##'     call a secondary hit "significant"
##' @export
##' @return data.table of coloc results, one row per pair of lead snps
##'     detected in each dataset
##' @author Chris Wallace
coloc.signals <- function(dataset1, dataset2,
                          MAF=NULL, LD=NULL, method="",
                          p1=1e-4, p2=1e-4, p12=NULL, maxhits=3, r2thr=0.01,
                          pthr = 1e-06) {
    if(is.null(p12))
        stop("default value for p12 has been removed. please read ... and choose a value appropriate for your study")
    ## error checking
    if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
        dataset1$MAF <- MAF
    if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
        dataset2$MAF <- MAF
    if(!("LD" %in% names(dataset1)) & !is.null(LD)) 
        dataset1$LD <- LD
    if(!("LD" %in% names(dataset2)) & !is.null(LD))
        dataset2$LD <- LD
    if(!("method" %in% names(dataset1)) & !is.null(method))
        dataset1$method <- method
    if(!("method" %in% names(dataset2)) & !is.null(method))
        dataset2$method <- method
    check.dataset(dataset1,1)
    check.dataset(dataset2,2)
    zthr = qnorm(pthr/2,lower.tail=FALSE)    

    ## fine map
    fm1 <- finemap.signals(dataset1,
                                 method=dataset1$method,
                                 maxhits=maxhits,r2thr=r2thr,pthr=pthr)
    fm2 <- finemap.signals(dataset2,
                                 method=dataset2$method,
                                 maxhits=maxhits,r2thr=r2thr,pthr=pthr)
    ## print(fm1)
    ## print(fm2)
    if(is.null(fm1)){
        ## warning("no signal with p <",pthr," in dataset1")
        fm1 <- find.best.signal(dataset1)
    }
    if(is.null(fm2)){
        ## warning("no signal with p <",pthr," in dataset1")
        fm2 <- find.best.signal(dataset2)
    }
    if(dataset1$method=="")
        fm1 <- fm1[1]
    if(dataset2$method=="")
        fm2 <- fm2[1]

    ## conditionals if needed
    if(!is.null(fm1) && dataset1$method=="cond") {
        cond1 <- est_all_cond(dataset1,fm1)
        X1 <- dataset1[ intersect(names(dataset1), c("N","sdY","type")) ]
    }
    if(!is.null(fm2) && dataset2$method=="cond")  {
        cond2 <- est_all_cond(dataset2,fm2)
        X2 <- dataset2[ intersect(names(dataset1), c("N","sdY","type")) ]
    }

    ## double mask
    if(dataset1$method=="mask" & dataset2$method=="mask") {
        col <- coloc.detail(dataset1,dataset2, p1=p1,p2=p2,p12=p12)
        res <- coloc.process(col, hits1=names(fm1), hits2=names(fm2),
                             LD1=dataset1$LD, LD2=dataset2$LD, r2thr=r2thr)
    } 

    ## double cond
    if(dataset1$method=="cond" & dataset2$method=="cond") {
        todo <- expand.grid(i=seq_along(cond1),j=seq_along(cond2))
        res <- vector("list",nrow(todo))
        for(k in 1:nrow(todo)) {
            col <- coloc.detail(dataset1=c(cond1[[ todo[k,"i"] ]],X1),
                                dataset2=c(cond2[[ todo[k,"j"] ]],X2),
                                p1=p1,p2=p2,p12=p12)
            res[[k]] <- coloc.process(col,
                                 hits1=names(fm1)[ todo[k,"i"] ],
                                 hits2=names(fm2)[ todo[k,"j"] ],
                                 LD1=dataset1$LD, LD2=dataset2$LD) 
        }
    }

    ## double -
    if(dataset1$method=="" && dataset2$method=="") {
        col <- coloc.detail(dataset1,
                            dataset2,
                            p1=p1,p2=p2,p12=p12)
        res <- coloc.process(col,
                             hits1=names(fm1)[ 1 ],
                             hits2=names(fm2)[1],
                             LD1=dataset1$LD, LD2=dataset2$LD,r2thr=r2thr)
  }
        
    ## cond-
    if(dataset1$method=="cond" && dataset2$method=="") {
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
    }

    ## cond-mask
    if(dataset1$method=="cond" && dataset2$method=="mask") {
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
    }

    ## -cond
    if(dataset1$method=="" && dataset2$method=="cond") {
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
    if(dataset1$method=="mask" && dataset2$method=="cond") {
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

    if(dataset1$method=="cond" || dataset2$method=="cond") {
        summary <- rbindlist(lapply(res,"[[","summary"))
        results <- c(list(res[[1]]$results),
                     lapply(res[-1], function(x) x$results[,grep("SNP.PP.H4|z.df1|z.df2",names(x$results)),drop=FALSE,with=FALSE]))
        results <- do.call("cbind",results)
        ppvars <- grep("SNP.PP.H4",names(results))
        names(results)[ppvars] <- paste0("SNP.PP.H4.row",seq_along(ppvars))
        ppvars <- grep("z.df1",names(results))
        names(results)[ppvars] <- paste0("z.df1.row",seq_along(ppvars))
        ppvars <- grep("z.df2",names(results))
        names(results)[ppvars] <- paste0("z.df2.row",seq_along(ppvars))
    } else {
        summary <- res[["summary"]]
        results <- res[["results"]]
    }

    ## ## if cond-cond, will have lost position
    ## if(dataset1$method=="cond" && dataset2$method=="cond" &&
    ##    ("position" %in% names(dataset1))) {
    ##     m <- match(results$snp,dataset1$position)
    ##     results <- cbind(snp=results$snp, results[,-1],
    ##                      position=dataset1$position[m])
    ## }

    d1 <- data.table(hit1=names(fm1),hit1.margz=c(fm1))
    res <- merge(summary,d1,by="hit1")
    d2 <- data.table(hit2=names(fm2),hit2.margz=c(fm2))
    res <- merge(res,d2,by="hit2")
    
    list(summary=res,
         results=results,
         priors=c(p1=p2,p2=p2,p12=p12))
}


    


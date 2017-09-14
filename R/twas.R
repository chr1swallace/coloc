#'Wrapper to use colocalization testing within a Bayesian model averaging
#'structure.
#'
#'Performs the colocalisation tests described in Plagnol et al (2009)
#' and Wallace et al (2012), adapted to the specific situation of
#' twas, where we are comparing eQTL (trait 1) and GWAS (trait 2)
#' data, and where the number of SNPs may be very large.
#'
#' The adaptations are first to split the SNPs according to LD, then
#' take through to colocalisation testing only those with some
#' evidence for association to the eQTL.  This is done by summing the
#' posterior support for H1, H3 and H4 given by the coloc.abf function
#' within each SNP group and using only SNPs in groups with posterior support of at least 0.05.
#'
#' Note that SNPs groups that show association only to the GWAS trait
#' are ignored, because it is possible for GWAS variants acting
#' through different pathways to be physically proximal, and our focus
#' is on whether the eQTL can explain the GWAS trait.  Note that there
#' *may* exist GWAS variants that are proximal in evolutionary
#' distance (ie in some LD) that act through different pathways.  In
#' that case, this test should reject colocalisation.  This is not
#' appropriate, but developing a test that copes with this better will
#' require further work.
#'
#' The second adaptation is to generate faster (though less accurate)
#' approximations to the Bayes Factors for multi SNP models that are
#' used to weight the colocalistaion test result for each pair of SNPs
#' by using the approximation that
#'
#' ABF = ( BIC1 - BIC0 ) / 2
#'
#' This is for reasons of speed, as more SNPs may enter this test.
#' 
#'@inheritParams coloc.bma
#' @param snps The SNPs to consider as potential explanatory variables
#' @param stratum1 optional column name of df1 that gives stratum information
#' @param stratum2 optional column name of df2 that gives stratum information
#' @param thr posterior probability threshold used to trim SNP list.
#'     Only SNPs with a marginal posterior probability of inclusion
#'     greater than this with one or other trait will be included in
#'     the full BMA analysis
#' @param nsnps number of SNPs required to model both traits.  The BMA
#'     analysis will average over all possible \code{nsnp} SNP models,
#'     subject to \code{thr} above.
#' @param n.approx number of values at which to numerically
#'     approximate the posterior
#' @param bayes.factor if true, compare specific models
#' @param plot.coeff deprecated
#' @param r2.trim for pairs SNPs with r2>\code{r2.trim}, only one SNP
#'     will be retained.  This avoids numerical instability problems
#'     caused by including two highly correlated SNPs in the model.
#' @param quiet suppress messages about how the model spaced is
#'     trimmed for BMA
#' @param bma if true (default), average over models
#' @param ... other parameters passed to \code{coloc.test}
#' @param df1,df2 Each is a dataframe, containing response and
#'     potential explanatory variables for two independent datasets.
#' @param response1,response2 The names of the response variables in
#'     \code{df1} and \code{df2} respectively
#' @param family1,family2 the error family for use in \code{glm}
#'@return a \code{coloc} or \code{colocBayes} object
#'@author Chris Wallace
#'@references Wallace et al (2012).  Statistical colocalisation of
#'     monocyte gene expression and genetic risk variants for type 1
#'     diabetes.  Hum Mol Genet 21:2815-2824.
#'     \url{http://europepmc.org/abstract/MED/22403184}
#'
#'Plagnol et al (2009).  Statistical independence of the colocalized
#'association signals for type 1 diabetes and RPS26 gene expression on
#'chromosome 12q13. Biostatistics 10:327-34.
#'\url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
#'@examples
#'
#'  ## simulate covariate matrix (X) and continuous response vector (Y)
#'  ## for two populations/triats Y1 and Y2 depend equally on f1 and f2
#'  ## within each population, although their distributions differ between
#'  ## populations.  They are compatible with a null hypothesis that they
#'  ## share a common causal variant
#'set.seed(1)
#'  X1 <- matrix(rbinom(2000,1,0.4),ncol=4)
#'  Y1 <- rnorm(500,rowSums(X1[,1:2]),2)
#'  X2 <- matrix(rbinom(2000,1,0.6),ncol=4)
#'  Y2 <- rnorm(500,rowSums(X2[,1:2]),5)
#'  
#'  boxplot(list(Y1,Y2),names=c("Y1","Y2"))
#'  
#'  ## fit and store linear model objects
#'  colnames(X1) <- colnames(X2) <- sprintf("f%s",1:ncol(X1))
#'  summary(lm1 <- lm(Y1~f1+f2+f3+f4,data=as.data.frame(X1)))
#'  summary(lm2 <- lm(Y2~f1+f2+f3+f4,data=as.data.frame(X2)))
#'  
#'
#'## test colocalisation using bma
#'df1=as.data.frame(cbind(Y1=Y1,X1))
#'df2=as.data.frame(cbind(Y2=Y2,X2))
#'
#'result <- coloc.bma( df1, df2, snps=colnames(X1), response1="Y1", response2="Y2",
#'family1="gaussian", family2="gaussian",
#'nsnps=2,bayes.factor=c(1,2,3))
#'result
#'plot(result)
#'
#' ## test colocalisation when one dataset contains a stratifying factor in column named "s"
#' df1$s <- rbinom(500,1,0.5)
#'result <- coloc.bma( df1, df2, snps=colnames(X1), response1="Y1", response2="Y2",
#' stratum1="s",
#'family1="gaussian", family2="gaussian",
#'nsnps=2,bayes.factor=c(1,2,3))
#'result
#'plot(result)
#'
#'@export
#' @importFrom flashClust hclust
#' @importFrom speedglm speedglm
#' @importFrom snpStats ld snp.rhs.estimates
coloc.twas <- function(df1,df2,snps=intersect(setdiff(colnames(df1),c(response1,stratum1)),
                                setdiff(colnames(df2),c(response2,stratum2))),
                      response1="Y", response2="Y",
                      family1="gaussian",family2="binomial",
                      stratum1=NULL, stratum2=NULL,
                      r2.window=200,
                      bayes=!is.null(bayes.factor),
                      thr=0.01,nsnps=2,n.approx=1001, bayes.factor=NULL,
                      plot.coeff=FALSE,r2.trim=0.95,quiet=FALSE,bma=FALSE,...) {
    snps <- unique(snps)
    n.orig <- length(setdiff(snps,c(response1,response2,stratum1,stratum2)))
    if(n.orig<nsnps)
        stop("require at least ",nsnps," SNPs to do colocalisation testing")
    df1 <- as.data.frame(df1[,c(response1,stratum1,snps)])
    df2 <- as.data.frame(df2[,c(response2,stratum2,snps)])
    
    ## remove missings
    message("remove any incomplete observations")
    rmmiss <- function(x){
        use <- complete.cases(x)
        if(!all(use))
            x <- x[which(use),]
        return(x)
    }
    df1 <- rmmiss(df1)
    df2 <- rmmiss(df2)

  ## remove rares (<0.05) and invariants
  message("remove rare variants (MAF < 0.01)")
  x <- rbind(df1[,snps],df2[,snps])
  v <- apply(x,2,var)
  m <- colMeans(x)
  drop <- which(v==0 | m < 0.01/2 | m > 2 - 0.01/2 )
  if(length(drop)) {
      snps <- snps[-drop]
      x <- rbind(df1[,snps],df2[,snps])
  }

  ## cluster SNPs, and calc posterior prob for association
  message("cluster SNPs to define groups for coloc testing")
 ## system.time(r2 <- WGCNA::cor(x)^2)
    x <- new("SnpMatrix",round(as.matrix(x)+1))
    r2 <- ld(x,depth=min(r2.window,ncol(x)),symmetric=TRUE,stat="R.squared")
    hf <- flashClust::hclust(as.dist(1-r2),method="single")
    hfc <- cutree(hf,h=0.9)
    g<-split(names(hfc),hfc)
    message("   ",length(g)," groups found")

    message("quantify support for association within each group")   
    RESULTS <- vector("list",length(g))
    sink("/dev/null")
    for(i in seq_along(g)) {
      d1 <- coloc:::getstats(g[[i]],df1,response1,stratum1)
      d2 <- coloc:::getstats(g[[i]],df2,response2,stratum2)
      RESULTS[[i]] <- coloc.abf(d1,d2)$summary
    }
  sink()
  results <- as.data.frame(do.call("rbind",RESULTS))
    ##results$id <- 1:nrow(results)
    ##results$sumH1 <- with(results,(PP.H1.abf + PP.H3.abf + PP.H4.abf))
    ## with(results,qplot(sumH1))
    use <- results$PP.H0.abf<0.1
    if(!(any(use)))
        stop("no evidence for association")
    plot.data <- vector("list",length(use))
    for(i in which(use)) {
        pcs <- pcs.prepare(df1[,snps], df2[,snps])
        pcs.1 <- pcs.model(pcs, group=1, Y=df1[,response1], threshold=0.8, family=family1, stratum=stratum1)
        pcs.2 <- pcs.model(pcs, group=2, Y=df2[,response2], threshold=0.8, family=family2, stratum=stratum2)
        ct.pcs <- coloc.test(pcs.1,pcs.2)
        if(!all(names(ct.pcs@result) %in% names(results)))
            results[,names(ct.pcs@result)] <- NA
        results[use,names(ct.pcs@result)] <- ct.pcs@result
        plot.data[[i]] <- ct.pcs@plot.data
    }
    results[,"p"] <- pchisq(results[,"chisquare"],df=results[,"n"]-1,lower.tail=FALSE)
    return(new("colocTWAS", result=results, plot.data=plot.data))
}



## coloc.twas.bad <- function(df1,df2,snps=intersect(setdiff(colnames(df1),c(response1,stratum1)),
##                                 setdiff(colnames(df2),c(response2,stratum2))),
##                       response1="Y", response2="Y",
##                       stratum1=NULL, stratum2=NULL,
##                       family1="binomial", family2="binomial",
##                       r2.window=200,
##                       bayes=!is.null(bayes.factor),
##                       thr=0.01,nsnps=2,n.approx=1001, bayes.factor=NULL,
##                       plot.coeff=FALSE,r2.trim=0.95,quiet=FALSE,bma=FALSE,...) {
##     snps <- unique(snps)
##   n.orig <- length(setdiff(snps,c(response1,response2,stratum1,stratum2)))
##   if(n.orig<nsnps)
##     stop("require at least ",nsnps," SNPs to do colocalisation testing")
##   df1 <- as.data.frame(df1[,c(response1,stratum1,snps)])
##   df2 <- as.data.frame(df2[,c(response2,stratum2,snps)])

##   ## remove missings
##   message("remove any incomplete observations")
##   rmmiss <- function(x){
##       use <- complete.cases(x)
##       if(!all(use))
##           x <- x[which(use),]
##     return(x)
##   }
##   df1 <- rmmiss(df1)
##   df2 <- rmmiss(df2)

##   ## remove rares (<0.05) and invariants
##   message("remove rare variants (MAF < 0.05)")
##   x <- rbind(df1[,snps],df2[,snps])
##   v <- apply(x,2,var)
##   m <- colMeans(x)
##   drop <- which(v==0 | m < 0.05/2 | m > 2 - 0.05/2 )
##   if(length(drop)) {
##       snps <- snps[-drop]
##       x <- rbind(df1[,snps],df2[,snps])
##   }

##   ## cluster SNPs, and calc posterior prob for association
##   message("cluster SNPs to define groups for coloc testing")
##  ## system.time(r2 <- WGCNA::cor(x)^2)
##     x <- new("SnpMatrix",round(as.matrix(x)+1))
##     r2 <- ld(x,depth=min(r2.window,ncol(x)),symmetric=TRUE,stat="R.squared")
##     hf <- flashClust::hclust(as.dist(1-r2),method="single")
##     hfc <- cutree(hf,h=0.9)
##     g<-split(names(hfc),hfc)
##     message("   ",length(g)," groups found")

##     message("quantify support for eQTL within each group")   
##     RESULTS <- vector("list",length(g))
##     sink("/dev/null")
##     for(i in seq_along(g)) {
##       d1 <- coloc:::getstats(g[[i]],df1,response1,stratum1)
##       d2 <- coloc:::getstats(g[[i]],df2,response2,stratum2)
##       RESULTS[[i]] <- coloc.abf(d1,d2)$summary
##     }
##   sink()
##   results <- as.data.frame(do.call("rbind",RESULTS))
##   results$id <- 1:nrow(results)
##   results$sumH1 <- with(results,(PP.H1.abf + PP.H3.abf + PP.H4.abf))
##   ## with(results,qplot(sumH1))
##   use <- results$sumH1>0.05
##   if(!(any(use)))
##       stop("no evidence for eQTL")
##   snps <- unlist(g[use])
##   message("   taking ",sum(use)," group(s) of ",length(snps)," SNPs (total) forward for coloc testing")
  
##   ## now tag
##   message("tag selected groups")
##   x <- x[,snps]
##   #summary(x)
##   r2 <- ld(x,x,stat="R.squared")
##   h <- flashClust::hclust(as.dist(1-r2))
##   hc <- cutree(h,h=0.5)
##   g <- split(names(hc),hc)
##   length(g)
##   snps <- sapply(g,"[[",1)
##   groupsize <- sapply(g,length)
##   message("   after tagging, ", length(snps), " remain")
  

## ## evaluate all pairs of interesting SNPs
##     cmb <- combn(1:length(snps),2)
##     message("evaluating ",ncol(cmb)," models")
## csize <- groupsize[cmb[1,]] * groupsize[cmb[2,]]
## cmb <- combn((snps),2)
## ln1 <- log(nrow(df1))
## ln2 <- log(nrow(df2))
## if(is.character(family2))
##     family2 <- if(family2=="binomial") { binomial() } else { gaussian() }
## if(is.character(family1))
##     family1 <- if(family1=="binomial") { binomial() } else { gaussian() }
## RESULTS <- vector("list",ncol(cmb))
##   drop1 <- drop2 <- 1
## if(!is.null(stratum1))
##     drop1 <- c(1,2)
## if(!is.null(stratum2))
##     drop2 <- c(1,2)
##   bic01 <- speedglm::speedglm(paste(response1," ~ ."), data=df1[,c(response1,stratum1),drop=FALSE], family=family1,
##                     k=ln1)$aic
## bic02 <- speedglm::speedglm(paste(response2," ~ ."), data=df2[,c(response2,stratum2),drop=FALSE], family=family2,
##                   k=ln2)$aic
## for(j in 1:ncol(cmb)) {
##     m1 <- speedglm::speedglm(paste(response1," ~ ."), data=df1[,c(response1,stratum1,cmb[,j])], family=family1,
##                    k=ln1)
##     b1 <- coefficients(m1)[-drop1]
##     v1 <- vcov(m1)[-drop1,-drop1]
##     m2 <- speedglm::speedglm(paste(response2," ~ ."), data=df2[,c(response2,stratum2,cmb[,j])], family=family2,
##                    k=ln2)
##     b2 <- coefficients(m2)[-drop2]
##     v2 <- vcov(m2)[-drop2,-drop2]
##     this.coloc <-  coloc.test.summary(b1,b2,v1,v2, plot.coeff=FALSE,bma=TRUE,n.approx=1001)
##    RESULTS[[j]] <- c(this.coloc@result, abf1=(m1$aic - bic01)/2, abf2= (m2$aic-bic01)/2)
## }

## results <- as.data.frame(do.call("rbind",RESULTS))
## results$tagsize <- csize
##     p1 <- results$abf1 + log(results$tagsize)
##     p1 <- p1 - logsum(p1)
##     p2 <- results$abf2 + log(results$tagsize)
##     p2 <- p2 - logsum(p2)
    
##     probs <- exp(p1) + exp(p2)
##     results$probs <- probs/sum(probs)

## ## average p value
## with(results,sum(pchisq(chisquare,n-1) * probs))
## with(results,sum(eta.hat * probs))
## stats <- with(results,
##               c(eta.hat=sum(eta.hat * probs),
##                 n=2,
##                 p=sum(pchisq(chisquare,n-1) * probs)))
## return(new("coloc",
##            result=c(stats["eta.hat"],chisquare=NA,stats["n"],stats["p"]),
##            method="BMA"))

## }

  ## ## require x matrices to be of full rank - but this is nuts when ncol(x) >> nrow(x)
  ## ranker <- function(x) {
  ##     Q <- qr(x)
  ##     colnames(x)[ Q$pivot[seq_len(Q$rank)] ]
  ## }
  ## sel1 <- ranker(df1[,-c(1:2)])


  

prepare.df <- function(df1,df2=NULL,drop.cols,r2.trim,dataset=1,quiet=FALSE) {
  df1 <- rmna(df1,1,quiet=quiet)
  if(!is.null(df2)) {
    df2 <- rmna(df2,2,quiet=quiet)
    x <- rbind(df1[,setdiff(colnames(df1),drop.cols)],df2[,setdiff(colnames(df2),drop.cols)])
  } else {
    x <- df1[,setdiff(colnames(df1),drop.cols)]
  }
  r2 <- cor(x)^2
  d <- as.dist(1-r2)
  h <- hclust(d,method="single")
  hc <- cutree(h,h=1-r2.trim)
  hgrp <- split(names(hc),hc)
  snps <- sapply(hgrp,"[[",1)
  return(list(df1=df1,df2=df2,snps=snps))
}
## prepare.df.old <- function(df1,df2=NULL,drop.cols,r2.trim,dataset=1,quiet=FALSE) {
##   df1 <- rmna(df1,1,quiet=quiet)
##   if(!is.null(df2)) {
##     df2 <- rmna(df2,2,quiet=quiet)
##     x <- rbind(df1[,setdiff(colnames(df1),drop.cols)],df2[,setdiff(colnames(df2),drop.cols)])
##   } else {
##     x <- df1[,setdiff(colnames(df1),drop.cols)]
##   }
##   r2 <- cor(x)^2
##   r2.high <- apply(r2, 1, function(x) which(x>r2.trim)[1])
##   snps <- colnames(x)[ r2.high == 1:length(r2.high) ]
##   return(list(df1=df1,df2=df2,snps=snps))
## }

getstats <- function(snps,df,response,stratum) {
    x <- df[,snps]
    y <- df[,response]
    s <- df[,stratum]
    x <- new("SnpMatrix",as.matrix(x+1))
    ty <- table(y)
    family <- if(length(ty)==2) { "binomial" } else { "gaussian" } # binomial
    f <- if(is.null(stratum)) {
             "y ~ 1"
         } else {
             "y ~ s"
         }
    ss <- snp.rhs.estimates(as.formula(f),snp.data=x,family=family)
    ret <- list(MAF=col.summary(x)[,"MAF"],
                beta=sapply(ss,"[[","beta"),
                varbeta=sapply(ss,"[[","Var.beta"),
                N=nrow(df),
                snp=snps)
    if(family=="binomial") {
        ret <- c(ret,
                 list(type="cc",
                      s=ty[[2]]/sum(ty)))
    } else {
        ret <- c(ret,
                 list(type="quant",
                      sdY=sd(y)))
    }
    ret
}

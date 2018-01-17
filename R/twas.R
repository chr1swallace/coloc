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
#' evidence for association to either trait.  This is done by summing
#' the posterior support for H1, H2, H3 and H4 given by the coloc.abf
#' function within each SNP group and using only SNPs in groups with
#' summed posterior support of at least 0.1.
#'
#' Note that there *may* exist GWAS variants that are proximal in
#' evolutionary distance (ie in some LD) that act through different
#' pathways.  In that case, this test should reject colocalisation.
#' This is not appropriate, but developing a test that copes with this
#' better will require further work.
#'
#'@inheritParams coloc.bma
#' @param maf.min SNPs with MAF < min.maf will be dropped before testing.
#' @param r2.window window (number of SNPs) within which to calculate
#'     rsq to define LD blocks
#' @param pca.thr proportion of variance explained threshold used to
#'     select principal components to include in the test
#' @param n.approx number of values at which to numerically
#'     approximate the posterior
#' @param bayes.factor if true, compare specific models
#' @param ... other parameters passed to \code{coloc.abf}
#' @param nsnps number of SNPs required to model both traits.  
#' @return a \code{colocTWAS}
#' @author Chris Wallace
#' @references Wallace et al (2012).  Statistical colocalisation of
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
#'result <- coloc.twas( df1, df2, snps=colnames(X1), response1="Y1", response2="Y2",
#'family1="gaussian", family2="gaussian",
#'nsnps=2,bayes.factor=c(1,2,3))
#'result
#'
#' ## test colocalisation when one dataset contains a stratifying factor in column named "s"
#' df1$s <- rbinom(500,1,0.5)
#'result <- coloc.twas( df1, df2, snps=colnames(X1), response1="Y1", response2="Y2",
#' stratum1="s",
#'family1="gaussian", family2="gaussian",
#'nsnps=2,bayes.factor=c(1,2,3))
#'result
#'
#'@export
#' @importFrom flashClust hclust
#' @importFrom speedglm speedglm
#' @importFrom snpStats ld snp.rhs.estimates
coloc.twas <- function(df1,df2,snps=intersect(colnames(df1), colnames(df2)),
                       maf.min=0.01,
                       response1="Y", response2="Y",
                      family1="gaussian",family2="binomial",
                      stratum1=NULL, stratum2=NULL,
                      r2.window=min(length(snps)-1,140),
                      pca.thr=0.8,
                      nsnps=2,
                      bayes=!is.null(bayes.factor),
                      thr=0.01,n.approx=1001, bayes.factor=NULL,
                      plot.coeff=FALSE,r2.trim=0.95,quiet=FALSE,bma=FALSE,...) {
    snps <- unique(setdiff(snps,c(response1,response2,stratum1,stratum2)))
    n.orig <- length(snps)
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

  ## remove rares (<0.01) and invariants
  message("remove rare variants : MAF < ",maf.min)
  x <- rbind(df1[,snps],df2[,snps])
  v <- apply(x,2,var)
  m <- colMeans(x)
  drop <- which(v==0 | m < maf.min/2 | m > 2 - maf.min/2 )
  if(length(drop)) {
      snps <- snps[-drop]
      x <- rbind(df1[,snps],df2[,snps])
  }

  ## cluster SNPs, and calc posterior prob for association
    message("cluster SNPs to define groups for coloc testing")
    ## system.time(r2 <- WGCNA::cor(x)^2)
    x <- new("SnpMatrix",round(as.matrix(x)+1))
    r2 <- snpStats::ld(x,depth=min(r2.window,ncol(x)),symmetric=TRUE,stat="R.squared")
    if(any(is.na(r2))) {
        message("NAs in r2 matrix.  setting to 0.  might be risky...")
        r2[is.na(r2)] <- 0
    }
    hf <- flashClust::hclust(as.dist(1-r2),method="single")
    for(hi in seq(0.9,1,by=0.01)) {
        hfc <- cutree(hf,h=hi)
        g<-split(names(hfc),hfc)
        sapply(g,length)
        minl <- min(sapply(g,length))
        if(minl>2)
            break
    }
    message("   ",length(g)," groups found, minimim size ", minl)

    message("quantify support for association within each group")   
    RESULTS <- vector("list",length(g))
    sink("/dev/null")
    for(i in seq_along(g)) {
      d1 <- getstats(g[[i]],df1,response1,stratum1)
      d2 <- getstats(g[[i]],df2,response2,stratum2)
      RESULTS[[i]] <- coloc.abf(d1,d2,...)$summary
    }
  sink()
  results <- as.data.frame(do.call("rbind",RESULTS))
    use <- results$PP.H0.abf<0.1
    if(!(any(use))) {
        message("no evidence for association")
        return(results)
    }
    plot.data <- vector("list",length(use))
    s1 <- if(is.null(stratum1)) { NULL } else { df1[,stratum1] }
    s2 <- if(is.null(stratum2)) { NULL } else { df2[,stratum2] }
    for(i in which(use)) {
        message(i,"\t",length(g[[i]]))
        pcs <- pcs.prepare(df1[,g[[i]]], df2[,g[[i]]])
        pcs.1 <- pcs.model(pcs, group=1, Y=df1[,response1], threshold=pca.thr, family=family1, stratum=s1)
        pcs.2 <- pcs.model(pcs, group=2, Y=df2[,response2], threshold=pca.thr, family=family2, stratum=s2)
        ct.pcs <- coloc.test(pcs.1,pcs.2)
        for(v in c("eta.hat","chisquare","n")) {
            if(!(v %in% names(results)))
                results[[v]] <- NA
            results[i,v] <- ct.pcs@result[v]
        }
        plot.data[[i]] <- ct.pcs@plot.data
    }
    results[,"p"] <- pchisq(results[,"chisquare"],df=results[,"n"]-1,lower.tail=FALSE)
    return(new("colocTWAS", result=results, plot.data=plot.data,snp.blocks=g))
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

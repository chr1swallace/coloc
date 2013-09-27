##' variance of MLE of beta for quantitative trait, assuming var(y)=0
##'
##' Internal function
##' @title Var.data
##' @param f minor allele freq
##' @param N sample number
##' @return variance of MLE beta
##' @author Claudia Giambartolomei
Var.data <- function(f, N) {
  1 / (2 * N * f * (1 - f))
}

##' variance of MLE of beta for case-control
##'
##' Internal function
##' @title Var.data
##' @inheritParams Var.data
##' @param s ???
##' @return variance of MLE beta
##' @author Claudia Giambartolomei
Var.data.cc <- function(f, N, s) {
  1 / (2 * N * f * (1 - f) * s * (1 - s))
}

##' Internal function, logsum
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' @title logsum
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
logsum <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}

##' Internal function, logdiff
##'
##' This function calculates the log of the difference of the exponentiated
##' logs taking out the max, i.e. insuring that the difference is not negative
##' @title logdiff
##' @param x numeric
##' @param y numeric
##' @return max(x) + log(exp(x - max(x,y)) - exp(y-max(x,y)))
##' @author Chris Wallace
logdiff <- function(x,y) {
  my.max <- max(x,y)                              ##take out the maximum value in log form
  my.res <- my.max + log(exp(x - my.max ) - exp(y-my.max))
  return(my.res)
}


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
##' Estimate trait standard deviation given vectors of variance of coefficients,  MAF and sample size
##'
##' Estimate is based on var(beta-hat) = var(Y) / (n * var(X))
##' var(X) = 2*maf*(1-maf)
##' so we can estimate var(Y) by regressing n*var(X) against 1/var(beta)
##' 
##' @title Estimate trait variance, internal function
##' @param vbeta vector of variance of coefficients
##' @param maf vector of MAF (same length as vbeta)
##' @param n sample size
##' @return estimated standard deviation of Y
##' 
##' @author Chris Wallace
sdY.est <- function(vbeta, maf, n) {
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  if(coef(m)[["oneover"]] < 0)
    stop("Trying to estimate trait variance from betas, and getting negative estimate.  Something is wrong.  You can 'fix' this by supplying an estimate of trait standard deviation yourself, as sdY=<value> in the dataset list.")
  return(sqrt(coef(m)[["oneover"]]))
}

##' Internal function, process each dataset list for coloc.abf
##'
##' @title process.dataset
##' @param d list
##' @param suffix "df1" or "df2"
##' @return data.frame with log(abf) or log(bf)
##' @author Chris Wallace
process.dataset <- function(d, suffix) {
  message('Processing dataset')

  nd <- names(d)
  if (! 'type' %in% nd)
    stop('The variable type must be set, otherwise the Bayes factors cannot be computed')

  if("beta" %in% nd && "varbeta" %in% nd && ("MAF" %in% nd || "sdY" %in% nd)) {
    if(length(d$beta) != length(d$varbeta))
      stop("Length of the beta vectors and variance vectors must match")
    if(!("snp" %in% nd))
      d$snp <- sprintf("SNP.%s",1:length(d$beta))
    if(length(d$snp) != length(d$beta))
      stop("Length of snp names and beta vectors must match")
 
    if(d$type == 'quant' & !('sdY' %in% nd)) 
      d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
    
    df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                              V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
    df$snp <- as.character(d$snp)
    return(df)
  }

  if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) {
    if (length(d$pvalues) != length(d$MAF))
      stop('Length of the P-value vectors and MAF vector must match')
    if(d$type=="cc" & !("s" %in% nd))
      stop("Must specify s if type=='cc' and you want to use approximate Bayes Factors")
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

  stop("Must give, as a minimum, either (beta, varbeta, type) or (pvalues, MAF, N, type)")
}

##' Bayesian colocalisation analysis
##'
##' This function calculates posterior probabilities of different
##' causal variant configurations under the assumption of a single
##' causal variant for each trait.
##'
##' If regression coefficients and variances are available, it
##' calculates Bayes factors for association at each SNP.  If only p
##' values are available, it uses an approximation that depends on the
##' SNP's MAF and ignores any uncertainty in imputation.  Regression
##' coefficients should be used if available.
##' 
##' @title Fully Bayesian colocalisation analysis using Bayes Factors
##' @param dataset1 a list with the following elements
##' \describe{
##' 
##'   \item{pvalues}{P-values for each SNP in dataset 1}
##'
##'   \item{N}{Number of samples in dataset 1}
##'
##'   \item{MAF}{minor allele frequency of the variants}
##'
##' \item{beta}{regression coefficient for each SNP from dataset 1}
##' 
##' \item{varbeta}{variance of beta}
##' 
##' \item{type}{the type of data in dataset 1 - either "quant" or "cc" to denote quantitative or case-control}
##'
##' \item{s}{the proportion of samples in dataset 1 that are cases (only relevant for case control samples)}
##'
##' \item{snp}{a character vector of snp ids, optional. If present, it will be used to merge dataset1 and dataset2.  Otherwise, the function assumes dataset1 and dataset2 contain results for the same SNPs in the same order.}
##'
##' }
##'
##' Some of these items may be missing, but you must give \code{type}
##' and then either \code{pvalues}, \code{N} and \code{s} (if
##' type="cc") or \code{beta} and \code{varbeta}.  If you use pvalues,
##' then the function needs to know minor allele frequencies, and will
##' either use the MAF given here or a global estimate of MAF supplied
##' separately.
##' @param dataset2 as above, for dataset 2
##' @param MAF Common minor allele frequency vector to be used for both dataset1 and dataset2
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.H4 of the SNP being causal for the shared signal
##' }
##' @author Claudia Giambartolomei, Chris Wallace
##' @export
coloc.abf <- function(dataset1, dataset2, MAF=NULL, 
                      p1=1e-4, p2=1e-4, p12=1e-5) {

  if(!is.list(dataset1) || !is.list(dataset2))
    stop("dataset1 and dataset2 must be lists.")
  if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
    dataset1$MAF <- MAF
  if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
    dataset2$MAF <- MAF
  
  df1 <- process.dataset(d=dataset1, suffix="df1")
  df2 <- process.dataset(d=dataset2, suffix="df2")
  merged.df <- merge(df1,df2)

   if(!nrow(merged.df))
    stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
  my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
  merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)
  
 
############################## 

  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
  
  output<-list(summary=results, results=merged.df)
  return(output)
}

##' Bayesian colocalisation analysis using data.frames
##'
##' Converts genetic data to snpStats objects, generates p values via score tests, then runs \code{\link{coloc.abf}}
##' 
##' @title Bayesian colocalisation analysis using data.frames
##' @param df1 dataset 1
##' @param df2 dataset 2
##' @param snps col.names for snps
##' @param response1 col.name for response in dataset 1
##' @param response2 col.name for response in dataset 2
##' 
##' @param ... parameters passed to \code{\link{coloc.abf.snpStats}}
##' @return output of \code{\link{coloc.abf}}
##' @export
##' @author Chris Wallace
coloc.abf.datasets <- function(df1,df2,
                               snps=intersect(setdiff(colnames(df1),response1),
                                 setdiff(colnames(df2),response2)),
                               response1="Y", response2="Y", ...) {
  if(length(snps)<2)
    stop("require at least two SNPs in common between df1 and df2 to do anything sensible")
  if(!(response1 %in% colnames(df1)))
    stop(paste("response1",response1,"not found in df1"))
  if(!(response2 %in% colnames(df2)))
    stop(paste("response2",response2,"not found in df2"))
  X1 <- new("SnpMatrix",as.matrix(df1[,snps]))
  X2 <- new("SnpMatrix",as.matrix(df2[,snps]))
  coloc.abf.snpStats(X1,X2,df1[,response1], df2[,response2], ...)
}
##' Bayesian colocalisation analysis using snpStats objects
##'
##' Generates p values via score tests, then runs \code{\link{coloc.abf}}
##' @title Bayesian colocalisation analysis using snpStats objects
##' @param X1 genetic data for dataset 1
##' @param X2 genetic data for dataset 2
##' @param Y1 response for dataset 1
##' @param Y2 response for dataset 2
##' @param snps optional subset of snps to use
##' @param type1 type of data in Y1, "quant" or "cc"
##' @param type2 type of data in Y2, "quant" or "cc"
##' @param s1 the proportion of samples in dataset 1 that are cases (only relevant for case control samples)
##' @param s2 the proportion of samples in dataset 2 that are cases (only relevant for case control samples)
##' @param ... parameters passed to \code{\link{coloc.abf}}
##' @return output of \code{\link{coloc.abf}}
##' @export
##' @author Chris Wallace
coloc.abf.snpStats <- function(X1,X2,Y1,Y2,snps=intersect(colnames(X1),colnames(X2)),
                               type1=c("quant","cc"),type2=c("quant","cc"),s1=NA,s2=NA,...) {
  type1 <- match.arg(type1)
  type2 <- match.arg(type2)
  if(!is(X1,"SnpMatrix") || !is(X2,"SnpMatrix"))
    stop("X1 and X2 must be SnpMatrices")
  if(length(Y1) != nrow(X1) || length(Y2) != nrow(X2))
    stop("length(Y1) != nrow(X1) || length(Y2) != nrow(X2)")
  if(length(snps)<2)
    stop("require at least two SNPs in common between X1 and X2 to do anything sensible")
  X1 <- X1[,snps]
  X2 <- X2[,snps]
  if(is.null(rownames(X1)))
    rownames(X1) <- paste("X1",1:nrow(X1),sep=".")
  if(is.null(rownames(X2)))
    rownames(X2) <- paste("X2",1:nrow(X2),sep=".")
  
  pval1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
  pval2 <- snpStats::p.value(single.snp.tests(phenotype=Y2, snp.data=X2),df=1)
  maf1 <- col.summary(X1)[,"MAF"]
  maf2 <- col.summary(X2)[,"MAF"]
  
  coloc.abf(dataset1=list(pvalues=pval1, N=nrow(X1), MAF=maf1, snp=snps, type=type1, s=s1),
            dataset2=list(pvalues=pval2, N=nrow(X2), MAF=maf2, snp=snps, type=type2, s=s2))
}

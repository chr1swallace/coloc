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



##' Internal function, approx.bf
##'
##' Calculate approximate Bayes Factors
##' @title Internal function, approx.bf
##' @param p p value
##' @param f MAF
##' @param type "quant" or "cc"
##' @param N sample size
##' @param s proportion of samples that are cases, ignored if type=="quant"
##' @param suffix suffix to append to column names of returned data.frame
##' @return data.frame containing lABF and intermediate calculations
##' @author Claudia Giambartolomei, Chris Wallace
approx.bf <- function(p,f,type, N, s, suffix) {
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
  colnames(ret) <- paste(colnames(ret), suffix, sep=".")
  return(ret)
  
}

##' Bayesian colocalisation analysis using summary p values
##'
##' This function takes a data frame obtained by merging p-values for
##' both eQTL and biomarker dataset and returns a list with [1]
##' summary df [2] original df with additional ABF and other values
##' Using MAF from eQTL dataset (column named "MAF.df2")
##' "pvalues.df1" and "pvalues.df2" : names of the colums with p-values
##' N.dataset1 and N.dataset2 number of indviduals used to get the p-values in each dataset
##' sd.prior = standard deviation of prior
##' @title Fully Bayesian colocalisation analysis
##' @param pvalues.dataset1 single variant P-values in dataset 1
##' @param pvalues.dataset2 single variant P-values in dataset 2
##' @param MAF minor allele frequency of the variants
##' @param N.dataset1 number of individuals in dataset 1
##' @param N.dataset2 number of individuals in dataset 2
##' @param type.dataset1 the type of data in dataset 1 - either "quant" or "cc" to denote quantitative or case-control
##' @param type.dataset2 the type of data in dataset 2 
##' @param p1 prior probability a SNP is associated with trait 1
##' @param p2 prior probability a SNP is associated with trait 2
##' @param p12 prior probability a SNP is associated with both traits
##' @param s.dataset1 the proportion of samples in dataset 1 that are cases (only relevant for case control samples)
##' @param s.dataset2 the proportion of samples in dataset 2 that are cases
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item results is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item merged.df is an annotated version of the input \code{data.frame}
##' }
##' @author Claudia Giambartolomei, Chris Wallace
##' @export
coloc.abf <- function(pvalues.dataset1, pvalues.dataset2, MAF , N.dataset1, N.dataset2,
                      type.dataset1="quant", type.dataset2="quant",
                      p1=1e-4, p2=1e-4, p12=1e-5,
                      s.dataset1=0.5, s.dataset2=0.5) {

  if (length(pvalues.dataset1) != length(pvalues.dataset2)) stop('Length of the P-value vectors must match')
  if (length(pvalues.dataset1) != length(MAF)) stop('Length of the P-value vectors and MAF vector must match')
  
  merged.df <- data.frame (pvalues.df1 = pvalues.dataset1,
                           pvalues.df2 = pvalues.dataset2,
                           MAF = MAF)
  
  merged.df <- subset(merged.df, apply(merged.df, 1, min)>0) # all p values and MAF > 0
  
  pvalues.df1 = merged.df[,1]
  pvalues.df2 = merged.df[,2]
  f = merged.df[,3]
  
####### Use different priors and different computation of variance of the mle for case/control vs. quantitative trait
  abf.df1 <- approx.bf(pvalues.df1, f, type.dataset1, N.dataset1, s.dataset1, suffix="df1")
  abf.df2 <- approx.bf(pvalues.df2, f, type.dataset2, N.dataset2, s.dataset2, suffix="df2")
  merged.df <- cbind(merged.df, abf.df1, abf.df2)  
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  
############################## 

  lH0.abf <-  0
  lH1.abf <- log(p1) + logsum(merged.df$lABF.df1)
  lH2.abf <- log(p2) + logsum(merged.df$lABF.df2)
  lH4.abf <- log(p12) + logsum(merged.df$internal.sum.lABF)
  lH3.abf <- log(p1) + log(p2) + logdiff(logsum(merged.df$lABF.df1) + logsum(merged.df$lABF.df2),
                                         logsum(merged.df$internal.sum.lABF))

  ## lH3new.abf = lH3new.f(lH3.abf, lH4.abf, p1, p2, p12)
  ## ## If x=(0.001*y), the difference (x-(0.001*y)) = 0, log(0) = -Inf, so fix this:
  ## ## If it is NaN and the difference between lH3.abf and (log(0.001) + lH4.abf) is very small (<0.001), then keep the old value of lH3.abf:
  ## if (is.na(lH3new.abf) & abs( lH3.abf - ( log(p1*p2/p12) + lH4.abf ) ) < p1*p2/p12)
  ##   lH3new.abf <- lH3.abf
  
#### Now we can compute the PP under each looping through all models as numerator: 
  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3.abf, lH4.abf)
  my.denom.log.abf <- logsum(all.abf)
  pp.abf <- exp(all.abf - my.denom.log.abf) 
  names(pp.abf) <- paste('PP.H', (1:length(pp.abf))-1, '.abf', sep='')
  
  ##  pp.abf <- signif(pp.abf*100,3)
  print(signif(pp.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.abf)
  
  output<-list(results, merged.df)
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
##' @param ... parameters passed to \code{\link{coloc.abf}}
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
##' @param ... parameters passed to \code{\link{coloc.abf}}
##' @return output of \code{\link{coloc.abf}}
##' @export
##' @author Chris Wallace
coloc.abf.snpStats <- function(X1,X2,Y1,Y2,snps=intersect(colnames(X1),colnames(X2)), ...) {
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
  
  p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
  p2 <- snpStats::p.value(single.snp.tests(phenotype=Y2, snp.data=X2),df=1)
  maf <- col.summary(X2)[,"MAF"]
  
  coloc.abf(pvalues.dataset1=p1,pvalues.dataset2=p2,MAF=maf,
            N.dataset1=nrow(X1), N.dataset2=nrow(X2), ...)
}

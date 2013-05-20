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

##' Internal function, logmean
##'
##' This function calculates the log of the sum of the exponentiated
##' logs taking out the max, i.e. insuring that the sum is not Inf
##' @title logmean
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
logmean <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}


##' Internal function, lH3new.f
##'
##'  This is to correct the posterior probability of H3 so that it
##' does not include H4 (same SNP)
##' @title Internal function, lH3new.f
##' @param x lH3
##' @param y lH4
##' @return log(exp(lH3)-0.001*exp(lH4))
##' @author Claudia Giambartolomei
lH3new.f <- function(x, y, prior1, prior2) {
  prop <- prior1^2/prior2
  my.max <- max(c(x,y))                             
  my.res <- my.max + log( exp(x - my.max) - prop*exp(y - my.max))     

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
##' N.df1 and N.df2 number of indviduals used to get the p-values in each dataset
##' sd.prior = standard deviation of prior
##' @title Fully Bayesian colocalisation analysis
##' @param merged.df data.frame containing (at least) three columns:
##' \itemize{ \item the p values in dataset 1 \item the p values in
##' dataset 2, and \item SNP minor allele frequencies, MAF.  } By
##' default, these are named pvalues.df1, pvalues.df2 and MAF, but you
##' may supply your own column names using the argument
##' \code{col.names}
##' @param N.df1 number of individuals in dataset 1
##' @param N.df2 number of individuals in dataset 2
##' @param type.df1, type.df2 the type of data in datasets 1 and 2 - either "quant" or "cc" to denote quantitative or case-control.
##' @param prior1 p_1 == p_2, prior probability a SNP is associated with one trait
##' @param prior2 p_12, prior probability a SNP is associated with both traits
##' @param col.names column names in merged.df that correspond to the
##' pvalues in datasets 1 and 2 and the SNP minor allele frequencies.
##' @param s.df1, s.df2  the proportion of samples in datasets 1 and 2 that are cases.  Only relevant for case control samples.
##' @param sd.prior standard deviation of prior. TODO! CHANGE THIS TO ALLOW SEPARATE PRIORS FOR EACH TRAIT?
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item results is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item merged.df is an annotated version of the input \code{data.frame}
##' }
##' @author Claudia Giambartolomei, Chris Wallace
coloc.config <- function(merged.df, N.df1, N.df2, type.df1="quant", type.df2="quant", prior1= log(10^(-4)), prior2= log(10^(-5)), s.df1=0.5, s.df2=0.5,
                   col.names=c("pvalues.df1","pvalues.df2","MAF")) {  # set s to 0 is it is not a case control study

  if(!is.data.frame(merged.df))
    stop("merged.df should be a data.frame")
  if(!all(col.names %in% colnames(merged.df)))
    stop("not all col.names found in colnames(merged.df)")
  merged.df <- subset(merged.df, apply(merged.df[,col.names], 1, min)>0) # all p values and MAF > 0
  
  pvalues.df1 = merged.df[,col.names[1]]
  pvalues.df2 = merged.df[,col.names[2]]
  f = merged.df[,col.names[3]]
  
####### Use different priors and different computation of variance of the mle for case/control vs. quantitative trait
  abf.df1 <- approx.bf(pvalues.df1, f, type.df1, N.df1, s.df1, suffix="df1")
  abf.df2 <- approx.bf(pvalues.df2, f, type.df2, N.df2, s.df2, suffix="df2")
  merged.df <- cbind(merged.df, abf.df1, abf.df2)  
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  
############################## 

  lH0.abf <- 0
  lH1.abf <- prior1 + logmean(merged.df$lABF.df1)
  lH2.abf <- prior1 + logmean(merged.df$lABF.df2)
  lH3.abf <- prior1 + logmean(merged.df$lABF.df1) + prior1 + logmean(merged.df$lABF.df2)
  lH4.abf <- prior2 + logmean(merged.df$internal.sum.lABF)

  lH3new.abf = lH3new.f(lH3.abf, lH4.abf, prior1, prior2)
  ## If x=(0.001*y), the difference (x-(0.001*y)) = 0, log(0) = -Inf, so fix this:
  ## If it is NaN and the difference between lH3.abf and (log(0.001) + lH4.abf) is very small (<0.001), then keep the old value of lH3.abf:
  if (is.na(lH3new.abf) & abs( lH3.abf -(log(0.001) + lH4.abf) ) <0.001)  (lH3new.abf = lH3.abf)
  
  my.denom.log.abf <- logmean(c(lH0.abf, lH1.abf, lH2.abf, lH3new.abf, lH4.abf))
  
#### Now we can compute the PP under each looping through all models as numerator: 
  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3new.abf, lH4.abf)
  pp.df.abf <- exp(all.abf - my.denom.log.abf) 
  names(pp.df.abf) <- paste('PP.H', (1:length(pp.df.abf))-1, '.abf', sep='')
  
  ##  pp.df.abf <- signif(pp.df.abf*100,3)
  print(signif(pp.df.abf,3))
  print(paste("PP abf for shared variant: ", signif(pp.df.abf["PP.H4.abf"],3)*100 , '%', sep=''))
  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.df.abf)
  
  output<-list(results, merged.df)
  return(output)
}
##' Bayesian colocalisation analysis using data.frames
##'
##' Converts genetic data to snpStats objects, generates p values via score tests, then runs \code{\link{coloc.config.summary}}
##' 
##' @title Bayesian colocalisation analysis using data.frames
##' @param df1, df2 datasets 1 and 2
##' @param snps col.names for snps
##' @param response1 col.name for response in dataset 1
##' @param response2 col.name for response in dataset 2
##' @param ... parameters passed to \code{\link{coloc.config.summary}}
##' @return output of \code{\link{coloc.config.summary}}
##' @author Chris Wallace
coloc.config.datasets <- function(df1,df2,snps=intersect(setdiff(colnames(df1),response1),
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
  coloc.config.snpStats(X1,X2,df1[,response1], df2[,response2], ...)
}
##' Bayesian colocalisation analysis using snpStats objects
##'
##' Generates p values via score tests, then runs \code{\link{coloc.config.summary}}
##' @title Bayesian colocalisation analysis using snpStats objects
##' @param X1 genetic data for dataset 1
##' @param X2 genetic data for dataset 2
##' @param Y1 response for dataset 1
##' @param Y2 response for dataset 2
##' @param snps optional subset of snps to use
##' @param ... parameters passed to \code{\link{coloc.config.summary}}
##' @return output of \code{\link{coloc.config.summary}}
##' @author Chris Wallace
coloc.config.snpStats <- function(X1,X2,Y1,Y2,snps=intersect(colnames(X1),colnames(X2)), ...) {
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
  
  coloc.config(merged.df=data.frame(pvalues.df1=p1,pvalues.df2=p2,MAF=maf),
                 N.df1=nrow(X1), N.df2=nrow(X2), ...)
}

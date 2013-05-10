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
lH3new.f <- function(x, y) {                                  
  my.max <- max(c(x,y))                             
  my.res <- my.max + log( exp(x - my.max) - 0.001*exp(y - my.max))     
  return(my.res)
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
##' @param prior1 p_1 == p_2, prior probability a SNP is associated with one trait
##' @param prior2 p_12, prior probability a SNP is associated with both traits
##' @param sd.prior standard deviation of prior. TODO! CHANGE THIS TO ALLOW SEPARATE PRIORS FOR EACH TRAIT?
##' @param col.names column names in merged.df that correspond to the
##' pvalues in datasets 1 and 2 and the SNP minor allele frequencies.
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item results is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item merged.df is an annotated version of the input \code{data.frame}
##' }
##' @author Claudia Giambartolomei, Chris Wallace
coloc.bayesian.summary <- function(merged.df, N.df1, N.df2, prior1= log(10^(-4)), prior2= log(10^(-5)), s=0,
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
  case_control <- s!=0 
  if (case_control) {     # s=0 if not a case-control studies
    sd.prior = 0.20
    merged.df$V.df1 <- Var.data.cc(f, N.df1, s) 
    merged.df$V.df2 <- Var.data.cc(f, N.df2, s) 
   } else {
      sd.prior = 0.15
      merged.df$V.df1 <- Var.data(f, N.df1) 
      merged.df$V.df2 <- Var.data(f, N.df2) 
    }
  
  merged.df$z.df1 <- qnorm(0.5 * (pvalues.df1), lower.tail = FALSE)  
  ## Shrinkage factor: ratio of the prior variance to the total variance
  merged.df$r.df1 <- sd.prior^2 / (sd.prior^2 + merged.df$V.df1)
  ## Approximate BF  # I want ln scale to compare in log natural scale with LR diff
  merged.df$lABF.df1 = 0.5 * (log(1-merged.df$r.df1) + (merged.df$r.df1 * merged.df$z.df1^2))    
    
  ## Do the same for other dataset:
  merged.df$z.df2 <- qnorm(0.5 * (pvalues.df2), lower.tail = FALSE)      
  ## Shrinkage factor: ratio of the prior variance to the total variance
  merged.df$r.df2 <- sd.prior^2 / (sd.prior^2 + merged.df$V.df2)
  ## Approximate BF
  merged.df$lABF.df2 = 0.5 * (log(1-merged.df$r.df2) + (merged.df$r.df2 * merged.df$z.df2^2))   
  
  merged.df$internal.sum.lABF <- (merged.df$lABF.df1 + merged.df$lABF.df2) 
   
############################## 

  lH0.abf <- 0
  lH1.abf <- prior1 + logmean(merged.df$lABF.df1)
  lH2.abf <- prior1 + logmean(merged.df$lABF.df2)
  lH3.abf <- prior1 + logmean(merged.df$lABF.df1) + prior1 + logmean(merged.df$lABF.df2)
  lH4.abf <- prior2 + logmean(merged.df$internal.sum.lABF)

  lH3new.abf = lH3new.f(lH3.abf, lH4.abf)
  ## If x=(0.001*y), the difference (x-(0.001*y)) = 0, log(0) = -Inf, so fix this:
  ## If it is NaN and the difference between lH3.abf and (log(0.001) + lH4.abf) is very small (<0.001), then keep the old value of lH3.abf:
  if (is.na(lH3new.abf) & abs( lH3.abf -(log(0.001) + lH4.abf) ) <0.001)  (lH3new.abf = lH3.abf)
  
  my.denom.log.abf <- logmean(c(lH0.abf, lH1.abf, lH2.abf, lH3new.abf, lH4.abf))
  
#### Now we can compute the PP under each looping through all models as numerator: 
  all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3new.abf, lH4.abf)
  pp.df.abf <- exp(all.abf - my.denom.log.abf) 
  names(pp.df.abf) <- paste('PP.H', 0:length(pp.df.abf), '.abf', sep='')
  
  ##  pp.df.abf <- signif(pp.df.abf*100,3)
  print(signif(pp.df.abf*100,3))
  print(paste("PP abf for shared variant: ", signif(pp.df.abf$PP.H4.abf,3) , '%', sep=''))
  
  common.snps <- nrow(merged.df)
  results <- c(nsnps=common.snps, pp.df.abf)
  
  output<-list(results, merged.df)
  ##output <- list(summary.df)
  return(output)
}

coloc.bayesian.datasets <- function(df1,df2,snps=intersect(setdiff(colnames(df1),response1),
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
  p1 <- p.value(single.snp.tests(phenotype=df1[,response1], snp.data=X1),df=1)
  p2 <- p.value(single.snp.tests(phenotype=df1[,response1], snp.data=X1),df=1)
  maf <- col.summary(rbind(X1,X2))[,"MAF"]
  
  coloc.bayesian(merged.df=data.frame(pvalues.df1=p1,pvalues.df2=p2,MAF=maf),
                 N.df1=nrow(X1), N.df2=nrow(X2), ...)
}

coloc.bayesian.snpStats <- function(X1,X2,Y1,Y2,snps=intersect(colnames(X1),colnames(X2)), ...) {
  if(!is(X1,"SnpMatrix") || !is(X2,"SnpMatrix"))
    stop("X1 and X2 must be SnpMatrices")
  if(length(Y1) != nrow(X1) || length(Y2) != nrow(X2))
    stop("length(Y1) != nrow(X1) || length(Y2) != nrow(X2)")
  if(length(snps)<2)
    stop("require at least two SNPs in common between X1 and X2 to do anything sensible")
  X1 <- X1[,snps]
  X2 <- X2[,snps]
  X2 <- new("SnpMatrix",as.matrix(df2[,snps]))
  p1 <- p.value(single.snp.tests(phenotype=df1[,response1], snp.data=X1),df=1)
  p2 <- p.value(single.snp.tests(phenotype=df1[,response1], snp.data=X1),df=1)
  maf <- col.summary(rbind(X1,X2))[,"MAF"]
  
  coloc.bayesian(merged.df=data.frame(pvalues.df1=p1,pvalues.df2=p2,MAF=maf),
                 N.df1=nrow(X1), N.df2=nrow(X2), ...)
}

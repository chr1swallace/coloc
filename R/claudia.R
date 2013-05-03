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

##' This function calculates the log of the sum of the exponentiated logs taking out the max, i.e. insuring that the sum is not Inf
##'
##' Internal function
##' @title logmean
##' @param x numeric vector
##' @return max(x) + log(sum(exp(x - max(x))))
##' @author Claudia Giambartolomei
logmean <- function(x) {
  my.max <- max(x)                              ##take out the maximum value in log form
  my.res <- my.max + log(sum(exp(x - my.max ))) 
  return(my.res)
}


# This is to correct the posterior probability of H3 so that it does not include H4 (same SNP)
# lH3new <- log(exp(lH3)-0.001*exp(lH4))
lH3new.f <- function(x, y) {                                  
  my.max <- max(c(x,y))                             
  my.res <- my.max + log( exp(x - my.max) - 0.001*exp(y - my.max))     
  return(my.res)
}



## This function takes a data frame obtained by merging p-values for both eQTL and biomarker dataset and returns a list with [1] summary df [2] original df with additional ABF and other values
# Using MAF from eQTL dataset (column named "MAF.df2")
# "pvalues.df1" and "pvalues.df2" : names of the colums with p-values
# N.df1 and N.df2 number of indviduals used to get the p-values in each dataset
# sd.prior = standard deviation of prior
abf.df <- function(merged.df, N.df1, N.df2, prior1= log(10^(-4)), prior2= log(10^(-5)), s=0 ) {  # set s to 0 is it is not a case control study

	merged.df <- subset(merged.df, merged.df$pvalues.df1 !=0)

	f = merged.df$MAF.df2
	pvalues.df1 = merged.df$pvalues.df1
	pvalues.df2 = merged.df$pvalues.df2

        ####### Use different priors and different computation of variance of the mle for case/control vs. quantitative trait
        if (s==0)  (case_control = FALSE) else (case_control = TRUE)
        if (case_control) {     # s=0 if not a case-control studies
            sd.prior = 0.20
 	    merged.df$V.df1 <- Var.data.cc(f, N.df1, s) 
 	    merged.df$V.df2 <- Var.data.cc(f, N.df2, s) } else {
            sd.prior = 0.15
 	    merged.df$V.df1 <- Var.data(f, N.df1) 
 	    merged.df$V.df2 <- Var.data(f, N.df2) 
        }

          
	merged.df$z.df1 <- qnorm(0.5 * (pvalues.df1), lower.tail = FALSE)  
	# Shrinkage factor: ratio of the prior variance to the total variance
	merged.df$r.df1 <- sd.prior^2 / (sd.prior^2 + merged.df$V.df1)
	# Approximate BF  # I want ln scale to compare in log natural scale with LR diff
	merged.df$lABF.df1 = 0.5 * (log(1-merged.df$r.df1) + (merged.df$r.df1 * merged.df$z.df1^2))    


	# Do the same for other dataset:
        merged.df$z.df2 <- qnorm(0.5 * (pvalues.df2), lower.tail = FALSE)      
	# Shrinkage factor: ratio of the prior variance to the total variance
        merged.df$r.df2 <- sd.prior^2 / (sd.prior^2 + merged.df$V.df2)
        # Approximate BF
        merged.df$lABF.df2 = 0.5 * (log(1-merged.df$r.df2) + (merged.df$r.df2 * merged.df$z.df2^2))   


        merged.df$internal.sum.lABF <- (merged.df$lABF.df1 + merged.df$lABF.df2) 


        ############################## 


	lH0.abf <- 0
        lH1.abf <- prior1 + logmean(merged.df$lABF.df1)
        lH2.abf <- prior1 + logmean(merged.df$lABF.df2)
        lH3.abf <- prior1 + logmean(merged.df$lABF.df1) + prior1 + logmean(merged.df$lABF.df2)
        lH4.abf <- prior2 + logmean(merged.df$internal.sum.lABF)


        lH3new.abf = lH3new.f(lH3.abf, lH4.abf)
        # If x=(0.001*y), the difference (x-(0.001*y)) = 0, log(0) = -Inf, so fix this:
        # If it is NaN and the difference between lH3.abf and (log(0.001) + lH4.abf) is very small (<0.001), then keep the old value of lH3.abf:
        if (is.na(lH3new.abf) & abs( lH3.abf -(log(0.001) + lH4.abf) ) <0.001)  (lH3new.abf = lH3.abf)
        
	my.denom.log.abf <- logmean(c(lH0.abf, lH1.abf, lH2.abf, lH3new.abf, lH4.abf))

	#### Now we can compute the PP under each looping through all models as numerator: 
	all.abf <- c(lH0.abf, lH1.abf, lH2.abf, lH3new.abf, lH4.abf)
	pp.df.abf <- data.frame(row.names=1)

	for (num in 1:length(all.abf)) {
             PP.abf <- exp(all.abf[num] - my.denom.log.abf) 
                     pp.df.abf <- cbind(pp.df.abf, PP.abf)
                     names(pp.df.abf)[names(pp.df.abf)=='PP.abf'] <- paste('PP.H', num-1, '.abf', sep='')
                }

	pp.df.abf <- signif(pp.df.abf*100,3)
	print(pp.df.abf)
	print(paste("PP abf for shared variant: ", signif(pp.df.abf$PP.H4.abf,3) , '%', sep=''))


        common.snps <- nrow(merged.df)
	summary.df<- cbind(common.snps, pp.df.abf)
        

#        output<-list(summary.df, merged.df)
        output <- list(summary.df)

        return(output)

        }

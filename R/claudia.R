##' variance of MLE of beta for quantitative trait, assuming var(y)=1
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
    my.res <- my.max + log(max(exp(x - my.max) - exp(y-my.max), 0))
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
##' @param quiet don't print posterior summary if TRUE. default=FALSE
##' @inheritParams coloc.abf
##' @return named numeric vector of posterior probabilities
##' @author Claudia Giambartolomei, Chris Wallace
combine.abf <- function(l1, l2, p1, p2, p12, quiet=FALSE) {
    stopifnot(length(l1)==length(l2))
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
    if(!quiet) {
        print(signif(pp.abf,3))
        print(paste("PP abf for shared variant: ", signif(pp.abf["PP.H4.abf"],3)*100 , '%', sep=''))
    }
    return(pp.abf)
}

combine_abf_weighted <- function(l1, l2, p1, p2, p12,
                                prior_weights1, prior_weights2,
                                quiet = FALSE) {

  stopifnot(length(l1) == length(l2))

  Q <- length(l1)

  if (is.null(prior_weights1)) {
    p1_vec <- rep(p1, Q)
  } else {
    prior_weights1 <- prior_weights1 / sum(prior_weights1)
    p1_vec <- Q * p1 * prior_weights1
  }

  if (is.null(prior_weights2)) {
    p2_vec <- rep(p2, Q)
  } else {
    prior_weights2 <- prior_weights2 / sum(prior_weights2)
    p2_vec <- Q * p2 * prior_weights2
  }

  stopifnot(length(p1_vec) == length(p2_vec))
  stopifnot(length(p1_vec) == length(l1))

  p12_vec <- p1_vec * p2_vec * (p12 / (p1 * p2))

  lsum <- l1 + l2
  lH0_abf <- 0
  lH1_abf <- logsum(log(p1_vec) + l1)
  lH2_abf <- logsum(log(p2_vec) + l2)
  lH3_abf <- logdiff(logsum(log(p1_vec) + l1) + logsum(log(p2_vec) + l2),
                     logsum(log(p1_vec) + log(p2_vec) + lsum))
  lH4_abf <- logsum(log(p12_vec) + lsum)

  all_abf <- c(lH0_abf, lH1_abf, lH2_abf, lH3_abf, lH4_abf)
  denom_log_abf <- logsum(all_abf)
  pp_abf <- exp(all_abf - denom_log_abf)
  names(pp_abf) <- paste0("PP.H", (1:length(pp_abf)) - 1, ".abf")

  if(!quiet) {
    print(signif(pp_abf,3))
    print(paste("PP abf for shared variant: ", signif(pp_abf["PP.H4.abf"],3) * 100 ,
                "%", sep = ""))
  }

  pp_abf
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
    warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
    oneover <- 1/vbeta
    nvx <- 2 * n * maf * (1-maf)
    m <- lm(nvx ~ oneover - 1)
    cf <- coef(m)[['oneover']]
    if(cf < 0)
        stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
    return(sqrt(cf))
}

##' Internal function, process each dataset list for coloc.abf.
##' 
##' Made public for another package to use, but not intended for users to use.
##'
##' @title process.dataset
##' @param d list
##' @param suffix "df1" or "df2"
##' @return data.frame with log(abf) or log(bf)
##' @export
##' @author Chris Wallace
process.dataset <- function(d, suffix) {
                                        #message('Processing dataset')

    nd <- names(d)
    ## if (! 'type' %in% nd)
    ##   stop("dataset ",suffix,": ",'The variable type must be set, otherwise the Bayes factors cannot be computed')

    ## if(!(d$type %in% c("quant","cc")))
    ##     stop("dataset ",suffix,": ","type must be quant or cc")
    
    ## if(d$type=="cc" & "pvalues" %in% nd) {
    ## if(!( "s" %in% nd))
    ##     stop("dataset ",suffix,": ","please give s, proportion of samples who are cases, if using p values")
    ## if(!("MAF" %in% nd))
    ##     stop("dataset ",suffix,": ","please give MAF if using p values")
    ## if(d$s<=0 || d$s>=1)
    ##     stop("dataset ",suffix,": ","s must be between 0 and 1")
    ## }
    
    ## if(d$type=="quant") {
    ##     if(!("sdY" %in% nd || ("MAF" %in% nd && "N" %in% nd )))
    ##         stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")
    ## }
    
    if("beta" %in% nd && "varbeta" %in% nd) {  ## use beta/varbeta.  sdY should be estimated by now for quant
        ## if(length(d$beta) != length(d$varbeta))
        ##   stop("dataset ",suffix,": ","Length of the beta vectors and variance vectors must match")
        ## if(!("snp" %in% nd))
        ##   d$snp <- sprintf("SNP.%s",1:length(d$beta))
        ## if(length(d$snp) != length(d$beta))
        ##   stop("dataset ",suffix,": ","Length of snp names and beta vectors must match")
        
        if(d$type=="quant" && !('sdY' %in% nd)) 
            d$sdY <- sdY.est(d$varbeta, d$MAF, d$N)
        df <- approx.bf.estimates(z=d$beta/sqrt(d$varbeta),
                                  V=d$varbeta, type=d$type, suffix=suffix, sdY=d$sdY)
        df$snp <- as.character(d$snp)
        if("position" %in% nd)
            df <- cbind(df,position=d$position)
        return(df)
    }

    if("pvalues" %in% nd & "MAF" %in% nd & "N" %in% nd) { ## no beta/varbeta: use p value / MAF approximation
        ## if (length(d$pvalues) != length(d$MAF))
        ##   stop('Length of the P-value vectors and MAF vector must match')
        ## if(!("snp" %in% nd))
        ##   d$snp <- sprintf("SNP.%s",1:length(d$pvalues))
        df <- data.frame(pvalues = d$pvalues,
                         MAF = d$MAF,
                         N=d$N,
                         snp=as.character(d$snp))    
        snp.index <- which(colnames(df)=="snp")
        colnames(df)[-snp.index] <- paste(colnames(df)[-snp.index], suffix, sep=".")
        ## keep <- which(df$MAF>0 & df$pvalues > 0) # all p values and MAF > 0
        ## df <- df[keep,]
        abf <- approx.bf.p(p=df$pvalues, f=df$MAF, type=d$type, N=df$N, s=d$s, suffix=suffix)
        df <- cbind(df, abf)
        if("position" %in% nd)
            df <- cbind(df,position=d$position)
        return(df)  
    }

    stop("Must give, as a minimum, one of:\n(beta, varbeta, type, sdY)\n(beta, varbeta, type, MAF)\n(pvalues, MAF, N, type)")
}

##' Bayesian finemapping analysis
##'
##' This function calculates posterior probabilities of different
##' causal variant for a single trait.
##'
##' If regression coefficients and variances are available, it
##' calculates Bayes factors for association at each SNP.  If only p
##' values are available, it uses an approximation that depends on the
##' SNP's MAF and ignores any uncertainty in imputation.  Regression
##' coefficients should be used if available.
##' 
##' @title Bayesian finemapping analysis
##' @param dataset a list with specifically named elements defining the dataset
##'   to be analysed. See \code{\link{check_dataset}} for details.
##'
##' @param p1 prior probability a SNP is associated with the trait 1, default 1e-4
##' @param prior_weights Non-negative weights for the prior probability a SNP is causal 
##' @return a \code{data.frame}:
##' \itemize{
##' \item an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability of the SNP being causal
##' }
##' @author Chris Wallace
##' @export
finemap.abf <- function(dataset, p1=1e-4, prior_weights = NULL) {

  check_dataset(dataset,"")

    df <- process.dataset(d=dataset, suffix="")
    nsnps <- nrow(df)
    p1=adjust_prior(p1,nsnps,"1")

    dfnull <- df[1,]
    for(nm in colnames(df))
        dfnull[,nm] <- NA
    dfnull[,"snp"] <- "null"
    dfnull[,"lABF."] <- 0
    df <- rbind(df,dfnull)
    ## data.frame("V."=NA,
    ##            z.=NA,
    ##            r.=NA,
    ##            lABF.=1,
    ##            snp="null"))
    if (!is.null(prior_weights)) {
      stopifnot(length(prior_weights) == nsnps)
      prior_vec <- p1 * nsnps * prior_weights / sum(prior_weights)
      df$prior <- c(prior_vec, 1 - nsnps * p1)
    } else {
      df$prior <- c(rep(p1,nsnps),1-nsnps*p1)
    }

    ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
    ## BUGFIX 16/5/19
    ## my.denom.log.abf <- logsum(df$lABF + df$prior)
    ## df$SNP.PP <- exp(df$lABF - my.denom.log.abf)
    my.denom.log.abf <- logsum(df$lABF + log(df$prior))
    df$SNP.PP <- exp(df$lABF + log(df$prior) - my.denom.log.abf)
    
    return(df)
}

adjust_prior=function(p,nsnps,suffix="") {
    if(nsnps * p >= 1) { ## for very large regions
        warning(paste0("p",suffix," * nsnps >= 1, setting p",suffix,"=1/(nsnps + 1)"))
        1/(nsnps + 1)
    } else {
        p
    }
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
##' @param dataset1 a list with specifically named elements defining the dataset
##'   to be analysed. See \code{\link{check_dataset}} for details.
##' @param dataset2 as above, for dataset 2
##' @param MAF Common minor allele frequency vector to be used for both dataset1 and dataset2, a shorthand for supplying the same vector as parts of both datasets
##' @param p1 prior probability a SNP is associated with trait 1, default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2, default 1e-4
##' @param p12 prior probability a SNP is associated with both traits, default 1e-5
##' @param prior_weights1 Non-negative weights for the prior probability a SNP is associated with trait 1 
##' @param prior_weights2 Non-negative weights for the prior probability a SNP is asscoiated with trait 2
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item summary is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item results is an annotated version of the input data containing log Approximate Bayes Factors and intermediate calculations, and the posterior probability SNP.PP.H4 of the SNP being causal for the shared signal *if* H4 is true. This is only relevant if the posterior support for H4 in summary is convincing.
##' }
##' @author Claudia Giambartolomei, Chris Wallace, Jeffrey Pullin
##' @export
coloc.abf <- function(dataset1, dataset2, MAF=NULL, 
                      p1=1e-4, p2=1e-4, p12=1e-5,
                      prior_weights1 = NULL, prior_weights2 = NULL) {

    if(!("MAF" %in% names(dataset1)) & !is.null(MAF))
        dataset1$MAF <- MAF
    if(!("MAF" %in% names(dataset2)) & !is.null(MAF))
        dataset2$MAF <- MAF
    check_dataset(d=dataset1,1)
    check_dataset(d=dataset2,2)
    
    df1 <- process.dataset(d=dataset1, suffix="df1")
    df2 <- process.dataset(d=dataset2, suffix="df2")
    p1=adjust_prior(p1,nrow(df1),"1")
    p2=adjust_prior(p2,nrow(df2),"2")

    merged.df <- merge(df1,df2)
    p12=adjust_prior(p12,nrow(merged.df),"12")

    prior_weights1 <- prior_weights1[which(df1$snp %in% merged.df$snp)]
    prior_weights2 <- prior_weights2[which(df2$snp %in% merged.df$snp)]

    if(!nrow(merged.df))
        stop("dataset1 and dataset2 should contain the same snps in the same order, or should contain snp names through which the common snps can be identified")

    merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
    ## add SNP.PP.H4 - post prob that each SNP is THE causal variant for a shared signal
    my.denom.log.abf <- logsum(merged.df$internal.sum.lABF)
    merged.df$SNP.PP.H4 <- exp(merged.df$internal.sum.lABF - my.denom.log.abf)

    is_weighted <- !is.null(prior_weights1) || !is.null(prior_weights2)
    if (is_weighted) {
      pp.abf <- combine_abf_weighted(merged.df$lABF.df1, merged.df$lABF.df2,
                                     p1, p2, p12, prior_weights1, prior_weights2)
    } else {
      pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12) 
    }
    
    common.snps <- nrow(merged.df)
    results <- c(nsnps=common.snps, pp.abf)
    output<-list(summary=results,
                 results=merged.df,
                 priors=c(p1=p1,p2=p2,p12=p12))

    if (is_weighted) {
      output$weights <- list(
        prior_weights1 = prior_weights1,
        prior_weights2 = prior_weights2
      )
    }

    class(output) <- c("coloc_abf",class(output))
    return(output)
}

##' Get credible sets from finemapping results
##'
##' @title credible.sets
##' @param datasets data.frame output of `finemap.abf()`
##' @param credible.size threshold of the credible set (Default: 0.95)
##' @return SNP ids of the credible set
##' @author Guillermo Reales, Chris Wallace
##' @export 
credible.sets <- function(dataset, credible.size = 0.95){
    if(!"SNP.PP" %in% names(dataset)) stop("Input must be finemap.abf() output and have a SNP.PP column.")
    t2 <- dataset[ order(dataset$SNP.PP, decreasing = TRUE),]
    t2$csum <- cumsum(t2$SNP.PP)
    w=which(cumsum(t2$SNP.PP)>=credible.size)[1]
    t2[ 1:w ,c("snp","SNP.PP")]
}


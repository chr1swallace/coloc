##' Internal function, approx.bf.imputed
##'
##' Calculate approximate Bayes Factors using supplied variance of the regression coefficients
##' @title Internal function, approx.bf.imputed
##' @param z normal deviate associated with regression coefficient and its variance
##' @param V its variance
##' @inheritParams approx.bf
##' @return data.frame containing lABF and intermediate calculations
##' @author Vincent Plagnol, Chris Wallace
approx.bf.imputed <- function (z, V, type, suffix) {
  if (type == "quant") {sd.prior <- 0.15} else {sd.prior <- 0.2}
  r <- sd.prior^2/(sd.prior^2 + V)
  lABF = 0.5 * (log(1 - r) + (r * z^2))
  ret <- data.frame(V, z, r, lABF)
  colnames(ret) <- paste(colnames(ret), suffix, sep = ".")
  return(ret)
}


##' Bayesian colocalisation analysis using summary coefficients and their variances
##'
##' This function makes a data frame obtained by merging p-values for
##' two traits over the same set of SNPs and calculates posterior
##' probabilities of different causal variant configurations under the
##' assumption of a single causal variant for each trait.
##'
##' It uses the variance of the regression coefficients to estimate the 
##' @title Fully Bayesian colocalisation analysis using regression coefficients
##' @inheritParams coloc.abf
##' @param beta.dataset1 coefficient from dataset 1
##' @param beta.dataset2  coefficient from dataset 2
##' @param varbeta.dataset1 variance of the coefficient from dataset 1
##' @param varbeta.dataset2  variance of the coefficient from dataset 2
##' @return a list of two \code{data.frame}s:
##' \itemize{
##' \item results is a vector giving the number of SNPs analysed, and the posterior probabilities of H0 (no causal variant), H1 (causal variant for trait 1 only), H2 (causal variant for trait 2 only), H3 (two distinct causal variants) and H4 (one common causal variant)
##' \item merged.df is an annotated version of the input \code{data.frame}
##' }
##' @author Vincent Plagnol, Chris Wallace
##' @export
coloc.abf.imputed <- function (beta.dataset1, beta.dataset2, varbeta.dataset1, varbeta.dataset2, type.dataset1, type.dataset2,
                                p1 = 1e-04, p2 = 1e-04, p12 = 1e-05) {

  if (length(beta.dataset1) != length(beta.dataset2)) stop("Length of the beta vectors must match")
  if (length(beta.dataset1) != length(varbeta.dataset1)) stop("Length of the beta vectors and variance vectors must match")
  if (length(beta.dataset1) != length(varbeta.dataset2)) stop("Length of the beta vectors and variance vectors must match")

  merged.df <- data.frame(z.df1 = beta.dataset1/sqrt(varbeta.dataset1),
                          z.df2 = beta.dataset2/sqrt(varbeta.dataset2),
                          V.df1 = varbeta.dataset1,
                          V.df2 = varbeta.dataset2)
  abf.df1 <- approx.bf.imputed(merged.df$z.df1, merged.df$V.df1, type.dataset1, suffix = "df1")
  abf.df2 <- approx.bf.imputed(merged.df$z.df2, merged.df$V.df2, type.dataset2, suffix = "df2")
  merged.df <- cbind(merged.df, abf.df1, abf.df2)
  merged.df$internal.sum.lABF <- with(merged.df, lABF.df1 + lABF.df2)
  pp.abf <- combine.abf(merged.df$lABF.df1, merged.df$lABF.df2, p1, p2, p12)
  common.snps <- nrow(merged.df)
  results <- c(nsnps = common.snps, pp.abf)
  output <- list(results, merged.df)
  return(output)
}

\name{coloc.abf.imputed}
\alias{coloc.abf.imputed}
\title{Fully Bayesian colocalisation analysis using regression coefficients}
\usage{
  coloc.abf.imputed(beta.dataset1, beta.dataset2,
    varbeta.dataset1, varbeta.dataset2, type.dataset1,
    type.dataset2, p1 = 1e-04, p2 = 1e-04, p12 = 1e-05)
}
\arguments{
  \item{beta.dataset1}{coefficient from dataset 1}

  \item{beta.dataset2}{coefficient from dataset 2}

  \item{varbeta.dataset1}{variance of the coefficient from
  dataset 1}

  \item{varbeta.dataset2}{variance of the coefficient from
  dataset 2}

  \item{p1}{prior probability a SNP is associated with
  trait 1}

  \item{p2}{prior probability a SNP is associated with
  trait 2}

  \item{p12}{prior probability a SNP is associated with
  both traits}
}
\value{
  a list of two \code{data.frame}s: \itemize{ \item results
  is a vector giving the number of SNPs analysed, and the
  posterior probabilities of H0 (no causal variant), H1
  (causal variant for trait 1 only), H2 (causal variant for
  trait 2 only), H3 (two distinct causal variants) and H4
  (one common causal variant) \item merged.df is an
  annotated version of the input \code{data.frame} }
}
\description{
  Bayesian colocalisation analysis using summary
  coefficients and their variances
}
\details{
  This function makes a data frame obtained by merging
  p-values for two traits over the same set of SNPs and
  calculates posterior probabilities of different causal
  variant configurations under the assumption of a single
  causal variant for each trait.

  It uses the variance of the regression coefficients to
  estimate the
}
\author{
  Vincent Plagnol, Chris Wallace
}


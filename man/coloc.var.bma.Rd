\name{coloc.var.bma}
\alias{coloc.var.bma}
\title{Wrapper to use colocalization testing within a Bayesian model averaging
structure for datasets with common controls.}
\usage{
  coloc.var.bma(df1,
    snps = setdiff(colnames(df1), response),
    response = "Y", bayes = !is.null(bayes.factor),
    thr = 0.01, nsnps = 2, n.approx = 1001,
    bayes.factor = NULL, plot.coeff = FALSE,
    r2.trim = 0.95, quiet = FALSE, ...)
}
\arguments{
  \item{df1}{A dataframe, containing response and potential
  explanatory variables for the dataset.}

  \item{snps}{The SNPs to consider as potential explanatory
  variables}

  \item{response}{The name of the response variable in
  \code{df1}}

  \item{thr}{posterior probability threshold used to trim
  SNP list.  Only SNPs with a marginal posterior
  probability of inclusion greater than this with one or
  other trait will be included in the full BMA analysis}

  \item{nsnps}{number of SNPs required to model both
  traits.  The BMA analysis will average over all possible
  \code{nsnp} SNP models, subject to \code{thr} above.}

  \item{n.approx}{number of values at which to numerically
  approximate the posterior}

  \item{r2.trim}{for pairs SNPs with r2> \code{r2.trim},
  only one SNP will be retained.  This avoids numerical
  instability problems caused by including two highly
  correlated SNPs in the model.}

  \item{quiet}{suppress messages about how the model spaced
  is trimmed for BMA}

  \item{...}{other parameters passed to \code{coloc.test}}

  \item{bayes}{Logical, indicating whether to perform
  Bayesian inference for the coefficient of
  proportionality, eta.  If \code{bayes.factor} is
  supplied, Bayes factors are additionally computed for the
  specificed values.  This can add a little time as it
  requires numerical integration, so can be set to FALSE to
  save time in simulations, for example.}

  \item{bayes.factor}{Calculate Bayes Factors to compare
  specific values of eta.  \code{bayes.factor} should
  either a numeric vector, giving single value(s) of
  \code{eta} or a list of numeric vectors, each of length
  two and specifying ranges of eta which should be compared
  to each other.  Thus, the vector or list needs to have
  length at least two.}

  \item{plot.coeff}{\code{TRUE} if you want to generate a
  plot showing the coefficients from the two regressions
  together with confidence regions.}
}
\value{
  a \code{coloc} or \code{colocBayes} object
}
\description{
  This is a test for proportionality of regression
  coefficients from a multinomial regression. Analysis is
  based on taking a hybrid-Bayesian approach and
  integrating the p value over the posterior distribution
  of \code{eta}, which gives a posterior predictive p
  value.  The Bayesian approach can also be used to give a
  credible interval for \code{eta}.
}
\author{
  Mary Fortune
}


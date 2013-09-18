\name{coloc.var.test}
\alias{coloc.var.test}
\title{Wrapper to use colocalization testing within a principle components
structure for datasets with common controls.}
\usage{
  coloc.var.test(X, vars.drop = NULL, ...)
}
\arguments{
  \item{X}{A glm object for the traits.  Any Intercept term
  is dropped, but other covariates should have distinct
  names or be listed in \code{vars.drop} to avoid them
  being included in the colocalisation test.}

  \item{vars.drop}{Character vector naming additional
  variables in either regression which are not SNPs and
  should not be used in the colocalisation test.  They
  should appear in \code{c(names(coefficients(X)))}}

  \item{...}{other arguments passed to
  \code{\link{coloc.test.summary}()}.}
}
\value{
  a numeric vector with 3 named elements:
  \item{eta.hat}{The estimated slope.} \item{chisquare}{The
  chisquared test statistic} \item{n}{The number of snps
  used in the test.  If eta were known, this would be the
  degrees of freedom of the test. Because eta has been
  replaced by its ML estimate, Plagnol et al suggest we
  expect the degrees of freedom to be n-1, but this
  requires the likelihood to be well behaved which is not
  always the case.  We prefer to consider the posterior
  predictive p value.} \item{ppp}{The posterior predictive
  p value}
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


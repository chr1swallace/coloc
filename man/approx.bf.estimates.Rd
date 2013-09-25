\name{approx.bf.estimates}
\alias{approx.bf.estimates}
\title{Internal function, approx.bf.estimates}
\usage{
  approx.bf.estimates(z, V, type, suffix = NULL, sdY = 1)
}
\arguments{
  \item{z}{normal deviate associated with regression
  coefficient and its variance}

  \item{V}{its variance}

  \item{sdY}{standard deviation of the trait. If not
  supplied, will be estimated.}

  \item{type}{"quant" or "cc"}

  \item{suffix}{suffix to append to column names of
  returned data.frame}
}
\value{
  data.frame containing lABF and intermediate calculations
}
\description{
  Internal function, approx.bf.estimates
}
\details{
  Calculate approximate Bayes Factors using supplied
  variance of the regression coefficients
}
\author{
  Vincent Plagnol, Chris Wallace
}


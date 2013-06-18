\name{approx.bf.imputed}
\alias{approx.bf.imputed}
\title{Internal function, approx.bf.imputed}
\usage{
  approx.bf.imputed(z, V, type, suffix = NULL)
}
\arguments{
  \item{z}{normal deviate associated with regression
  coefficient and its variance}

  \item{V}{its variance}

  \item{type}{"quant" or "cc"}

  \item{suffix}{suffix to append to column names of
  returned data.frame}
}
\value{
  data.frame containing lABF and intermediate calculations
}
\description{
  Internal function, approx.bf.imputed
}
\details{
  Calculate approximate Bayes Factors using supplied
  variance of the regression coefficients
}
\author{
  Vincent Plagnol, Chris Wallace
}


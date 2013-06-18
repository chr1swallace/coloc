\name{approx.bf}
\alias{approx.bf}
\title{Internal function, approx.bf}
\usage{
  approx.bf(p, f, type, N, s, suffix)
}
\arguments{
  \item{p}{p value}

  \item{f}{MAF}

  \item{type}{"quant" or "cc"}

  \item{N}{sample size}

  \item{s}{proportion of samples that are cases, ignored if
  type=="quant"}

  \item{suffix}{suffix to append to column names of
  returned data.frame}
}
\value{
  data.frame containing lABF and intermediate calculations
}
\description{
  Internal function, approx.bf
}
\details{
  Calculate approximate Bayes Factors
}
\author{
  Claudia Giambartolomei, Chris Wallace
}


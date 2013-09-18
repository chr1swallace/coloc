\name{my.logsum}
\alias{my.logsum}
\title{my.logsum}
\usage{
  my.logsum(x)
}
\arguments{
  \item{x}{numeric vector}
}
\value{
  my.max + log(sum(exp(x - my.max )))
}
\description{
  Internal function, my.logsum
}
\details{
  This function calculates the log of the sum of the
  exponentiated logs taking out the max, i.e. insuring that
  the sum is not Inf This sum is weighted by some constants
  w
}
\author{
  Mary Fortune
}


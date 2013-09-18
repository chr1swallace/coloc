\name{wlogsum}
\alias{wlogsum}
\title{wlogsum}
\usage{
  wlogsum(x, w = NULL)

  wlogsum(x, w = NULL)
}
\arguments{
  \item{x}{numeric vector}

  \item{w}{numeric vector}

  \item{x}{numeric vector}

  \item{w}{numeric vector}
}
\value{
  my.max + log(sum(exp(x - my.max )*w))

  my.max + log(sum(exp(x - my.max )*w))
}
\description{
  Internal function, wlogsum

  Internal function, wlogsum
}
\details{
  This function calculates the log of the sum of the
  exponentiated logs taking out the max, i.e. insuring that
  the sum is not Inf This sum is weighted by some constants
  w

  This function calculates the log of the sum of the
  exponentiated logs taking out the max, i.e. insuring that
  the sum is not Inf This sum is weighted by some constants
  w
}
\author{
  Chris Wallace

  Chris Wallace
}


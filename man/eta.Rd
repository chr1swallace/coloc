\name{eta}
\alias{bf}
\alias{bf,colocBayes-method}
\alias{chisquare}
\alias{chisquare,coloc-method}
\alias{ci}
\alias{ci,colocBayes-method}
\alias{df}
\alias{df,coloc-method}
\alias{eta}
\alias{eta,coloc-method}
\alias{p.value}
\alias{p.value,coloc-method}
\alias{plot,colocPCs,missing-method}
\alias{ppp.value}
\alias{ppp.value,colocBayes-method}
\alias{show,coloc-method}
\alias{theta}
\alias{theta,coloc-method}
\title{Functions to extract information from a \code{coloc} or \code{colocBayes}
object}
\arguments{
  \item{object}{Object returned by \code{coloc.test()}
  function.}
}
\description{
  Extract information from a \code{coloc} object.
}
\details{
  \code{eta()} returns eta.hat, the maximum likelihood
  value of eta.

  \code{theta()} returns theta.hat, the maximum likelihood
  value of eta.

  \code{lhood()} returns -2 times the log-likelihood
  ratio..

  \code{chisquare()} returns the value of the chisquare
  statistic derived by Wallace et al (in preparation).

  \code{df()} returns the associated degrees of freedom.

  \code{p.value()} returns the associated p value.

  \code{ppp.value()} returns the posterior predicted p
  value, or \code{NA} if not calculated.

  \code{ci()} returns the credible interval, or \code{NA}
  if not calculated.
}
\author{
  Chris Wallace.
}
\seealso{
  \code{\link{coloc.test}}, \code{\link{pcs.prepare}}
}
\keyword{methods}


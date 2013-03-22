\name{Methods for coloc}
\alias{eta}
\alias{theta}
\alias{df}
\alias{ppp.value}
\alias{p.value}
\alias{chisquare}
\alias{df,coloc-method}
\alias{n,coloc-method}
\alias{eta,coloc-method}
\alias{ppp.value,coloc-method}
\alias{p.value,coloc-method}
\alias{show,coloc-method}
\alias{chisquare,coloc-method}
\alias{theta,coloc-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to extract information from a \code{coloc} object
}
\usage{
\S4method{eta}{coloc}(object)
}
\description{
  Extract information from a \code{coloc} object.
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{ Object returned by \code{coloc.test()} function.}
}
\details{
  eta() returns eta.hat, the maximum likelihood value of eta.
  
  theta() returns theta.hat, the maximum likelihood value of eta.
  
  lhood() returns -2 times the log-likelihood ratio..
  
  chisquare() returns the value of the chisquare statistic derived by
  Wallace et al (in preparation).
  
  df() returns the associated degrees of freedom.
  
  p.value() returns the associated p value.
  
  ppp.value() returns the posterior predicted p value.
  
}

\author{
  Chris Wallace.
}

\seealso{\code{\link{coloc.test}}}
%% ~Make other sections like Warning with \section{Warning }{....} ~

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

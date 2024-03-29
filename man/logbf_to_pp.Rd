% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/susie.R
\name{logbf_to_pp}
\alias{logbf_to_pp}
\title{logbf 2 pp}
\usage{
logbf_to_pp(bf, pi, last_is_null)
}
\arguments{
\item{bf}{an L by p or p+1 matrix of log Bayes factors}

\item{pi}{\emph{either} a scalar representing the prior probability for any snp
to be causal, \emph{or} a full vector of per snp / null prior probabilities}

\item{last_is_null}{TRUE if last value of the bf vector or last column of a
bf matrix relates to the null hypothesis of no association. This is
standard for SuSiE results, but may not be for BF constructed in other
ways.}
}
\value{
matrix of posterior probabilities, same dimensions as bf
}
\description{
generic convenience function to convert logbf matrix to PP matrix
}
\author{
Chris Wallace
}

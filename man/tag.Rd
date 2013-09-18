\name{tag}
\alias{tag}
\title{tag}
\usage{
  tag(X, tag.threshold = 0.99, snps = NULL, samples = NULL)
}
\arguments{
  \item{X}{a SnpMatrix}

  \item{snps}{colnames of the SnpMatrix object to be used}

  \item{tag.threshold}{threshold to cut tree, default=0.99}

  \item{samples}{optional, subset of samples to use}
}
\value{
  character vector, names are \code{snps}, values are the
  tag for each SNP
}
\description{
  Derive tag SNPs for a SnpMatrix object using heirarchical
  clustering
}
\details{
  Uses complete linkage and the \code{\link{hclust}}
  function to define clusters, then cuts the tree at
  1-tag.threshold
}
\author{
  Chris Wallace
}


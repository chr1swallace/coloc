\name{pcs.var.prepare}
\alias{pcs.var.prepare}
\title{Prepares principal components of dataset with common control for
colocalisation testing.}
\usage{
  pcs.var.prepare(X)
}
\arguments{
  \item{X}{Either a SnpMatrix or numeric matrix of genetic
  data. Columns index SNPs, rows index samples.}
}
\value{
  a \code{colocPCs} object.
}
\description{
  If \code{X1} is a \code{SnpMatrix} object, it are checked
  for missing data, and any missing values imputed by
  repeated use of \code{impute.snps} from the
  \code{snpStats} package.
}
\details{
  Principal components are calculated using \code{prcomp}.

  \code{pcs.var.model} can then be invoked to create
  \code{glm} objects.
}
\author{
  Mary Fortune
}


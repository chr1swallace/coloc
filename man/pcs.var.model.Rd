\name{pcs.var.model}
\alias{pcs.var.model}
\title{Prepares models of response based on principal components of two datasets for
colocalisation testing.}
\usage{
  pcs.var.model(object, Y, threshold = 0.8)
}
\arguments{
  \item{object}{A colocPCs object, result of
  \code{pcs.var.prepare()}.}

  \item{Y}{Numeric phenotype vector}

  \item{threshold}{The minimum number of principal
  components which captures at least threshold proportion
  of the variance will be selected.  Simulations suggest
  \code{threshold=0.8} is a good default value.}
}
\value{
  \code{pcs.prepare} returns a \code{colocPCs} object.
  \code{pcs.model} returns a \code{glm} object.
}
\description{
  Prepares models of response based on principal components
  of two datasets for colocalisation testing.
}
\author{
  Mary Fortune
}


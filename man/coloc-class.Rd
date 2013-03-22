\docType{class}
\name{coloc-class}
\alias{coloc-class}
\alias{colocBayes-class}
\title{Classes \code{"coloc"} and \code{"colocBayes"}}
\description{
  %% ~~ A concise (1-5 lines) description of what the class
  is. ~~ Classes designed to hold objects returned by
  function \code{\link{coloc.test}} which performs a test
  of the null hypothesis that two genetic traits colocalise
  - that they share a common causal variant.
}
\section{Objects from the Class}{
  Objects can be created by calls to the function
  \code{\link{coloc.test}()}.  Class \code{colocBayes}
  extends class \code{coloc}. %% ~~ describe objects here
  ~~
}
\examples{
showClass("coloc")
showClass("colocBayes")
}
\author{
  Chris Wallace.
}
\references{
  Wallace et al (2012).  Statistical colocalisation of
  monocyte gene expression and genetic risk variants for
  type 1 diabetes.  Hum Mol Genet 21:2815-2824.
  \url{http://europepmc.org/abstract/MED/22403184}

  Plagnol et al (2009).  Statistical independence of the
  colocalized association signals for type 1 diabetes and
  RPS26 gene expression on chromosome 12q13. Biostatistics
  10:327-34.
  \url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
}
\seealso{
  \code{\link{coloc.test}}
}
\keyword{classes}


validColoc <- function(object) {
  if(length(object@result) != 3 ||
     !all(c("eta.hat","chisquare","n") %in% names(object@result))) {
    return("result slot in coloc objects should be a named numeric vector with 3 elements: eta.hat, chisquare, n")
  }
  ## if(length(object@u) != nrow(object@V))
  ##   return("dim of V should equal length of u")
}

### colocABF class
#' Class \code{"colocABF"} holds objects returned by the \code{coloc.abf} function
#'
#' Objects can be created by calls to the
#'function \code{\link{coloc.abf}()}.  
#'@name colocABF-class
#'@aliases colocABF-class
#'@docType class
#'@author Chris Wallace.
#'@seealso \code{\link{coloc.abf}}
#'@keywords classes
#'@examples
#'showClass("colocABF")
#'@exportClass colocABF
setClass("colocABF",
         representation(summary="data.frame", results="data.frame"))

#'Classes \code{"coloc"} and \code{"colocBayes"}
#'
#'Classes 
#'designed to hold objects returned by function \code{\link{coloc.test}} which
#'performs a test of the null hypothesis that two genetic traits colocalise -
#'that they share a common causal variant.
#'
#'
#'@name coloc-class
#'@aliases coloc-class colocBayes-class colocBayesBMA colocBMA
#'@docType class
#'@section Objects from the Class:
#' Objects can be created by calls to the
#'function \code{\link{coloc.test}()}.  Class \code{colocBayes} extends class
#'\code{coloc}. 
#'@author Chris Wallace.
#'@seealso \code{\link{coloc.test}}, \code{\link{coloc.test.summary}}, \code{\link{coloc.bma}}
#'@references Wallace et al (2012).  Statistical colocalisation of monocyte
#'gene expression and genetic risk variants for type 1 diabetes.  Hum Mol Genet
#'21:2815-2824.  \url{http://europepmc.org/abstract/MED/22403184}
#'
#'Plagnol et al (2009).  Statistical independence of the colocalized
#'association signals for type 1 diabetes and RPS26 gene expression on
#'chromosome 12q13. Biostatistics 10:327-34.
#'\url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
#'@keywords classes
#'@examples
#'
#'showClass("coloc")
#'showClass("colocBayes")
#'@exportClass coloc
#'@exportClass colocBayes
setClass("coloc",
         representation(result="numeric", method="character", plot.data="list"))
setClass("colocTWAS",
         representation(result="data.frame", method="character", plot.data="list"))
setClass("colocBayes",
         representation(ppp="numeric",
                        credible.interval="list",
                        bayes.factor="numeric"),
         contains="coloc")
setClass("colocBMA",
         representation(bma="numeric"),
         contains="coloc")
setClass("colocBayesBMA",
         representation(bma="numeric"),
         contains="colocBayes")

show.coloc <- function(object) {
  if(object@method=="single") {
     res <- c(object@result,  p.value=p.value(object))
  } else {
    res <- c(object@result[c("eta.hat","n","p")])
    names(res)[3] <- "ppp.value"
  }
  print(res)
}
setMethod("summary","coloc",show.coloc)
setMethod("show","coloc",show.coloc)
setMethod("summary","colocBayes",show.coloc)
setMethod("show","colocBayes",show.coloc)

## eta <- function(object) {
##   object@result["eta.hat"]
## }
## chisquare <- function(object) {
##   object@result["chisquare"]
## }
## df <- function(object) {
##   object@result["df"]
## }
## df <- function(object) {
##   object@result["n"]
## }
## ppp.value <- function(object) {
##   object@result["ppp"]
## }
## p.value <- function(object) {
##   pchisq(object@result["chisquare"],df=object@result["df"],lower.tail=FALSE)
## }




validColocPCs <- function(object) {
  if(ncol(object@pcs) != length(object@vars))
    return("ncol(pc matrix) != length(cvars)")
  if(nrow(object@pcs) != length(object@group[object@use]))
     return("nrow(pc matrix) != length(group)")
}
#'Class \code{"colocPCs"}
#'
#'%% ~~ A concise (1-5 lines) description of what the class is. ~~ Class
#'designed to hold objects returned by function \code{\link{pcs.prepare}} which
#'generates a principal component summary of two genotype matrices in a form
#'suitable for use in the function \code{\link{pcs.model}}.
#'
#'
#'@name colocPCs-class
#'@docType class
#'@section Objects from the Class: Objects can be created by calls to the
#'function \code{\link{pcs.prepare}()}. %% ~~ describe objects here ~~
#'@author Chris Wallace.
#'@seealso \code{\link{pcs.prepare}}, \code{\link{pcs.model}}
#'@references Wallace et al (2012).  Statistical colocalisation of monocyte
#'gene expression and genetic risk variants for type 1 diabetes.  Hum Mol Genet
#'21:2815-2824.  \url{http://europepmc.org/abstract/MED/22403184}
#'
#'Plagnol et al (2009).  Statistical independence of the colocalized
#'association signals for type 1 diabetes and RPS26 gene expression on
#'chromosome 12q13. Biostatistics 10:327-34.
#'\url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
#'@keywords classes
#'@export
#'@examples
#'
#'showClass("colocPCs")
#'
setClass("colocPCs",
         representation(pcs="matrix",
                        group="numeric",
                        use="logical",
                        vars="numeric"),
         validity=validColocPCs)

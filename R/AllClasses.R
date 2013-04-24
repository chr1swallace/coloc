### coloc class: simple class to hold results of coloc.test()
validColoc <- function(object) {
  if(length(object@result) != 3 ||
     !all(c("eta.hat","chisquare","n") %in% names(object@result))) {
    return("result slot in coloc objects should be a named numeric vector with 3 elements: eta.hat, chisquare, n")
  }
  ## if(length(object@u) != nrow(object@V))
  ##   return("dim of V should equal length of u")
}
#'Classes \code{"coloc"} and \code{"colocBayes"}
#'
#'Classes
#'designed to hold objects returned by function \code{\link{coloc.test}} which
#'performs a test of the null hypothesis that two genetic traits colocalise -
#'that they share a common causal variant.
#'
#'
#'@name coloc-class
#'@aliases coloc-class colocBayes-class
#'@docType class
#'@section Objects from the Class: Objects can be created by calls to the
#'function \code{\link{coloc.test}()}.  Class \code{colocBayes} extends class
#'\code{coloc}. 
#'@author Chris Wallace.
#'@seealso \code{\link{coloc.test}}
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
#'
setClass("coloc",
         representation(result="numeric"))
setClass("colocBayes",
         representation(ppp="numeric",
                        credible.interval="list",
                        bayes.factor="numeric"),
         contains="coloc")
setClass("colocBMA",
         representation(ppp="numeric",
                        bma="numeric",
                        bayes.factor="numeric"),
         contains="coloc")

show.coloc <- function(object) {
  if(!is.na(object@result["chisquare"])) {
    res <- c(object@result,  p.value=p.value(object))
  } else {
    res <- c(object@result[c("eta.hat","n")],ppp.value=ppp.value(object))
  }
  print(res)
}
setMethod("show","coloc",show.coloc)
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

#'Functions to extract information from a \code{coloc} or \code{colocBayes}
#'object
#'
#'Extract information from a \code{coloc} object.
#'
#'\code{eta()} returns eta.hat, the maximum likelihood value of eta.
#'
#'\code{theta()} returns theta.hat, the maximum likelihood value of eta.
#'
#'\code{lhood()} returns -2 times the log-likelihood ratio..
#'
#'\code{chisquare()} returns the value of the chisquare statistic derived by
#'Wallace et al (in preparation).
#'
#'\code{df()} returns the associated degrees of freedom.
#'
#'\code{p.value()} returns the associated p value.
#'
#'\code{ppp.value()} returns the posterior predicted p value, or \code{NA} if
#'not calculated.
#'
#'\code{ci()} returns the credible interval, or \code{NA} if not calculated.
#'
#'@aliases eta theta df ppp.value p.value chisquare ci bf eta,coloc-method
#'p.value,coloc-method show,coloc-method chisquare,coloc-method
#'theta,coloc-method df,coloc-method ci,colocBayes-method bf,colocBayes-method
#'ppp.value,colocBayes-method plot,colocPCs,missing-method
#'@param object Object returned by \code{coloc.test()} function.
#'@author Chris Wallace.
#'@seealso \code{\link{coloc.test}}, \code{\link{pcs.prepare}}
#'@keywords methods
setGeneric("eta",function(object) standardGeneric("eta"))
setMethod("eta","coloc",function(object) object@result["eta.hat"])
setGeneric("theta",function(object) standardGeneric("theta"))
setMethod("theta","coloc",function(object) {
          theta <- atan(object@result["eta.hat"])
          names(theta) <- "theta.hat"
          theta
          })
setGeneric("chisquare",function(object) standardGeneric("chisquare"))
setMethod("chisquare","coloc",function(object) object@result["chisquare"])
setGeneric("ci",function(object) standardGeneric("ci"))
setMethod("ci","colocBayes",function(object) {
  if(is.na(object@credible.interval$level.observed)) {
    return(object@credible.interval[c("eta.mean","eta.mode","lower","upper","level","interior")])
  } else {
    return(object@credible.interval[c("eta.mean","eta.mode","lower","upper","level.observed","interior")])
  }})
setGeneric("df",function(object) standardGeneric("df"))
setMethod("df","coloc",function(object) object@result["n"]-1)
## setGeneric("n",function(object) standardGeneric("n"))
## setMethod("n","coloc",function(object) object@result["n"])
setGeneric("ppp.value",function(object) standardGeneric("ppp.value"))
setMethod("ppp.value","colocBayes",function(object) object@ppp)
setMethod("ppp.value","colocBMA",function(object) object@ppp)
setGeneric("p.value",function(object) standardGeneric("p.value"))
setMethod("p.value","coloc",function(object)
          pchisq(object@result["chisquare"],df=object@result["n"]-1,lower.tail=FALSE))

setGeneric("bf",function(object) standardGeneric("bf"))
setMethod("bf","colocBayes",function(object) {
  if(!length(object@bayes.factor))
    stop("No Bayes factor calculations stored.\n")
  if(length(object@bayes.factor)==1)
    warning("Comparing a single value or interval for eta to itself.  Probably not what you meant to do.  bayes.factor should have length at least two.")
  bayes.factor.table <- outer(object@bayes.factor, object@bayes.factor, "/")
  dimnames(bayes.factor.table) <- list(values.for=names(object@bayes.factor),
                                       values.against=names(object@bayes.factor))
  return(bayes.factor.table)
})

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
setMethod("show","colocPCs",function(object) {
  tt <- table(object@group)
  cat("Principal component matrix with",ncol(object@pcs),"components\ncovering",nrow(object@pcs),"individuals from",length(tt),"groups,\nwith representation:")
  print(tt)
})
setGeneric("plot")
setMethod("plot", signature(x="colocPCs",y="missing"),
          function(x, threshold=0.8, ggplot2=FALSE) {
            n <- length(object@vars)
            npc <- which(object@vars>threshold)[1]
                                        #  title <- paste("Number of components required to capture >=",threshold,"of the variance")
            if(ggplot2) {
              require(ggplot2)
              ggplot(data.frame(n=1:n,v=object@vars),
                     aes_string(x = 'n', y = 'v') # to get around R CMD check complaining about next line
                                        # aes(x=n,y=v)
                     ) + xlab("Number of components") +
                       ylab("Proportion of variance explained") + geom_vline(xintercept=npc,lty=2,col="grey") +
                         geom_line()
            } else {
              plot(1:n, object@vars, xlab="Number of components", ylab="Proportion of variance explained",
                   type="l")
              abline(v=npc,col="grey",lty=2)
            }
          })


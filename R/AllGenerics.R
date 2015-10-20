#' Plotting functions for the coloc package
#'
#' You can plot objects of class \code{coloc}, \code{colocBayes} and \code{colocABF}
#'
#' @param x object to be plotted
#' @param y ignored
#' @param ... other arguments
#' @return no return value
#' @export
#' @docType methods
#' @rdname plot-methods
setGeneric("plot")

#'Methods to extract information from a \code{coloc} or \code{colocBayes}
#'object
#'
#'Extract information from a \code{coloc} object.
#'
#'\code{eta()} returns eta.hat, the maximum likelihood value of eta.
#'
#'\code{theta()} returns theta.hat, the maximum likelihood value of eta.
#'
#'\code{summary()} returns a summary, giving eta, chisquare statistic, number of SNPs/PCs, p value and, if a \code{colocBayes} object, the ppp.value
#'
#'\code{ci()} returns the credible interval, or \code{NA} if not calculated.
#'
#'@aliases eta theta summary ci
#'eta,coloc-method theta,coloc-method summary,coloc-method ci,coloc-method
#'eta,colocBayes-method theta,colocBayes-method summary,colocBayes-method ci,colocBayes-method
#'@param object Object returned by \code{coloc.test()} or \code{coloc.bma()} functions.
#'@author Chris Wallace.
#'@seealso \code{\link{coloc.test}}, \code{\link{pcs.prepare}}
#'@exportMethod eta
#'@exportMethod theta
#'@exportMethod ci
#'@exportMethod summary
#'@keywords methods
setGeneric("eta",function(object) standardGeneric("eta"))
setGeneric("theta",function(object) standardGeneric("theta"))
setGeneric("ci",function(object) standardGeneric("ci"))
setGeneric("chisquare",function(object) standardGeneric("chisquare"))
setGeneric("df",function(object) standardGeneric("df"))
setGeneric("ppp.value",function(object) standardGeneric("ppp.value"))
setGeneric("p.value",function(object) standardGeneric("p.value"))

##' Summarise the evidence for/against specific values or ranges of eta using bayes factors
##'
##' Only available for \code{colocBayes} objects, and you need to specify the specific values of interest using the \code{bayes.factor} argument when doing the proportional coloc analysis
##' @title Bayes factors to compare specific values of eta
##' @param object of class \code{colocBayes}
##' @return a matrix of Bayes factors
##' @rdname bf-methods
##' @aliases bf
#' @docType methods
##' @exportMethod  bf
##' @author Chris Wallace
##' @keywords methods
setGeneric("bf",function(object) standardGeneric("bf"))

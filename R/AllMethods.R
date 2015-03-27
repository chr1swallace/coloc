#' @rdname plot-methods
#' @aliases plot,coloc,missing-method
setMethod("plot", signature(x="coloc",y="missing"),
          function(x) {
            coeff.plot(b1=x@plot.data$coef1,
                     b2=x@plot.data$coef2,
                     s1=x@plot.data$var1,
                     s2=x@plot.data$var2,
                     eta=x@result[["eta.hat"]],
                     main="Coefficients",
                     xlab=expression(b[1]),ylab=expression(b[2]))})
#' @rdname plot-methods
#' @aliases plot,colocABF,missing-method
setMethod("plot", signature(x="colocABF",y="missing"),
          function(x,...) {
            abf.plot(x, ...)
          })
  
#' @rdname plot-methods
#' @aliases plot,coloc,missing-method
setMethod("plot", signature(x="coloc",y="missing"),
          function(x,...) {
            if(!("model.prob" %in% names(x@plot.data)))
              x@plot.data$model.prob <- NULL
            coeff.plot(b1=x@plot.data$coef1,
                       b2=x@plot.data$coef2,
                       s1=sqrt(x@plot.data$var1),
                       s2=sqrt(x@plot.data$var2),
                       alpha=x@plot.data$model.prob,
                       slope=x@result[["eta.hat"]],
                       annot=paste("p =",format.pval(p.value(x))),
                       ...)
          })
setMethod("plot", signature(x="colocPCs",y="missing"),
          function(x, threshold=0.8, ggplot2=FALSE) {
            n <- length(object@vars)
            npc <- which(object@vars>threshold)[1]
                                                    #  title <- paste("Number of components required to capture >=",threshold,"of the variance")
##             if(ggplot2) {
##               require(ggplot2)
##               ggplot(data.frame(n=1:n,v=object@vars),
##                      aes_string(x = 'n', y = 'v') # to get around R CMD check complaining about next line
##                                         # aes(x=n,y=v)
##                      ) + xlab("Number of components") +
##                        ylab("Proportion of variance explained") + geom_vline(xintercept=npc,lty=2,col="grey") +
##                          geom_line()
##             } else {
              plot(1:n, object@vars, xlab="Number of components", ylab="Proportion of variance explained",
                   type="l")
              abline(v=npc,col="grey",lty=2)
##             }
          })


setMethod("show","colocPCs",function(object) {
  tt <- table(object@group)
  cat("Principal component matrix with",ncol(object@pcs),"components\ncovering",nrow(object@pcs),"individuals from",length(tt),"groups,\nwith representation:")
  print(tt)
})

setMethod("eta","coloc",function(object) object@result["eta.hat"])
setMethod("theta","coloc",function(object) {
          theta <- atan(object@result["eta.hat"])
          names(theta) <- "theta.hat"
          theta
          })
setMethod("ci","colocBayes",function(object) {
  if(is.null(object@credible.interval$level.observed)) {
    return(object@credible.interval[c("eta.mode","lower","upper","level","interior")])
  } else {
    return(object@credible.interval[c("eta.mode","lower","upper","level.observed","interior")])
  }})

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


## these are not exported
setMethod("chisquare","coloc",function(object) object@result["chisquare"])
setMethod("df","coloc",function(object) object@result["n"]-1)
## setGeneric("n",function(object) standardGeneric("n"))
## setMethod("n","coloc",function(object) object@result["n"])
setMethod("ppp.value","colocBayes",function(object) object@ppp)
setMethod("ppp.value","colocBMA",function(object) object@ppp)
setMethod("p.value","coloc",function(object) pchisq(object@result["chisquare"],df=object@result["n"]-1,lower.tail=FALSE))

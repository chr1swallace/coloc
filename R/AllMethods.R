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
            coeff.plot(b1=x@plot.data$coef1,
                       b2=x@plot.data$coef2,
                       s1=sqrt(x@plot.data$var1),
                       s2=sqrt(x@plot.data$var2),
                       alpha=x@plot.data$model.prob,
                       slope=x@result[["eta.hat"]],
                       annot=paste("p =",format.pval(x@result[["p"]])),
                       ...)
          })
  

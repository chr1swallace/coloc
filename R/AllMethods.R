#' @rdname plot-methods
#' @aliases plot,coloc,missing-method
setMethod("plot", signature(x="coloc",y="missing"),
          coeff.plot(b1=x@plot.data$coef1,
                     b2=x@plot.data$coef2,
                     s1=x@plot.data$var1,
                     s2=x@plot.data$var2,
                     eta=x@result[["eta.hat"]],
                     main="Coefficients",
                     xlab=expression(b[1]),ylab=expression(b[2])))
  

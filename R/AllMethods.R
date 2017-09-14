#' @rdname plot-methods
#' @aliases plot,coloc,missing-method
setMethod("plot", signature(x="colocTWAS",y="missing"),
          function(x) {
              use <- which(!sapply(x@plot.data,is.null))
              if(length(use)==0)
                  return(NULL)
              if(length(use)>1) 
                  opar <- par(ask=TRUE)
              for(i in use) {
                  last.plot <- coeff.plot(b1=x@plot.data[[i]]$coef1,
                             b2=x@plot.data[[i]]$coef2,
                             s1=x@plot.data[[i]]$var1,
                             s2=x@plot.data[[i]]$var2,
                             eta=x@result[i,"eta.hat"],
                             annot=paste("p =",format.pval(pchisq(x@result[i,"chisquare"],x@result[i,"n"]-1,lower.tail=FALSE))),
                             xlab=expression(b[1]),ylab=expression(b[2]))
                  print(last.plot)
              }
              if(length(use)>1)
                  par(opar)
              })
#' @rdname plot-methods
#' @aliases plot,coloc,missing-method
setMethod("plot", signature(x="coloc",y="missing"),
          function(x) {
            coeff.plot(b1=x@plot.data$coef1,
                     b2=x@plot.data$coef2,
                     s1=x@plot.data$var1,
                     s2=x@plot.data$var2,
                     eta=x@result[["eta.hat"]],
                     annot=paste("p =",format.pval(pchisq(x@result[["chisquare"]],x@result[["n"]]-1,lower.tail=FALSE))),
                     xlab=expression(b[1]),ylab=expression(b[2]))})
#' @rdname plot-methods
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
#' @rdname plot-methods
#' @aliases plot,colocPCs,missing-method
setMethod("plot", signature(x="colocPCs",y="missing"),
          function(x) {
              threshold=0.8
              n <- length(x@vars)
              npc <- which(x@vars>threshold)[1]
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
              plot(1:n, x@vars, xlab="Number of components", ylab="Proportion of variance explained",
                   type="l")
              abline(v=npc,col="grey",lty=2)
##             }
          })

setMethod("show","colocTWAS", function(object) {
    print(object@result)
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

##' @rdname bf-methods
##' @aliases bf,colocBayes-method
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
setMethod("p.value","coloc",function(object) {
    if("p" %in% names(object@result)) {
        return(object@result["p"])
    }
    return(pchisq(object@result["chisquare"],df=object@result["n"]-1,lower.tail=FALSE))
})


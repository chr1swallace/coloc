## bootp <- function(b1,b2,s1,s2,eta) {
## }

extra.plot <- function(plot.data, plots.extra, theta.hat, eta.hat) {
  if(!is.list(plots.extra) || length(plots.extra)!=2 ||
     !("x" %in% names(plots.extra)) || !("y" %in% names(plots.extra)) ||
     !all(plots.extra$x %in% names(plot.data)) ||
     !all(plots.extra$y %in% names(plot.data))) {
    warning("plots.extra must be of the format list(x=..., y=...) where x and y are equal length character vectors with elements from theta, eta, lhood, chisq, post.theta.  Skipping plots.extra.")
    return(NULL)
  }
  
  labels <- list(theta="theta",
                 eta="eta",
                 chisq="chisquare",
                 post.theta="posterior for theta",
                 lhood="-2 log likelihood")
  if(length(plots.extra$y)>length(plots.extra$x)) {
    plots.extra$x <- rep(plots.extra$x,length=length(plots.extra$y))
  }
  if(length(plots.extra$x)>length(plots.extra$y)) {
    plots.extra$y <- rep(plots.extra$y,length=length(plots.extra$x))
  }
  for(i in 1:length(plots.extra$x)) {
    plot(plot.data[[ plots.extra$x[i] ]],
         plot.data[[ plots.extra$y[i] ]],
         type="l",axes=FALSE,
         xlab=labels[[ plots.extra$x[i] ]],ylab=labels[[ plots.extra$y[i] ]],
         main=paste(labels[[ plots.extra$y[i] ]],"vs",labels[[ plots.extra$x[i] ]]))
    if(plots.extra$x[i]=="theta") {
      axis(1,at=seq(0,pi,length=5),
           labels=c("0",expression(pi / 4),expression(pi / 2), expression(3 * pi / 4), expression(pi)))
    } else {
      axis(1)
    }
    axis(2); box()
    if(plots.extra$x[i]=="theta")
      abline(v=theta.hat,col="blue",lty=3)
    if(plots.extra$x[i]=="eta")
      abline(v=eta.hat,col="blue",lty=3)

  }
}

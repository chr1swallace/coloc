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


coeff.plot <- function(b1,b2,s1,s2,eta,add=NULL, ...) {
  c1 <- cbind(b1,sqrt(s1),b1+1.96*sqrt(s1),b1-1.96*sqrt(s1))
  c2 <- cbind(b2,sqrt(s2),b2+1.96*sqrt(s2),b2-1.96*sqrt(s2))
  if(!is.null(add)) {
    c1 <- cbind(c1,add[[1]])
    c2 <- cbind(c2,add[[2]])
  }
  xr <- range(c1)
  yr <- range(c2)
  if(xr[1]>0)
    xr[1] <- -0.005
  if(xr[2]<0)
    xr[2] <- 0.005
  if(yr[1]>0)
    yr[1] <- -0.005
  if(yr[2]<0)
    yr[2] <- 0.005
  plot(c1[,1],c2[,1],pch=".",xlim=xr,ylim=yr, ...)
  T <- seq(0,2,by=0.05)*pi
  for(j in 1:nrow(c1)) {
    X <- c1[j,1] + 1.96*c1[j,2]*cos(T)
    Y <- c2[j,1] + 1.96*c2[j,2]*sin(T)
    polygon(X,Y,lty=0,col="grey90")
  }
  for(j in 1:nrow(c1)) {
    X <- c1[j,1] + 1.96*c1[j,2]*cos(T)
    Y <- c2[j,1] + 1.96*c2[j,2]*sin(T)
    polygon(X,Y,lty=3)
  }
  abline(0,eta,col="green",lwd=1.5)
  abline(v=0,lty=2,col="green4",lwd=1.5)
  abline(h=0,lty=2,col="green4",lwd=1.5)
  arrows(c1[,3,drop=FALSE],c2[,1,drop=FALSE],c1[,4,drop=FALSE],c2[,1,drop=FALSE],angle=90,length=0.02,code=3)
  arrows(c1[,1,drop=FALSE],c2[,3,drop=FALSE],c1[,1,drop=FALSE],c2[,4,drop=FALSE],angle=90,length=0.02,col="red",code=3)
  points(c1[,1,drop=FALSE],c2[,1,drop=FALSE],pch="*",cex=1)
  if(!is.null(add)) {
    points(add[[1]],add[[2]],pch="*",cex=2,col="blue")
    segments(c1[,5],c2[,5],c1[,1],c2[,1],col="blue",lty=3)
  }
}



                                        # ds1: dataset1
# ds2: dataset2
# both ds1 and ds2 should contain the same snps in the same order

# Pos: positions of all snps in ds1 or in ds2
# chr: Chromosome
# pos.start: lower bound of positions
# pos.end: upper bound of positions
# trait1: name of trait 1
# trait2: name of trait 2   
ymin <- NULL
ymax <- NULL


##' Print summary of a coloc.abf run
##'
##' @title print.coloc_abf
##' @param x object of class \code{coloc_abf} returned by coloc.abf()
##'   or coloc.signals()
##' @param ... optional arguments: "trait1" name of trait 1, "trait2"
##'   name of trait 2
##' @return x, invisibly
##' @author Chris Wallace
##' @export
##' @docType methods
##' @rdname print-methods
print.coloc_abf <- function(x,...) {
  trait1="trait 1"
  trait2="trait 2"
  args <- list(...)
  if("trait1" %in% names(args))
    trait1=args$trait1
  if("trait2" %in% names(args))
    trait2=args$trait2
  
    message("Coloc analysis of ",trait1,", ",trait2)
    message("\nSNP Priors")
    print(x$priors)
    message("\nHypothesis Priors")
    ns <- if(is.data.table(x$summary)) { x$summary$nsnps[1] } else { x$summary["nsnps"] }
    hprior <- prior.snp2hyp(nsnp=ns,
                        p1=x$priors["p1"],
                        p2=x$priors["p2"],
                        p12=x$priors["p12"])
    rownames(hprior) <- ""
    print(hprior)
    message("\nPosterior")
    if(is.data.table(x$summary)) {
        summ <- copy(x$summary)
        setnames(summ, gsub("PP.|.abf","",names(summ)))
        print(summ[,list(nsnps,hit1,hit2,H0,H1,H2,H3,H4)])
    } else {
        summ <- x$summary
        names(summ) <- gsub("PP.|.abf","",names(summ))
        print(summ)
    }
    invisible(x)
}
    
##' Plot results of a coloc.abf run
##'
##' If obj is missing, it will be created as obj=coloc.abf(ds1,ds2).  Both ds1 and ds2 should contain the same snps in the same order
##' @title plot.coloc_abf
##' @param obj object of class \code{colocABF} returned by coloc.abf()
##' @param Pos positions of all snps in ds1 or in ds2
##' @param chr Chromosome
##' @param pos.start lower bound of positions
##' @param pos.end upper bound of positions
##' @param trait1 name of trait 1
##' @param trait2 name of trait 2   
##' @return a ggplot object
##' @author Hui Guo, Chris Wallace
##' @docType methods
##' @rdname plot-methods
plot.coloc_abf <- function(obj, Pos=1:nrow(obj$results),
                     chr=NULL, pos.start=min(Pos), pos.end=max(Pos),
                     trait1="trait 1", trait2="trait 2") {

  
  d.pp1 = signif(obj$summary["PP.H1.abf"], 3)
  d.pp2 = signif(obj$summary["PP.H2.abf"], 3)
  d.pp4 = signif(obj$summary["PP.H4.abf"], 3)
  df = obj$results
  
    df$pp1 <- exp(df$lABF.df1 - logsum(df$lABF.df1))
    df$pp2 <- exp(df$lABF.df2 - logsum(df$lABF.df2))
    if(!("position" %in% colnames(df)))
       df$position <- Pos

  df <- melt(df[,c("snp","position","pp1","pp2","SNP.PP.H4")], id.vars=c("snp","position"))
  df$variable <- sub("pp1", paste0("H1 (", trait1, "), PP = ", d.pp1) ,df$variable)
  df$variable <- sub("pp2", paste0("H2 (", trait2, ") PP = ", d.pp2) ,df$variable)
  df$variable <- sub("SNP.PP.H4", paste0("H4 (Both) = PP ", d.pp4) ,df$variable)


  ## identify and label the top 3 SNPs that have highest pp1, pp2 or pp4
  
  df.ord = df[order(df$value, decreasing=TRUE), ]
  snps = unique(df.ord$snp)[1:3]

    label <- NULL # avoid R CMD check NOTE
  df$label <- ifelse(df$snp %in% snps, df$snp,"")
  ttl <- paste0(trait1, ' & ', trait2, ' (chr', chr, ': ', pos.start, '-', pos.end, ')')
 
  ggplot(df, aes_string(x="position",y="value")) +
    geom_point(data=subset(df,label==""),size=1.5) +
    geom_point(data=subset(df,label!=""),col="red",size=1.5) +
    geom_text(aes_string(label="label"),hjust=-0.1,vjust=0.5,size=2.5,col="red") +
    facet_grid(variable ~ .) +
    theme(legend.position="none") + xlab(paste("Chromosome", chr, sep=' ')) + ylab("Posterior probability") + 
    ggtitle(ttl)
    
}



coeff.plot <- function(b1,b2,s1,s2,eta,lower=NULL,upper=NULL,add=NULL,alpha=NULL,slope=NULL,annot=NULL, ...) {
##   c1 <- cbind(b1,sqrt(s1),b1+1.96*sqrt(s1),b1-1.96*sqrt(s1))
##   c2 <- cbind(b2,sqrt(s2),b2+1.96*sqrt(s2),b2-1.96*sqrt(s2))
##   if(!is.null(add)) {
##     c1 <- cbind(c1,add[[1]])
##     c2 <- cbind(c2,add[[2]])
##   }
##   xr <- range(c1)
##   yr <- range(c2)
##   if(xr[1]>0)
##     xr[1] <- -0.005
##   if(xr[2]<0)
##     xr[2] <- 0.005
##   if(yr[1]>0)
##     yr[1] <- -0.005
##   if(yr[2]<0)
##     yr[2] <- 0.005

  ## ggplot version
  df <- data.frame(x=b1, y=b2, x.se=s1, y.se=s2,id=1:length(b1))
  if(!is.null(alpha)) {
    nr <- length(b1)/length(alpha)
    df$alpha <- rep(alpha,each=nr)
  } else {
    df$alpha <- 1
  }
  ##  T <- seq(0,2,by=0.05)*pi
##   df.path <- do.call("rbind",lapply(1:nrow(c1), function(j) {
##    tmp <- data.frame(x=c1[j,1] + 1.96*c1[j,2]*cos(T),
##                      y=c2[j,1] + 1.96*c2[j,2]*sin(T))
##    tmp$id <- j
##    return(tmp)
##   }))
  x <- y <- x.se <- y.se <- NULL
    p <- ggplot(df) +
        geom_hline(yintercept=0,linetype="dotted") +
        geom_vline(xintercept=0,linetype="dotted")
     if(!is.null(lower) && !is.null(upper)) {
        rib <- with(df,data.frame(x=seq(1.1*min(x-1.96*x.se),1.1*max(x+1.96*x.se),length.out=100)))
        rib$ymin <- rib$x * lower
        rib$ymax <- rib$x * upper
        rib$xmin <- rib$xmax <- rib$x; rib$alpha <- 1
        p <- p +
            geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax),data=rib,fill="steelblue4",alpha=0.1) +
            scale_x_continuous(limits=c(min(rib$xmin),max(rib$xmax)))
    }
    p <- p + geom_point(aes(x=x,y=y,alpha=alpha),
                   col="steelblue3") +
        geom_errorbar(aes(x=x,
                       ymin=y-1.96*y.se,ymax=y+1.96*y.se,alpha=alpha),
                   col="steelblue4") + 
        geom_errorbarh(aes(y=y,
                       xmin=x-1.96*x.se,xmax=x+1.96*x.se, alpha=alpha),
                   col="steelblue4") +
        theme(legend.position="none") +
        geom_abline(slope=eta,colour="steelblue3",size=1) +
        labs(x="beta.1",y="beta.2")
    if(!is.null(annot))
        p <- p + annotate("text", y= max(df$y+1.96*df$y.se), x =max(df$x+1.96*df$x.se),label=annot,hjust=1)
  return(p)
  
##   plot(c1[,1],c2[,1],pch=".",xlim=xr,ylim=yr, ...)

##   for(j in 1:nrow(c1)) {
##     X <- c1[j,1] + 1.96*c1[j,2]*cos(T)
##     Y <- c2[j,1] + 1.96*c2[j,2]*sin(T)
##     polygon(X,Y,lty=0,col="grey90")
##   }
##   for(j in 1:nrow(c1)) {
##     X <- c1[j,1] + 1.96*c1[j,2]*cos(T)
##     Y <- c2[j,1] + 1.96*c2[j,2]*sin(T)
##     polygon(X,Y,lty=3)
##   }
##   abline(0,eta,col="green",lwd=1.5)
##   abline(v=0,lty=2,col="green4",lwd=1.5)
##   abline(h=0,lty=2,col="green4",lwd=1.5)
##   arrows(c1[,3,drop=FALSE],c2[,1,drop=FALSE],c1[,4,drop=FALSE],c2[,1,drop=FALSE],angle=90,length=0.02,code=3)
##   arrows(c1[,1,drop=FALSE],c2[,3,drop=FALSE],c1[,1,drop=FALSE],c2[,4,drop=FALSE],angle=90,length=0.02,col="red",code=3)
##   points(c1[,1,drop=FALSE],c2[,1,drop=FALSE],pch="*",cex=1)
##   if(!is.null(add)) {
##     points(add[[1]],add[[2]],pch="*",cex=2,col="blue")
##     segments(c1[,5],c2[,5],c1[,1],c2[,1],col="blue",lty=3)
##   }
}

#' @importFrom graphics points
globalVariables(c("variable","pp","position","value"))

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

##' Plot a coloc structured dataset
##'
##' @title plot a coloc dataset
##' @param d a coloc dataset
##' @param susie_obj optional, the output of a call to runsusie()
##' @param credset  optional, a list of vectors, with each vector = one credible set
##' @param alty default is to plot a standard manhattan. If you wish to plot a
##'   different y value, pass it here. You may also want to change ylab to
##'   describe what you are plotting.
##' @param ylab label for y axis, default is -log10(p) and assumes you are
##'   plotting a manhattan
##' @param show_legend optional, show the legend or not. default is TRUE
##' @param color optional, specify the colours to use for each credible set when
##'   susie_obj is supplied. Default is shamelessly copied from
##'   susieR::susie_plot() so that colours will match
##' @param ... other arguments passed to the base graphics plot() function
##' @author Chris Wallace
##' @export
plot_dataset <- function(d,
                         susie_obj=NULL,
                         credset=NULL,
                         alty=NULL,ylab="-log10(p)",
                         show_legend=TRUE,
                         color = c("dodgerblue2", "green4", "#6A3D9A", "#FF7F00",
                                   "gold1", "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6",
                                   "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1", "deeppink1",
                                   "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
                                   "yellow3", "darkorange4", "brown"),
                         ...) {
  if(!("position" %in% names(d)))
    stop("no position element given")
  if(is.null(alty)) {
    y= - ( pnorm(-abs(d$beta)/sqrt(d$varbeta),log.p=TRUE) + log(2) ) / log(10)
  } else {
    y=alty
  }
  plot(d$position,y,xlab="Position",ylab=ylab,pch=16,col="grey",...)
  if(!is.null(susie_obj)) {
    cs=susie_obj$sets$cs
    for(i in 1:length(cs)) {
      w=which(d$snp %in% names(cs[[i]]))
      ## print(w)
      ## print(d$position[w])
      ## print(p[w])
      points(d$position[w], y[w], col=color[i], cex=2)
    }
    if(show_legend)
      legend("topright",col=color[1:length(cs)],pch=rep(1,length(cs)),legend=1:length(cs))
  }
  if(!is.null(credset)) {
    for(i in 1:length(credset)) {
      w=which(d$snp %in% credset[[i]])
      points(d$position[w], y[w], col=color[i], cex=2)
    }
    if(show_legend)
      legend("topright",col=color[1:length(credset)],pch=rep(1,length(credset)),legend=1:length(credset))
  }
}

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

##' plot a coloc_abf object
##'
##' @title plot a coloc_abf object
##' @param x coloc_abf object to be plotted
##' @param ... other arguments
##' @return ggplot object
##' @docType methods
##' @export
##' @rdname plot-methods
##' @author Chris Wallace
plot.coloc_abf <- function(x,...) {
  x=x$results
  if(!("position" %in% names(x)))
    x$position <- 1:nrow(x)
  m=melt(x,id.vars=c("snp","position"), measure.vars=grep("z.df",names(x),value=TRUE))
  m2=melt(x,id.vars=c("snp","position"), measure.vars=grep("PP",names(x),value=TRUE))
  setnames(m2,"value","pp")
  if(length(grep("row",m$variable))) {
    m[,c("z","df","row"):=tstrsplit(variable,"\\.")]
    m2[,row:=sub(".*\\.","",variable)]
    m <- merge(m,m2[,.(snp,row,pp)],by=c("snp","row"))
    ggplot(m, aes(x=position,y=abs(value),col=pp,size=pp)) +
      geom_point() +
      facet_grid(df ~ row)
  } else {
    m[,c("z","df"):=tstrsplit(variable,"\\.")]
    m <- merge(m,m2[,.(snp,pp)],by=c("snp"))
    ggplot(m, aes(x=position,y=abs(value),col=pp,size=pp)) +
      geom_point() +
      facet_grid(df ~ .)
  }
}


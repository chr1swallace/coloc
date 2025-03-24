##' Check coloc dataset inputs for errors
##'
##' A coloc dataset is a list, containing a mixture of vectors
##' capturing quantities that vary between snps (these vectors must
##' all have equal length) and scalars capturing quantities that
##' describe the dataset.
##'
##' Coloc is flexible, requiring perhaps only p values, or z scores, or effect
##' estimates and standard errors, but with this flexibility, also comes
##' difficulties describing exactly the combinations of items required.
##'
##' Required vectors are some subset of
##'
##' \describe{
##'   \item{beta}{regression coefficient for each SNP from dataset 1}
##'   \item{varbeta}{variance of beta}
##'   \item{pvalues}{P-values for each SNP in dataset 1}
##'   \item{MAF}{minor allele frequency of the variants}
##'   \item{snp}{a character vector of snp ids, optional. It will be used to merge dataset1 and dataset2 and will be retained in the results.}
##' }
##'
##' Preferably, give \code{beta} and \code{varbeta}.  But if these are not available, sufficient statistics can be approximated from \code{pvalues} and \code{MAF}.
##'
##' Required scalars are some subset of
##'
##' \describe{
##'   \item{N}{Number of samples in dataset 1}
##'   \item{type}{the type of data in dataset 1 - either "quant" or "cc" to denote quantitative or case-control}
##'   \item{s}{for a case control dataset, the proportion of samples in dataset 1 that are cases}
##'   \item{sdY}{for a quantitative trait, the population standard deviation of the trait.  if not given, it can be estimated from the vectors of varbeta and MAF}
##' }
##'
##' You must always give {\code{type}}.  Then,
##' \describe{
##' \item{if \code{type}=="cc"}{\code{s}}
##' \item{if \code{type}=="quant" and \code{sdY} known}{\code{sdY}}
##' \item{if beta, varbeta not known}{\code{N}}
##' }
##' If \code{sdY} is unknown, it will be approximated, and this will require
##' \describe{
##' \item{summary data to estimate \code{sdY}}{\code{beta}, \code{varbeta}, \code{N}, \code{MAF}}
##' }
##'
##' Optional vectors are
##'
##' \describe{
##'   \item{position}{a vector of snp positions, required for \code{plot_dataset}}
##' }
##'
##' \code{check_dataset} calls stop() unless a series of expectations on dataset
##' input format are met
##'
##' This is a helper function for use by other coloc functions, but
##' you can use it directly to check the format of a dataset to be
##' supplied to coloc.abf(), coloc.signals(), finemap.abf(), or
##' finemap.signals().
##' @title check_dataset
##' @param d dataset to check
##' @param suffix string to identify which dataset (1 or 2)
##' @param req names of elements that must be present
##' @param warn.minp print warning if no p value < warn.minp
##' @return NULL if no errors found
##' @export
##' @author Chris Wallace
check_dataset <- function(d,
                          suffix="",
                          req=c("type","snp"),
                          warn.minp=1e-6) {
    if(!is.list(d) )
        stop("dataset ",suffix,": is not a list")
    recognised_items=c("beta","varbeta","pvalues","MAF","snp","position","N","type","s","sdY","LD")
    nd <- intersect(names(d), recognised_items)

    ## no missing values - make people clean their own data rather
    ## than make assumptions here for datasets I don't know
    n <- 0
    
    if(length(setdiff(req,nd)))
        stop("dataset ",suffix,": missing required element(s) ",paste(setdiff(req,nd),collapse=", "))

    for(v in nd) {
        if(any(is.na(d[[v]])))
            stop("dataset ",suffix,": ",v," contains missing values")
    }

    if (!(d$type %in% c("quant", "cc")))
        stop("dataset ", suffix, ": ", "type must be quant or cc")

    ## snps should be unique
    if("snp" %in% nd && any(duplicated(d$snp)))
        stop("dataset ",suffix,": duplicated snps found")
    if("snp" %in% nd && is.factor(d$snp))
        stop("dataset ",suffix,": snp should be a character vector but is a factor")

    ## MAF should be > 0, < 1
    if("MAF" %in% nd && (!is.numeric(d$MAF) || any(is.na(d$MAF)) ||
                         any(d$MAF<=0) || any(d$MAF>=1)))
        stop("dataset ",suffix,": MAF should be a numeric, strictly >0 & <1")

    ## lengths of vector arguments should match
    l <- -1 # impossible length
    shouldmatch <- intersect(nd, c("pvalues","MAF","beta","varbeta","snp","position"))
    for(v in shouldmatch)
        if(l<0) { ## update
            l <- length(d[[v]])
        } else { ## check
            if(length(d[[v]])!=l) {
                stop("dataset ",suffix,": lengths of inputs don't match: ")
                print(shouldmatch)
            }
        }

    ## type of data
    if (! ('type' %in% nd))
        stop("dataset ",suffix,": variable type not set")
    if(!(d$type %in% c("quant","cc")))
        stop("dataset ",suffix,": ","type must be quant or cc")

    ## no beta/varbeta
    if(("s" %in% nd) && (!is.numeric(d$s) || d$s<=0 || d$s>=1))
        stop("dataset ",suffix,": ","s must be between 0 and 1")
    if(!("beta" %in% nd) || !("varbeta" %in% nd)) { # need to estimate var (Y)
        if(!("pvalues" %in% nd) || !( "MAF" %in% nd))
            stop("dataset ",suffix,": ","require p values and MAF if beta, varbeta are unavailable")
        if(any(d$pvalues<=0))
            stop("pvalues should not be negative or exactly 0")
        if(d$type=="cc" && !("s" %in% nd))
            stop("dataset ",suffix,": ","require, s, proportion of samples who are cases, if beta, varbeta are unavailable")
        if (!('N' %in% nd) || is.null(d$N) || any(d$N<=0) )
            stop("dataset ",suffix,": sample size N <=0 or not set")
        p=d$pvalues
    } else {
        p=pnorm( -abs( d$beta/sqrt(d$varbeta) ) ) * 2
    }

    # Infinite values.
    if (any(is.infinite(d$beta) | is.infinite(d$varbeta))) {
      stop("dataset ", suffix,": Infinite values in beta and/or varbeta")
    }

    ## minp
    if(min(p) > warn.minp)
        warning("minimum p value is: ",format.pval(min(p)),"\nIf this is what you expected, this is not a problem.\nIf this is not as small as you expected, please check you supplied var(beta) and not sd(beta) for the varbeta argument. If that's not the explanation, please check the 02_data vignette.")

    ## sdY
    if(d$type=="quant" && !("sdY" %in% nd))
        if(!("MAF" %in% nd && "N" %in% nd ))
            stop("dataset ",suffix,": ","must give sdY for type quant, or, if sdY unknown, MAF and N so it can be estimated")

    if("LD" %in% nd) {
        if(nrow(d$LD)!=ncol(d$LD))
            stop("LD not square")
        if(!identical(colnames(d$LD),rownames(d$LD)))
            stop("LD rownames != colnames")
        if(length(setdiff(d$snp,colnames(d$LD))))
            stop("colnames in LD do not contain all SNPs")
    }

    ## if we reach here, no badness detected
    NULL
}

#'@rdname check_dataset
#' @param ... arguments passed to check_dataset()
#'@export
check.dataset=function(...) {
  warning("Deprecated, use check_dataset() in future")
  check_dataset(...)
}


check_ld <- function(D,LD) {
    if(is.null(LD))
      stop("LD required")
    if(!is.matrix(LD))
      stop("LD must be of class matrix")
    ## if(any(LD[upper.tri(LD)]==1))
    ##   stop("LD includes SNPs in perfect LD")
    if(nrow(LD)!=ncol(LD))
        stop("LD not square")
    if(!identical(colnames(LD),rownames(LD)))
      stop("LD rownames != colnames")
    if(length(setdiff(D$snp,colnames(LD))))
        stop("colnames in LD do not contain all SNPs")
}

##' check alignment between beta and LD
##'
##' @title check alignment
##' @param D a coloc dataset
##' @param thr plot SNP pairs in absolute LD > thr
##' @param do_plot if TRUE (default) plot the diagnostic
##' @export
##' @return proportion of pairs that are positive
##' @author Chris Wallace
check_alignment <- function(D,thr=0.2,do_plot=TRUE) {
  check_dataset(D)
  bprod=outer(D$beta/sqrt(D$varbeta),D$beta/sqrt(D$varbeta),"*")
  ## plot(bprod,D$LD[D$snp,D$snp],xlab="product of z scores",ylab="LD")
  tmp=(bprod/D$LD)[abs(D$LD) > 0.2]
  if(do_plot) {
  hist(tmp,
       xlab="ratio of product of Z scores to LD",
       main="alignment check plot",
       sub="expect most values to be positive\nsymmetry is a warning sign of potentially poor alignment")
  legend("topright",legend=paste0("% positive = ",100*round(mean( tmp > 0 ),3)))
  abline(v=0,col="red")
  }
  return(mean(tmp > 0))
}
#'@rdname check_alignment
#'@export
#' @param ... arguments passed to check_alignment()
check.alignment=function(...) {
  warning("Deprecated, use check_alignment() in future")
  check_alignment(...)
}

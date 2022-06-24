##' tries to be smart about detecting the interesting subregion to finemap/coloc.
##'
##' @title trim a dataset to central peak(s)
##' @param d a coloc dataset
##' @param maxz keep all snps between the leftmost and rightmost snp with |z| >
##'   maxz
##' @param maxr2 expand window to keep all snps between snps with r2 > maxr2
##'   with the left/rightmost snps defined by the maxz threshold
##' @param do.plot if TRUE, plot dataset + boundaries
##' @return logical vector of length d$position indicating which snps to keep
##' @author Chris Wallace
##' @seealso findpeaks
##' @export
findends=function(d, maxz=4, maxr2=0.1, do.plot=FALSE) {
  if(!("position" %in% names(d)))
    stop("need position")
  ord=order(d$position)
  find_bound=function(o) {
    ## find left boundary
    left=cummax(abs(d$z[o]))
    left5=which(left>maxz)
    if(length(left5))
      left5=left5[1]
    else
      stop("no z > supplied maxz ", maxz)
    all5=which(abs(d$z)[o]>maxz)
    if(maxr2>0) {
      ld5=apply(d$LD[d$snp,d$snp][o,o][1:left5,all5,drop=FALSE], 1, max)
      leftr2=which(ld5^2 > maxr2)
      if(length(leftr2))
        leftr2=min(leftr2)
    } else {
      leftr2=left5
    }
    left_bound=d$position[o][leftr2]
  }
  left_bound=find_bound(ord)
  right_bound=find_bound(rev(ord))
  if(do.plot) {
    plot_dataset(d); abline(v=c(left_bound, right_bound))
  }
  keep=d$position >= left_bound & d$position <= right_bound
}

##' tries to be smart about detecting the interesting subregion to finemap/coloc.
##'
##' Differs from findends by finding multiple separate regions if there are multiple peaks
##'
##' @title trim a dataset to only peak(s)
##' @inheritParams findends
##' @return logical vector of length d$position indicating which snps to keep
##' @author Chris Wallace
##' @seealso findends
findpeaks=function(d, maxz=4, maxr2=0.1, do.plot=FALSE) {
  if(!("position" %in% names(d)))
    stop("need position")
  ord=order(d$position)
  find_peak=function(o) {
    ## find left boundary
    peak=which.max(abs(d$z[o]))
    if(abs(d$z[o][peak])< maxz)
      return(NULL)
    ld5=d$LD[d$snp,d$snp][o,o][,peak]
    limr2=which(ld5^2 > maxr2)
    bounds=d$position[o][c(min(limr2),max(limr2))]
    keep=d$position >= bounds[1] & d$position <= bounds[2]
  }
  keep=find_peak(ord)
  while(any(abs(d$z)[!keep] > maxz))
    keep = keep | find_peak(ord[!keep])
  if(do.plot)
    plot(d$position,-log10(pnorm(-abs(d$z))*2),pch=16,col=ifelse(keep,"red","grey"))
  keep
}

##' Subset a coloc dataset
##'
##' @title subset_dataset
##' @param dataset coloc dataset
##' @param index vector of indices of snps to KEEP
##' @return a copy of dataset, with only the data relating to snps in index remaining
##' @author Chris Wallace
##' @export
subset_dataset=function(dataset, index) {
  if(!length(index))
    return(dataset)
  if(!is.numeric(dataset))
    stop("index must be a numeric vector of indexing which snps to keep")
  d=dataset;
  n=length(d$snp) # how many we started with
  if(n<=1) {
    stop("not trimming length 1 dataset. check snp element exists")
  }
  if(max(index)>n)
    stop("cannot subset to more than the number of snps in the dataset. check index")
  message("trimming dataset from ",n," snps to ",length(index))
  for(v in names(dataset)) {
    if(is.matrix(d[[v]]) && ncol(d[[v]])==n && nrow(d[[v]])) {
      d[[v]]=d[[v]][index,index]
      next
    }
    if(is.vector(d[[v]]) && length(d[[v]])==n) {
      d[[v]]=d[[v]][index]
      next
    }
  }
  return(d)
}

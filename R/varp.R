
#' estimate prior probability of snp causality for gene expression based on
#' distance from TSS
#' @importFrom mgcv PredictMat
#' @importFrom scam Predict.matrix.mpd.smooth
#' @param dist distance in bp from TSS for an possible eqtl
#' @param prior_at_mb prior value at dist=1e+6 (1 mega base). default is 1e-4,
#'   which is the default value for p1 in coloc.abf
#' @author Chrs Wallace
#' @examples
#' dist=seq(0,1e+6,by=1000)
#' p=eqtl_prior(dist)
#' plot(dist,p,type="l")
.calc_eqtl_prior=function(dist, prior_at_mb=1e-4) {
  np <- length(dist)
  ## dist[ dist==0 ]=0.01
  predict_base=function(dist) {
    data <- data.frame(dist_bin=abs(dist))
    X <- PredictMat(eqtl_prior_data$smooth[[1]], data)
    if (!is.matrix(X))
      X<- matrix(X, nrow = nrow(data))
    fit <- X %*% eqtl_prior_data$coefficients.t
    linkinv <- eqtl_prior_data$family$linkinv
    as.vector(linkinv(fit))
  }
  sc=predict_base(1e+6)
  pre=predict_base(dist)
  pre*prior_at_mb/sc
}

.pmean=function(x,y) {
  apply(cbind(x,y),1,mean)
}

eqtl_priors=function(position, tss1=NULL, tss2=NULL, p1_at_mb=1e-4, p2_at_mb=1e-4, p12_at_mb=5e-6) {
  if(is.null(tss1) && is.null(tss2)) {
    warning("no TSS given, please supply tss1 or tss2 for variable priors. using values given")
    return(list(p1=p1_at_mb, p2=p2_at_mb, p12=p12_at_mb))
  }
  if(!is.numeric(position))
    stop("position must be numeric")
  ## one variable or both
  if(is.null(tss1) || is.null(tss2)) { # one variable prior
    if(is.null(tss1)) {
      p1=p1_at_mb
      p2=.calc_eqtl_prior(abs(position-tss2), p2_at_mb)
      p12=p12_at_mb * p2 / p2_at_mb
    } else {
      p2=p2_at_mb
      p1=.calc_eqtl_prior(abs(position-tss1), p2_at_mb)
      p12=p12_at_mb * p1 / p1_at_mb
    }
  } else { # both variable
    p1=.calc_eqtl_prior(abs(position-tss1), p2_at_mb)
    p2=.calc_eqtl_prior(abs(position-tss2), p2_at_mb)
    if(tss1==tss2) { # assume proportional
      p12 = p12_at_mb * .pmean(p1,p2) / .pmean(p1_at_mb,p2_at_mb)
    } else { # differnt
      p12 = p12_at_mb * p1 * p2 / p1_at_mb / p2_at_mb
    }
  }
  list(p1=p1,p2=p2,p12=p12)
}

## define an S4 class for priors
.diffequal=function(x,y)
  !( length(x)==1 || length(y)==1 || length(x)==length(y) )

check_prior <- function(object) {
  errors <- character()
  if(.difflength(p1,p2) || .difflength(p1,p12) || .difflength(p2,p12))
    errors=c(errors, "p1, p2, p12 don't have equal length")
}
setClass("coloc_prior",
         representation(p1="numeric",
                        p2="numeric",
                        p12="numeric",
                        ## snp="character",
                        ## position="numeric",
                        tss1="numeric",
                        tss2="numeric"),
         prototype(p1=1e-4,
                   p2=1e-4,
                   p12=5e-6,
                   ## snp=NA_character_,
                   ## position= NA_real_,
                   tss1=NA_real_,
                   tss2=NA_real_),
         validity = check_prior)

##' calculate p1, p2, p12
##'
##' @title get_priors
##' @param priorobj object of class coloc_prior
##' @param position numeric vector of snp positions
##' @return list of p1, p2, p12
##' @export
##' @author Chris Wallace
setGeneric("get_priors", function(priors, position) standardGeneric("get_priors"))
setMethod("get_priors", signature("coloc_prior", "numeric"),
          function(priors, position) {
            eqtl_priors(position, tss1=priors@tss1, tss2=priors@tss2,
                        p1_at_mb=priors@p1, p2_at_mb=priors@p2, p12_at_mb=priors@p12)
          })
setMethod("get_priors", signature("list", "numeric"),
          function(priors, position) {
            priors
          })


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
eqtl_prior=function(dist, prior_at_mb=1e-4) {
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

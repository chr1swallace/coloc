#' estimate prior probability of snp causality for gene expression based on
#' distance from TSS
#' @param dist distance in bp from TSS for an possible eqtl
#' @author Chrs Wallace
#' @examples
#' dist=seq(0,1e+6,by=1000)
#' p=eqtl_prior(dist)
#' plot(dist,p,type="l")
#' @importFrom @mgcv PredictMat
eqtl_prior=function(dist) {
  np <- length(dist)
  data <- data.frame(dist_bin=dist)
  X <- PredictMat(eqtl_prior_data$smooth[[1]], data)
  if (!is.matrix(X))
    X<- matrix(X, nrow = nrow(data))
  fit <- X %*% eqtl_prior_data$coefficients.t
  linkinv <- eqtl_prior_data$family$linkinv
  as.vector(linkinv(fit))
}

p1_p2_to_p12=function(d1, d2) {
  p12
}

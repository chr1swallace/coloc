#' @import data.table
## helper functions
est.sxx <- function(b,vb,fx,v,N) {
    Sx <- 2*N*fx
    ((vb * N + b ^ 2) * Sx ^ 2 + v * N ^ 2) / (vb * N + b ^ 2) / N
}
est.Sxy <- function(b1,Sx,Sxx,N) {
    (N * Sxx * b1 - Sx ^ 2 * b1) / N
} 
##' reconstruct joint model ioutput for single snp. only vb1b2 really needed
##'
##' The SNP is X, the traits are 1 and 2
##' 
##' note in this case, the coefficients, and their variances, are
##' equal in single or joint models.  it is just the var(b1,b2) term
##' we need.
##'
##' This provides cor(b1,b2), and must be multiplied by
##' sqrt(v(b1)*v(b2)) to get var(b1,b2)
##'
##' @title qq.onesnp
##' @param b1 coef for snp X, trait 1
##' @param b2 coef for snp X, trait 2
##' @param vb1 var(coef) for snp X, trait 1
##' @param vb2 var(coef) for snp X, trait 2
##' @param fx MAF of X
##' @param v1 variance of trait 1
##' @param v2 variance of trait 2
##' @param eta cor(trait1, trait2)
##' @param N number of samples
##' @export
##' @return var(b1,b2)
##' @author Chris Wallace
qq.onesnp <- function(b1,b2,vb1,vb2,fx,v1,v2,eta,N) {
    Sx <- 2*N*fx
    Sxx <- est.sxx(b1,vb1,fx,v1,N) #2*N*fx*(1+fx)
    Sy1 <- Sy2 <- 0
    Sxy1 <- est.Sxy(b1,Sx,Sxx,N)
    Sxy2 <- est.Sxy(b2,Sx,Sxx,N)
    vbr=(4 * N ^ 2 * sqrt(v1 * v2) * eta * fx ^ 2 - N * sqrt(v1 * v2) * Sxx * eta + Sxy2 * Sxy1) / (4 * N * fx ^ 2 - Sxx) * ((4 * N ^ 2 * fx ^ 2 * v1 - N * Sxx * v1 + Sxy1 ^ 2) * (4 * N ^ 2 * fx ^ 2 * v2 - N * Sxx * v2 + Sxy2 ^ 2) / N ^ 2 / (4 * N * fx ^ 2 - Sxx) ^ 2) ^ (-0.5) / N
}
##' reconstruct joint model ioutput for two snp models
##'
##' The SNPs are X and Z, the traits are 1 and 2
##'
##' @title qq.twosnp
##' @inheritParams qq.onesnp
##' @param g1 coef for snp Z, trait 1
##' @param g2 coef for snp Z, trait 2
##' @param vg1 var(coef) for snp Z, trait 1
##' @param vg2 var( coef ) for snp Z, trait 2
##' @param fz MAF of Z
##' @param rho cor(X,Z)
##' @param eta cor(Y1,Y2)
##' @export
##' @return list of coefficients and their (co-)variances that would
##'     be expected from a joint model
##' @author Chris Wallace
qq.twosnp <- function(b1,b2,g1,g2,vb1,vb2,vg1,vg2,fx,fz,rho,v1,v2,eta,N) {
    Sz <- 2*N*fz
    Sx <- 2*N*fx
    Sxx <- est.sxx(b1,vb1,fx,v1,N) #2*N*fx*(1+fx)
    Szz <- est.sxx(g1,vg1,fz,v1,N) #2*N*fx*(1+fx)
    ## Sxz <- N*(2*rho*sqrt(fx*(1-fx)*fz*(1-fz)) + 4*fx*fz)
    Sxz <- (rho * N ^ 2 * sqrt((N * Sxx - Sx ^ 2) * (Szz * N - Sz ^ 2) / N ^ 4) + Sx * Sz) / N
    Sy1 <- Sy2 <- 0
    Sy1y1 <- N*v1
    Sy2y2 <- N*v2
    Sy1y2 <- N * eta * sqrt(v1 * v2)
    Sxy1 <- est.Sxy(b1,Sx,Sxx,N)
    Sxy2 <- est.Sxy(b2,Sx,Sxx,N)
    Szy1 <- est.Sxy(g1,Sz,Szz,N)
    Szy2 <- est.Sxy(g2,Sz,Szz,N)
    
    b1star <-  (N * Sxy1 * Szz - N * Sxz * Szy1 + Sx * Sz * Szy1 - Sxy1 * Sz ^ 2) / (N * Szz * Sxx - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sx * Sxz * Sz - Sxx * Sz ^ 2)
    g1star <- (N * Sxx * Szy1 - N * Sxy1 * Sxz - Sx ^ 2 * Szy1 + Sx * Sxy1 * Sz) / (N * Szz * Sxx - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sx * Sxz * Sz - Sxx * Sz ^ 2)
    b2star <-  (N * Sxy2 * Szz - N * Sxz * Szy2 + Sx * Sz * Szy2 - Sxy2 * Sz ^ 2) / (N * Szz * Sxx - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sx * Sxz * Sz - Sxx * Sz ^ 2)
    g2star <- (N * Sxx * Szy2 - N * Sxy2 * Sxz - Sx ^ 2 * Szy2 + Sx * Sxy2 * Sz) / (N * Szz * Sxx - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sx * Sxz * Sz - Sxx * Sz ^ 2)

    GG1xx <-  (Szz * N - Sz ^ 2) / (N * Sxx * Szz - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sz * Sx * Sxz - Sz ^ 2 * Sxx)
    GG1zz <-(N * Sxx - Sx ^ 2) / (N * Sxx * Szz - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sz * Sx * Sxz - Sz ^ 2 * Sxx) 
    GG1xz <-  -(N * Sxz - Sx * Sz) / (N * Sxx * Szz - N * Sxz ^ 2 - Sx ^ 2 * Szz + 2 * Sz * Sx * Sxz - Sz ^ 2 * Sxx)


    S11 <- (Sy1y1 -  ((Sxx * Szy1 ^ 2 + Sxy1 ^ 2 * Szz - 2 * Sxy1 * Sxz * Szy1) * N - (Sx * Szy1 - Sxy1 * Sz) ^ 2) / ((Sxx * Szz - Sxz ^ 2) * N - Sx ^ 2 * Szz + 2 * Sz * Sx * Sxz - Sz ^ 2 * Sxx))/N
    S22 <- (Sy2y2 - ((Sxx * Szy2 ^ 2 + Sxy2 ^ 2 * Szz - 2 * Sxy2 * Sxz * Szy2) * N - (Sx * Szy2 - Sxy2 * Sz) ^ 2) / ((Sxx * Szz - Sxz ^ 2) * N - Sx ^ 2 * Szz + 2 * Sz * Sx * Sxz - Sz ^ 2 * Sxx))/N
    S12 <- (Sy1y2 - (((Sxy2 * Szz - Sxz * Szy2) * Sxy1 + Szy1 * (Sxx * Szy2 - Sxy2 * Sxz)) * N - (Sx * Szy1 - Sxy1 * Sz) * (Sx * Szy2 - Sxy2 * Sz)) / ((Sxx * Szz - Sxz ^ 2) * N - Sx ^ 2 * Szz + 2 * Sz * Sx * Sxz - Sz ^ 2 * Sxx))/N
    
    list(b1=b1star,b2=b2star,g1=g1star,g2=g2star,
         vb1=GG1xx * S11,
         vb2=GG1xx * S22,
         vb1b2=GG1xx * S12,
         vg1=GG1zz * S11,
         vg2=GG1zz * S22,
         vg1g2=GG1zz*S12,
         vb1g1=GG1xz*S11,
         vb2g2=GG1xz*S22,
         vb1g2=GG1xz*S12)
}
##' Calculate the relative likelihood of a vector of regression
##' coefficients given their var-covariance matrix under a prior that
##' the true effect is 0 or is unknown, sampled from a normal with
##' known variance (W)
##'
##' This is for single snp models. Assume traits are 1 and 2
##' @title Bayes factor for one snp model with t distribution
##' @param b1 coefficient snp A trait 1
##' @param b2 coefficient snp A trait 2
##' @param v1 var(b1)
##' @param v2 var(b2)
##' @param rho cor(b1,b2)
##' @param w1 prior variance for trait 1
##' @param w2 prior variance for trait 2
##' @param N number of samples
##' @return log ABF
##' @author Chris Wallace
tbf2 <- function(b1,b2,v1,v2,rho,w1=0.04,w2=0.02,N) {
    v12 <- sqrt(v1)*sqrt(v2)*rho
    innerprod.V <- (b1 ^ 2 * v2 - 2 * b1 * b2 * v12 + b2 ^ 2 * v1) / (v1 * v2 - v12 ^ 2) 
    innerprod.VW <-((v2 + w2) * b1 ^ 2 - 2 * b1 * b2 * v12 + b2 ^ 2 * (v1 + w1)) / ((v2 + w2) * v1 + w1 * v2 + w1 * w2 - v12 ^ 2)
    p <- 2
    nu <- N-p
    det.V <- v1 * v2 - v12 ^ 2
    det.VW <- v1 * v2 + v1 * w2 - v12 ^ 2 + w1 * v2 + w1 * w2
    ll.V <- -(nu+p) * log(1 + innerprod.V/nu)/2 - log(det.V)/2
    ll.VW <- -(nu+p) * log(1 + innerprod.VW/nu)/2 - log(det.VW)/2
    ll.VW - ll.V
}

##' Calculate the relative likelihood of a vector of regression
##' coefficients given their var-covariance matrix under a prior that
##' the true effect is 0 or is unknown, sampled from a normal with
##' known variance (W)
##'
##' Assume SNPs are A and B, traits are 1 and 2
##' @title Bayes factor for two snp model with t distribution
##' @param bA1 coefficient snp A trait 1
##' @param bB1 coefficient snp B trait 1
##' @param bA2 coefficient snp A trait 2
##' @param bB2 coefficient snp B trait 2
##' @param vA1 var(bA1)
##' @param vA2 var(bA2)
##' @param vA12 var(bA1,bA2)
##' @param vB1 var(bB1)
##' @param vB2 var(bB2)
##' @param vB12 var(bB1,bB2)
##' @param vAB1 var(bA1,bB1)
##' @param vAB2 var(bA2,bB2)
##' @param vAB12 var(bA1,bB2) = var(bB1,bA2)
##' @param W1 prior variance for trait 1
##' @param W2 prior variance for trait 1
##' @param N number of samples
##' @return log ABF
##' @author Chris Wallace
tbf4 <- function(bA1,bB1,bA2,bB2,vA1,vA2,vA12,vB1,vB2,vB12,vAB1,vAB2,vAB12,W1,W2,N) {
    innerprod.V <-  ((-bA1 ^ 2 * vB2 + 2 * bA1 * bB2 * vAB12 - bB2 ^ 2 * vA1) * vAB12 ^ 2 + (2 * bA2 * bB1 * vAB12 ^ 2 + (-2 * bB1 * bB2 * vA12 - 2 * bA2 * bB2 * vAB1 - 2 * bA1 * (bA2 * vB12 + bB1 * vAB2)) * vAB12 + (2 * bB2 ^ 2 * vAB1 - 2 * bA1 * (-bB1 * vB2 + bB2 * vB12)) * vA12 - 2 * bA1 * (-bA2 * vB2 + bB2 * vAB2) * vAB1 + (2 * bA1 ^ 2 * vB12 + 2 * bB1 * bB2 * vA1) * vAB2 + 2 * bA2 * vA1 * (-bB1 * vB2 + bB2 * vB12)) * vAB12 + (-bA2 ^ 2 * vB1 - bB1 ^ 2 * vA2) * vAB12 ^ 2 + ((2 * vAB2 * bB1 ^ 2 - 2 * bA2 * (bB1 * vB12 - bB2 * vB1)) * vA12 + (2 * bA2 ^ 2 * vB12 - 2 * bA2 * bB1 * vAB2 + 2 * bB1 * bB2 * vA2) * vAB1 + 2 * (vAB2 * bA2 * vB1 + vA2 * (bB1 * vB12 - bB2 * vB1)) * bA1) * vAB12 + (-bB1 ^ 2 * vB2 + 2 * bB1 * bB2 * vB12 - bB2 ^ 2 * vB1) * vA12 ^ 2 + ((-2 * vAB2 * bB1 * bB2 - 2 * bA2 * (-bB1 * vB2 + bB2 * vB12)) * vAB1 - 2 * bA1 * ((bB1 * vB12 - bB2 * vB1) * vAB2 - bA2 * (-vB1 * vB2 + vB12 ^ 2))) * vA12 + (-bA2 ^ 2 * vB2 + 2 * bA2 * bB2 * vAB2 - bB2 ^ 2 * vA2) * vAB1 ^ 2 + 2 * bA1 * (vAB2 ^ 2 * bB1 - vB12 * bA2 * vAB2 + vA2 * (-bB1 * vB2 + bB2 * vB12)) * vAB1 + (-bA1 ^ 2 * vB1 - bB1 ^ 2 * vA1) * vAB2 ^ 2 + 2 * bA2 * vA1 * (bB1 * vB12 - bB2 * vB1) * vAB2 + (-bA1 ^ 2 * vA2 - bA2 ^ 2 * vA1) * vB12 ^ 2 - 2 * bB1 * bB2 * vA1 * vA2 * vB12 + bA1 ^ 2 * vA2 * vB1 * vB2 + vA1 * (bA2 ^ 2 * vB1 * vB2 + vA2 * (bB1 ^ 2 * vB2 + bB2 ^ 2 * vB1))) / ((-vA1 * vB2 + vAB12 ^ 2) * vAB12 ^ 2 + ((-2 * vA12 * vB12 - 2 * vAB1 * vAB2) * vAB12 + 2 * vA12 * vAB1 * vB2 + 2 * vA1 * vAB2 * vB12) * vAB12 - vA2 * vAB12 ^ 2 * vB1 + (2 * vA12 * vAB2 * vB1 + 2 * vA2 * vAB1 * vB12) * vAB12 + (-vB1 * vB2 + vB12 ^ 2) * vA12 ^ 2 - 2 * vA12 * vAB1 * vAB2 * vB12 + (-vA2 * vB2 + vAB2 ^ 2) * vAB1 ^ 2 - (vAB2 ^ 2 * vB1 + vA2 * (-vB1 * vB2 + vB12 ^ 2)) * vA1)
    innerprod.VW <- (((vB2 + W2) * bA1 ^ 2 - 2 * bA1 * bB2 * vAB12 + bB2 ^ 2 * (vA1 + W1)) * vAB12 ^ 2 + (((-2 * vAB12 ^ 2 + 2 * (vB2 + W2) * (vA1 + W1)) * bB1 + ((-2 * vB2 - 2 * W2) * bA1 + 2 * vAB12 * bB2) * vAB1 - 2 * vB12 * (-vAB12 * bA1 + bB2 * (vA1 + W1))) * bA2 + (((-2 * vB2 - 2 * W2) * bA1 + 2 * vAB12 * bB2) * vA12 - 2 * (-vAB12 * bA1 + bB2 * (vA1 + W1)) * vAB2) * bB1 - 2 * (-bA1 * vB12 + bB2 * vAB1) * (-bA1 * vAB2 + bB2 * vA12)) * vAB12 + ((vB2 + W2) * vAB1 ^ 2 - 2 * vAB1 * vAB12 * vB12 + (vA1 + W1) * vB12 ^ 2 - (-vAB12 ^ 2 + (vB2 + W2) * (vA1 + W1)) * vB1) * bA2 ^ 2 + ((((-2 * vB2 - 2 * W2) * vAB1 + 2 * vB12 * vAB12) * vA12 - 2 * (vB12 * (vA1 + W1) - vAB1 * vAB12) * vAB2) * bB1 + (2 * vB12 * bB2 * vAB1 - 2 * bA1 * vB12 ^ 2 + 2 * ((vB2 + W2) * bA1 - vAB12 * bB2) * vB1) * vA12 + 2 * (-bB2 * vAB1 ^ 2 + vB12 * bA1 * vAB1 + (-vAB12 * bA1 + bB2 * (vA1 + W1)) * vB1) * vAB2) * bA2 + ((vB2 + W2) * vA12 ^ 2 - 2 * vA12 * vAB12 * vAB2 + (vA1 + W1) * vAB2 ^ 2 - (-vAB12 ^ 2 + (vB2 + W2) * (vA1 + W1)) * vA2) * bB1 ^ 2 + (-2 * vB12 * bB2 * vA12 ^ 2 + 2 * vAB2 * (bA1 * vB12 + bB2 * vAB1) * vA12 + (-2 * bA1 * vAB2 ^ 2 + 2 * ((vB2 + W2) * bA1 - vAB12 * bB2) * vA2) * vAB1 + 2 * vB12 * (-vAB12 * bA1 + bB2 * (vA1 + W1)) * vA2) * bB1 + bB2 ^ 2 * vA12 ^ 2 * vB1 - 2 * bA1 * bB2 * vA12 * vAB2 * vB1 + bB2 ^ 2 * vA2 * vAB1 ^ 2 - 2 * bA1 * bB2 * vA2 * vAB1 * vB12 + bA1 ^ 2 * vAB2 ^ 2 * vB1 - (-bA1 ^ 2 * vB12 ^ 2 + ((vB2 + W2) * bA1 ^ 2 - 2 * bA1 * bB2 * vAB12 + bB2 ^ 2 * (vA1 + W1)) * vB1) * vA2) / ((-vAB12 ^ 2 + (vB2 + W2) * (vA1 + W1)) * vAB12 ^ 2 + (((-2 * vB2 - 2 * W2) * vAB1 + 2 * vB12 * vAB12) * vA12 - 2 * (vB12 * (vA1 + W1) - vAB1 * vAB12) * vAB2) * vAB12 + (-vB12 ^ 2 + vB1 * (vB2 + W2)) * vA12 ^ 2 + 2 * vAB2 * (vAB1 * vB12 - vAB12 * vB1) * vA12 + (-vAB2 ^ 2 + vA2 * (vB2 + W2)) * vAB1 ^ 2 - 2 * vA2 * vAB1 * vAB12 * vB12 + vB1 * (vA1 + W1) * vAB2 ^ 2 - ((-W1 - vA1) * vB12 ^ 2 + (-vAB12 ^ 2 + (vB2 + W2) * (vA1 + W1)) * vB1) * vA2)
 det.V <- (-vA1 * vB2 + vAB12 ^ 2) * vAB12 ^ 2 + ((-2 * vA12 * vB12 - 2 * vAB1 * vAB2) * vAB12 + 2 * vA12 * vAB1 * vB2 + 2 * vA1 * vAB2 * vB12) * vAB12 - vA2 * vAB12 ^ 2 * vB1 + (2 * vA12 * vAB2 * vB1 + 2 * vA2 * vAB1 * vB12) * vAB12 + (-vB1 * vB2 + vB12 ^ 2) * vA12 ^ 2 - 2 * vA12 * vAB1 * vAB2 * vB12 + (-vA2 * vB2 + vAB2 ^ 2) * vAB1 ^ 2 + (-vAB2 ^ 2 * vB1 + (vB1 * vB2 - vB12 ^ 2) * vA2) * vA1
det.VW <- ((-W2 - vB2) * vA1 - W1 * W2 - W1 * vB2 + vAB12 ^ 2) * vAB12 ^ 2 + (((2 * W2 + 2 * vB2) * vAB1 - 2 * vB12 * vAB12) * vA12 + 2 * (vB12 * (vA1 + W1) - vAB1 * vAB12) * vAB2) * vAB12 + (((vB2 + W2) * vA1 + W1 * W2 + W1 * vB2 - vAB12 ^ 2) * vB1 + (-W2 - vB2) * vAB1 ^ 2 + 2 * vAB1 * vAB12 * vB12 - (vA1 + W1) * vB12 ^ 2) * vA2 + ((-W2 - vB2) * vA12 ^ 2 + 2 * vA12 * vAB12 * vAB2 - (vA1 + W1) * vAB2 ^ 2) * vB1 + (vA12 * vB12 - vAB1 * vAB2) ^ 2
    p <- 4
    nu <- N-p
    ll.V <- -(nu+p) * log(1 + innerprod.V/nu)/2 - log(det.V)/2
    ll.VW <- -(nu+p) * log(1 + innerprod.VW/nu)/2 - log(det.VW)/2
    ll.VW - ll.V
}


##' Bayesian colocalisation analysis for two quantitative traits, measured on the same individuals
##'
##' This function calculates posterior probabilities of different
##' causal variant configurations under the assumption of a single
##' causal variant for each trait.
##'
##' Unlike coloc.abf, this specifically requires coefficients and
##' their standard errors - ie p values and MAF are not enough,
##' because the direction of effect matters
##'
##' @param N number of samples
##' @param VY1 variance of trait 1
##' @param VY2 variance of trait 2
##' @param corY cor(Y1,Y2)
##' @param MAF minor allele frequency of the variants, named vector
##' @param LD matrix giving correlation between SNPs.  NB - this is r,
##'     not rsquared!  Expect a mix of positive and negative values.
##'     This should be estimated from reference data.
##' @param p1 prior probability a SNP is associated with trait 1,
##'     default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2,
##'     default 1e-4
##' @param p12 prior probability a SNP is associated with both traits,
##'     default 1e-5
##' @param W1 the effect size of any variant for trait 1 under the
##'     assumption it is causal is assumed to follow a normal with
##'     mean 0 and variance W1
##' @param W2 as above, for trait 2
##' @return a list, with main entry "results" a vector giving the
##'     number of SNPs analysed, and the posterior probabilities of H0
##'     (no causal variant), H1 (causal variant for trait 1 only), H2
##'     (causal variant for trait 2 only), H3 (two distinct causal
##'     variants) and H4 (one common causal variant)
##' @export
##' @author Chris Wallace
coloc.qq <- function(N, VY1, VY2, corY,
                     MAF, LD, 
                     p1=1e-4, p2=1e-4, p12=1e-5,
                     W1=0.15^2*VY1, W2=0.15^2*VY2) {
    erho <- b1 <- b2 <- v1 <- v2 <- lbf1 <- lbf2 <- lbf3 <- lbf4 <- f0 <- snpA <- snpB <- isnpA <- isnpB <- r <- bA1 <- bB1 <- bA <- bB2 <- vA1 <- vB1 <- vA2 <- vB2 <- fa <- fB <- vA12 <- vB12 <- vAB1 <- vAB2 <- vAB12 <- .SD  <- NULL

    snps <- intersect(dataset1$snp,dataset2$snp)
    df <- data.table(snp=c(snps))
    df <- merge(df,as.data.table(as.data.frame(dataset1[c("snp","beta","varbeta")])),
                by.x="snp",by.y="snp",all.x=TRUE)
    setnames(df,c("beta","varbeta"),c("b1","v1"))
    df <- merge(df,as.data.table(as.data.frame(dataset2[c("beta","varbeta","snp")])),
                by.x="snp",by.y="snp",all.x=TRUE) 
    setnames(df,c("beta","varbeta"),c("b2","v2"))
    fq <- data.table(snp=snps)
    fq[,f0:=MAF[snp]]
    df <- merge(df,fq[,.(snp,f0)],by.x="snp",by.y="snp",all.x=TRUE)
    
    ## H1, H2, H4 all two beta models
    ## estimate cor(b1,b2)
    df[,erho:=qq.onesnp(b1,b2,v1,v2,f0,VY1,VY2,corY,N)]

    ## t distribution
    df[, lbf1:=tbf2(b1,b2,v1,v2,erho,w1=W1,w2=0,N=N)]
    df[, lbf2:=tbf2(b1,b2,v1,v2,erho,w1=0,W2,N=N)]
    df[, lbf4:=tbf2(b1,b2,v1,v2,erho,w1=W1,w2=W2,N=N)]
    
    ## normal approximation to t
    ## df[, lbf1:=bf2(b1,b2,v1,v2,erho,w1=W1,w2=0)]
    ## df[, lbf2:=bf2(b1,b2,v1,v2,erho,w1=0,W2)]
    ## df[, lbf4:=bf2(b1,b2,v1,v2,erho,w1=W1,w2=W2)]
    
    ## for h3, we need pairs of snps
    df3 <- as.data.table(as.data.frame(expand.grid(snps,snps)))
    setnames(df3,c("snpA","snpB"))
    df3 <- df3[snpA!=snpB,]
    df3[,isnpA:=match(as.character(snpA),rownames(LD))]
    df3[,isnpB:=match(as.character(snpB),rownames(LD))]
    df3[,r:=LD[cbind(isnpA,isnpB)]]
    df3 <- merge(df3,df[,.(snp,f0,b1,v1)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 1, snp A
    setnames(df3,c("f0","b1","v1"),c("fA","bA1","vA1"))
    df3 <- merge(df3,df[,.(snp,f0,b1,v1)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 1, snp B
    setnames(df3,c("f0","b1","v1"),c("fB","bB1","vB1"))
    df3 <- merge(df3, df[,.(snp,b2,v2)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 2, snp A
    setnames(df3,c("b2","v2"),c("bA2","vA2"))
    df3 <- merge(df3, df[,.(snp,b2,v2)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 2, snp B
    setnames(df3,c("b2","v2"),c("bB2","vB2"))

    ## special case: r==1. This won't allow a model to be fitted at all.  Instead, add in the bf from h4 for one of the two SNPs
    df3.r1 <- df3[abs(r)>=0.99,]
    m <- match(df3.r1$snpA, df$snp)
    df3.r1[, lbf3:=df$lbf4[m]]

    ## special case: r=0. Can fit two bivariate normals, which means we can vectorize
    df3.r0 <- df3[abs(r)<0.001,]
    mA <- match(df3.r0$snpA, df$snp)
    mB <- match(df3.r0$snpB, df$snp)
    df3.r0[, lbf3:=df$lbf1[mA] + df$lbf2[mB] ]

    ## what's left requires 4 dim Gaussians
    df3 <- df3[abs(r)<0.99 & abs(r)>=0.001,]
    df3[,W1:=W1]
    df3[,W2:=W2]
    df3[,c("bA1old","bB1old","bA2old","bB2old"):=list(bA1,bB1,bA2,bB2)]
    df3[,c("bA1","bA2","bB1","bB2","vA1","vA2","vA12","vB1","vB2","vB12","vAB1","vAB2","vAB12"):=
           qq.twosnp(b1=bA1,b2=bA2,g1=bB1,g2=bB2,vb1=vA1,vb2=vA2,vg1=vB1,vg2=vB2,fx=fA,fz=fB,rho=r,v1=VY1,v2=VY2,eta=corY,N=N)]

    ## t distribution
    df3[,lbf3:=tbf4(bA1=bA1,bA2=bA2,bB1=bB1,bB2=bB2,
                            vA1=vA1,vA2=vA2,vA12=vA12,vB1=vB1,vB2=vB2,vB12=vB12,
                            vAB1=vAB1,vAB2=vAB2,vAB12=vAB12,
                            W1=W1,W2=W2,N=N)]
    ## Normal approximation to t
    ## df3[,lbf3:=lbfh3.qq.vec(bA1=bA1,bA2=bA2,bB1=bB1,bB2=bB2,
    ##                         vA1=vA1,vA2=vA2,vA12=vA12,vB1=vB1,vB2=vB2,vB12=vB12,
    ##                         vAB1=vAB1,vAB2=vAB2,vAB12=vAB12,
    ##                         W1=W1,W2=W2)]
    ## df3[,erho:=alt.twosnp(bA1,bB2,fA,fB,rho=r,eta=corY*sqrt(VY1*VY2),N=N)]
    ## df3[, lbf3:=bf2(bA1,bB2,vA1,vB2,erho,w1=W1,w2=W2)]
 
    

    ## ## reconstruct coefficients from bivariate models
    ## df3[,c("bA1old","bB1old","bA2old","bB2old"):=list(bA1,bB1,bA2,bB2)]
    ## df3[,c("bA1","bB1"):=bstar_from_b(bA1,bB1,vA1,vB1,r)]
    ## df3[,c("bA2","bB2"):=bstar_from_b(bA2,bB2,vA2,vB2,r)]
    ## df3[,VA:=2*fA*(1-fA)]
    ## df3[,VB:=2*fB*(1-fB)]
    ## df3[,VAB := r * sqrt(VA) * sqrt(VB)]
    ## df3[,sumA2:=2*n*fA*(1+fA)] #sum(bvec^2 * rowSums(geno))
    ## df3[,sumB2:=2*n*fB*(1+fB)] #sum(bvec^2 * colSums(geno))
    ## df3[,sumAB:=n*( sqrt(VA)*sqrt(VB)*r+4*fA*fB)] # sum( (bvec %*% t(bvec)) * geno )
    ## df3[,det:=sumA2*sumB2-sumAB^2]
    ## df3[,EE:=covY - bA1*bA2*VA - bB1*bB2*VB - (bA1*bB2 + bB1*bA2)*VAB]
    ## df3[,sumA2:=EE*sumA2/det]
    ## df3[,sumB2:=EE*sumB2/det]
    ## df3[,sumAB:=EE*sumAB/det]
    ## df3[,W1:=VY1*0.15^2]
    ## df3[,W2:=VY2*0.15^2]

    ## wh <- head(which(df3$r < 0.99 & df3$snpA < df3$snpB),50)
    ## microbenchmark({for(i in wh) {
    ##                     tmp[[i]] <- with(df3[i,],
    ##                                 lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
    ##                                            vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
    ##                                            W1=W1,W2=W2,sumA2,sumAB,sumB2,r))})

     ## for(i in 1:nrow(df3)) { #& df3$snpA < df3$snpB)) {
        ## df3[i,lbf3:=lbfh3.qq(bA1=bA1,bA2=bA2,bB1=bB1,bB2=bB2,
        ##                        vA1=vA1,vA2=vA2,vA12=vA12,vB1=vB1,vB2=vB2,vB12=vB12,
        ##                        vAB1=vAB1,vAB2=vAB2,vAB12=vAB12,
        ##                        W1=W1,W2=W2)]
## }
        ## microbenchmark({tmp <- with(df3[i,],
        ##                             lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
        ##                                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
        ##                                        W1=W1,W2=W2,sumA2,sumAB,sumB2,r))})
        ## microbenchmark({df3[i,lbf3:=lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
        ##                                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
        ##                                        W1=W1,W2=W2,sumA2,sumAB,sumB2,r)]})
        ## microbenchmark({df3[i,lbf3:= lbf.h3.qq2(bA1=df3$bA1[i],bB1=df3$bB1[i],bA2=df3$bA2[i],bB2=df3$bB2[i],
        ##                                        vA1=df3$vA1[i],vB1=df3$vB1[i],vA2=df3$vA2[i],vB2=df3$vB2[i],
        ##                                        W1=df3$W1[i],W2=df3$W2[i],df3$sumA2[i],df3$sumAB[i],df3$sumB2[i],df3$r[i])]})
        ## j <- which(df3$snpA==df3$snpB[i] & df3$snpB==df3$snpA[i])
        ## microbenchmark({tmpj <- with(df3[j,],
        ##              lbf.h3.qq(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
        ##                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
        ##                        n,
        ##                        VY1,VY2,covY,fA,fB,r,W=0.04))})
        ## df3[i,lbf3:=tmp$asis]
        ## df3[j,lbf3:=tmp$reversed]

    ## put it all together
    o<-intersect(names(df3),names(df3.r1))
    df3 <- rbind(df3.r1[,o,with=FALSE],df3.r0[,o,with=FALSE],df3[,o,with=FALSE])
    sdf <- df[, lapply(.SD, logsum), .SDcols=grep("lbf",names(df),value=TRUE)]
    sdf3 <- df3[, lapply(.SD, logsum), .SDcols=grep("lbf",names(df3),value=TRUE)]
    sdf <- rbind(melt(sdf,measure.vars=names(sdf)),
                 melt(sdf3,measure.vars=names(sdf3)))
    ## add h0
    sdf0 <- melt(sdf3,measure.vars=names(sdf3))
    sdf0[,value:=0]
    sdf0[,variable:=sub("3","0",variable)]
    sdf <- rbind(sdf,sdf0)
    
    sdf[,h:=as.integer(sub(".*lbf","",variable))]
    sdf[,method:=sub("[0-4]","",variable)]
    sdf[,abf:=0]
    sdf[h==1,abf:=log(p1) + value]
    sdf[h==2,abf:=log(p2) + value]
    sdf[h==3,abf:=log(p1) + log(p2) + value]
    sdf[h==4,abf:=log(p12) + value]
    sdf[,pp:=exp(abf - logsum(abf)),by="method"]
    sdf <- sdf[order(method,h),]
    ## while checking, return this
    ## return(sdf)
    
    ret <- c(nsnps=length(snps), sdf[method=="lbf",]$pp)
    names(ret) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", 
                    "PP.H4.abf")
    print(ret)
    return(list(summary=ret,results=df,df3=df3,sdf=sdf))
}
 




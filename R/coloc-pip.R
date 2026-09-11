##' @title convert CLPP to coloc posteriors
##' @param clPP CLPP value
##' @inheritParams coloc.abf
##' @param alpha1 size of credible set 1 (the sum of PP associated with variants in the credible set)
##' @param alpha2 size of credible set 2
##' @param h4_only set TRUE if you only want PP.H4 returned
##' @return a named vector containing the posterior probabilities of H3 and H4
##' @export 
##' @author Jeffrey Pullin, Chris Wallace
clPP_to_coloc <- function(clPP, alpha1, alpha2, p1=1e-4, p2=1e-4, p12=5e-6, h4_only=FALSE) {
    stopifnot("p1, p2, p12 must be probabilities"=is_prob(p1) && is_prob(p2) && is_prob(p12))
    stopifnot("clPP must be in [0,1]"=is_prob(clPP,strict=FALSE))
    stopifnot("alpha1, alpha2 must be in [0,1]"=is_prob(clPP,strict=FALSE))
    pr_h4 <- p12 * clPP / ((p12 - p1*p2) * clPP + p1*p2 * (alpha1 * alpha2))
    if(!h4_only)
        return(c(PP.H3=1-pr_h4, PP.H4=pr_h4))
    pr_h4
}

##' An aPProximation to coloc using only information from fine-maPPed credible sets.
##'
##' A credible set should be passed in the form of a list, data.frame,
##' or equivalent tabular object containing at least two equal-length
##' vectors.  The object should contain at a minimim
##'
##' - snp: a character vector of snp identifiers
##' - PP: a numeric vector of posterior probabilities
##'
##' @title coloc using credible sets
##' @param cs1 credible set 1
##' @param cs2 credible set 2
##' @inheritParams coloc.abf
##' @return a named vector containing the number of overlaPPing SNPs and the posterior probabilities of H3 and H4
##' @export 
##' @author Jeffrey Pullin, Chris Wallace
coloc_pip <- function(cs1, cs2, p1=1e-4, p2=1e-4, p12=5e-6) {
    check_PP(cs1)
    check_PP(cs2)
    stopifnot("p1, p2, p12 must be probabilities"=is_prob(p1) && is_prob(p2) && is_prob(p12))
    alpha1 <- sum(cs1$PP)
    alpha2 <- sum(cs2$PP)

    cs1 <- as.data.frame(cs1[c("snp","PP")])
    cs2 <- as.data.frame(cs2[c("snp","PP")])
    
    m <- merge(cs1,cs2,by="snp")
    if (!nrow(m)) {
        pr_h4 <- 0
    } else {
        clPP <- sum(m$PP.x * m$PP.y)
        pr_h4 <- clPP_to_coloc(clPP=clPP, p1=p1, p2=p2, p12=p12, alpha1=alpha1, alpha2=alpha2, h4_only=TRUE)
    }
    c(nsnps=nrow(m), PP.H3=1-pr_h4, PP.H4=pr_h4)
}


check_PP <- function(cs) {
    nm <- names(cs)
    if(!all(c("snp","PP") %in% nm))
        stop("credible set must include snp and PP")
    if(!is.character(cs$snp))
        stop("snp must be a character vector")
    if(!is.numeric(cs$PP))
        stop("PP must be a numeric vector")
    if(length(cs$PP) != length(cs$snp))
        stop("PP and snp should be equal length")
    NULL
}

is_prob <- function(x,strict=TRUE) {
    if(strict) {
        is.numeric(x) && x > 0 && x < 1 && length(x)==1
    } else {
        is.numeric(x) && x >= 0 && x <= 1 && length(x)==1
    }
}

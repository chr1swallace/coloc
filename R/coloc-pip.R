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
    stopifnot("clPP must be in [0,1]"=is_prob(clPP))
    stopifnot("alpha1, alpha2 must be in [0,1]"=is_prob(clPP))
    pr_h4 <- p12 * clPP / ((p12 - p1*p2) * clPP + p1*p2 * (alpha1 * alpha2))
    if(!h4_only)
        return(c(PP.H3=1-pr_h4, PP.H4=pr_h4))
    pr_h4
}

##' An approximation to coloc using only information from fine-mapped credible sets.
##'
##' A credible set should be passed in the form of a list, data.frame,
##' or equivalent tabular object containing at least two equal-length
##' vectors.  The object should contain at a minimim
##'
##' - snp: a character vector of snp identifiers
##' - PP: a numeric vector of posterior probabilities
##'
##' Each cs may additionally contain a column
##'
##' - beta: the estimated per all snp effect on the trait
##'
##' If they both do, then the returned object will also include
##' "freq_sign_agree", which is the proportion of snps in both
##' credible sets for which the difrection of effects (sign of beta)
##' agrees. This is a useful quantity to be able to say
##' "higher gene expression associates with higher (lower) disease risk"
##' when freq_sign_agree=1 (0). If 0 < freq_sign_agree < 1 then there
##' is some disagreement between SNPs and users are invited to
##' investigate potential causes of this further.
##'
##' @param cs1 credible set 1
##' @param cs2 credible set 2
##' @inheritParams coloc.abf
##' @return a named vector containing the number of SNPs found in each
##'     credible set (nsnps1, nsnps2), the number in both credible
##'     sets (nsnps_int) and the posterior probabilities of H3 and H4
##'     (PP.H3 and PP.H4)
##' @export
##' @author Jeffrey Pullin, Chris Wallace
##' @examples
##' cs1 <- data.frame(snp=c("s1","s2"), PP=c(0.05,0.9))
##' cs2 <- data.frame(snp=c("s1","s2","s3"), PP=c(0.1,0.8,0.05))
##' coloc_pip(cs1,cs2)
coloc_pip <- function(cs1, cs2, p1=1e-4, p2=1e-4, p12=5e-6) {
    check_PP(cs1)
    check_PP(cs2)
    stopifnot("p1, p2, p12 must be probabilities"=is_prob(p1) && is_prob(p2) && is_prob(p12))
    alpha1 <- sum(cs1$PP)
    alpha2 <- sum(cs2$PP)

    if("beta" %in% names(cs1) && "beta" %in% names(cs2)) {
        cs1 <- as.data.frame(cs1)[c("snp","PP","beta")]
        cs2 <- as.data.frame(cs2)[c("snp","PP","beta")]
    } else {
        cs1 <- as.data.frame(cs1)[c("snp","PP")]
        cs2 <- as.data.frame(cs2)[c("snp","PP")]
    }
    
    m <- merge(cs1,cs2,by="snp")
    if (!nrow(m)) {
        pr_h4 <- 0
    } else {
        clPP <- sum(m$PP.x * m$PP.y)
        pr_h4 <- clPP_to_coloc(clPP=clPP, p1=p1, p2=p2, p12=p12, alpha1=alpha1, alpha2=alpha2, h4_only=TRUE)
    }
    ret <- c(nsnps1=nrow(cs1), nsnps2=nrow(cs2), nsnps_int=nrow(m), PP.H3=1-pr_h4, PP.H4=pr_h4)
    ret["freq_sign_agree"] <-    if(all(c("beta.x","beta.y") %in% names(m)) && nrow(m)) {
                                     mean(sign(m$beta.x)==sign(m$beta.y))
                                 } else {
                                     NA
                                 }
    ret
}


check_PP <- function(cs) {
    nm <- names(cs)
    if(!all(c("snp","PP") %in% nm))
        stop("credible set must include snp and PP")
    if(!is.character(cs$snp))
        stop("snp must be a character vector")
    if(any(duplicated(cs$snp)))
        stop("duplicated snps")
    if(!is.numeric(cs$PP))
        stop("PP must be a numeric vector")
    if(length(cs$PP) != length(cs$snp))
        stop("PP and snp should be equal length")
    if(!is_prob(sum(cs$PP)))
        stop("sum of cs$PP not in [0,1]")
    NULL
}

is_prob <- function(x) {
    length(x)==1 && !is.na(x) && is.numeric(x) &&
        (x >= 0 || all.equal(x, 0)) && # allows for numerical rounding errors
        (x <= 1 || all.equal(x, 1))
}

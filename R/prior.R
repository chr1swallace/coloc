## ---------------------------------------------------------------------------
## The prior and control objects.
##
## This is where prior drift stops. In 6.0.3 the same three probabilities are
## re-declared in every entry point signature and then diverge: coloc.abf()
## defaults p12 to 1e-5, coloc.bf_bf() to 5e-6, and coloc.susie_bf() declares
## p1/p2/p12 and never forwards them. The bounds check and the region-size
## guard exist in only some of those places.
##
## One object, validated once, passed down.
## ---------------------------------------------------------------------------

new_coloc_prior <- function(p1, p2, p12, weights1, weights2) {
    structure(list(p1 = p1, p2 = p2, p12 = p12,
                   weights1 = weights1, weights2 = weights2),
              class = "coloc_prior")
}

validate_coloc_prior <- function(x) {
    ## First, because the bounds check below is vector-tolerant (any()/pmin())
    ## and would otherwise accept a vector silently. p1 is the prior that any
    ## ONE snp is causal for trait 1 only: a single number applying to every
    ## snp. Varying the prior across snps is what weights1/weights2 are for,
    ## and combine_abf_weighted() builds Q * p * w/sum(w) from them. A vector
    ## p1 is a second spelling of the same thing, and the confusing one.
    .check_scalar_prior(x)
    if (any(x$p12 < x$p1 * x$p2) || any(x$p12 > pmin(x$p1, x$p2)))
        stop("p12 must lie in [p1*p2, min(p1, p2)]")
    ## The region-size guard needs the snp counts, which only exist once the
    ## inputs are adapted, so coloc() applies it then, as coloc.abf() does.
    for (w in list(x$weights1, x$weights2))
        if (!is.null(w) && any(is.na(w) | w <= 0))
            stop("prior weights must be non-NA and positive")
    if ((!is.null(x$weights1) || !is.null(x$weights2)) &&
        (length(x$p1) != 1 || length(x$p2) != 1 || length(x$p12) != 1))
        stop("prior weights can only be used with scalar p1, p2, p12")
    x
}

## The engine takes scalar priors only. coloc_prior() checks that through
## validate_coloc_prior(), but coloc() can also be handed a plain list, which
## is what the old entry points pass it, and coloc.bf_bf() calls the engine
## directly. Both check here, so the message is the same everywhere.
.check_scalar_prior <- function(x) {
    for (nm in c("p1", "p2", "p12"))
        if (length(x[[nm]]) != 1)
            stop("p1, p2 and p12 must each be a single number; ",
                 "to vary the prior across SNPs, use weights1/weights2")
    invisible(x)
}

##' Prior probabilities for a colocalisation analysis
##'
##' One object, validated once, in place of the same three probabilities being
##' re-declared in every entry point.
##'
##' `p1`, `p2` and `p12` are each the prior probability for any one SNP: that
##' it is causal for trait 1 only, trait 2 only, or both. Each is a single
##' number. To let the prior vary across SNPs, give weights, which scale it
##' SNP by SNP.
##'
##' @param p1 the prior probability for any one SNP that it is associated with
##'   trait 1 only
##' @param p2 the prior probability for any one SNP that it is associated with
##'   trait 2 only
##' @param p12 the prior probability for any one SNP that it is associated with
##'   both traits
##' @param weights1 optional: one weight per SNP of trait 1, positive, no `NA`
##' @param weights2 optional: one weight per SNP of trait 2, positive, no `NA`
##' @return an object of class coloc_prior
##' @export
coloc_prior <- function(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6,
                        weights1 = NULL, weights2 = NULL) {
    prior <- validate_coloc_prior(new_coloc_prior(p1, p2, p12, weights1, weights2))
    ## coloc.abf()'s advice (R/claudia.R:439-443), given here so it outlives
    ## coloc.abf() at step 7
    if (!is.null(weights1) && !is.null(weights2))
        warning("Prior weights specified for both traits.\n",
                "  The two weight vectors should be derived from independent sources of\n",
                "  information, i.e. not computed using the same method.")
    prior
}

##' Control parameters for a colocalisation analysis
##'
##' @param overlap.min minimum proportion of a signal's posterior that must
##'   fall on snps shared with the other trait before the pair is considered
##' @param trim_by_posterior whether to drop signal pairs failing overlap.min
##' @return an object of class coloc_control
##' @export
coloc_control <- function(overlap.min = 0.5, trim_by_posterior = TRUE) {
    stopifnot(is.numeric(overlap.min), length(overlap.min) == 1,
              overlap.min > 0, overlap.min <= 1,
              is.logical(trim_by_posterior), length(trim_by_posterior) == 1)
    structure(list(overlap.min = overlap.min, trim_by_posterior = trim_by_posterior),
              class = "coloc_control")
}

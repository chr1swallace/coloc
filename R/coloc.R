## ---------------------------------------------------------------------------
## The front door.
##
## Stages 1 and 2 and the three adapters, wired together. Additive: every
## existing entry point still works and still returns what it always did.
## ---------------------------------------------------------------------------

##' Colocalisation of two genetic traits
##'
##' One entry point for every combination of input coloc accepts. `x1` and
##' `x2` may each be a coloc dataset, a [runsusie()] fit, or a named vector or
##' matrix of log Bayes factors, in any combination.
##'
##' `method` applies PER TRAIT, and only to traits supplied as datasets. Give
##' it length 2 to fine-map one trait and not the other, which is the common
##' case when colocalising your own fine-mapped GWAS against eQTL Catalogue
##' SuSiE output, or whenever one trait has no LD panel: SuSiE cannot run
##' without LD, and a single global method would force the other trait onto
##' `"single"` as well, collapsing its credible sets into one assumed causal
##' variant. For a trait already supplied as log Bayes factors, `method` is
##' inapplicable and supplying it warns.
##'
##' @param x1 trait 1: a coloc dataset, a runsusie() fit, or log Bayes factors
##' @param x2 trait 2, same options
##' @param method `"single"` (one causal variant assumed, no fine-mapping) or
##'   `"susie"`; length 1 for both traits, length 2 for per-trait
##' @param prior a [coloc_prior()] object
##' @param control a [coloc_control()] object
##' @param susie.args passed to [runsusie()] for traits fine-mapped here
##' @return an object of class coloc_result
##' @export
coloc <- function(x1, x2,
                  method = c("single", "susie"),
                  prior = coloc_prior(),
                  control = coloc_control(),
                  susie.args = list()) {

    ## `method` is a match.arg default VECTOR, so an unsupplied method arrives
    ## here as c("single", "susie") and is indistinguishable from a deliberate
    ## per-trait request. Per-trait scope and match.arg-style defaults do not
    ## compose, so the default is collapsed explicitly.
    supplied <- !missing(method)
    if (!supplied) method <- "single"

    kinds  <- c(.classify_input(x1), .classify_input(x2))
    method <- .resolve_method(method, kinds, supplied = supplied)

    adapt <- function(x, kind, m) switch(kind,
        susie_fit = .adapt_susie(x, susie.args = susie.args),
        bf        = .adapt_bf(x),
        dataset   = if (identical(m, "susie")) .adapt_susie(x, susie.args = susie.args)
                    else .adapt_single(x))
    a1 <- adapt(x1, kinds[1], method[1])
    a2 <- adapt(x2, kinds[2], method[2])

    res <- .coloc_engine(a1$bf, a2$bf, prior = prior, control = control)

    ## The engine reports idx1/idx2 as row numbers of the TRIMMED matrices,
    ## which is not what a caller wants back, so they are mapped through
    ## meta$idx exactly as coloc.susie() does at R/susie.R:96-98. meta$idx is
    ## what makes that work the same way whichever adapter ran.
    if (nrow(res$summary)) {
        res$summary$idx1 <- a1$meta$idx[res$summary$idx1]
        res$summary$idx2 <- a2$meta$idx[res$summary$idx2]
    }
    ## $meta fills in two stages: the engine writes what it knows, and the
    ## front door adds what only it can know, since the engine never learns
    ## which adapter ran.
    res$meta$method <- c(a1$meta$method, a2$meta$method)
    res$meta$idx1 <- a1$meta$idx
    res$meta$idx2 <- a2$meta$idx
    validate_coloc_result(res)
}

## ---- accessors, so nobody reaches into internals to ask -------------------

##' Posterior probabilities from a coloc_result
##'
##' Returns the named vector that `coloc.abf()`'s `$summary` used to be, for
##' one row of the result.
##'
##' @param x a coloc_result
##' @param row which signal pair
##' @return named numeric, length 6
##' @export
pp <- function(x, row = 1) {
    stopifnot(inherits(x, "coloc_result"))
    if (!nrow(x$summary)) return(stats::setNames(rep(NA_real_, 6L),
                                                 c("nsnps", paste0("PP.H", 0:4, ".abf"))))
    unlist(x$summary[row, c("nsnps", paste0("PP.H", 0:4, ".abf")), with = FALSE])
}

##' Did a colocalisation analysis produce usable posteriors?
##'
##' `res$summary$PP.H4.abf > 0.9` yields `logical(0)` on an empty result, so
##' an `if()` around it errors. Ask this instead.
##'
##' @param x a coloc_result
##' @return logical(1)
##' @export
is_ok <- function(x) {
    stopifnot(inherits(x, "coloc_result"))
    identical(x$status, "ok")
}

##' How many signal pairs a coloc_result carries
##'
##' @param x a coloc_result
##' @return integer(1)
##' @export
n_signals <- function(x) {
    stopifnot(inherits(x, "coloc_result"))
    nrow(x$summary)
}

##' @export
print.coloc_result <- function(x, ...) {
    cat("coloc_result: status", x$status, "with", nrow(x$summary), "signal pair(s)\n")
    if (!is.null(x$meta$method))
        cat("  method:", paste(ifelse(is.na(x$meta$method), "bf", x$meta$method),
                               collapse = ", "), "\n")
    cat("  priors:", paste(names(x$priors), signif(x$priors, 3),
                           sep = " = ", collapse = ", "), "\n")
    if (nrow(x$summary)) print(x$summary)
    for (m in x$messages) cat("  note:", m, "\n")
    invisible(x)
}

## ---- adapting the old entry points ----------------------------------------

## The old doors take priors and control as loose arguments in `...`. Collect
## them here so one place knows the defaults. Deliberately plain lists rather
## than coloc_prior()/coloc_control(): those validate, and the old doors keep
## the permissiveness they have always had. coloc() is the validating door.
.prior_from_dots <- function(dots) {
    pick <- function(nm, default) if (is.null(dots[[nm]])) default else dots[[nm]]
    out <- list(p1 = pick("p1", 1e-4), p2 = pick("p2", 1e-4), p12 = pick("p12", 5e-6),
                weights1 = dots$prior_weights1, weights2 = dots$prior_weights2)

    .check_prior_weights(out)
    out
}

## The BF-route weight rules, in one place. coloc.susie() and
## coloc.susie_bf() used to inherit these by calling coloc.bf_bf(); they now
## call coloc(), so without this a negative weight would reach log() and
## produce NaN rather than a clean error.
##
## coloc.abf() deliberately does NOT use this: its own check permits a zero
## weight, where the BF route rejects one. That asymmetry is 6.0.3's, not a
## new one, and it disappears when the old doors do.
.check_prior_weights <- function(prior) {
    if ((length(prior$p1) != 1 || length(prior$p2) != 1 || length(prior$p12) != 1) &&
        (!is.null(prior$weights1) || !is.null(prior$weights2)))
        stop("Prior weights can only be used with single values of the prior parameters.")
    if (!is.null(prior$weights1) && any(is.na(prior$weights1) | (prior$weights1 <= 0)))
        stop("prior_weights1 must contain non-NA/non-negative weights.")
    if (!is.null(prior$weights2) && any(is.na(prior$weights2) | (prior$weights2 <= 0)))
        stop("prior_weights2 must contain non-NA/non-negative weights.")
    invisible(prior)
}

.control_from_dots <- function(dots) {
    pick <- function(nm, default) if (is.null(dots[[nm]])) default else dots[[nm]]
    list(overlap.min = pick("overlap.min", 0.5),
         trim_by_posterior = pick("trim_by_posterior", TRUE))
}

## coloc() warns that `method` is inapplicable to a precomputed BF input, which
## is right for someone calling coloc() directly. coloc.susie_bf()'s entire
## signature is "one dataset, one BF matrix", so there the warning is noise
## about a choice the caller was never offered.
.without_method_warning <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
        if (grepl("method is inapplicable", conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
    })
}

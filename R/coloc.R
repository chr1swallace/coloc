## ---------------------------------------------------------------------------
## The front door.
##
## Stages 1 and 2 and the three adapters, wired together. Additive: every
## existing entry point still works and still returns what it always did.
## ---------------------------------------------------------------------------

##' Colocalisation of two genetic traits
##'
##' Asks whether two traits share a causal variant in a region. `coloc()` is
##' the one entry point: `x1` and `x2` may each be a coloc dataset, a SuSiE
##' fit, or log Bayes factors, in any combination.
##'
##' @section Inputs:
##' Each of `x1` and `x2` is one of three things.
##'
##' * **A coloc dataset**: a list holding one trait's summary statistics for
##'   the region. It has vectors with one value per SNP (`snp`, and `beta`
##'   with `varbeta`, or `pvalues` with `MAF`) and single values that describe
##'   the study (`type`, which is `"quant"` or `"cc"`; `N`; and `sdY` or `s`).
##'   Fine-mapping with SuSiE also needs `LD`, the correlation matrix of the
##'   SNPs. [check_dataset()] lists every element and which combinations are
##'   needed, and checks a dataset for you. [coloc_test_data] holds four
##'   examples.
##' * **A SuSiE fit**: the result of [runsusie()] on a coloc dataset, or of
##'   `susieR::susie_rss()` passed through [annotate_susie()]. `coloc()` uses
##'   the log Bayes factors of the credible sets the fit found, one signal per
##'   credible set.
##' * **Log Bayes factors**: a numeric vector with one value per SNP, named by
##'   SNP id, for a single signal; or a matrix with one row per signal and one
##'   column per SNP, with SNP ids as its column names. A column named
##'   `"null"`, if there is one, holds the log Bayes factor for no
##'   association.
##'
##' The two traits are matched by SNP id, and only the SNPs both have are
##' compared.
##'
##' @section Choosing `method`:
##' `method` says how a dataset is turned into log Bayes factors, and so how
##' many signals a trait can have.
##'
##' * `"single"` assumes the trait has at most one causal variant in the
##'   region, and computes one Bayes factor per SNP from `beta` and `varbeta`
##'   (or `pvalues` and `MAF`), as [coloc.abf()] did. This is the default.
##' * `"susie"` fine-maps the dataset with SuSiE through [runsusie()], which
##'   needs `LD` and `N`, and gives one signal per credible set, so a trait
##'   can have several.
##'
##' One value applies to both traits. Two values apply to trait 1 and trait 2
##' in turn, so the traits can be treated differently: `method = c("single",
##' "susie")` fine-maps trait 2 only, for example when only trait 2 has an LD
##' matrix. `method` applies only to datasets, because a SuSiE fit is already
##' fine-mapped and log Bayes factors are already computed. Asking for
##' `"single"` for a SuSiE fit is an error, and naming a method for a trait
##' given as log Bayes factors gives a warning.
##'
##' @param x1 trait 1: a coloc dataset, a SuSiE fit or log Bayes factors (see
##'   Inputs)
##' @param x2 trait 2, with the same choices
##' @param method `"single"` or `"susie"`: one value for both traits, or one
##'   for each trait (see Choosing `method`)
##' @param prior the prior probabilities, made with [coloc_prior()]
##' @param control when to leave a pair of signals out, made with
##'   [coloc_control()]
##' @param susie.args other arguments for [runsusie()], for datasets
##'   fine-mapped here
##' @return A `coloc_result`. Its `$summary` has one row per pair of signals
##'   compared, with the posterior probabilities of the five hypotheses,
##'   `PP.H0.abf` to `PP.H4.abf`; `$results` gives, for each SNP, the
##'   posterior probability that it is the shared causal variant; and
##'   `$status` says what happened. [coloc_result] describes every field.
##' @seealso [coloc_prior()], [coloc_control()], [pp()], [is_ok()],
##'   [n_signals()], [log_bf()], [coloc_result]
##' @export
##' @examples
##' data(coloc_test_data)
##' D1 <- coloc_test_data$D1
##' D2 <- coloc_test_data$D2
##'
##' ## two datasets, with at most one causal variant assumed in each trait
##' res <- coloc(D1, D2)
##' res
##'
##' ## reading the result
##' is_ok(res)              # the posteriors were computed
##' n_signals(res)          # pairs of signals compared
##' pp(res)                 # number of snps, and PP.H0 to PP.H4 for the first pair
##' bf1 <- log_bf(res, 1)   # trait 1's log Bayes factors, to fine-map
##' finemap.bf(bf1)$summary
##'
##' if (requireNamespace("susieR", quietly = TRUE)) {
##'     D3 <- coloc_test_data$D3
##'     D4 <- coloc_test_data$D4
##'
##'     ## fine-map both traits with SuSiE: one row per pair of credible sets
##'     res2 <- coloc(D3, D4, method = "susie")
##'     res2$summary
##'
##'     ## fine-map trait 2 only: D1 stays one signal, D3 gives two
##'     res3 <- coloc(D1, D3, method = c("single", "susie"))
##'     res3$summary
##'
##'     ## two matrices of log Bayes factors, one row per signal: here, the
##'     ## credible sets of two SuSiE fits
##'     fit3 <- runsusie(D3)
##'     fit4 <- runsusie(D4)
##'     bf3 <- fit3$lbf_variable[fit3$sets$cs_index, , drop = FALSE]
##'     bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
##'     res4 <- coloc(bf3, bf4)
##'     res4$summary
##' }
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

    ## a coloc_prior has been checked already, but a plain list has not
    .check_scalar_prior(prior)

    kinds  <- c(.classify_input(x1), .classify_input(x2))
    method <- .resolve_method(method, kinds, supplied = supplied)

    ## each adapter checks its own input, and labels its messages with the
    ## trait's number
    adapt <- function(x, kind, m, trait) switch(kind,
        susie_fit = .adapt_susie(x, trait, susie.args = susie.args),
        bf        = .adapt_bf(x, trait),
        dataset   = if (identical(m, "susie")) .adapt_susie(x, trait, susie.args = susie.args)
                    else .adapt_single(x, trait))
    a1 <- adapt(x1, kinds[1], method[1], 1L)
    a2 <- adapt(x2, kinds[2], method[2], 2L)

    ## The region-size guard, coloc.abf()'s own, applied once the snp counts
    ## are known. The prior for any one snp times the number of snps must stay
    ## below 1, or nothing is left for "no causal variant". adjust_prior()
    ## warns and caps it at 1/(nsnps + 1), so coloc() gives coloc.abf()'s
    ## numbers on a large region. p12 counts the snps the traits share, as
    ## coloc.abf() does.
    snps1 <- setdiff(colnames(a1$bf), "null")
    snps2 <- setdiff(colnames(a2$bf), "null")
    prior$p1  <- adjust_prior(prior$p1,  length(snps1), "1")
    prior$p2  <- adjust_prior(prior$p2,  length(snps2), "2")
    prior$p12 <- adjust_prior(prior$p12, length(intersect(snps1, snps2)), "12")

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
    ## Each trait's log BFs as the adapter handed them to the engine, before
    ## the snps were intersected or any pair trimmed, so they can be taken on
    ## to fine-mapping without being computed again. Set here rather than in
    ## the engine, which has five exits and would need it at each.
    res$bf <- list(bf1 = a1$bf, bf2 = a2$bf)
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
##' @examples
##' data(coloc_test_data)
##' res <- coloc(coloc_test_data$D1, coloc_test_data$D2)
##' pp(res)
pp <- function(x, row = 1) {
    stopifnot(inherits(x, "coloc_result"))
    ## No rows means no signal pair was summarised, so the count is 0 and
    ## agrees with n_signals(). The posteriors stay NA under the governing
    ## rule: they could not be computed.
    if (!nrow(x$summary))
        return(stats::setNames(c(0, rep(NA_real_, 5L)),
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
##' @examples
##' data(coloc_test_data)
##' res <- coloc(coloc_test_data$D1, coloc_test_data$D2)
##' if (is_ok(res) && pp(res)["PP.H4.abf"] > 0.8)
##'     message("the traits share a causal variant")
is_ok <- function(x) {
    stopifnot(inherits(x, "coloc_result"))
    identical(x$status, "ok")
}

##' How many signal pairs a coloc_result carries
##'
##' @param x a coloc_result
##' @return integer(1)
##' @export
##' @examples
##' data(coloc_test_data)
##' res <- coloc(coloc_test_data$D1, coloc_test_data$D2)
##' n_signals(res)
n_signals <- function(x) {
    stopifnot(inherits(x, "coloc_result"))
    nrow(x$summary)
}

##' Log Bayes factors a coloc_result was computed from
##'
##' Each trait's log Bayes factors, as they went into the colocalisation:
##' signals in rows, snps in columns, covering the trait's own snps before
##' the two traits were intersected or any signal pair trimmed. Rows are named
##' as in `x$meta$idx1` and `x$meta$idx2`, unless log Bayes factors given to
##' [coloc()] had row names of their own. Use them to fine-map a trait from
##' the same numbers [coloc()] used.
##'
##' @param x a coloc_result
##' @param trait 1 or 2
##' @return numeric matrix, signals x snps
##' @export
##' @examples
##' data(coloc_test_data)
##' res <- coloc(coloc_test_data$D1, coloc_test_data$D2)
##' bf1 <- log_bf(res, 1)
##' dim(bf1)                         # one signal, 500 SNPs
##' ## fine-map trait 1 from the same numbers coloc() used
##' finemap.bf(bf1)$summary
log_bf <- function(x, trait = 1) {
    stopifnot(inherits(x, "coloc_result"))
    if (length(trait) != 1 || !trait %in% 1:2)
        stop("trait must be 1 or 2")
    if (is.null(x$bf)) {
        if (isTRUE(x$meta$schema < 2))
            stop("this coloc_result was made before coloc() kept the log Bayes ",
                 "factors; rerun coloc() to get them")
        stop("this coloc_result carries no log Bayes factors")
    }
    x$bf[[paste0("bf", trait)]]
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
## them here, starting from the canonical defaults so the numbers stay in one
## place. The RESULT is an unclassed list, and values the caller supplied are
## not put through validate_coloc_prior(): the old doors keep the
## permissiveness they have always had, and coloc() is the validating door.
.prior_from_dots <- function(dots) {
    out <- .default_prior()
    for (nm in c("p1", "p2", "p12"))
        if (!is.null(dots[[nm]])) out[[nm]] <- dots[[nm]]
    ## x[["nm"]] <- NULL DELETES the element; x["nm"] <- list(NULL) sets it.
    ## The engine reads either as NULL, but keeping the shape identical to
    ## .default_prior() is what stops the two drifting apart.
    out["weights1"] <- list(dots$prior_weights1)
    out["weights2"] <- list(dots$prior_weights2)
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
    out <- .default_control()
    for (nm in c("overlap.min", "trim_by_posterior"))
        if (!is.null(dots[[nm]])) out[[nm]] <- dots[[nm]]
    out
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

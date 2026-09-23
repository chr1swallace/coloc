## ---------------------------------------------------------------------------
## Adapters: any input becomes a log Bayes factor matrix.
##
## The contract is one function per method:
##
##   adapter(x, ...) -> list(bf   = <signals x snps, dimnames(signal_id, snp)>,
##                           meta = list(method, idx, has_null))
##
## That is the extensibility story. Adding a fine-mapping method means writing
## one function returning list(bf, meta) and adding one `method` value, and
## nothing downstream changes: the engine sees only the matrices.
##
## Dispatch is two-stage, and `method` is the SECOND stage. Stage 1 classifies
## each input from its class alone; stage 2 arbitrates only the dataset case.
## ---------------------------------------------------------------------------

## ---- stage 1: what did the caller hand us? --------------------------------

## Order is load-bearing. A susie fit is also a list, and lbf_variable is also
## a matrix, so a later test would claim an input an earlier one should have
## caught. The package writes this check as "susie" %in% class(x)
## (R/susie.R:75); inherits() is the same test.
.classify_input <- function(x) {
    if (inherits(x, "susie"))                 return("susie_fit")
    if (is.matrix(x))                         return("bf")
    if (is.numeric(x) && !is.null(names(x)))  return("bf")
    if (is.list(x) && !is.null(x$snp))        return("dataset")
    stop("input is not a susie fit, a named log-BF vector or matrix, or a coloc dataset")
}

## ---- stage 2: which method, per trait -------------------------------------

.resolve_method <- function(method, kinds, supplied = TRUE) {
    method <- unname(rep_len(method, 2L))        ## length 1 recycles to both traits
    ## USE.NAMES = FALSE is functional, not cosmetic: vapply would otherwise
    ## name the result from the character values it iterated over, making
    ## method[1] return c(single = "single"). Because names are attributes,
    ## every downstream identical(m, "susie") would then be silently FALSE and
    ## every dataset would route to the single adapter.
    method <- vapply(method, match.arg, character(1L),
                     choices = c("single", "susie"), USE.NAMES = FALSE)
    for (i in seq_len(2L)) {
        if (kinds[i] == "susie_fit" && method[i] == "single")
            stop("trait ", i, ': method = "single" on a susie fit; a fit carries no ',
                 "beta/varbeta to run the Wakefield calculation on")
        if (kinds[i] == "bf") {
            if (supplied)
                warning("trait ", i, ": method is inapplicable to precomputed log BFs",
                        call. = FALSE)
            method[i] <- NA_character_
        }
    }
    method
}

## ---- the three adapters ----------------------------------------------------

## coloc.abf()'s first step, lifted out. One row: the single causal variant
## assumed for the whole trait. No null column, because coloc.abf() appends
## none; only finemap.abf() does (R/claudia.R:329-333).
.adapt_single <- function(x, suffix = "df1", ...) {
    pd <- process.dataset(d = x, suffix = suffix)
    bf <- matrix(pd[[paste0("lABF.", suffix)]], nrow = 1L,
                 dimnames = list("sig1", pd$snp))
    list(bf = bf, meta = list(method = "single", idx = 1L, has_null = FALSE))
}

## A ROW-SELECTION POLICY, not an identity adapter with extra steps. A susie
## fit carries $sets$cs_index naming which of the L components resolved into a
## credible set; a bare matrix carries no such information. That is the whole
## reason .adapt_susie() and .adapt_bf() cannot be collapsed into one.
.adapt_susie <- function(x, susie.args = list(), ...) {
    fit <- if (inherits(x, "susie")) x
           else do.call(runsusie, c(list(d = x), susie.args))
    idx <- fit$sets$cs_index
    if (is.null(idx) || !length(idx)) {          ## -> engine status "no_signal"
        snps <- colnames(fit$lbf_variable)
        return(list(bf = matrix(numeric(0), nrow = 0L, ncol = length(snps),
                                dimnames = list(NULL, snps)),
                    meta = list(method = "susie", idx = integer(0), has_null = FALSE)))
    }
    bf <- fit$lbf_variable[idx, , drop = FALSE]
    rownames(bf) <- paste0("sig", idx)           ## lbf_variable carries no rownames
    list(bf = bf, meta = list(method = "susie", idx = idx,
                              has_null = "null" %in% colnames(bf)))
}

## Validation only. Every row becomes a signal, because nothing says otherwise.
.adapt_bf <- function(x, ...) {
    if (is.null(dim(x))) {
        if (is.null(names(x))) stop("a log-BF vector must be named by snp")
        x <- matrix(x, nrow = 1L, dimnames = list("sig1", names(x)))
    }
    if (is.null(colnames(x))) stop("a log-BF matrix needs snp ids as colnames")
    if (is.null(rownames(x))) rownames(x) <- paste0("sig", seq_len(nrow(x)))
    list(bf = x, meta = list(method = NA_character_, idx = seq_len(nrow(x)),
                             has_null = "null" %in% colnames(x)))
}

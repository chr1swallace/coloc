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

## Stage 1 answers what KIND of thing this is, and nothing else. Whether the
## contents are usable, snp ids present, values numeric, is the adapter's job:
## it is the layer that knows what a log BF has to look like, and it can say so
## in a message that names the actual problem. Asking here instead means an
## input that is obviously a vector gets rejected as "not any of these".
##
## Order is load-bearing. A susie fit is also a list, and lbf_variable is also
## a matrix, so a later test would claim an input an earlier one should have
## caught. The package writes this check as "susie" %in% class(x)
## (R/susie.R:75); inherits() is the same test.
.classify_input <- function(x) {
    if (inherits(x, "susie"))          return("susie_fit")
    if (is.matrix(x))                  return("bf")
    if (is.atomic(x))                  return("bf")   ## type checked by .adapt_bf()
    if (is.list(x) && !is.null(x$snp)) return("dataset")
    stop("input is not a susie fit, a log-BF vector or matrix, or a coloc dataset")
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
        ## a fit implies "susie", so only a "single" the caller asked for is
        ## an error; the unsupplied default must not reach it
        if (kinds[i] == "susie_fit" && method[i] == "single") {
            if (supplied)
                stop("trait ", i, ': method = "single" on a susie fit; a fit carries no ',
                     "beta/varbeta to run the Wakefield calculation on")
            method[i] <- "susie"
        }
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

## coloc.abf()'s first step, lifted out whole: check_dataset() and then
## process.dataset(), as coloc.abf() runs them (R/claudia.R:418-422).
## process.dataset() checks almost nothing itself, so without the check one NA
## in beta reaches the engine and comes back as NA posteriors with a named hit.
## Validation is the adapter's job because the adapter is what knows what its
## input has to look like; runsusie() does the same for .adapt_susie().
##
## One row: the single causal variant assumed for the whole trait. No null
## column, because coloc.abf() appends none; only finemap.abf() does
## (R/claudia.R:329-333).
.adapt_single <- function(x, trait = 1L, ...) {
    check_dataset(x, trait)
    suffix <- paste0("df", trait)
    pd <- process.dataset(d = x, suffix = suffix)
    bf <- matrix(pd[[paste0("lABF.", suffix)]], nrow = 1L,
                 dimnames = list("sig1", pd$snp))
    ## FALSE by construction, not by inspection: process.dataset() never
    ## appends a null pseudo-snp. finemap.abf() does (R/claudia.R:329-333),
    ## coloc.abf() does not, and it is coloc.abf() this reproduces.
    list(bf = bf, meta = list(method = "single", idx = 1L, has_null = FALSE))
}

## A ROW-SELECTION POLICY, not an identity adapter with extra steps. A susie
## fit carries $sets$cs_index naming which of the L components resolved into a
## credible set; a bare matrix carries no such information. That is the whole
## reason .adapt_susie() and .adapt_bf() cannot be collapsed into one.
.adapt_susie <- function(x, trait = 1L, susie.args = list(), ...) {
    if (inherits(x, "susie")) {
        fit <- x
    } else {
        ## runsusie() checks the dataset and labels its messages with suffix
        if ("suffix" %in% names(susie.args))
            stop("susie.args cannot set suffix: coloc() sets it to each trait's number")
        fit <- do.call(runsusie, c(list(d = x, suffix = trait), susie.args))
    }
    if (is.null(fit$lbf_variable))
        stop("trait ", trait, ": the susie fit has no lbf_variable")
    if (is.null(colnames(fit$lbf_variable)))
        stop("trait ", trait, ": the susie fit's lbf_variable needs snp ids as ",
             "colnames, which runsusie() adds")
    idx <- fit$sets$cs_index
    if (is.null(idx) || !length(idx)) {          ## -> engine status "no_signal"
        snps <- colnames(fit$lbf_variable)
        return(list(bf = matrix(numeric(0), nrow = 0L, ncol = length(snps),
                                dimnames = list(NULL, snps)),
                    meta = list(method = "susie", idx = integer(0),
                                ## the fit's colnames are kept above, so this
                                ## branch can carry a "null" column too
                                has_null = "null" %in% snps)))
    }
    bf <- fit$lbf_variable[idx, , drop = FALSE]
    rownames(bf) <- paste0("sig", idx)           ## lbf_variable carries no rownames
    list(bf = bf, meta = list(method = "susie", idx = idx,
                              has_null = "null" %in% colnames(bf)))
}

## Validation only. Every row becomes a signal, because nothing says otherwise.
.adapt_bf <- function(x, trait = 1L, ...) {
    ## First: is.matrix() is true of a character matrix too, and without this
    ## one it passes every remaining guard and fails three layers down inside
    ## logbf_to_pp() with "non-numeric argument to binary operator". That is
    ## 6.0.3's behaviour through coloc.bf_bf(), which is untouched; coloc() can
    ## do better.
    if (!is.numeric(x))
        stop("trait ", trait, ": log Bayes factors must be numeric, not ", typeof(x))
    if (is.null(dim(x))) {
        if (is.null(names(x)))
            stop("trait ", trait, ": a log-BF vector must be named by snp")
        x <- matrix(x, nrow = 1L, dimnames = list("sig1", names(x)))
    }
    if (is.null(colnames(x)))
        stop("trait ", trait, ": a log-BF matrix needs snp ids as colnames")
    ## check_dataset()'s rule for a dataset, applied to log BFs: no missing
    ## values, so the caller cleans their own data. One NA makes every pair
    ## using its row NA, and 6.0.3's coloc.bf_bf() reported an all-NA row as a
    ## certainty (allowlist AL12).
    if (any(!is.finite(x)))
        stop("trait ", trait, ": log Bayes factors contain missing or infinite values")
    ## the engine matches snps by name, and match() would take the first of two
    if (anyDuplicated(colnames(x)))
        stop("trait ", trait, ": duplicated snps found in the log Bayes factors")
    ## a matrix with no rows gets no names: assigning character(0) is an error
    if (is.null(rownames(x)) && nrow(x)) rownames(x) <- paste0("sig", seq_len(nrow(x)))
    list(bf = x, meta = list(method = NA_character_, idx = seq_len(nrow(x)),
                             has_null = "null" %in% colnames(x)))
}

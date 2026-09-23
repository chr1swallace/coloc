## ---------------------------------------------------------------------------
## Golden comparison against 6.0.3.
##
## D6 makes 6.0.3 the specification: anything that differs and is not on the
## allowlist is a regression, and anything on the allowlist that has STOPPED
## differing means the allowlist has gone stale. expect_matches_6_0_3() fails
## in both directions, which is the check every harness forgets.
##
## Deliberately not a tolerance map. A tolerance says "this number may wobble",
## which is the wrong statement for every allowlist entry: p12 5e-6 to 1e-5 is
## a different model, not a rounding difference, and a tolerance wide enough to
## absorb it absorbs real regressions too. One global knob, defaulting to zero.
##
## CITATIONS. A file:line reference in the harness always means 6.0.3, at the
## pinned commit 8f20f0b. That commit never moves, so the line is stable
## forever; live code in this tree is named by symbol instead, never by line.
## ---------------------------------------------------------------------------

golden_tolerance <- 0

## ---- running things and recording what they emitted -----------------------

## Returns the value plus every condition signalled, because AL03 is only ever
## observable as a MISSING warning.
capture_run <- function(expr) {
    conds <- list()
    value <- withCallingHandlers(
        tryCatch(expr, error = function(e) {
            conds[[length(conds) + 1L]] <<- c(class = "error", message = conditionMessage(e))
            structure(list(), class = "coloc_case_error")
        }),
        warning = function(w) {
            conds[[length(conds) + 1L]] <<- c(class = "warning", message = conditionMessage(w))
            invokeRestart("muffleWarning")
        },
        message = function(m) {
            conds[[length(conds) + 1L]] <<- c(class = "message",
                                              message = trimws(conditionMessage(m)))
            invokeRestart("muffleMessage")
        })
    list(value = value, conditions = conds)
}

## ---- canonicalisation ------------------------------------------------------

## .internal.selfref on data.tables and row.names on data.frames are the
## classic false-positive generators. Strip them before anything is compared.
.flatten_table <- function(x) {
    if (is.null(x)) return(NULL)
    ## coloc.abf()'s $summary is a NAMED VECTOR, not a table. as.data.frame()
    ## would move its names into rownames, which are stripped below, making a
    ## renamed PP column invisible. sensitivity() requires those names to be
    ## exact, so an atomic input is kept atomic.
    if (is.atomic(x) && is.null(dim(x))) return(x)
    out <- as.data.frame(x, stringsAsFactors = FALSE)
    rownames(out) <- NULL
    attr(out, ".internal.selfref") <- NULL
    out
}

canonicalise <- function(res) {
    if (inherits(res, "coloc_case_error")) return(res)
    if (!is.list(res)) return(res)
    ## a table is sorted into a canonical row order; an atomic summary (what
    ## coloc.abf() returns) has no rows to sort and is left exactly as it is
    .tidy <- function(x, by) {
        x <- .flatten_table(x)
        if (is.null(x) || is.null(dim(x))) return(x)
        if (nrow(x) && all(by %in% names(x)))
            x <- x[do.call(order, unname(as.list(x[, by, drop = FALSE]))), , drop = FALSE]
        rownames(x) <- NULL
        x
    }
    out <- list()
    out$summary <- .tidy(res$summary, c("idx1", "idx2"))
    out$results <- .tidy(res$results, "snp")

    out$priors  <- res$priors
    out$weights <- res$weights
    ## shape is part of the contract, so record which fields were present at all
    out$.fields <- names(res)
    out$.class  <- class(res)
    out
}

## ---- field-level differences ----------------------------------------------

.same <- function(a, b) isTRUE(all.equal(a, b, tolerance = golden_tolerance))

## Returns the names of the fields that differ. Column-level granularity on
## $summary, because that is what the allowlist needs to name.
diff_fields <- function(new, old) {
    new <- canonicalise(new)
    old <- canonicalise(old)
    d <- character(0)

    if (!.same(old$.fields, new$.fields)) d <- c(d, "shape")
    if (!.same(old$.class, new$.class))   d <- c(d, "class")

    ns <- new$summary; os <- old$summary
    if (is.null(ns) != is.null(os)) {
        d <- c(d, "summary")
    } else if (!is.null(ns)) {
        atomic <- is.atomic(os) && is.null(dim(os))
        if (atomic != (is.atomic(ns) && is.null(dim(ns)))) {
            d <- c(d, "summary")                 ## a table became a vector, or vice versa
        } else if (atomic) {
            ## names are part of the contract: prior.adjust() checks them with
            ## identical() at R/sensitivity.R:5
            if (!identical(names(os), names(ns))) d <- c(d, "summary$names")
            if (!.same(unname(os), unname(ns)))   d <- c(d, "summary$values")
        } else if (!identical(names(ns), names(os)) || nrow(ns) != nrow(os)) {
            d <- c(d, "summary")
        } else {
            for (nm in names(os)) if (!.same(os[[nm]], ns[[nm]]))
                d <- c(d, paste0("summary$", nm))
        }
    }
    if (!.same(old$results, new$results)) d <- c(d, "results")
    if (!.same(old$priors,  new$priors))  d <- c(d, "priors")
    if (!.same(old$weights, new$weights)) d <- c(d, "weights")
    unique(d)
}

diff_conditions <- function(new_conds, old_conds) {
    norm <- function(cs) sort(vapply(cs, function(c) paste(c[["class"]], c[["message"]]),
                                     character(1)))
    if (.same(norm(old_conds), norm(new_conds))) character(0) else "conditions"
}

## ---- locating the goldens --------------------------------------------------

golden_dir <- function() {
    p <- system.file("testdata", "golden", package = "coloc")
    if (nzchar(p) && dir.exists(p)) return(p)
    ## running from source without an install
    for (cand in c("inst/testdata/golden", "../../inst/testdata/golden"))
        if (dir.exists(cand)) return(normalizePath(cand))
    stop("golden directory not found; run data-raw/make_golden.R")
}

read_golden <- function(id) {
    f <- file.path(golden_dir(), paste0(id, ".rds"))
    if (!file.exists(f)) stop("no golden for case '", id, "'; run data-raw/make_golden.R")
    readRDS(f)
}

read_manifest <- function() readRDS(file.path(golden_dir(), "MANIFEST.rds"))

## ---- the assertion ---------------------------------------------------------

## layer is "B" (engine output) or "C" (end to end through an old entry point).
expect_matches_6_0_3 <- function(new, case_id, layer = "B") {
    g <- read_golden(case_id)
    old <- g[[paste0("layer", layer)]]
    testthat::expect_false(is.null(old),
                           label = paste0("golden layer ", layer, " for ", case_id))

    differing <- c(diff_fields(new$value, old$value),
                   diff_conditions(new$conditions, old$conditions))
    allowed <- allowlist_fields(case_id, layer)

    unexpected <- setdiff(differing, allowed)
    stale      <- setdiff(allowed, differing)

    testthat::expect(
        length(unexpected) == 0,
        sprintf("%s [layer %s]: unexpected difference from 6.0.3 in %s",
                case_id, layer, paste(unexpected, collapse = ", ")))
    testthat::expect(
        length(stale) == 0,
        sprintf("%s [layer %s]: allowlist is stale, no longer differs in %s",
                case_id, layer, paste(stale, collapse = ", ")))
    invisible(differing)
}

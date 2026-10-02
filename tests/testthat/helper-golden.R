## ---------------------------------------------------------------------------
## Golden comparison against 6.0.3.
##
## D6 makes 6.0.3 the specification: anything that differs and is not on the
## allowlist is a regression, and anything on the allowlist that has STOPPED
## differing means the allowlist has gone stale. test-regression.R fails in
## both directions, which is the check every harness forgets.
##
## Deliberately not a tolerance map. A tolerance says "this number may wobble",
## which is the wrong statement for every allowlist entry: p12 5e-6 to 1e-5 is
## a different model, not a rounding difference, and a tolerance wide enough to
## absorb it absorbs real regressions too. The smallest allowlisted difference
## is 9.9e-2, millions of times above the tolerance below.
##
## Comparisons use all.equal()'s default tolerance, about 1.5e-8 in mean
## relative difference, as any numeric test does. That allows for machine
## arithmetic and nothing else. The goldens were recorded on one machine
## (Ubuntu 24.04, glibc 2.39, reference BLAS). A Rocky Linux 8 machine with
## glibc 2.28 and OpenBLAS differs in the last bit, and SuSiE's prior-variance
## optimiser magnifies that to about 1e-8 in the log BFs.
##
## A different susieR version can move a fit by more than any tolerance, so
## the golden checks that refit SuSiE also need the version the goldens were
## recorded with: see skip_if_susieR_differs().
##
## CITATIONS. A file:line reference in the harness always means 6.0.3, at the
## pinned commit 8f20f0b. That commit never moves, so the line is stable
## forever; live code in this tree is named by symbol instead, never by line.
## ---------------------------------------------------------------------------

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

.same <- function(a, b) isTRUE(all.equal(a, b))

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

## A run's warnings and messages as sorted "class: message" lines, so two
## runs can be compared and a difference read.
conditions_text <- function(cs)
    sort(vapply(cs, function(c) paste0(c[["class"]], ": ", c[["message"]]), character(1)))

diff_conditions <- function(new_conds, old_conds) {
    if (.same(conditions_text(old_conds), conditions_text(new_conds))) character(0)
    else "conditions"
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

## Skips the calling test unless the installed susieR is the one the goldens
## were recorded with. Only the checks that refit SuSiE need it: layer B
## replays frozen matrices and touches no susieR.
skip_if_susieR_differs <- function() {
    want <- read_manifest()$susieR
    have <- if (requireNamespace("susieR", quietly = TRUE))
                as.character(utils::packageVersion("susieR")) else "none"
    testthat::skip_if(!identical(have, want),
                      paste0("goldens recorded with susieR ", want, ", you have ", have))
}

## ---- what 6.0.3 recorded, and coloc()'s answer in the same shape ----------

## What 6.0.3 produced for one case: the two log BF matrices it fed its engine
## (bf1, bf2), and for layers B and C the object its call returned (value)
## and the warnings and messages it emitted (conditions).
get_match_6_0_3 <- function(case_id, layer = "B") {
    layer <- match.arg(layer, c("A", "B", "C"))
    g <- read_golden(case_id)
    out <- list(bf1 = g$layerA$bf1, bf2 = g$layerA$bf2)
    if (layer == "A") return(out)
    run <- g[[paste0("layer", layer)]]
    if (is.null(run)) stop("case '", case_id, "' has no layer ", layer, " golden")
    c(out, list(value = run$value, conditions = run$conditions))
}

## coloc() returns a coloc_result, and a golden holds whatever 6.0.3's call
## returned, so the first is put into the second's shape before they are
## compared. coloc.bf_bf(), coloc.susie() and coloc.susie_bf() returned one
## shape, which .coloc_result_as_6_0_3() builds; the old entry points use it
## too, and it moves here when step 7 deletes them. coloc.abf() returned
## another, built here.
as_6_0_3 <- function(res, like = c("coloc.bf_bf", "coloc.abf")) {
    like <- match.arg(like)
    if (like == "coloc.bf_bf") return(.coloc_result_as_6_0_3(res))
    s <- res$summary
    out <- list(summary = c(nsnps = s$nsnps,
                            unlist(s[, paste0("PP.H", 0:4, ".abf"), with = FALSE])),
                results = data.frame(snp = res$results$snp,
                                     SNP.PP.H4 = res$results$SNP.PP.H4.abf),
                priors  = res$priors)
    if (!is.null(res$weights))
        out$weights <- list(prior_weights1 = res$weights$w1,
                            prior_weights2 = res$weights$w2)
    structure(out, class = c("coloc_abf", "list"))
}

## coloc.abf()'s $results columns that coloc() does not produce: each
## dataset's per-snp working from process.dataset(), the snp position, and
## the sum of the two lABFs. They are dropped from the golden before it is
## compared, which leaves snp and SNP.PP.H4.
abf_only_results <- c("position", "internal.sum.lABF",
                      paste0(rep(c("pvalues", "MAF", "N", "V", "z", "r", "lABF"),
                                 each = 2), ".df", 1:2))

## testthat prints numbers to 7 significant digits, so two values that differ
## beyond that fail with "actual != expected but don't know how to show the
## difference". This lists the largest differences at 17 digits, each with the
## snp (or row) and column it belongs to, for expect_equal()'s info argument.
## NULL when there is nothing numeric to compare element by element.
largest_differences <- function(got, expected, n = 5) {
    if (is.data.frame(expected) && is.data.frame(got) &&
        identical(dim(got), dim(expected)) && identical(names(got), names(expected))) {
        at <- if ("snp" %in% names(expected)) expected$snp
              else paste("row", seq_len(nrow(expected)))
        cols <- Filter(function(cl) is.numeric(expected[[cl]]) && is.numeric(got[[cl]]),
                       names(expected))
        d <- do.call(rbind, lapply(cols, function(cl)
            data.frame(where = paste(at, cl), expected = expected[[cl]], got = got[[cl]])))
    } else if (is.numeric(expected) && is.numeric(got) &&
               identical(length(got), length(expected))) {
        where <- if (is.matrix(expected)) {
                     r <- rownames(expected)
                     if (is.null(r)) r <- paste("row", seq_len(nrow(expected)))
                     cl <- colnames(expected)
                     if (is.null(cl)) cl <- paste("column", seq_len(ncol(expected)))
                     paste(r[row(expected)], cl[col(expected)])
                 } else if (!is.null(names(expected))) names(expected)
                 else paste("row", seq_along(expected))
        d <- data.frame(where = where, expected = as.vector(expected), got = as.vector(got))
    } else return(NULL)
    if (is.null(d) || !nrow(d)) return(NULL)

    ## a number against NA ranks first; NA against NA is no difference
    na_g <- is.na(d$got); na_e <- is.na(d$expected)
    d$diff <- ifelse(na_g | na_e, ifelse(na_g & na_e, 0, Inf), abs(d$got - d$expected))
    d <- d[d$diff > 0, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    d <- d[order(-d$diff), , drop = FALSE][seq_len(min(n, nrow(d))), , drop = FALSE]
    paste0("largest differences, to 17 significant digits:\n",
           paste0("  ", d$where, ": expected ", sprintf("%.17g", d$expected),
                  ", got ", sprintf("%.17g", d$got), collapse = "\n"))
}

## A log BF matrix compared row by row, each row a vector named by snp, so a
## difference prints with the snps around it rather than as the whole matrix.
expect_bf_equal <- function(got, expected, what) {
    if (is.null(dim(expected)) || !identical(dim(got), dim(expected)))
        return(testthat::expect_equal(got, expected, label = what,
                                      expected.label = paste("6.0.3", what),
                                      info = largest_differences(got, expected)))
    testthat::expect_equal(dimnames(got), dimnames(expected),
                           label = paste("dimnames of", what),
                           expected.label = paste("6.0.3 dimnames of", what))
    for (i in seq_len(nrow(expected)))
        testthat::expect_equal(got[i, ], expected[i, ], label = paste(what, "row", i),
                               expected.label = paste("6.0.3", what, "row", i),
                               info = largest_differences(got[i, ], expected[i, ]))
}

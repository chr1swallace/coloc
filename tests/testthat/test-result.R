## The coloc_result: its contract, and the functions that read it

bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
bf_x <- bf
colnames(bf_x) <- paste0("x", 1:4)            # no snps in common with bf

test_that("validate_coloc_result() accepts a result of the right shape", {
    ok <- .new_coloc_result(prior = .default_prior(), status = "no_signal")
    expect_s3_class(ok, "coloc_result")
    expect_silent(validate_coloc_result(ok))
    expect_equal(ok$meta$schema, 2L)
    expect_equal(ok$priors, c(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6))
    ## $bf is optional; when present it is two numeric matrices named by snp
    ok$bf <- list(bf1 = bf, bf2 = bf)
    expect_silent(validate_coloc_result(ok))
})

test_that("validate_coloc_result() says what is wrong, and that it is a bug", {
    ok <- .new_coloc_result(prior = .default_prior(), status = "no_signal")
    expect_bug <- function(x, found, should) {
        msg <- tryCatch({ validate_coloc_result(x); "" }, error = conditionMessage)
        expect_match(msg, "^internal error in coloc: ")
        expect_match(msg, found, fixed = TRUE)
        expect_match(msg, should, fixed = TRUE)
        expect_match(msg, "This is a bug in coloc; please report it$")
    }
    expect_bug(unclass(ok), "was given an object of class 'list'", "not a coloc_result")
    bad <- ok; bad$priors <- NULL
    expect_bug(bad, "the result is missing [priors]", "should have at least [summary")
    ## $summary's columns: the difference is named, whichever it is
    bad <- ok; bad$summary$extra <- integer(0)
    expect_bug(bad, "has extra column(s) [extra]", "it should have [nsnps, hit1")
    bad <- ok; bad$summary$idx2 <- NULL
    expect_bug(bad, "is missing column(s) [idx2]", "it should have [nsnps, hit1")
    bad <- ok; bad$summary <- data.table::copy(ok$summary)
    data.table::setcolorder(bad$summary, c("hit1", "nsnps"))
    expect_bug(bad, "has its columns in the wrong order", "it should have [nsnps, hit1")
    bad <- ok; bad$status <- "fine"
    expect_bug(bad, "$status is 'fine'", "should be one of [ok, no_overlap")
    bad <- ok; bad$meta$schema <- NULL
    expect_bug(bad, "$meta$schema is missing", "should be a number")
    bad <- ok; bad$meta$schema <- "2"
    expect_bug(bad, "$meta$schema is of class 'character'", "should be a number")
    ## rows against status
    bad <- coloc(bf, bf); bad$status <- "no_overlap"
    expect_bug(bad, "status 'no_overlap' has 1 row(s) in $summary", "should have none")
    bad <- ok; bad$status <- "all_trimmed"
    expect_bug(bad, "has no rows in $summary", "should keep one for each signal pair")
    bad <- coloc(bf, bf); bad$summary$PP.H4.abf <- NA_real_
    expect_bug(bad, "status 'ok' has NA in $summary$PP.H4.abf, row(s) [1]", "should have none")
    ## $bf, when present
    bad <- ok; bad$bf <- list(bf1 = "x", bf2 = "y")
    expect_bug(bad, "$bf is list(bf1 = character, bf2 = character)",
               "should be list(bf1, bf2) of numeric matrices")
    bad <- ok; bad$bf <- list(bf1 = bf, bf2 = unname(bf))
    expect_bug(bad, "bf2 = matrix without snp ids", "with snp ids as column names")
    bad <- ok; bad$bf <- list(a = bf, b = bf)
    expect_bug(bad, "$bf is list(a = numeric matrix, b = numeric matrix)", "list(bf1, bf2)")
    bad <- ok; bad$bf <- bf
    expect_bug(bad, "$bf is of class 'matrix'", "should be list(bf1, bf2)")
})

test_that("pp() gives the number of snps and the five posteriors of one pair", {
    res <- coloc(bf, bf)
    p <- pp(res)
    expect_named(p, c("nsnps", paste0("PP.H", 0:4, ".abf")))
    expect_equal(unname(p), unname(unlist(res$summary[1, names(p), with = FALSE])))
})

test_that("pp() stops on a row that does not exist", {
    res <- coloc(bf, bf)
    expect_equal(pp(res, row = 1), pp(res))
    expect_error(pp(res, row = 2), "between 1 and 1")
    expect_error(pp(res, row = 0), "between 1 and 1")
})

test_that("pp(), is_ok() and n_signals() answer for a result with no rows", {
    ## bf and bf_x share no snps, so the result has no rows: reading the column
    ## directly, empty$summary$PP.H4.abf > 0.9 is logical(0), which if() cannot
    ## use. These three functions still answer: not ok, 0 pairs, NA posteriors
    empty <- suppressWarnings(coloc(bf, bf_x))
    expect_false(is_ok(empty))
    expect_equal(n_signals(empty), 0L)
    expect_equal(unname(pp(empty)["nsnps"]), 0)
    expect_true(all(is.na(pp(empty)[-1])))
})

test_that("is_ok() is TRUE only when the posteriors were computed", {
    expect_true(is_ok(coloc(bf, bf)))
    ## every pair trimmed: rows kept, PP.H1 to PP.H4 NA
    far <- matrix(c(10, 0, 0, 0, 0, 0), nrow = 1, dimnames = list("sig1", paste0("s", 5:10)))
    expect_false(is_ok(suppressWarnings(coloc(cbind(bf, far), bf))))
    expect_equal(n_signals(coloc(rbind(bf, bf), bf)), 2L)
})

test_that("log_bf() returns each trait's log BFs as the adapter handed them over", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D1 <- coloc_test_data$D1
    D2 <- coloc_test_data$D2
    D3 <- coloc_test_data$D3
    D4 <- coloc_test_data$D4
    fit3 <- suppressMessages(runsusie(D3))
    fit4 <- suppressMessages(runsusie(D4))
    bf3 <- fit3$lbf_variable[fit3$sets$cs_index, , drop = FALSE]
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    runs <- list(
        single = list(res = coloc(D1, D2),
                      a1 = .adapt_single(D1, 1), a2 = .adapt_single(D2, 2)),
        susie = list(res = suppressMessages(coloc(D3, D4, method = "susie")),
                     a1 = suppressMessages(.adapt_susie(D3, 1)),
                     a2 = suppressMessages(.adapt_susie(D4, 2))),
        fit = list(res = coloc(fit3, fit4),
                   a1 = .adapt_susie(fit3, 1), a2 = .adapt_susie(fit4, 2)),
        matrix = list(res = coloc(bf3, bf4), a1 = .adapt_bf(bf3, 1), a2 = .adapt_bf(bf4, 2)),
        vector = list(res = coloc(bf3[1, ], bf4[1, ]),
                      a1 = .adapt_bf(bf3[1, ], 1), a2 = .adapt_bf(bf4[1, ], 2)))
    for (kind in names(runs)) {
        r <- runs[[kind]]
        expect_equal(log_bf(r$res, 1), r$a1$bf, label = paste(kind, "trait 1"))
        expect_equal(log_bf(r$res, 2), r$a2$bf, label = paste(kind, "trait 2"))
        expect_equal(rownames(log_bf(r$res, 1)), paste0("sig", r$res$meta$idx1))
        expect_equal(rownames(log_bf(r$res, 2)), paste0("sig", r$res$meta$idx2))
    }
    ## each trait's own snps, not the shared ones
    sub <- coloc(D1, suppressMessages(subset_dataset(D2, 1:300)))
    expect_equal(ncol(log_bf(sub, 1)), 500L)
    expect_equal(ncol(log_bf(sub, 2)), 300L)
})

test_that("a result with no rows still keeps its log BFs", {
    none <- suppressWarnings(coloc(bf, bf_x))
    expect_equal(none$status, "no_overlap")
    expect_equal(log_bf(none, 1), bf)
    expect_equal(log_bf(none, 2), bf_x)
})

test_that("log_bf() asks for trait 1 or 2", {
    res <- coloc(bf, bf)
    expect_error(log_bf(res, 3), "trait must be 1 or 2")
})

test_that("log_bf() says to rerun coloc() for a result saved before $bf existed", {
    ## from before coloc() kept each trait's log BFs in $bf: schema 1, and no $bf
    res <- coloc(bf, bf)
    res$bf <- NULL
    res$meta$schema <- 1L
    expect_error(log_bf(res, 1), "rerun coloc()", fixed = TRUE)
})

test_that("pp(), is_ok(), n_signals() and log_bf() say to run coloc() first", {
    ## a coloc.abf() result (class coloc_abf) is the likeliest mistake. It has
    ## no $status, $meta or $bf, and its $summary is a named vector, so without
    ## the check is_ok() would say FALSE, n_signals() NULL, and pp() and
    ## log_bf() would stop with misleading errors
    invisible(utils::capture.output(old <- coloc.abf(coloc_test_data$D1, coloc_test_data$D2)))
    fns <- list(pp = pp, is_ok = is_ok, n_signals = n_signals, log_bf = log_bf)
    for (nm in names(fns)) {
        expect_error(fns[[nm]](old),
                     "x must be the result of coloc(), but is a coloc_abf; run coloc() first",
                     fixed = TRUE, label = nm)
        expect_error(fns[[nm]](list()), "but is a list", fixed = TRUE, label = nm)
    }
    ## the error names the call made, not the helper behind it
    expect_equal(conditionCall(tryCatch(pp(old), error = identity)), quote(pp(old)))
})

test_that("print() shows the status, the methods, the priors and the summary", {
    out <- utils::capture.output(print(coloc(bf, bf)))
    expect_match(out[1], "status ok with 1 signal pair")
    expect_true(any(grepl("method: bf, bf", out)))
    expect_true(any(grepl("priors: p1 = 1e-04, p2 = 1e-04, p12 = 5e-06", out)))
    expect_true(any(grepl("PP.H4.abf", out)))
    ## a result's own warnings are printed as notes; its log BFs are not shown
    out <- utils::capture.output(print(suppressWarnings(coloc(bf, bf_x))))
    expect_match(out[1], "status no_overlap with 0 signal pair")
    expect_true(any(grepl("note: no snps in common", out)))
    expect_false(any(grepl("bf1|sig1", out)))
})

## ---------------------------------------------------------------------------
## Step 5: the old entry points now run on the one engine.
##
## The regression grid already proves they return what they always returned.
## What is asserted here is the part the grid cannot see: the contracts other
## code depends on, and the two behaviours that deliberately changed.
## ---------------------------------------------------------------------------

test_that("coloc.abf() keeps the exact summary contract sensitivity() requires", {
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    invisible(utils::capture.output(res <- coloc.abf(D1, D2)))

    ## prior.adjust() checks these names with identical() at R/sensitivity.R:5,
    ## so the named vector and its order are load-bearing, not cosmetic
    expect_true(is.numeric(res$summary))
    expect_null(dim(res$summary))
    expect_identical(names(res$summary),
                     c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                       "PP.H3.abf", "PP.H4.abf"))
    expect_s3_class(res, "coloc_abf")

    ## and $results keeps every per-trait column, which is why coloc.abf()
    ## still does its own merge rather than taking the engine's two columns
    expect_true(all(c("snp", "lABF.df1", "lABF.df2", "internal.sum.lABF",
                      "SNP.PP.H4", "z.df1", "V.df1", "r.df1") %in% names(res$results)))
    expect_equal(sum(res$results$SNP.PP.H4), 1, tolerance = 1e-10)
})

test_that("coloc.abf() still prints its posteriors", {
    ## combine.abf() printed unless told to be quiet, and the engine always
    ## tells it to be quiet. Dropping the print would be a silent change to
    ## what users see at the console.
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    out <- utils::capture.output(coloc.abf(D1, D2))
    expect_true(any(grepl("PP.H4.abf", out)))
    expect_true(any(grepl("PP abf for shared variant", out)))
})

test_that("sensitivity() still accepts a coloc.abf() result", {
    ## sensitivity() has no tests of its own in 6.0.3, and it is the main
    ## consumer of the shape above. Step 6 rewrites it; this guards the
    ## handover.
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    invisible(utils::capture.output(res <- coloc.abf(D1, D2)))
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_error(suppressMessages(
        sensitivity(res, rule = "H4 > 0.5", plot.manhattans = FALSE)))
})

test_that("AL02: coloc.susie_bf() now forwards the priors it declares", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    fit4 <- suppressMessages(runsusie(D4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]

    ## 6.0.3 declared p1/p2/p12 in the signature and forwarded none of them,
    ## so these two calls used to give identical answers
    a <- suppressMessages(coloc.susie_bf(D3, bf4))
    b <- suppressMessages(coloc.susie_bf(D3, bf4, p12 = 1e-4))
    expect_false(isTRUE(all.equal(a$summary$PP.H4.abf, b$summary$PP.H4.abf,
                                  tolerance = 0)))

    ## and the default still reproduces what it always did
    expect_equal(a$priors[["p12"]], 5e-6)
})

test_that("the susie doors keep their old return shapes", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    res <- suppressMessages(coloc.susie(D3, D4))

    expect_true(data.table::is.data.table(res$summary))
    expect_true(all(c("nsnps", "hit1", "hit2", "idx1", "idx2") %in% names(res$summary)))
    ## idx must be susie component ids, not row numbers of the trimmed matrices
    expect_true(all(res$summary$idx1 %in% seq_len(nrow(
        suppressMessages(runsusie(D3))$lbf_variable))))
    ## no coloc_result leaks out of the old doors
    expect_false(inherits(res, "coloc_result"))
})

test_that("an empty credible set still returns nsnps = NA, with a new message", {
    ## The wording moved when the check moved into the engine: the engine
    ## speaks of traits and signals because it never sees a dataset. The
    ## SHAPE, which is what callers branch on, is unchanged.
    fit <- structure(list(lbf_variable = matrix(0, 3, 4,
                                                dimnames = list(NULL, paste0("s", 1:4))),
                          sets = list(cs = list(), cs_index = integer(0))),
                     class = "susie")
    expect_warning(res <- coloc.susie(fit, fit), "no signals")
    expect_true(data.table::is.data.table(res))
    expect_true(is.na(res$nsnps))
})

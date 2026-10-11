## The old entry points, coloc.abf(), coloc.bf_bf(), coloc.susie() and
## coloc.susie_bf(), which now run on the new engine until step 7 (the hard
## break, in a later stage) deletes them. test-v603.R checks they return the
## same values as 6.0.3's. These tests check what comparing values cannot:
## - that other code can still use their results: coloc.abf()'s $summary keeps
##   the names and order sensitivity() needs, and sensitivity() still accepts
##   it; the susie functions keep their old shapes;
## - that coloc.abf() still prints its posteriors;
## - that coloc() takes the same route as coloc.abf() and coloc.susie();
## - how they take their arguments: priors and control as loose arguments,
##   with 6.0.3's defaults kept;
## - the errors they give: a vector prior is refused, and a character matrix
##   still fails coloc.bf_bf() as it did in 6.0.3;
## - one deliberate change: coloc.susie_bf() now uses the priors it is given
##   (AL02 in test-v603.R).

D1 <- coloc_test_data$D1
D2 <- coloc_test_data$D2
D3 <- coloc_test_data$D3
D4 <- coloc_test_data$D4
bf <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0), nrow = 1, dimnames = list("sig1", paste0("s", 1:8)))

test_that("coloc.abf() keeps the summary sensitivity() requires", {
    invisible(utils::capture.output(res <- coloc.abf(D1, D2)))
    ## prior.adjust() checks these names with identical(), so their order matters
    expect_identical(names(res$summary),
                     c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"))
    expect_null(dim(res$summary))
    expect_s3_class(res, "coloc_abf")
    ## and each dataset's per-snp working, which coloc() does not keep
    expect_true(all(c("snp", "lABF.df1", "lABF.df2", "internal.sum.lABF", "SNP.PP.H4",
                      "z.df1", "V.df1", "r.df1") %in% names(res$results)))
})

test_that("coloc.abf() still prints its posteriors", {
    out <- utils::capture.output(coloc.abf(D1, D2))
    expect_true(any(grepl("PP.H4.abf", out)))
    expect_true(any(grepl("PP abf for shared variant", out)))
})

test_that("sensitivity() still accepts a coloc.abf() result", {
    ## sensitivity() has no tests of its own; step 6 rewrites it
    invisible(utils::capture.output(res <- coloc.abf(D1, D2)))
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    expect_no_error(suppressMessages(sensitivity(res, rule = "H4 > 0.5",
                                                 plot.manhattans = FALSE)))
})

test_that("coloc() takes the same route as coloc.abf() and coloc.susie()", {
    ## both sides run on the same engine, so this checks routing: the p12
    ## default, the adapters and the renumbering of signals
    invisible(utils::capture.output(old <- coloc.abf(D1, D2)))
    new <- coloc(D1, D2, prior = coloc_prior(p12 = 1e-5))
    expect_equal(as.list(pp(new)), as.list(old$summary))
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    old <- suppressMessages(coloc.susie(D3, D4))
    new <- suppressMessages(coloc(D3, D4, method = "susie"))
    expect_equal(new$summary, old$summary)
})

test_that("the susie doors keep their old return shapes", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    res <- suppressMessages(coloc.susie(D3, D4))
    expect_true(data.table::is.data.table(res$summary))
    expect_named(res, c("summary", "results", "priors"))
    expect_false(inherits(res, "coloc_result"))
    ## no credible sets: a bare data.table, as in 6.0.3, with the engine's wording
    fit <- structure(list(lbf_variable = matrix(0, 3, 4, dimnames = list(NULL, paste0("s", 1:4))),
                          sets = list(cs = list(), cs_index = integer(0))),
                     class = "susie")
    expect_warning(res <- coloc.susie(fit, fit), "no signals")
    expect_true(is.na(res$nsnps))
})

test_that("AL02: coloc.susie_bf() now uses the priors it declares", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    fit4 <- suppressMessages(runsusie(D4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    a <- suppressMessages(coloc.susie_bf(D3, bf4))
    b <- suppressMessages(coloc.susie_bf(D3, bf4, p12 = 1e-4))
    expect_equal(a$priors[["p12"]], 5e-6)
    expect_equal(b$priors[["p12"]], 1e-4)
    expect_false(isTRUE(all.equal(a$summary$PP.H4.abf, b$summary$PP.H4.abf)))
})

test_that("the old doors take priors and control as loose arguments", {
    expect_equal(.prior_from_dots(list()), .default_prior())
    expect_equal(.control_from_dots(list()), .default_control())
    expect_equal(.prior_from_dots(list(p12 = 1e-7))$p12, 1e-7)
    expect_equal(.control_from_dots(list(overlap.min = 0.9))$overlap.min, 0.9)
    ## and keep 6.0.3's own defaults in their signatures (AL01)
    expect_equal(formals(coloc.abf)$p12, 1e-5)
    expect_equal(formals(coloc.bf_bf)$p12, 5e-6)
})

test_that("the old doors reject a vector prior, as coloc_prior() does", {
    ## 6.0.3 accepted a vector p1 or p2 (AL10 in test-v603.R).
    ## coloc.susie() and coloc.susie_bf() go through coloc(), which checks
    vec <- rep(1e-4, 8)
    msg <- tryCatch(coloc_prior(p1 = vec), error = conditionMessage)
    expect_error(coloc.bf_bf(bf, bf, p1 = vec), msg, fixed = TRUE)
    expect_error(coloc.bf_bf(bf, bf, p2 = vec), msg, fixed = TRUE)
    fit <- structure(list(lbf_variable = bf, sets = list(cs = list(3L), cs_index = 1L)),
                     class = "susie")
    expect_error(coloc.susie(fit, fit, p1 = vec), msg, fixed = TRUE)
    expect_error(coloc.susie_bf(fit, bf, p2 = vec), msg, fixed = TRUE)
    ## with weights too, 6.0.3's own message comes first
    expect_error(coloc.bf_bf(bf, bf, p1 = vec, prior_weights1 = rep(1, 8)),
                 "only be used with single values")
    ## coloc.abf() stops earlier, in adjust_prior(), with R's own message
    abf <- tryCatch(coloc.abf(D1, D2, p1 = rep(1e-4, 500)), error = conditionMessage)
    expect_false(grepl("single number", abf))
})

test_that("coloc.bf_bf() still fails as 6.0.3 did on a character matrix", {
    ## it never calls .adapt_bf(), whose type check is coloc()'s
    chr <- matrix("a", 2, 3, dimnames = list(NULL, paste0("s", 1:3)))
    num <- matrix(0, 2, 3, dimnames = list(NULL, paste0("s", 1:3)))
    expect_error(coloc.bf_bf(chr, num), "non-numeric argument")
})

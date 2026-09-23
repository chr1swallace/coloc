## ---------------------------------------------------------------------------
## Step 4: the coloc() front door, additive.
##
## Nothing is removed, so the gate here is that coloc() reproduces the old
## entry points when asked for the same model, and that the per-trait method
## does what decision D2 says it does.
## ---------------------------------------------------------------------------

test_that("coloc() reproduces coloc.abf() on the single-variant route", {
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    ## p12 = 1e-5 is coloc.abf()'s default; coloc_prior() defaults to 5e-6.
    ## Two defaults for one parameter is exactly the drift the prior object
    ## exists to stop, so it has to be asked for explicitly here.
    new <- coloc(D1, D2, method = "single", prior = coloc_prior(p12 = 1e-5))
    old <- suppressWarnings(utils::capture.output(res <- coloc.abf(D1, D2)))
    old <- res

    expect_s3_class(new, "coloc_result")
    expect_equal(new$summary$PP.H4.abf, unname(old$summary["PP.H4.abf"]), tolerance = 0)
    for (h in 0:4)
        expect_equal(new$summary[[paste0("PP.H", h, ".abf")]],
                     unname(old$summary[[paste0("PP.H", h, ".abf")]]), tolerance = 0)
    expect_equal(new$summary$nsnps, unname(old$summary["nsnps"]))
})

test_that("coloc() reproduces coloc.susie(), including the index renumbering", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    new <- suppressMessages(coloc(D3, D4, method = "susie"))
    old <- suppressMessages(coloc.susie(D3, D4))

    expect_equal(new$summary$PP.H4.abf, old$summary$PP.H4.abf, tolerance = 0)
    ## idx must be the susie component ids, not row numbers of the trimmed
    ## matrices, which is what coloc.susie() does at R/susie.R:96-98
    expect_equal(new$summary$idx1, old$summary$idx1)
    expect_equal(new$summary$idx2, old$summary$idx2)
    expect_equal(new$summary$hit1, old$summary$hit1)
})

test_that("AL01: coloc()'s default p12 differs from coloc.abf()'s, on purpose", {
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    expect_equal(coloc_prior()$p12, 5e-6)
    expect_equal(formals(coloc.abf)$p12, 1e-5)

    default <- coloc(D1, D2, method = "single")
    matched <- coloc(D1, D2, method = "single", prior = coloc_prior(p12 = 1e-5))
    expect_false(isTRUE(all.equal(default$summary$PP.H4.abf,
                                  matched$summary$PP.H4.abf, tolerance = 0)))
})

test_that("method is per trait: one LD-less trait does not drag the other down", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D1_noLD <- coloc_test_data$D1; D1_noLD$LD <- NULL     ## no LD panel
    D3 <- coloc_test_data$D3                               ## good LD, 2 credible sets

    ## a global "susie" cannot run at all: runsusie() requires LD
    expect_error(suppressMessages(coloc(D1_noLD, D3, method = "susie")),
                 "missing required element")

    ## a global "single" runs, but collapses D3's credible sets into one
    forced <- suppressMessages(coloc(D1_noLD, D3, method = "single"))
    mixed  <- suppressMessages(coloc(D1_noLD, D3, method = c("single", "susie")))

    expect_equal(n_signals(forced), 1L)
    expect_equal(n_signals(mixed), 2L)
    expect_equal(mixed$meta$method, c("single", "susie"))
})

test_that("an unsupplied method is not read as a per-trait request", {
    ## method = c("single", "susie") is a match.arg default VECTOR. If the
    ## front door recycled it blindly, trait 2 would be fine-mapped whenever
    ## the caller said nothing at all.
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    silent <- coloc(D1, D2)
    explicit <- coloc(D1, D2, method = "single")
    expect_equal(silent$meta$method, c("single", "single"))
    expect_equal(silent$summary$PP.H4.abf, explicit$summary$PP.H4.abf, tolerance = 0)
})

test_that("every input combination reaches the right adapter", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    fit4 <- suppressMessages(runsusie(D4))
    bf4  <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]

    ## dataset against an already-fitted susie object
    a <- suppressMessages(coloc(D3, fit4, method = "susie"))
    ## dataset against a raw BF matrix: method is inapplicable to trait 2
    expect_warning(b <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   "trait 2: method is inapplicable")
    expect_equal(a$summary$PP.H4.abf, b$summary$PP.H4.abf, tolerance = 0)
    expect_equal(b$meta$method, c("susie", NA_character_))

    ## a susie fit cannot be re-fitted under the single-variant model
    expect_error(coloc(fit4, D3, method = "single"), "susie fit")
    ## supplying method for a BF input warns, once per trait
    w <- capture_warnings(coloc(bf4, bf4, method = "single"))
    expect_length(w, 2L)
    expect_true(all(grepl("method is inapplicable", w)))
})

test_that("the accessors answer without reaching into internals", {
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    res <- coloc(D1, D2, method = "single")
    expect_true(is_ok(res))
    expect_equal(n_signals(res), 1L)
    p <- pp(res)
    expect_named(p, c("nsnps", paste0("PP.H", 0:4, ".abf")))
    expect_equal(unname(p["PP.H4.abf"]), res$summary$PP.H4.abf)

    ## an empty result answers rather than returning logical(0)
    bf <- matrix(0, 1, 3, dimnames = list("sig1", paste0("s", 1:3)))
    dis <- bf; colnames(dis) <- paste0("x", 1:3)
    empty <- suppressWarnings(coloc(bf, dis))
    expect_false(is_ok(empty))
    expect_equal(n_signals(empty), 0L)
    expect_true(all(is.na(pp(empty))))
})

test_that("coloc_prior() and coloc_control() validate once, in one place", {
    expect_error(coloc_prior(p12 = 2e-4), "p12 must lie")
    expect_error(validate_coloc_prior(coloc_prior(), nsnps = 20000), "exceeds 1")
    expect_error(coloc_prior(weights1 = c(1, -1)), "non-NA and positive")
    expect_error(coloc_prior(p1 = c(1e-4, 1e-4), weights1 = c(1, 1)), "scalar")
    expect_error(coloc_control(overlap.min = 2))
    expect_s3_class(coloc_prior(), "coloc_prior")
    expect_s3_class(coloc_control(), "coloc_control")
})

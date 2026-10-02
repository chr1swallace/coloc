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
    ## nsnps is a count, so it agrees with n_signals() rather than saying
    ## "unknown"; the posteriors stay NA because they could not be computed
    expect_equal(unname(pp(empty)["nsnps"]), 0)
    expect_true(all(is.na(pp(empty)[paste0("PP.H", 0:4, ".abf")])))
})

test_that("the prior and control defaults are declared in exactly one place", {
    ## 6.0.3 has p12 = 1e-5 in coloc.abf() and 5e-6 in coloc.bf_bf() (AL01)
    ## precisely because the same defaults were re-declared per entry point.
    ## Every new path must derive from the constructors, so that changing a
    ## default in R/prior.R changes it everywhere at once.
    expect_equal(.default_prior(), unclass(coloc_prior()))
    expect_equal(.default_control(), unclass(coloc_control()))
    expect_equal(.prior_from_dots(list()), .default_prior())
    expect_equal(.control_from_dots(list()), .default_control())

    ## and a supplied value still wins over the default
    expect_equal(.prior_from_dots(list(p12 = 1e-7))$p12, 1e-7)
    expect_equal(.control_from_dots(list(overlap.min = 0.9))$overlap.min, 0.9)

    ## the OLD doors keep their own literals on purpose: that is 6.0.3's
    ## published API, frozen until the hard break
    expect_equal(formals(coloc.abf)$p12, 1e-5)
    expect_equal(formals(coloc.bf_bf)$p12, 5e-6)
})

test_that("p1, p2 and p12 must be scalar; per-snp variation goes through weights", {
    ## p1 is the prior that any ONE snp is causal for trait 1 only, so it is a
    ## single number. Varying it across snps is what weights1/weights2 express,
    ## and combine_abf_weighted() builds Q * p * w/sum(w) from them. A vector p1
    ## is a second spelling of the same thing.
    for (arg in c("p1", "p2", "p12")) {
        a <- list(rep(1e-4, 8)); names(a) <- arg
        expect_error(do.call(coloc_prior, a), "must each be a single number",
                     label = paste0("vector ", arg))
        expect_error(do.call(coloc_prior, a), "weights1/weights2",
                     label = paste0(arg, " error names the alternative"))
    }
    ## the guard fires before the bounds check, so the message is the useful one
    expect_error(coloc_prior(p1 = c(1e-4, 1e-4), weights1 = c(1, 1)),
                 "must each be a single number")
    ## and the supported way still works
    expect_s3_class(coloc_prior(p1 = 1e-4, weights1 = c(1, 2, 3)), "coloc_prior")
})

test_that("no door hands the engine a vector prior", {
    ## The engine takes p1, p2 and p12 as single numbers. 6.0.3's
    ## coloc.bf_bf() accepted a vector and recycled it into 33 hypotheses
    ## (AL10), and until now the engine re-expressed one as prior weights.
    ## Every door that can reach the engine without coloc_prior() refuses with
    ## coloc_prior()'s message; this pins which door checks where.
    bf  <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0), nrow = 1,
                  dimnames = list("sig1", paste0("s", 1:8)))
    vec <- rep(1e-4, 8)
    msg <- tryCatch(coloc_prior(p1 = vec), error = conditionMessage)
    expect_match(msg, "must each be a single number")

    ## coloc(): a coloc_prior refuses when it is made, and a plain list, which
    ## never went through coloc_prior(), is refused by coloc() itself
    expect_error(coloc(bf, bf, prior = list(p1 = vec, p2 = 1e-4, p12 = 5e-6)),
                 msg, fixed = TRUE)

    ## coloc.bf_bf() calls the engine directly, so it checks for itself. With
    ## weights as well, 6.0.3's own message still comes first.
    expect_error(coloc.bf_bf(bf, bf, p1 = vec), msg, fixed = TRUE)
    expect_error(coloc.bf_bf(bf, bf, p2 = vec), msg, fixed = TRUE)
    expect_error(coloc.bf_bf(bf, bf, p1 = vec, prior_weights1 = rep(1, 8)),
                 "only be used with single values")

    ## coloc.susie() and coloc.susie_bf() go through coloc(), which refuses
    fit <- structure(list(lbf_variable = bf, sets = list(cs = list(3L), cs_index = 1L)),
                     class = "susie")
    expect_error(coloc.susie(fit, fit, p1 = vec), msg, fixed = TRUE)
    expect_error(coloc.susie_bf(fit, bf, p2 = vec), msg, fixed = TRUE)

    ## coloc.abf() builds a plain list too, but a vector p1 dies earlier in
    ## adjust_prior()'s if (nsnps * p >= 1), which R 4.2 made an error on a
    ## length > 1 condition. That is an accident, not a guard, so assert only
    ## that it fails and that it is NOT the deliberate message; the base R
    ## wording is not ours to pin.
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    abf <- tryCatch(coloc.abf(D1, D2, p1 = rep(1e-4, length(D1$snp))),
                    error = conditionMessage)
    expect_type(abf, "character")
    expect_false(grepl("must each be a single number", abf))
})

test_that("$priors is three named numbers", {
    bf <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0), nrow = 1,
                 dimnames = list("sig1", paste0("s", 1:8)))
    res <- .coloc_engine(bf, bf)
    expect_true(is.numeric(res$priors))
    expect_named(res$priors, c("p1", "p2", "p12"))
})

test_that("coloc_prior() and coloc_control() validate once, in one place", {
    expect_error(coloc_prior(p12 = 2e-4), "p12 must lie")
    expect_error(coloc_prior(weights1 = c(1, -1)), "non-NA and positive")
    expect_error(coloc_control(overlap.min = 2))
    expect_s3_class(coloc_prior(), "coloc_prior")
    expect_s3_class(coloc_control(), "coloc_control")
})

test_that("log_bf() returns each trait's log BFs as the adapter handed them over", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    fit3 <- suppressMessages(runsusie(D3)); fit4 <- suppressMessages(runsusie(D4))
    bf3 <- fit3$lbf_variable[fit3$sets$cs_index, , drop = FALSE]
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]

    ## each input kind, through coloc() and through its adapter on its own
    runs <- list(
        single = list(res = coloc(D1, D2, method = "single"),
                      a1 = .adapt_single(D1, 1), a2 = .adapt_single(D2, 2)),
        susie  = list(res = suppressMessages(coloc(D3, D4, method = "susie")),
                      a1 = suppressMessages(.adapt_susie(D3, 1)),
                      a2 = suppressMessages(.adapt_susie(D4, 2))),
        fit    = list(res = coloc(fit3, fit4),
                      a1 = .adapt_susie(fit3, 1), a2 = .adapt_susie(fit4, 2)),
        matrix = list(res = coloc(bf3, bf4), a1 = .adapt_bf(bf3, 1), a2 = .adapt_bf(bf4, 2)),
        vector = list(res = coloc(bf3[1, ], bf4[1, ]),
                      a1 = .adapt_bf(bf3[1, ], 1), a2 = .adapt_bf(bf4[1, ], 2)))
    for (kind in names(runs)) {
        r <- runs[[kind]]
        expect_equal(log_bf(r$res, 1), r$a1$bf, label = paste(kind, "trait 1"))
        expect_equal(log_bf(r$res, 2), r$a2$bf, label = paste(kind, "trait 2"))
        ## rows line up with the indices the result reports
        expect_equal(rownames(log_bf(r$res, 1)), paste0("sig", r$res$meta$idx1))
        expect_equal(rownames(log_bf(r$res, 2)), paste0("sig", r$res$meta$idx2))
    }
    ## the trait's own snps, not the intersection: D1 against a subset of D2
    sub <- coloc(D1, subset_dataset(D2, 1:300), method = "single")
    expect_equal(ncol(log_bf(sub, 1)), length(D1$snp))
    expect_equal(ncol(log_bf(sub, 2)), 300L)
})

test_that("a result with no rows still keeps its log BFs", {
    a <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    b <- a; colnames(b) <- paste0("x", 1:4)
    none <- suppressWarnings(coloc(a, b))
    expect_equal(none$status, "no_overlap")
    expect_equal(log_bf(none, 1), a)
    expect_equal(colnames(log_bf(none, 2)), colnames(b))

    empty_fit <- structure(list(lbf_variable = matrix(0, 3, 4,
                                                      dimnames = list(NULL, paste0("s", 1:4))),
                                sets = list(cs = list(), cs_index = integer(0))),
                           class = "susie")
    nosig <- suppressWarnings(coloc(empty_fit, a))
    expect_equal(nosig$status, "no_signal")
    expect_equal(nrow(log_bf(nosig, 1)), 0L)
    expect_equal(log_bf(nosig, 2), a)
})

test_that("log_bf() on a result saved before it existed says to rerun coloc()", {
    a <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    old <- coloc(a, a)
    old$bf <- NULL; old$meta$schema <- 1L
    expect_error(log_bf(old, 1), "rerun coloc()", fixed = TRUE)
    expect_error(log_bf(coloc(a, a), 3), "trait must be 1 or 2")
    ## and print() leaves it out
    out <- utils::capture.output(print(coloc(a, a)))
    expect_false(any(grepl("bf1|sig1", out)))
})

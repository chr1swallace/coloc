## coloc(), the entry point, tested through its own steps:
## - which method each trait gets: "single" when none is given, and one
##   method per trait when two are given;
## - that each kind of input (a dataset, a SuSiE fit, log BFs) reaches the
##   adapter for it;
## - the checks made before the engine runs: a prior given as a plain list
##   must hold scalar p1, p2 and p12, and any weights1 and weights2 must be
##   positive and not NA (vectors, one per snp, as in coloc_prior()); a bad
##   dataset or bad log BFs stop coloc() with the trait named;
## - the region-size guard, which caps p1, p2 and p12 when a prior times its
##   number of snps reaches 1, as coloc.abf() did.
## Some tests compare with 6.0.3: get_match_6_0_3(), in helper-v603.R, runs a
## call as 6.0.3 made it.

D1 <- coloc_test_data$D1
D2 <- coloc_test_data$D2
D3 <- coloc_test_data$D3
D4 <- coloc_test_data$D4

## the messages of every warning an expression gives, in order
warnings_of <- function(expr) {
    w <- character(0)
    withCallingHandlers(expr, warning = function(cond) {
        w <<- c(w, conditionMessage(cond))
        invokeRestart("muffleWarning")
    })
    w
}

test_that("an unsupplied method means single for datasets", {
    ## method's default is a vector, c("single", "susie"), not a request
    res <- coloc(D1, D2)
    expect_equal(res$meta$method, c("single", "single"))
    expect_equal(res$summary, coloc(D1, D2, method = "single")$summary)
})

test_that("method is per trait: a trait without LD can stay single", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D1_noLD <- D1
    D1_noLD$LD <- NULL
    ## SuSiE needs LD, so "susie" for both cannot run
    expect_error(suppressMessages(coloc(D1_noLD, D3, method = "susie")),
                 "missing required element")
    mixed <- suppressMessages(coloc(D1_noLD, D3, method = c("single", "susie")))
    expect_equal(mixed$meta$method, c("single", "susie"))
    expect_equal(n_signals(mixed), 2L)
    ## "single" for both runs, but takes D3's two signals as one
    expect_equal(n_signals(coloc(D1_noLD, D3)), 1L)
})

test_that("every kind of input reaches its adapter", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    fit4 <- suppressMessages(runsusie(D4))
    bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]
    from_fit <- suppressMessages(coloc(D3, fit4, method = "susie"))
    ## AL14: the method meant for D3 reaches the log BFs too, and warns
    expect_warning(from_bf <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   "trait 2: method is inapplicable")
    expect_equal(from_bf$summary$PP.H4.abf, from_fit$summary$PP.H4.abf)
    expect_equal(from_bf$meta$method, c("susie", NA))
    ## a fit needs no method: with none named, it is taken as "susie". One
    ## method applies to both traits, so "single" here reaches the fit too,
    ## and stops coloc(): a fit has no beta or varbeta for "single" to use
    expect_equal(coloc(fit4, fit4)$meta$method, c("susie", "susie"))
    expect_error(coloc(fit4, D3, method = "single"), "susie fit")
    ## a method named for log BFs warns once for each trait
    expect_length(warnings_of(coloc(bf4, bf4, method = "single")), 2L)
})

test_that("a prior given as a plain list must still hold scalars", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    expect_error(coloc(bf, bf, prior = list(p1 = rep(1e-4, 4), p2 = 1e-4, p12 = 5e-6)),
                 "must each be a scalar")
})

## one mistake each, made to a copy of a good dataset
mistakes <- list(
    "no type" = function(d) { d$type <- NULL; d },
    "a wrong type" = function(d) { d$type <- "quantitative"; d },
    "an NA beta" = function(d) { d$beta[5] <- NA; d },
    "a duplicated snp" = function(d) { d$snp[2] <- d$snp[1]; d },
    "mismatched lengths" = function(d) { d$beta <- d$beta[-1]; d },
    "a zero varbeta" = function(d) { d$varbeta[5] <- 0; d })

test_that("a bad dataset stops coloc() with 6.0.3's coloc.abf() message, for either trait", {
    for (what in names(mistakes)) for (trait in 1:2) {
        d1 <- if (trait == 1) mistakes[[what]](D1) else D1
        d2 <- if (trait == 2) mistakes[[what]](D2) else D2
        expected <- tryCatch(get_match_6_0_3(coloc.abf(d1, d2)), error = conditionMessage)
        expect_match(expected, paste("dataset", trait), fixed = TRUE,
                     label = paste(what, "in trait", trait))
        expect_error(coloc(d1, d2), expected, fixed = TRUE,
                     label = paste(what, "in trait", trait))
    }
})

test_that("on the SuSiE route a bad dataset is reported with its own trait number", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    fit3 <- suppressMessages(runsusie(D3))
    fit4 <- suppressMessages(runsusie(D4))
    problem <- c("no type" = "variable type not set",
                 "a wrong type" = "type must be quant or cc",
                 "an NA beta" = "beta contains missing values",
                 "a duplicated snp" = "duplicated snps found",
                 "mismatched lengths" = "lengths of inputs don't match",
                 "a zero varbeta" = "Zero values in varbeta")
    for (what in names(mistakes)) {
        expect_error(suppressMessages(coloc(mistakes[[what]](D3), fit4, method = "susie")),
                     paste0("dataset 1: ", problem[[what]]), fixed = TRUE,
                     label = paste(what, "in trait 1"))
        expect_error(suppressMessages(coloc(fit3, mistakes[[what]](D4), method = "susie")),
                     paste0("dataset 2: ", problem[[what]]), fixed = TRUE,
                     label = paste(what, "in trait 2"))
    }
})

test_that("bad log BFs stop coloc() with the trait named", {
    good <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    bad <- good
    bad[1, 2] <- NA
    expect_error(coloc(bad, good), "trait 1: log Bayes factors contain missing")
    expect_error(coloc(good, bad), "trait 2: log Bayes factors contain missing")
    expect_error(coloc(c(1, 2, 3), good), "trait 1: a log-BF vector must be named by snp")
})

test_that("the region-size guard caps a prior as 6.0.3's coloc.abf() did", {
    ## 500 snps at p1 = p2 = 0.01 imply five causal variants per trait
    old_w <- warnings_of(old <- get_match_6_0_3(
        coloc.abf(D1, D2, p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    new_w <- warnings_of(new <- coloc(D1, D2,
                                      prior = coloc_prior(p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    expect_equal(new_w, old_w)
    expect_equal(new$priors, old$priors)
    expect_equal(new$priors[["p1"]], 1 / 501)
    expect_equal(as.list(pp(new)), as.list(old$summary))
})

test_that("the region-size guard starts where nsnps * p reaches 1", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    ## 4 snps at 0.25 is exactly 1: capped at 1 / (4 + 1)
    w <- warnings_of(res <- coloc(bf, bf, prior = coloc_prior(p1 = 0.25, p2 = 0.25,
                                                               p12 = 0.0625)))
    expect_equal(w, c("p1 * nsnps >= 1, setting p1=1/(nsnps + 1)",
                      "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)"))
    expect_equal(res$priors[["p1"]], 1 / 5)
    ## 4 snps at 0.2 is 0.8: nothing is capped
    expect_silent(res <- coloc(bf, bf, prior = coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.04)))
    expect_equal(res$priors[["p1"]], 0.2)
})

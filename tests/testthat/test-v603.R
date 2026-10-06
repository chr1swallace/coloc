## coloc() against 6.0.3. get_match_6_0_3(), in helper-v603.R, runs a call with
## 6.0.3's own functions: get_match_6_0_3(coloc.abf(D1, D2)) is what 6.0.3
## returned for coloc.abf(D1, D2).

data(coloc_test_data)
D1 <- coloc_test_data$D1
D2 <- coloc_test_data$D2
D3 <- coloc_test_data$D3
D4 <- coloc_test_data$D4

fit3 <- suppressMessages(runsusie(D3))
fit4 <- suppressMessages(runsusie(D4))
bf4 <- fit4$lbf_variable[fit4$sets$cs_index, , drop = FALSE]

## small log BF matrices over 8 snps: two signals in trait 1, one in trait 2
snps <- paste0("s", 1:8)
bf_a <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 7, 1, 0), nrow = 2, byrow = TRUE,
               dimnames = list(c("sig1", "sig2"), snps))
bf_b <- matrix(c(0, 0, 7, 3, 0, 0, 0, 0), nrow = 1, dimnames = list("sig1", snps))
with_null <- function(m) cbind(m, null = 0)

## the same datasets in other forms
as_cc <- function(d) {
    d$type <- "cc"
    d$s <- 0.3
    d$sdY <- NULL
    d
}
to_pvalues <- function(d) {
    d$pvalues <- 2 * pnorm(-abs(d$beta / sqrt(d$varbeta)))
    d[c("beta", "varbeta", "sdY")] <- NULL
    d
}
reverse <- function(d) {
    for (v in c("beta", "varbeta", "MAF", "snp", "position"))
        d[[v]] <- rev(d[[v]])
    d
}

## coloc.abf()'s SNP.PP.H4, in the snp order of a coloc() result
pp4_of <- function(abf, snps) abf$results$SNP.PP.H4[match(snps, abf$results$snp)]

## P(no causal variant) for each row of a log BF matrix, with the null at log
## BF 0 and prior p for each snp
p0 <- function(bf, p = 1e-4) {
    bf <- bf[, setdiff(colnames(bf), "null"), drop = FALSE]
    n <- ncol(bf)
    unname((1 - n * p) / ((1 - n * p) + p * rowSums(exp(bf))))
}

## ---- coloc() gives 6.0.3's answers ---------------------------------------------

test_that("coloc() matches 6.0.3's coloc.bf_bf(): two signals against one", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b))
    got <- coloc(bf_a, bf_b)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): named vectors", {
    v <- bf_a[1, ]
    expected <- get_match_6_0_3(coloc.bf_bf(v, v))
    got <- coloc(v, v)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): a null column", {
    for (pair in list(both = list(with_null(bf_a), with_null(bf_b)),
                      one = list(with_null(bf_a), bf_b),
                      ## nulllast() has to move it
                      not_last = list(with_null(bf_a)[, c("null", snps)], bf_b),
                      ## one row with a null column, the shape #210 broke
                      one_row = list(with_null(bf_b), with_null(bf_b)))) {
        expected <- get_match_6_0_3(coloc.bf_bf(pair[[1]], pair[[2]]))
        got <- coloc(pair[[1]], pair[[2]])
        expect_equal(got$summary, expected$summary)
        expect_equal(got$results, expected$results)
    }
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): partial snp overlap", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a[, 1:6], bf_b))
    got <- coloc(bf_a[, 1:6], bf_b)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): trimming off", {
    far <- bf_a[, 5:8]
    expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b, trim_by_posterior = FALSE))
    got <- coloc(far, bf_b, control = coloc_control(trim_by_posterior = FALSE))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.bf_bf(): p12 and prior weights", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p12 = 1e-5))
    got <- coloc(bf_a, bf_b, prior = coloc_prior(p12 = 1e-5))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$priors, expected$priors)

    w <- c(1, 1, 1, 1, 1, 1, 1, 100)
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, prior_weights1 = w, prior_weights2 = w))
    prior <- suppressWarnings(coloc_prior(weights1 = w, weights2 = w))   # AL15
    got <- coloc(bf_a, bf_b, prior = prior)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.abf()", {
    D1_sub <- suppressMessages(subset_dataset(D1, 1:400))
    D2_sub <- suppressMessages(subset_dataset(D2, 101:500))
    pairs <- list(quant = list(D1, D2),
                  cc = list(as_cc(D1), as_cc(D2)),
                  pvalues = list(to_pvalues(D1), to_pvalues(D2)),
                  mixed_type = list(D1, as_cc(D2)),
                  partial_overlap = list(D1_sub, D2_sub),
                  reversed = list(D1, reverse(D2)))
    for (nm in names(pairs)) {
        d1 <- pairs[[nm]][[1]]
        d2 <- pairs[[nm]][[2]]
        ## p12 = 1e-5 is coloc.abf()'s default (AL01)
        expected <- get_match_6_0_3(coloc.abf(d1, d2))
        got <- coloc(d1, d2, prior = coloc_prior(p12 = 1e-5))
        ## as lists, so each posterior meets the tolerance on its own scale
        expect_equal(as.list(pp(got)), as.list(expected$summary), label = nm)
        expect_equal(got$results$SNP.PP.H4.abf, pp4_of(expected, got$results$snp), label = nm)
        expect_equal(got$priors, expected$priors, label = nm)
    }
})

test_that("coloc() matches 6.0.3's coloc.abf(): sdY estimated", {
    d1 <- D1
    d2 <- D2
    d1$sdY <- d2$sdY <- NULL
    ## both estimate sdY, and say so
    expected <- suppressWarnings(get_match_6_0_3(coloc.abf(d1, d2)))
    got <- suppressWarnings(coloc(d1, d2, prior = coloc_prior(p12 = 1e-5)))
    expect_equal(as.list(pp(got)), as.list(expected$summary))
    expect_equal(got$results$SNP.PP.H4.abf, pp4_of(expected, got$results$snp))
})

test_that("coloc() stops on a bad dataset with 6.0.3's coloc.abf() message", {
    bad <- D2
    bad$beta[1] <- NA
    expected <- tryCatch(get_match_6_0_3(coloc.abf(D1, bad)), error = conditionMessage)
    expect_equal(expected, "dataset 2: beta contains missing values")
    expect_error(coloc(D1, bad), expected, fixed = TRUE)
})

test_that("no fixture reaches susieR's random subsampling", {
    ## the SuSiE comparisons fit each dataset once for each side, and agree only
    ## if fitting is deterministic. Its one random step, in susie_get_cs(),
    ## samples only from credible sets of more than n_purity snps
    threshold <- eval(formals(susieR::susie_get_cs)$n_purity)
    for (nm in c("D1", "D2", "D3", "D4")) {
        fit <- suppressMessages(runsusie(coloc_test_data[[nm]]))
        expect_lt(max(lengths(fit$sets$cs)), threshold, label = paste(nm, "largest credible set"))
    }
})

test_that("coloc() matches 6.0.3's coloc.susie(): two datasets", {
    expected <- suppressMessages(get_match_6_0_3(coloc.susie(D3, D4)))
    got <- suppressMessages(coloc(D3, D4, method = "susie"))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

test_that("coloc() matches 6.0.3's coloc.susie(): fits", {
    ## D4 has one credible set, the shape #210 broke when it came first
    for (pair in list(list(fit4, fit3), list(fit3, fit3))) {
        expected <- get_match_6_0_3(coloc.susie(pair[[1]], pair[[2]]))
        got <- coloc(pair[[1]], pair[[2]])
        expect_equal(got$summary, expected$summary)
        expect_equal(got$results, expected$results)
    }
})

test_that("coloc() matches 6.0.3's coloc.susie_bf()", {
    expected <- get_match_6_0_3(coloc.susie_bf(fit3, bf4))
    got <- coloc(fit3, bf4)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

## ---- today's versions of the changed functions, before and after ---------------
## These go at step 7, with the functions.

test_that("coloc.abf() returns what 6.0.3's did", {
    expected <- get_match_6_0_3(coloc.abf(D1, D2))
    invisible(utils::capture.output(got <- coloc.abf(D1, D2)))
    expect_equal(got, expected)
    expect_equal(as.list(got$summary), as.list(expected$summary))
})

test_that("coloc.bf_bf() returns what 6.0.3's did", {
    expect_equal(coloc.bf_bf(bf_a, bf_b), get_match_6_0_3(coloc.bf_bf(bf_a, bf_b)))
})

test_that("coloc.susie() returns what 6.0.3's did", {
    expect_equal(coloc.susie(fit3, fit4), get_match_6_0_3(coloc.susie(fit3, fit4)))
})

test_that("coloc.susie_bf() returns what 6.0.3's did", {
    expect_equal(coloc.susie_bf(fit3, bf4), get_match_6_0_3(coloc.susie_bf(fit3, bf4)))
})

test_that("check_dataset() names a missing type, where 6.0.3 failed obscurely", {
    d <- D3
    d$type <- NULL
    req <- c("beta", "varbeta", "LD", "snp", "N")   # runsusie()'s, without type
    expect_error(get_match_6_0_3(check_dataset(d, 1, req = req)), "argument is of length zero")
    expect_error(check_dataset(d, 1, req = req), "dataset 1: variable type not set")
    ## with type required, as coloc.abf() asks, the two agree
    expect_error(get_match_6_0_3(check_dataset(d, 1)), "missing required element")
    expect_error(check_dataset(d, 1), "missing required element")
})

## ---- the deliberate differences, one test each -------------------------------
## Each test shows 6.0.3's answer and the new one. The ids are those of the
## allowlist in the design notes.

test_that("AL01: coloc() defaults p12 to 5e-6, where coloc.abf() used 1e-5", {
    expected <- get_match_6_0_3(coloc.abf(D1, D2))
    expect_equal(expected$priors[["p12"]], 1e-5)
    expect_equal(coloc(D1, D2)$priors[["p12"]], 5e-6)
    expect_equal(as.list(pp(coloc(D1, D2, prior = coloc_prior(p12 = 1e-5)))),
                 as.list(expected$summary))
})

test_that("AL02: the priors given are used, which 6.0.3's coloc.susie_bf() ignored", {
    expected <- get_match_6_0_3(coloc.susie_bf(fit3, bf4, p12 = 1e-4))
    expect_equal(expected$priors[["p12"]], 5e-6)
    got <- coloc(fit3, bf4, prior = coloc_prior(p12 = 1e-4))
    expect_equal(got$priors[["p12"]], 1e-4)
})

test_that("AL08: PP.H0 in a trimmed result is a probability, where 6.0.3 left rounding noise", {
    far <- bf_a[, 5:8]
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b)), "overlap too small")
    expect_warning(got <- coloc(far, bf_b), "overlap too small")
    expect_true(all(abs(expected$summary$PP.H0.abf) < 1e-15))
    expect_equal(got$summary$PP.H0.abf, p0(far) * p0(bf_b))
    ## and nothing else moves
    cols <- c("nsnps", "hit1", "hit2", "idx1", "idx2")
    expect_equal(got$summary[, cols, with = FALSE], expected$summary[, cols, with = FALSE])
})

test_that("AL10: a vector p1 or p2 is rejected, where 6.0.3 reported 33 hypotheses", {
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8), p2 = rep(1e-4, 8)))
    expect_equal(sum(startsWith(names(expected$summary), "PP.H")), 33)
    expect_error(coloc(bf_a, bf_b, prior = coloc_prior(p1 = rep(1e-4, 8))), "single number")
    expect_error(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8)), "single number")
})

test_that("AL11: traits with no snps in common warn, where 6.0.3 was silent", {
    bf_x <- bf_b
    colnames(bf_x) <- paste0("x", 1:8)
    expect_silent(expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_x)))
    expect_true(is.na(expected$nsnps))
    expect_warning(got <- coloc(bf_a, bf_x), "no snps in common")
    expect_equal(got$status, "no_overlap")
})

test_that("AL12: an all-NA signal is rejected, where 6.0.3 reported PP.H2 = 1", {
    bf_na <- bf_a
    bf_na[1, ] <- NA
    expected <- get_match_6_0_3(coloc.bf_bf(bf_na, bf_b, trim_by_posterior = FALSE))
    expect_equal(expected$summary$PP.H2.abf[1], 1)
    expect_equal(expected$summary$hit1[1], "-")
    expect_error(coloc(bf_na, bf_b), "missing or infinite")
})

test_that("AL13: trimmed PP.H0 is the joint of the traits, where 6.0.3 took the smaller", {
    far <- with_null(bf_a[, 5:8])
    b <- with_null(bf_b)
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, b)), "overlap too small")
    expect_warning(got <- coloc(far, b), "overlap too small")
    expect_equal(expected$summary$PP.H0.abf, pmin(p0(far), p0(b)))
    expect_equal(got$summary$PP.H0.abf, p0(far) * p0(b))
    cols <- c("nsnps", "hit1", "hit2", "idx1", "idx2")
    expect_equal(got$summary[, cols, with = FALSE], expected$summary[, cols, with = FALSE])
})

test_that("AL14: a dataset with a BF matrix warns that method cannot apply to the matrix", {
    expected <- suppressMessages(get_match_6_0_3(coloc.susie_bf(D3, bf4)))
    expect_warning(got <- suppressMessages(coloc(D3, bf4, method = "susie")),
                   "method is inapplicable")
    expect_equal(got$summary, expected$summary)
})

test_that("AL15: coloc_prior() warns when both traits have weights, as coloc.abf() did", {
    w <- c(1, 1, 1, 1, 1, 1, 1, 100)
    expect_silent(expected <- get_match_6_0_3(
        coloc.bf_bf(bf_a, bf_b, prior_weights1 = w, prior_weights2 = w)))
    expect_warning(prior <- coloc_prior(weights1 = w, weights2 = w), "both traits")
    expect_equal(coloc(bf_a, bf_b, prior = prior)$summary, expected$summary)
})

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
## the same snps listed in the opposite order, as two studies sorted differently would be
reverse <- function(d) {
    for (v in c("beta", "varbeta", "MAF", "snp", "position"))
        d[[v]] <- rev(d[[v]])
    d
}

## coloc.abf()'s SNP.PP.H4, in the snp order of a coloc() result
pp4_of <- function(abf, snps) abf$results$SNP.PP.H4[match(snps, abf$results$snp)]

## Seven pairs of datasets for the two coloc.abf() comparisons below: coloc()
## against 6.0.3's coloc.abf(), and today's coloc.abf() against 6.0.3's.
## The first five cover the kinds of dataset coloc.abf() takes: quant or
## case-control traits, given as beta and varbeta or as p-values, and one
## trait of each type. In the last two, the traits share only some snps (300
## of 500), or list the same snps in a different order, so coloc.abf() must
## match them up by snp name.
D1_sub <- suppressMessages(subset_dataset(D1, 1:400))
D2_sub <- suppressMessages(subset_dataset(D2, 101:500))
abf_pairs <- list(quant = list(D1, D2),                            # beta and varbeta, sdY given
                  cc = list(as_cc(D1), as_cc(D2)),                 # both case-control
                  pvalues = list(to_pvalues(D1), to_pvalues(D2)),  # p-values and MAF only
                  cc_pvalues = list(to_pvalues(as_cc(D1)),         # case-control and p-values
                                    to_pvalues(as_cc(D2))),
                  mixed_type = list(D1, as_cc(D2)),                # one of each type
                  partial_overlap = list(D1_sub, D2_sub),          # 300 of 500 SNPs shared
                  reversed = list(D1, reverse(D2)))                # D2's SNPs in reverse order

## P(no causal variant) for each row of a log BF matrix, with the null at log
## BF 0 and prior p for each snp
p0 <- function(bf, p = 1e-4) {
    bf <- bf[, setdiff(colnames(bf), "null"), drop = FALSE]
    n <- ncol(bf)
    unname((1 - n * p) / ((1 - n * p) + p * rowSums(exp(bf))))
}

## ---- coloc() gives 6.0.3's answers ---------------------------------------------
## Both sides run with their default priors unless a test says otherwise:
## p1 = p2 = 1e-4 and p12 = 5e-6, the same in coloc() and in 6.0.3's
## coloc.bf_bf(), coloc.susie() and coloc.susie_bf(), so each test also checks
## that the defaults agree. coloc.abf()'s p12 default is 1e-5 (AL01), so its
## tests give p12 = 1e-5.

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
                      ## one row with a null column: a subset of one row
                      ## becomes a vector unless the code says drop = FALSE
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

test_that("coloc() matches 6.0.3's coloc.bf_bf(): trim_by_posterior = FALSE", {
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
    ## coloc_prior() warns when both traits have weights (AL15, below)
    prior <- suppressWarnings(coloc_prior(weights1 = w, weights2 = w))
    got <- coloc(bf_a, bf_b, prior = prior)
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
})

test_that("coloc() matches 6.0.3's coloc.abf()", {
    for (nm in names(abf_pairs)) {
        d1 <- abf_pairs[[nm]][[1]]
        d2 <- abf_pairs[[nm]][[2]]
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

test_that("SuSiE fits each test dataset the same way every time", {
    ## The SuSiE comparisons fit each dataset twice, once for each side, so they
    ## agree only if SuSiE gives the same fit both times. It has one random step:
    ## for a credible set of more than n_purity (100) snps, it checks a random 100
    ## of them. So every test dataset's credible sets must stay smaller than that.
    threshold <- eval(formals(susieR::susie_get_cs)$n_purity)
    for (nm in c("D1", "D2", "D3", "D4")) {
        fit <- suppressMessages(runsusie(coloc_test_data[[nm]]))
        expect_lt(max(lengths(fit$sets$cs)), threshold, label = paste(nm, "largest credible set"))
    }
})

test_that("coloc() matches 6.0.3's coloc.susie(): two coloc datasets", {
    expected <- suppressMessages(get_match_6_0_3(coloc.susie(D3, D4)))
    got <- suppressMessages(coloc(D3, D4, method = "susie"))
    expect_equal(got$summary, expected$summary)
    expect_equal(got$results, expected$results)
    expect_equal(got$priors, expected$priors)
})

test_that("coloc() matches 6.0.3's coloc.susie(): two SuSiE fits", {
    ## D4 has one credible set, so as trait 1 it is a one-row matrix: a subset
    ## of one row becomes a vector unless the code says drop = FALSE
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
## The tests of coloc.abf(), coloc.bf_bf(), coloc.susie() and coloc.susie_bf()
## go at step 7, the hard break, with those four functions.

test_that("coloc.abf() returns what 6.0.3's did", {
    for (nm in names(abf_pairs)) {
        d1 <- abf_pairs[[nm]][[1]]
        d2 <- abf_pairs[[nm]][[2]]
        expected <- get_match_6_0_3(coloc.abf(d1, d2))
        invisible(utils::capture.output(got <- coloc.abf(d1, d2)))
        expect_equal(got, expected, label = nm)
        expect_equal(as.list(got$summary), as.list(expected$summary), label = nm)
    }
})

test_that("coloc.bf_bf() returns what 6.0.3's did", {
    expect_equal(coloc.bf_bf(bf_a, bf_b), get_match_6_0_3(coloc.bf_bf(bf_a, bf_b)))
})

test_that("coloc.bf_bf() with every pair trimmed returns 6.0.3's result, apart from PP.H0", {
    ## PP.H0 is the one intended change here (AL08 and AL13, below)
    far <- bf_a[, 5:8]
    expect_warning(got <- coloc.bf_bf(far, bf_b), "overlap too small")
    expect_warning(expected <- get_match_6_0_3(coloc.bf_bf(far, bf_b)), "overlap too small")
    expect_equal(names(got), names(expected))
    expect_equal(sapply(got$summary, class), sapply(expected$summary, class))
    cols <- setdiff(names(expected$summary), "PP.H0.abf")
    expect_equal(got$summary[, cols, with = FALSE], expected$summary[, cols, with = FALSE])
})

test_that("coloc.susie() returns what 6.0.3's did", {
    expect_equal(coloc.susie(fit3, fit4), get_match_6_0_3(coloc.susie(fit3, fit4)))
})

test_that("coloc.susie_bf() returns what 6.0.3's did", {
    expect_equal(coloc.susie_bf(fit3, bf4), get_match_6_0_3(coloc.susie_bf(fit3, bf4)))
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

test_that("AL03: coloc.susie() now applies the region-size guard, which 6.0.3's never did", {
    ## Today's coloc.susie() runs through coloc(), which applies the region-size
    ## guard; 6.0.3's ran through coloc.bf_bf(), which does not. With 500 SNPs,
    ## p1 = p2 = 0.01 implies five causal variants in each trait, so coloc() caps
    ## both at 1/501. coloc.susie_bf() also runs through coloc() now, but cannot
    ## show this: 6.0.3's ignored the priors it was given (AL02).
    w <- capture_warnings(expected <- get_match_6_0_3(
        coloc.susie(fit3, fit4, p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    expect_length(w, 0)
    expect_equal(expected$priors, c(p1 = 0.01, p2 = 0.01, p12 = 1e-4))
    w <- capture_warnings(got <- coloc.susie(fit3, fit4, p1 = 0.01, p2 = 0.01, p12 = 1e-4))
    expect_equal(w, c("p1 * nsnps >= 1, setting p1=1/(nsnps + 1)",
                      "p2 * nsnps >= 1, setting p2=1/(nsnps + 1)"))
    expect_equal(got$priors, c(p1 = 1 / 501, p2 = 1 / 501, p12 = 1e-4))
    expect_gt(got$summary$PP.H4.abf[1], expected$summary$PP.H4.abf[1])
    ## and nothing else moves: the same as 6.0.3 given the capped priors
    capped <- get_match_6_0_3(coloc.susie(fit3, fit4, p1 = 1 / 501, p2 = 1 / 501, p12 = 1e-4))
    expect_equal(got$summary, capped$summary)
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

test_that("AL10: a vector p1 or p2 is rejected, where 6.0.3 allowed it, with extra hypotheses", {
    ## p1 and p2 of length 8, one per snp: 6.0.3 accepted them and reported
    ## 4 * 8 + 1 = 33 hypotheses in place of five; the count grows with the
    ## length of the vectors
    expected <- get_match_6_0_3(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8), p2 = rep(1e-4, 8)))
    expect_equal(sum(startsWith(names(expected$summary), "PP.H")), 33)
    expect_error(coloc(bf_a, bf_b, prior = coloc_prior(p1 = rep(1e-4, 8))),
                 "must each be a scalar")
    expect_error(coloc.bf_bf(bf_a, bf_b, p1 = rep(1e-4, 8)), "must each be a scalar")
})

test_that("AL11: traits with no snps in common warn, where 6.0.3 returned nsnps = NA silently", {
    ## 6.0.3's coloc.bf_bf() returned early with data.table(nsnps = NA) and no warning;
    ## coloc() warns "no snps in common" and returns status "no_overlap"
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

test_that("AL16: a p-value above 1 is refused, which 6.0.3 accepted", {
    msg <- "dataset 1: pvalues must be above 0 and at most 1"
    above <- to_pvalues(D1)
    above$pvalues[5] <- 1.5
    expect_null(get_match_6_0_3(check_dataset(above, 1)))
    expect_error(check_dataset(above, 1), msg, fixed = TRUE)
    expect_error(coloc(above, D2), msg, fixed = TRUE)
    expect_error(finemap.abf(above), "pvalues must be above 0 and at most 1", fixed = TRUE)
    ## 0 and below were refused already; only the message moved
    zero <- to_pvalues(D1)
    zero$pvalues[5] <- 0
    expect_error(get_match_6_0_3(check_dataset(zero, 1)),
                 "pvalues should not be negative or exactly 0")
    expect_error(check_dataset(zero, 1), msg, fixed = TRUE)
    ## exactly 1 is a p-value
    one <- to_pvalues(D1)
    one$pvalues[5] <- 1
    expect_null(check_dataset(one, 1))
})

test_that("AL17: check_dataset() reports a dataset with no type, where 6.0.3 gave R's error", {
    ## runsusie() leaves type out of req, so check_dataset() has to notice a
    ## dataset without one by itself
    d <- D3
    d$type <- NULL
    req <- c("beta", "varbeta", "LD", "snp", "N")   # runsusie()'s, without type
    expect_error(get_match_6_0_3(check_dataset(d, 1, req = req)), "argument is of length zero")
    expect_error(check_dataset(d, 1, req = req), "dataset 1: variable type not set")
    ## with type required, as coloc.abf() asks, the two agree
    expect_error(get_match_6_0_3(check_dataset(d, 1)), "missing required element")
    expect_error(check_dataset(d, 1), "missing required element")
})

test_that("AL18: a SuSiE fit with a duplicated snp id is refused, where 6.0.3 used the first", {
    ## today's coloc.susie() runs through coloc(), which checks a fit's log
    ## BFs as it checks log BFs given directly; 6.0.3's went to coloc.bf_bf()
    dup <- fit4
    colnames(dup$lbf_variable)[2] <- colnames(dup$lbf_variable)[1]
    expected <- get_match_6_0_3(coloc.susie(fit3, dup))
    expect_equal(nrow(expected$summary), 2L)
    expect_error(coloc.susie(fit3, dup),
                 "trait 2: duplicated snps found in the log Bayes factors")
    ## an NA in a credible set's log BFs stopped 6.0.3 too, with R's own message
    na <- fit4
    na$lbf_variable[na$sets$cs_index[1], 10] <- NA
    expect_error(get_match_6_0_3(coloc.susie(fit3, na)), "missing value where TRUE/FALSE needed")
    expect_error(coloc.susie(fit3, na),
                 "trait 2: log Bayes factors contain missing or infinite values")
})

## ---- what both sides share ---------------------------------------------------
## The 6.0.3 copies and coloc() both use the package's process.dataset(), with
## the approx.bf.estimates(), approx.bf.p() and sdY.est() it calls, and
## runsusie(). They are not copied and are unchanged since 6.0.3, so an edit to
## one would move both sides together, unseen above. These values, worked out
## with them as they are, pin them. logsum(), logdiff() and logbf_to_pp() are
## copied into helper-v603.R with the 6.0.3 functions, so the comparisons above
## already see an edit to the package's versions.

test_that("the shared code still gives 6.0.3's log ABFs", {
    at_s105 <- function(d, suffix) {
        pd <- process.dataset(d, suffix = suffix)
        pd[[paste0("lABF.", suffix)]][pd$snp == "s105"]
    }
    ## each value is written to 12 digits, so compared to 10
    expect_equal(at_s105(D1, "df1"), 28.343423232, tolerance = 1e-10)
    expect_equal(at_s105(to_pvalues(D1), "df1"), 29.9818369694, tolerance = 1e-10)
    expect_equal(at_s105(as_cc(D2), "df2"), 25.698511433, tolerance = 1e-10)
    expect_equal(suppressWarnings(sdY.est(D1$varbeta, D1$MAF, D1$N)), 1.28815854798,
                 tolerance = 1e-10)
})

test_that("the shared code still finds D3's credible sets", {
    rows <- fit3$lbf_variable[fit3$sets$cs_index, ]
    expect_equal(fit3$sets$cs_index, c(1, 2))
    expect_equal(unname(lengths(fit3$sets$cs)), c(3L, 13L))
    expect_equal(colnames(rows)[apply(rows, 1, which.max)], c("s105", "s89"))
    ## a SuSiE fit can differ across machines by about 1e-8, so 6 digits
    expect_equal(apply(rows, 1, max), c(25.82456833, 13.32711173), tolerance = 1e-6)
})

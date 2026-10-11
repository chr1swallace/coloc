## The engine, .coloc_engine(), tested on its own, with small log BF matrices
## written into each test: each status it can return, what every result must
## satisfy (the five posteriors sum to 1), and values the test works out from
## a formula itself, such as the softmax over the shared snps.
## test-v603.R compares its answers, through coloc(), with 6.0.3's;
## get_match_6_0_3(), in helper-v603.R, runs a call as 6.0.3 made it.

## P(no causal variant) for each row of a log BF matrix, with the null at log
## BF 0 and prior p for each snp. Written out, not taken from logbf_to_pp(),
## so the tests below do not lean on the code they check.
p0 <- function(bf, p) {
    bf <- bf[, setdiff(colnames(bf), "null"), drop = FALSE]
    n <- ncol(bf)
    unname((1 - n * p) / ((1 - n * p) + p * rowSums(exp(bf))))
}

snps <- paste0("s", 1:8)
a <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0), nrow = 1, dimnames = list("sig1", snps))
b <- matrix(c(0, 0, 7, 3, 0, 0, 0, 0), nrow = 1, dimnames = list("sigA", snps))

test_that("every status is reachable, and each result validates", {
    dis <- b
    colnames(dis) <- paste0("x", snps)
    na1 <- a
    na1[1, ] <- NA
    got <- list(ok = .coloc_engine(a, b),
                no_signal = suppressWarnings(.coloc_engine(a[0, , drop = FALSE], b)),
                no_overlap = suppressWarnings(.coloc_engine(a, dis)),
                all_trimmed = suppressWarnings(.coloc_engine(a[, 5:8, drop = FALSE], b)),
                degenerate = suppressWarnings(.coloc_engine(na1, b)))
    expect_equal(vapply(got, `[[`, character(1), "status"),
                 c(ok = "ok", no_signal = "no_signal", no_overlap = "no_overlap",
                   all_trimmed = "all_trimmed", degenerate = "degenerate"))
    for (nm in names(got)) expect_silent(validate_coloc_result(got[[nm]]))
    ## rows only where a pair could be formed; NA only where it could not be computed
    expect_equal(nrow(got$no_signal$summary), 0L)
    expect_equal(nrow(got$no_overlap$summary), 0L)
    expect_true(all(is.na(got$all_trimmed$summary$PP.H4.abf)))
    expect_false(anyNA(got$ok$summary$PP.H4.abf))
    ## $messages keeps what warning() said
    expect_match(got$no_overlap$messages, "no snps in common")
})

test_that("the five posteriors sum to 1, and H4 beats H3 only when the signals coincide", {
    peak <- function(at) {
        v <- rep(0, 6)
        v[at] <- 8
        matrix(v, nrow = 1, dimnames = list("sig1", paste0("s", 1:6)))
    }
    same <- .coloc_engine(peak(2), peak(2))
    diff <- .coloc_engine(peak(2), peak(5))
    for (r in list(same, diff))
        expect_equal(sum(unlist(r$summary[1, paste0("PP.H", 0:4, ".abf"), with = FALSE])), 1)
    expect_gt(same$summary$PP.H4.abf, same$summary$PP.H3.abf)
    expect_gt(diff$summary$PP.H3.abf, diff$summary$PP.H4.abf)
})

test_that("SNP.PP.H4 is the softmax of bf1 + bf2 over the shared snps", {
    x <- matrix(c(0, 4, 1), nrow = 1, dimnames = list("sig1", c("s1", "s2", "s3")))
    y <- matrix(c(0, 3, 2), nrow = 1, dimnames = list("sigA", c("s1", "s2", "s3")))
    res <- .coloc_engine(x, y)
    lsum <- as.numeric(x) + as.numeric(y)
    expect_equal(res$results$SNP.PP.H4.abf, exp(lsum) / sum(exp(lsum)))
})

test_that("one shared snp is enough", {
    x <- matrix(c(0, 0, 0, 5), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    y <- matrix(c(5, 0, 0, 0), nrow = 1, dimnames = list("sig1", paste0("s", 4:7)))
    res <- .coloc_engine(x, y, control = coloc_control(trim_by_posterior = FALSE))
    expect_equal(res$summary$nsnps, 1L)
    expect_equal(res$results$snp, "s4")
    expect_equal(res$results$SNP.PP.H4.abf, 1)
    expect_equal(sum(unlist(res$summary[1, paste0("PP.H", 0:4, ".abf"), with = FALSE])), 1)
})

test_that("prior weights move the posteriors, and are kept for the shared snps", {
    with_weights <- function(w1, w2)
        modifyList(.default_prior(), list(weights1 = w1, weights2 = w2))
    ## equal weights are no weights: the weighted sum is the unweighted one
    flat <- .coloc_engine(a, b, prior = with_weights(rep(1, 8), rep(1, 8)))
    expect_equal(flat$summary$PP.H4.abf, .coloc_engine(a, b)$summary$PP.H4.abf)
    ## weight on the snp both traits peak at makes a shared variant likelier
    w <- c(1, 1, 10, 1, 1, 1, 1, 1)
    up <- .coloc_engine(a, b, prior = with_weights(w, w))
    expect_gt(up$summary$PP.H4.abf, flat$summary$PP.H4.abf)
    ## with partial overlap, $weights holds those of the shared snps
    part <- .coloc_engine(a[, 1:6, drop = FALSE], b,
                          prior = with_weights(1:6, rep(1, 8)))
    expect_equal(part$weights, list(w1 = 1:6, w2 = rep(1, 6)))
    expect_error(.coloc_engine(a, b, prior = with_weights(1:3, NULL)),
                 "Length of prior_weights1")
})

test_that("a trimmed pair with no evidence reports the prior, not a certainty", {
    ## Neither trait has any evidence (every log BF 0). 3 of 8 snps shared puts
    ## 0.375 of each posterior on them, under overlap.min, so the pair is
    ## trimmed. 6.0.3 then reported PP.H0 as rounding noise (AL08) and took the
    ## smaller of the traits' values (AL13); both are in test-v603.R
    x <- matrix(0, 1, 8, dimnames = list("sig1", paste0("s", 1:8)))
    y <- matrix(0, 1, 8, dimnames = list("sigA", paste0("s", 6:13)))
    res <- suppressWarnings(.coloc_engine(x, y))
    expect_equal(res$status, "all_trimmed")
    ## each trait's P(no causal variant) is its prior, 1 - 8p; H0 is the joint
    expect_equal(p0(x, 1e-4), 1 - 8 * 1e-4)
    expect_equal(res$summary$PP.H0.abf, (1 - 8 * 1e-4)^2)
    ## 6.0.3's 1 - rowSums() of a row that already sums to 1 left noise
    expect_lt(abs(1 - sum(logbf_to_pp(x, 1e-4, last_is_null = FALSE))), 1e-12)
})

test_that("trimmed pairs report PP.H0 as the product of each trait's", {
    res <- suppressWarnings(.coloc_engine(a[, 5:8, drop = FALSE], b))
    expect_equal(res$summary$PP.H0.abf, p0(a[, 5:8, drop = FALSE], 1e-4) * p0(b, 1e-4))
})

test_that("an all-NA signal changes only its own pair, which 6.0.3 rewrote to a certainty", {
    ## 6.0.3 reported an all-NA signal as PP.H2 = 1 (AL12 in test-v603.R).
    ## coloc() rejects the NA, but today's coloc.bf_bf() still brings it to the
    ## engine, until step 7 (the hard break, in a later stage)
    x <- rbind(a, matrix(c(0, 0, 0, 0, 0, 7, 1, 0), nrow = 1, dimnames = list("sig2", snps)))
    x[1, ] <- NA
    expected <- get_match_6_0_3(coloc.bf_bf(x, b, trim_by_posterior = FALSE))
    got <- suppressWarnings(.coloc_engine(x, b,
                                          control = coloc_control(trim_by_posterior = FALSE)))
    expect_equal(expected$summary$PP.H2.abf[1], 1)
    expect_true(all(is.na(unlist(got$summary[1, paste0("PP.H", 0:4, ".abf"), with = FALSE]))))
    expect_equal(got$summary[2], expected$summary[2])
})

test_that("one NA in a row marks that row's pairs, without trimming", {
    ## 6.0.3 stopped on any NA in drop; the engine keeps the pair, marked NA
    bad <- a
    bad[1, 2] <- NA
    x <- rbind(bad, a)
    rownames(x) <- c("bad", "good")
    res <- suppressWarnings(.coloc_engine(x, a))
    expect_equal(res$status, "degenerate")
    expect_equal(nrow(res$summary), 2L)
    expect_true(is.na(res$summary$PP.H4.abf[res$summary$idx1 == 1]))
    expect_false(is.na(res$summary$PP.H4.abf[res$summary$idx1 == 2]))
})

test_that("the engine counts the signal pairs in $meta, and leaves the method to coloc()", {
    x <- a[, 1:4, drop = FALSE]
    y <- rbind(x, x)
    rownames(y) <- c("sigA", "sigB")
    res <- .coloc_engine(x, y)
    expect_equal(res$meta$n_pairs_in, 2L)
    expect_equal(res$meta$n_pairs_kept, 2L)
    expect_equal(res$meta$schema, 2L)
    ## the engine sees only log BF matrices, not which adapter made them, so
    ## coloc() records the method after it returns (R/coloc.R:171)
    expect_null(res$meta$method)
})

## ---------------------------------------------------------------------------
## Step 1: the engine. Private, and since steps 2 and 5 it has three callers:
## coloc(), coloc.bf_bf() and coloc.abf().
##
## The engine's own shape is new, so bit-identity THROUGH those callers is
## step 2's gate and lives in test-regression.R. What is asserted here is that
## the contract holds, that every status is reachable, and that the NUMBERS
## match 6.0.3 wherever 6.0.3 produced any.
## ---------------------------------------------------------------------------

## P(no causal variant) for each row of a bf matrix, from first principles:
## with the null at log BF 0 and a scalar prior p per snp,
##   P(null) = (1 - n*p) / ((1 - n*p) + p * sum(exp(bf)))
## Written out rather than calling logbf_to_pp() so that the expectations below
## are independent of the code they check. Handles both shapes: a "null" column
## is excluded from the snp count and from the sum, being the 1 in the formula.
.p0 <- function(bf, p) {
    sn <- setdiff(colnames(bf), "null")
    n  <- length(sn)
    unname(apply(bf[, sn, drop = FALSE], 1, function(r)
        (1 - n * p) / ((1 - n * p) + p * sum(exp(r)))))
}

## 6.0.3's $summary, or NULL for the shapes that carry none
.old_summary <- function(v) {
    if (inherits(v, "coloc_case_error") || !is.list(v)) return(NULL)
    if (!"summary" %in% names(v)) return(NULL)          ## bare data.table(nsnps=NA)
    as.data.frame(v$summary)
}

test_that("validate_coloc_result() rejects every breach of the frozen contract", {
    ok <- .new_coloc_result(prior = .default_prior(), status = "no_signal")
    expect_s3_class(ok, "coloc_result")
    expect_silent(validate_coloc_result(ok))

    ## a drifted column set
    bad <- ok; bad$summary$extra <- integer(0)
    expect_error(validate_coloc_result(bad), "drifted")

    ## a missing field
    bad <- ok; bad$priors <- NULL
    expect_error(validate_coloc_result(bad), "missing required field")

    ## an unknown status
    bad <- ok; bad$status <- "fine"
    expect_error(validate_coloc_result(bad), "unknown status")

    ## the governing rule, both halves
    bad <- ok; bad$status <- "all_trimmed"
    expect_error(validate_coloc_result(bad), "keeps its rows")

    bad <- ok; bad$meta$schema <- NULL
    expect_error(validate_coloc_result(bad), "schema")
})

test_that("all five statuses are reachable and each one validates", {
    a <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0), nrow = 1,
                dimnames = list("sig1", paste0("s", 1:8)))
    b <- matrix(c(0, 0, 7, 3, 0, 0, 0, 0), nrow = 1,
                dimnames = list("sigA", paste0("s", 1:8)))
    dis <- b; colnames(dis) <- paste0("x", colnames(b))
    na1 <- a; na1[1, ] <- NA_real_

    got <- list(
        ok          = .coloc_engine(a, b),
        no_signal   = suppressWarnings(.coloc_engine(a[0, , drop = FALSE], b)),
        no_overlap  = suppressWarnings(.coloc_engine(a, dis)),
        all_trimmed = suppressWarnings(.coloc_engine(a[, 5:8, drop = FALSE], b)),
        degenerate  = suppressWarnings(.coloc_engine(na1, b)))

    expect_equal(vapply(got, `[[`, character(1), "status"),
                 c(ok = "ok", no_signal = "no_signal", no_overlap = "no_overlap",
                   all_trimmed = "all_trimmed", degenerate = "degenerate"))
    for (nm in names(got)) expect_silent(validate_coloc_result(got[[nm]]))

    ## rows exist iff a pair was considered
    expect_equal(nrow(got$no_signal$summary), 0L)
    expect_equal(nrow(got$no_overlap$summary), 0L)
    expect_gt(nrow(got$all_trimmed$summary), 0L)

    ## posteriors are NA iff they could not be computed. PP.H0 is the one
    ## the trimmed branch does compute, and since AL08 it is a real number
    ## rather than the 1 - 1 that 1 - rowSums() used to leave there.
    expect_equal(got$all_trimmed$summary$PP.H0.abf,
                 .p0(a[, 5:8, drop = FALSE], 1e-4) * .p0(b, 1e-4),
                 tolerance = 1e-10)
    expect_true(all(is.na(got$all_trimmed$summary$PP.H4.abf)))
    expect_true(all(is.na(got$degenerate$summary$PP.H4.abf)))
    expect_false(anyNA(got$ok$summary$PP.H4.abf))

    ## and $messages carries what warning() emitted
    expect_match(got$no_overlap$messages, "no snps in common")
})

test_that("the engine reproduces 6.0.3's numbers on every golden case", {
    for (case in coloc_cases) {
        g <- read_golden(case$id)
        os <- .old_summary(g$layerB$value)
        if (is.null(os)) next                       ## 6.0.3 produced no summary

        prior <- modifyList(.default_prior(),
                            case$args[intersect(names(case$args),
                                                c("p1", "p2", "p12"))])
        if (!is.null(case$args$prior_weights1)) prior$weights1 <- case$args$prior_weights1
        if (!is.null(case$args$prior_weights2)) prior$weights2 <- case$args$prior_weights2
        control <- modifyList(.default_control(),
                              case$args[intersect(names(case$args),
                                                  c("overlap.min", "trim_by_posterior"))])

        ## the engine takes matrices; promoting a named vector is .adapt_bf()'s
        ## job from step 3, and coloc.bf_bf()'s until then
        promote <- function(x) if (is.null(dim(x)))
            matrix(x, nrow = 1L, dimnames = list("sig1", names(x))) else x

        new <- suppressWarnings(.coloc_engine(promote(g$layerA$bf1),
                                              promote(g$layerA$bf2),
                                              prior = prior, control = control))
        ns <- as.data.frame(new$summary)

        ## The three cases that reach the all-trimmed branch, where PP.H0 has
        ## moved twice: AL08 gave ph0 a null to be the probability of, and
        ## AL13 combines the two traits' values as a joint rather than a
        ## minimum. Assert the new value, and assert the change is CONFINED to
        ## that one column.
        if (case$id %in% c("t0-all-trimmed", "t0-overlap-min",
                           "t0-all-trimmed-null")) {         ## AL08, AL13
            a0 <- .p0(g$layerA$bf1, 1e-4); b0 <- .p0(g$layerA$bf2, 1e-4)
            expect_equal(ns$PP.H0.abf, a0 * b0, tolerance = 1e-10)   ## AL13

            if (case$id == "t0-all-trimmed-null") {
                ## This case is AL08's control: 6.0.3 has a null column here,
                ## so it computed ph0 properly and the engine's ph0 is its own
                ## expression unchanged. AL13 moves the OUTPUT column, so the
                ## control can no longer be "differs in nothing"; it is made
                ## one level down instead, on the inputs AL13 combines.
                expect_equal(pmin(a0, b0), as.numeric(os$PP.H0.abf),
                             tolerance = 1e-10)
            } else {
                ## AL08: no null column, so 6.0.3's 1 - rowSums() of a row
                ## that already sums to 1 left noise, of either sign
                expect_true(all(abs(os$PP.H0.abf) < 1e-15))
            }

            for (col in c("nsnps", "hit1", "hit2", "idx1", "idx2"))
                expect_equal(ns[[col]], os[[col]], tolerance = 0,
                             label = paste0(case$id, " $", col))
            for (col in paste0("PP.H", 1:4, ".abf"))
                expect_true(all(is.na(ns[[col]])), label = paste0("ns$", col))
            next
        }

        ## AL10, pending in helper-allowlist.R: 6.0.3 recycles a per-snp
        ## prior through combine.abf() and reports 33 hypotheses. The engine
        ## collapses it to five, so there is nothing to compare column-wise.
        if (case$id == "t0-vector-priors") {          ## AL10
            expect_equal(ncol(os), 38L)                 ## 6.0.3: PP.H0..PP.H32
            expect_equal(ncol(ns), 10L)                 ## engine: the frozen shape
            next
        }

        ## AL12: 6.0.3 rewrites an all-NA pair to a certainty. Assert the
        ## divergence rather than skipping past it, and assert that it is
        ## CONFINED to that pair: the ordinary row beside it still matches to
        ## the last bit, which is what makes this a targeted change and not a
        ## shift in the arithmetic.
        if (case$id == "t0-degenerate") {             ## AL12
            hh <- paste0("PP.H", 0:4, ".abf")
            expect_equal(os$hit1, c("-", "s6"))       ## 6.0.3 names a hit it did not find
            expect_equal(as.numeric(unlist(os[1, hh])), c(0, 0, 1, 0, 0))
            expect_true(all(is.na(unlist(ns[1, hh]))))
            expect_true(is.na(ns$hit1[1]))
            for (col in c("nsnps", "hit1", "hit2", hh))
                expect_equal(ns[[col]][2], os[[col]][2], tolerance = 0,
                             label = paste0("t0-degenerate row 2 $", col))
            next
        }

        expect_equal(nrow(ns), nrow(os), label = paste0(case$id, " n rows"))
        for (col in c("nsnps", "hit1", "hit2", "idx1", "idx2"))
            expect_equal(ns[[col]], os[[col]], tolerance = 0,
                         label = paste0(case$id, " $summary$", col))
        ## 6.0.3 writes a bare NA into the all-trimmed branch, which makes the
        ## PP columns logical there while the frozen contract says numeric.
        ## Compare the values; coloc.bf_bf()'s wrapper restores the old type
        ## at step 2.
        for (col in paste0("PP.H", 0:4, ".abf"))
            expect_equal(as.numeric(ns[[col]]), as.numeric(os[[col]]), tolerance = 0,
                         label = paste0(case$id, " $summary$", col))
    }
})

test_that("a mixed drop, part NA part FALSE, does not trim and marks only the bad pair", {
    ## prop is a ratio of rowSums, so ONE NA anywhere in a row makes coverage
    ## unassessable for every pair using that row. With two rows, one poisoned
    ## and one clean, drop is c(NA, FALSE).
    ##
    ## This is the case the engine's drop[is.na(drop)] <- FALSE was written for
    ## and, until this test, the case no test covered: the existing degenerate
    ## case uses a single row, so drop is a lone NA and only exercises 6.0.3's
    ## if(all(drop)) guard. all() short-circuits on FALSE and any() on TRUE, so
    ## c(NA, FALSE) gets past if(all(drop)) at R/susie.R:285 and is caught by
    ## if(any(drop)) at :297 instead.
    sn  <- paste0("s", 1:8)
    ok  <- matrix(c(0, 0, 8, 2, 0, 0, 0, 0), nrow = 1, dimnames = list("sigA", sn))
    bad <- ok; bad[1, 2] <- NA_real_              ## a single NA is enough
    mixed <- rbind(bad, ok); rownames(mixed) <- c("bad", "good")

    res <- suppressWarnings(.coloc_engine(mixed, ok))
    expect_equal(res$status, "degenerate")
    expect_equal(nrow(res$summary), 2L)           ## nothing was trimmed
    expect_silent(validate_coloc_result(res))

    ## posteriors are NA iff they could not be computed, applied per row
    expect_true(is.na(res$summary$PP.H4.abf[res$summary$idx1 == 1]))
    expect_false(is.na(res$summary$PP.H4.abf[res$summary$idx1 == 2]))
})

test_that("6.0.3 errors on both NA shapes, which is what the guard replaces", {
    ## Pinning the behaviour being departed from. all() and any() short-circuit
    ## on opposite values, so between them 6.0.3 rejects every NA-bearing drop.
    expect_true(is.na(all(NA)))                    ## lone NA  -> if(all(...)) errors
    expect_false(all(c(NA, FALSE)))                ## mixed    -> all() is FALSE...
    expect_true(is.na(any(c(NA, FALSE))))          ##          ...but any() is NA
    expect_error(if (all(NA)) TRUE, "missing value")
    expect_error(if (any(c(NA, FALSE))) TRUE, "missing value")
})

test_that("a trimmed pair with no evidence reports the prior, not a certainty", {
    ## Every log BF is 0, so the data say nothing and P(no causal variant) is
    ## exactly the prior, 1 - n*p1. This is the case AL08 was written for: with
    ## no null column 6.0.3 reports ~0 for it, a certainty that both traits are
    ## associated, drawn from two traits that showed no evidence at all.
    ##
    ## Three snps in common out of eight makes prop 0.375 for both traits,
    ## under the 0.5 default, so every pair trims and the branch is reached.
    a <- matrix(0, 1, 8, dimnames = list("sig1", paste0("s", 1:8)))
    b <- matrix(0, 1, 8, dimnames = list("sigA", paste0("s", 6:13)))
    res <- suppressWarnings(.coloc_engine(a, b))

    expect_equal(res$status, "all_trimmed")
    ## each trait's own marginal is the prior, 1 - n*p1; the reported column
    ## is the JOINT, so it is the product. Under 6.0.3's pmin it would be one
    ## trait's 0.9992 rather than 0.99840064. AL13.
    expect_equal(res$summary$PP.H0.abf, (1 - 8 * 1e-4)^2, tolerance = 1e-12)
    expect_equal(.p0(a, 1e-4), 1 - 8 * 1e-4, tolerance = 1e-12)
    ## and the expression it replaces, so the test says why, not just what:
    ## the row already sums to 1, so 1 - rowSums() has nothing left to report
    expect_lt(abs(1 - sum(logbf_to_pp(a, 1e-4, last_is_null = FALSE))), 1e-12)
})

test_that("SNP.PP.H4 is the softmax of bf1 + bf2 over the shared snps", {
    a <- matrix(c(0, 4, 1), nrow = 1, dimnames = list("sig1", c("s1", "s2", "s3")))
    b <- matrix(c(0, 3, 2), nrow = 1, dimnames = list("sigA", c("s1", "s2", "s3")))
    res <- .coloc_engine(a, b)

    lsum <- as.numeric(a) + as.numeric(b)
    by_hand <- exp(lsum) / sum(exp(lsum))          ## no logsum(), computed directly
    expect_equal(res$results$SNP.PP.H4.abf, by_hand, tolerance = 1e-12)
    expect_equal(sum(res$results$SNP.PP.H4.abf), 1, tolerance = 1e-12)
})

test_that("the five posteriors sum to 1, and H4 beats H3 only when signals coincide", {
    snps <- paste0("s", 1:6)
    peak <- function(at) {
        v <- rep(0, 6); v[at] <- 8
        matrix(v, nrow = 1, dimnames = list("sig1", snps))
    }
    same <- .coloc_engine(peak(2), peak(2))
    diff <- .coloc_engine(peak(2), peak(5))

    for (r in list(same, diff))
        expect_equal(sum(unlist(r$summary[1, paste0("PP.H", 0:4, ".abf"), with = FALSE])),
                     1, tolerance = 1e-12)
    expect_gt(same$summary$PP.H4.abf, same$summary$PP.H3.abf)
    expect_gt(diff$summary$PP.H3.abf, diff$summary$PP.H4.abf)
})

test_that("$meta records what the engine knows, and not what it cannot", {
    a <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    b <- rbind(a, a); rownames(b) <- c("sigA", "sigB")
    res <- .coloc_engine(a, b)
    expect_equal(res$meta$n_pairs_in, 2L)
    expect_equal(res$meta$n_pairs_kept, 2L)
    expect_equal(res$meta$schema, 1L)
    ## the engine cannot know which adapter ran: coloc() fills that in later
    expect_null(res$meta$method)
})

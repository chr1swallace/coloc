## coloc_prior() and coloc_control(), and the checks behind them

test_that("coloc_prior() holds the defaults, declared once", {
    p <- coloc_prior()
    expect_s3_class(p, "coloc_prior")
    expect_equal(c(p$p1, p$p2, p$p12), c(1e-4, 1e-4, 5e-6))
    expect_null(p$weights1)
    ## every path takes its defaults from here; 6.0.3 declared them once per
    ## entry point, which is how p12 came to differ (AL01)
    expect_equal(.default_prior(), unclass(coloc_prior()))
    expect_equal(.default_control(), unclass(coloc_control()))
})

test_that("p1, p2 and p12 must each be a single number", {
    for (arg in c("p1", "p2", "p12")) {
        a <- list(rep(1e-4, 8))
        names(a) <- arg
        expect_error(do.call(coloc_prior, a), "single number.*weights1/weights2", label = arg)
    }
    ## checked before the bounds, which would accept a vector
    expect_error(coloc_prior(p1 = c(1e-4, 1e-4), weights1 = c(1, 1)), "single number")
    expect_error(.check_scalar_prior(list(p1 = 1e-4, p2 = c(1e-4, 1e-4), p12 = 5e-6)),
                 "single number")
    expect_silent(.check_scalar_prior(.default_prior()))
})

test_that("p12 may sit exactly on its bounds", {
    ## 0.2 * 0.2 rounds to 0.04000000000000001, just above a p12 of 0.04
    expect_s3_class(coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.04), "coloc_prior")
    expect_s3_class(coloc_prior(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4), "coloc_prior")
    expect_error(coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.039), "p12 must lie")
    expect_error(coloc_prior(p12 = 2e-4), "p12 must lie")
})

test_that("prior weights must be positive and not NA", {
    expect_equal(coloc_prior(weights1 = c(1, 2, 3))$weights1, c(1, 2, 3))
    expect_error(coloc_prior(weights1 = c(1, -1)), "non-NA and positive")
    expect_error(coloc_prior(weights1 = c(1, 0)), "non-NA and positive")
    expect_error(coloc_prior(weights2 = c(1, NA)), "non-NA and positive")
})

test_that("weights for both traits warn, as coloc.abf() did", {
    ## coloc.abf() warned that the two traits' weights should come from
    ## independent sources; coloc_prior() now gives that warning, so coloc() and
    ## every old entry point do (AL15 in test-v603.R)
    expect_warning(coloc_prior(weights1 = c(1, 2), weights2 = c(1, 2)), "both traits")
    expect_silent(coloc_prior(weights1 = c(1, 2)))
})

test_that("coloc_control() holds the defaults", {
    ctl <- coloc_control()
    expect_s3_class(ctl, "coloc_control")
    expect_equal(ctl$overlap.min, 0.5)
    expect_true(ctl$trim_by_posterior)
})

test_that("coloc_control() takes overlap.min in (0, 1]", {
    expect_equal(coloc_control(overlap.min = 1)$overlap.min, 1)
    expect_error(coloc_control(overlap.min = 0))
    expect_error(coloc_control(overlap.min = 1.5))
    expect_error(coloc_control(overlap.min = NA_real_))
    expect_error(coloc_control(overlap.min = c(0.5, 0.6)))
    expect_error(coloc_control(overlap.min = "0.5"))
})

test_that("coloc_control() takes TRUE or FALSE for trim_by_posterior", {
    expect_false(coloc_control(trim_by_posterior = FALSE)$trim_by_posterior)
    expect_error(coloc_control(trim_by_posterior = NA))
    expect_error(coloc_control(trim_by_posterior = "yes"))
})

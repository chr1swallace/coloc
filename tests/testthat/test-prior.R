## coloc_prior() and coloc_control(), and the checks behind them

test_that("p12 may sit exactly on its bounds", {
    ## 0.2 * 0.2 rounds to 0.04000000000000001, just above a p12 of 0.04
    expect_s3_class(coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.04), "coloc_prior")
    expect_s3_class(coloc_prior(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4), "coloc_prior")
    expect_error(coloc_prior(p1 = 0.2, p2 = 0.2, p12 = 0.039), "p12 must lie")
    expect_error(coloc_prior(p12 = 2e-4), "p12 must lie")
})

test_that("coloc_control() takes TRUE or FALSE for trim_by_posterior", {
    expect_false(coloc_control(trim_by_posterior = FALSE)$trim_by_posterior)
    expect_error(coloc_control(trim_by_posterior = NA))
    expect_error(coloc_control(trim_by_posterior = "yes"))
})


data(coloc_test_data)
attach(coloc_test_data)

# finemap.abf()

test_that("Prior weight lenght consistency is checked (finemap.abf())", {
  expect_error(finemap.abf(D1, prior_weights = 1))
})

# coloc.abf()

test_that("Prior weights are saved in output (coloc.abf())", {

  expect_warning(coloc_out <- coloc.abf(D1, D2, prior_weights1 = rep(1, 500), prior_weights2 = rep(1, 500)))
  expect_named(coloc_out, c("summary", "results", "priors", "weights"))
  expect_equal(coloc_out$weights$prior_weights1, rep(1, 500))
  expect_equal(coloc_out$weights$prior_weights2, rep(1, 500))

})

## The (coloc()) twins give the weights through coloc_prior(), and read them
## back from $weights$w1 and $w2. coloc_prior() gives coloc.abf()'s warning
## when both are given, once it has checked them.

test_that("Prior weights are saved in output (coloc.abf()) (coloc())", {

  expect_warning(coloc_out <- coloc(D1, D2, prior = coloc_prior(weights1 = rep(1, 500),
                                                                weights2 = rep(1, 500))))
  expect_named(coloc_out, c("summary", "results", "priors", "status", "messages", "meta",
                            "weights", "bf"))
  expect_equal(coloc_out$weights$w1, rep(1, 500))
  expect_equal(coloc_out$weights$w2, rep(1, 500))

})

test_that("Prior weight length consistency is checked (coloc.abf())", {

  expect_error(
    coloc.abf(D1, D2, prior_weights1 = 1, prior_weights2 = rep(1, 500))
  )
  expect_error(
    coloc.abf(D1, D2, prior_weights1 = rep(1, 500), prior_weights2 = 1)
  )

})

test_that("Prior weight length consistency is checked (coloc.abf()) (coloc())", {

  ## the length is checked once the snps are known, after coloc_prior() warns
  expect_warning(prior <- coloc_prior(weights1 = 1, weights2 = rep(1, 500)), "both traits")
  expect_error(
    coloc(D1, D2, prior = prior)
  )
  expect_warning(prior <- coloc_prior(weights1 = rep(1, 500), weights2 = 1), "both traits")
  expect_error(
    coloc(D1, D2, prior = prior)
  )

})

test_that("Negative prior weights cause errors (coloc.abf())", {

  expect_error(
    coloc.abf(D1, D2, prior_weights1 = c(-1, rep(1, 499)), prior_weights2 = rep(1, 500))
  )
  expect_error(
    coloc.abf(D1, D2, prior_weights1 = rep(1, 500), prior_weights2 = c(-1, rep(1, 499)))
  )

})

test_that("Negative prior weights cause errors (coloc.abf()) (coloc())", {

  expect_error(
    coloc(D1, D2, prior = coloc_prior(weights1 = c(-1, rep(1, 499)), weights2 = rep(1, 500)))
  )
  expect_error(
    coloc(D1, D2, prior = coloc_prior(weights1 = rep(1, 500), weights2 = c(-1, rep(1, 499))))
  )

})

test_that("NA prior weights cause errors (coloc.abf())", {
  expect_error(
    coloc.abf(D1, D2, prior_weights1 = rep(NA, 500))
  )
})

test_that("NA prior weights cause errors (coloc.abf()) (coloc())", {
  expect_error(
    coloc(D1, D2, prior = coloc_prior(weights1 = rep(NA, 500)))
  )
})

susie_D1 <- runsusie(D1)
susie_D2 <- runsusie(D2)

# coloc.bf_bf()

test_that("Prior weight length consistency is checked (coloc.susie())", {

  expect_error(
    coloc.susie(susie_D1, susie_D2, prior_weights1 = 1, prior_weights2 = rep(1, 500))
  )
  expect_error(
    coloc.susie(susie_D1, susie_D2, prior_weights1 = rep(1, 500), prior_weights2 = 1)
  )

})

test_that("Prior weight length consistency is checked (coloc.susie()) (coloc())", {

  expect_warning(prior <- coloc_prior(weights1 = 1, weights2 = rep(1, 500)), "both traits")
  expect_error(
    coloc(susie_D1, susie_D2, prior = prior)
  )
  expect_warning(prior <- coloc_prior(weights1 = rep(1, 500), weights2 = 1), "both traits")
  expect_error(
    coloc(susie_D1, susie_D2, prior = prior)
  )

})

test_that("NA prior weights cause errors (coloc.susie())", {
  expect_error(
    coloc.susie(susie_D1, susie_D2, prior_weights1 = rep(NA, 500))
  )
})

test_that("NA prior weights cause errors (coloc.susie()) (coloc())", {
  expect_error(
    coloc(susie_D1, susie_D2, prior = coloc_prior(weights1 = rep(NA, 500)))
  )
})

test_that("Negative prior weights cause errors (coloc.susie())", {

  expect_error(
    coloc.susie(susie_D1, susie_D2, prior_weights1 = c(-1, rep(1, 499)))
  )

})

test_that("Negative prior weights cause errors (coloc.susie()) (coloc())", {

  expect_error(
    coloc(susie_D1, susie_D2, prior = coloc_prior(weights1 = c(-1, rep(1, 499))))
  )

})

# combine_abf_weighted()

test_that("combine_abf_weighted() errors appropriately", {
  expect_error(combine_abf_weighted(1:10, 1:5, 0.1, 0.1, 0.1, 1:10, 1:10))
  expect_error(combine_abf_weighted(1:10, 1:10, 0.1, 0.1, 0.1, 1:10, 1:5))
  expect_error(combine_abf_weighted(1:5, 1:5, 0.1, 0.1, 0.1, 1:10, 1:10))
})

test_that("combine_abf_weighted() and combine.abf() are consistent", {

  l1 <- rnorm(10)
  l2 <- rnorm(10)

  expect_equal(
    combine.abf(l1, l2, 1e-4, 1e-4, 1e-5),
    combine_abf_weighted(l1, l2, 1e-4, 1e-4, 1e-5, NULL, NULL)
  )
})

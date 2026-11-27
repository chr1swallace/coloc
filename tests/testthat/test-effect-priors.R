test_that("coloc.abf and process.dataset forward effect_priors to approx.bf.estimates", {
  data(coloc_test_data)
  d1 <- coloc_test_data$D1
  d2 <- coloc_test_data$D2
  custom_priors <- c(quant = 0.5, cc = 0.5)

  pd_default <- process.dataset(d1, suffix = "df1")
  pd_custom <- process.dataset(d1, suffix = "df1", effect_priors = custom_priors)
  expected_abf <- coloc:::approx.bf.estimates(
    z = d1$beta / sqrt(d1$varbeta),
    V = d1$varbeta,
    type = d1$type,
    suffix = "df1",
    sdY = d1$sdY,
    effect_priors = custom_priors
  )
  expect_false(isTRUE(all.equal(pd_default$lABF.df1, pd_custom$lABF.df1)))
  expect_equal(pd_custom$lABF.df1, expected_abf$lABF.df1)

  res_default <- suppressMessages(coloc.abf(d1, d2))
  res_custom <- suppressMessages(coloc.abf(d1, d2, effect_priors = custom_priors))
  expect_false(isTRUE(all.equal(res_default$results$lABF.df1, res_custom$results$lABF.df1)))
  expect_false(isTRUE(all.equal(res_default$results$lABF.df2, res_custom$results$lABF.df2)))
  expect_equal(
    res_custom$results$lABF.df1,
    expected_abf$lABF.df1[match(res_custom$results$snp, d1$snp)]
  )
})

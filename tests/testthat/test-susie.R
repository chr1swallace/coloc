library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

test_that("passing null_weight works", {
  result.nonull=runsusie(D1,nref=500,null_weight=1e-64)
  expect_equal(result.nonull$pi[ length(D1$beta)+1 ], 1e-64)
})

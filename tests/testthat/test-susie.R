library(coloc)
data(coloc_test_data)
attach(coloc_test_data)
  result.nonull=runsusie(D1)
  pi0=1e-4
  result.null=runsusie(D1,null_weight=pi0)

test_that("passing null_weight works", {
  expect_equal(result.null$pi[ length(D1$beta)+1 ], pi0)
  expect_equal(ncol(result.null$lbf_variable), ncol(result.nonull$lbf_variable)+1)
  expect_equal(length(result.null$pip), length(result.nonull$pip))
})

test_that("adding dimnames", {
  expect_equal(names(result.null$pip), D1$snp)
  expect_equal(colnames(result.null$lbf_variable), c(D1$snp,"null"))
  expect_equal(colnames(result.nonull$lbf_variable), c(D1$snp))
  expect_equal(colnames(result.null$alpha), c(D1$snp,"null"))
  expect_equal(colnames(result.nonull$alpha), c(D1$snp))
})

## p now removed
## test_that("runsusie argument p is deprecated", {
##   expect_warning(runsusie(D1,p=pi0), "Argument p is deprecated")
## })

test_that("null_weight=0 and null_weight=NULL are equivalent", {
  # Confirm that upstream susieR behavior hasn't changed
  expect_equal(runsusie(D1,null_weight=NULL), runsusie(D1,null_weight=0))
})

library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

test_that("missing required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  expect_error(check_dataset(D1, req = "test"), "missing required element test")
})

test_that("LD matrix must have dimnames", {
  expect_null(check_ld(D3, D3$LD))
  ld_no_dimnames <- D3$LD
  attr(ld_no_dimnames, "dimnames") <- NULL
  expect_error(check_ld(D3, ld_no_dimnames), "LD required to have row and column names")
})

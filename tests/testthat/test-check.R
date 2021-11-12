library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

test_that("missing required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  expect_error(check_dataset(D1, req = "test"), "missing required element test")
})

## The coloc_result: its contract, and the functions that read it

test_that("pp() stops on a row that does not exist", {
    bf <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    res <- coloc(bf, bf)
    expect_equal(pp(res, row = 1), pp(res))
    expect_error(pp(res, row = 2), "between 1 and 1")
    expect_error(pp(res, row = 0), "between 1 and 1")
})

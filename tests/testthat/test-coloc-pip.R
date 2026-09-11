library(testthat)
context("coloc_pip and related functions")

## ---- is_prob ----

test_that("is_prob accepts valid probabilities", {
  expect_true(is_prob(0.5))
  expect_true(is_prob(0.001))
  expect_true(is_prob(0.999))
})

test_that("is_prob rejects boundary values in strict mode", {
  expect_false(is_prob(0))
  expect_false(is_prob(1))
})

test_that("is_prob accepts boundary values in non-strict mode", {
  expect_true(is_prob(0, strict = FALSE))
  expect_true(is_prob(1, strict = FALSE))
})

test_that("is_prob rejects invalid inputs", {
  expect_false(is_prob(-0.1))
  expect_false(is_prob(1.1))
  expect_false(is_prob("a"))
  expect_false(is_prob(c(0.1, 0.2))) # length > 1
  expect_false(is_prob(NA_real_))
})

## ---- check_pip ----

test_that("check_pip passes for valid credible set", {
  cs <- data.frame(snp = c("rs1", "rs2"), pip = c(0.6, 0.4))
  expect_null(check_pip(cs))
})

test_that("check_pip errors when snp or pip missing", {
  expect_error(check_pip(data.frame(snp = "rs1")), "must include snp and pip")
  expect_error(check_pip(data.frame(pip = 0.5)), "must include snp and pip")
})

test_that("check_pip errors for wrong types", {
  expect_error(check_pip(data.frame(snp = 1, pip = 0.5)), "snp must be a character")
  expect_error(check_pip(data.frame(snp = "rs1", pip = "a")), "pip must be a numeric")
})

test_that("check_pip errors for mismatched lengths", {
  cs <- list(snp = c("rs1", "rs2"), pip = 0.5)
  expect_error(check_pip(cs), "equal length")
})

## ---- clpp_to_coloc ----

test_that("clpp_to_coloc returns named vector with H3 and H4", {
  res <- clpp_to_coloc(clpp = 0.5, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, alpha1 = 0.95, alpha2 = 0.95)
  expect_named(res, c("PP.H3", "PP.H4"))
  expect_true(all(res >= 0 & res <= 1))
  expect_equal(sum(res), 1)
})

test_that("clpp_to_coloc h4_only returns scalar", {
  res <- clpp_to_coloc(clpp = 0.5, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, alpha1 = 0.95, alpha2 = 0.95, h4_only = TRUE)
  expect_length(res, 1)
  expect_null(names(res))
})

test_that("clpp_to_coloc: clpp=0 gives PP.H4=0", {
  res <- clpp_to_coloc(clpp = 0, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, alpha1 = 0.95, alpha2 = 0.95)
  expect_equal(res[["PP.H4"]], 0)
  expect_equal(res[["PP.H3"]], 1)
})

test_that("clpp_to_coloc validates probability inputs", {
  expect_error(clpp_to_coloc(clpp = 0.5, p1 = -1, p2 = 1e-4, p12 = 1e-5, alpha1 = 0.95, alpha2 = 0.95), "p1, p2, p12")
  expect_error(clpp_to_coloc(clpp = 2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5, alpha1 = 0.95, alpha2 = 0.95), "clpp")
})

## ---- coloc_pip ----

test_that("coloc_pip returns correct structure", {
  cs1 <- data.frame(snp = c("rs1", "rs2", "rs3"), pip = c(0.5, 0.3, 0.2))
  cs2 <- data.frame(snp = c("rs1", "rs2", "rs4"), pip = c(0.6, 0.3, 0.1))
  res <- coloc_pip(cs1, cs2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  expect_named(res, c("nsnps", "PP.H3", "PP.H4"))
})

test_that("coloc_pip with no overlap gives PP.H4=0", {
  cs1 <- data.frame(snp = c("rs1", "rs2"), pip = c(0.6, 0.4))
  cs2 <- data.frame(snp = c("rs3", "rs4"), pip = c(0.7, 0.3))
  res <- coloc_pip(cs1, cs2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  expect_equal(res[["nsnps"]], 0)
  expect_equal(res[["PP.H4"]], 0)
  expect_equal(res[["PP.H3"]], 1)
})

test_that("coloc_pip with full overlap and high pip gives high PP.H4", {
  cs1 <- data.frame(snp = "rs1", pip = 0.99)
  cs2 <- data.frame(snp = "rs1", pip = 0.99)
  res <- coloc_pip(cs1, cs2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  expect_equal(res[["nsnps"]], 1)
  expect_gt(res[["PP.H4"]], 0.5) # should be high when single shared SNP dominates
})

test_that("coloc_pip PP.H3 + PP.H4 sums to 1", {
  cs1 <- data.frame(snp = c("rs1", "rs2"), pip = c(0.6, 0.4))
  cs2 <- data.frame(snp = c("rs1", "rs3"), pip = c(0.7, 0.3))
  res <- coloc_pip(cs1, cs2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
  expect_equal(res[["PP.H3"]] + res[["PP.H4"]], 1)
})

test_that("coloc_pip validates probability inputs", {
  cs1 <- data.frame(snp = "rs1", pip = 0.9)
  cs2 <- data.frame(snp = "rs1", pip = 0.9)
  expect_error(coloc_pip(cs1, cs2, p1 = -1, p2 = 1e-4, p12 = 1e-5), "p1, p2, p12")
})

test_that("coloc_pip validates credible set format", {
  bad_cs <- data.frame(id = "rs1", prob = 0.5)
  good_cs <- data.frame(snp = "rs1", pip = 0.5)
  expect_error(coloc_pip(bad_cs, good_cs, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5))
})

## library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

## D3 should have two signals
s <- runsusie(D3)
a <- s$alpha[ s$sets$cs_index[1], s$sets$cs[[1]] ]
cs.D3a <- data.frame(snp=names(a),PP=unname(a))
b <- s$alpha[ s$sets$cs_index[2], s$sets$cs[[2]] ]
cs.D3b <- data.frame(snp=names(b),PP=unname(b))

## others have one signal
cs.D1 <- finemap.abf(D1) |> credible.sets()
cs.D2 <- finemap.abf(D2) |> credible.sets()
cs.D4 <- finemap.abf(D4) |> credible.sets()

## everything should colocalise
testthat("colocalises when it should", {
    coloc_pip(cs.D1,cs.D3a)["PP.H4"] > 0.99
    coloc_pip(cs.D2,cs.D3a)["PP.H4"] > 0.99
    coloc_pip(cs.D4,cs.D3a)["PP.H4"] > 0.99
})
testthat("doesn't colocalise when it shouldn't", {
    coloc_pip(cs.D1,cs.D3b)["PP.H4"] < 0.1
    coloc_pip(cs.D2,cs.D3b)["PP.H4"] < 0.1
    coloc_pip(cs.D4,cs.D3b)["PP.H4"] < 0.1
})
detach(coloc_test_data)

library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

B1=C1=D1
C1[c("type","s")]=list(type="cc",s=0.5)
B1[c("type","s")]=list(type="cc",s=2)

test_that("missing or bad required elements throws error", {
  expect_null(check_dataset(D1, req = names(D1)))
  expect_error(check_dataset(D1, req = "test"))
  ## type="cc" and s
  expect_null(check_dataset(C1))
  expect_error(check_dataset(B1))
  expect_that(check_dataset(list(), ""), throws_error())
  expect_that(check_dataset(list(beta=1,p=2,type="blah"), ""), throws_error())
})

test_that("LD matrix must have dimnames", {
  expect_null(check_ld(D3, D3$LD))
  ld_no_dimnames <- D3$LD
  attr(ld_no_dimnames, "dimnames") <- NULL
  expect_error(check_ld(D3, ld_no_dimnames))
})

test_that("issue 79", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="quant",
          MAF=rep(0.5,5))
  expect_error(coloc.abf(d1,d2), "must give sdY")
})

test_that("issue 160", {
  d1=list(snp=letters[1:5],
          position=1:5,
          N=200000,
          MAF=runif(5)/2,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          type="cc")
  d2=list(snp=letters[1:5],
          position=1:5,
          beta=rnorm(5),
          varbeta=rep(0.01,5),
          N=rep(150000,5),
          type="quant",
          MAF=rep(0.5,5))
  expect_error(coloc.abf(d1,d2), "per-snp sample sizes")
})

test_that("Infinite values in beta/varbeta triggers specific error", {
  D1_beta_inf <- D1
  D1_beta_inf$beta[[1]] <- Inf
  expect_error(check_dataset(D1_beta_inf), "Infinite")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- Inf
  expect_error(check_dataset(D1_varbeta_inf), "Infinite")

  D1_varbeta_inf <- D1
  D1_varbeta_inf$varbeta[[1]] <- 0
  expect_error(check_dataset(D1_varbeta_inf), "Zero")
})

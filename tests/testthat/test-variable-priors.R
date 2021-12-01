if(interactive()) {
  devtools::load_all("~/RP/coloc")
} else {
  library(coloc)
}

data(coloc_test_data)
attach(coloc_test_data)
expect_not_equal=function(x,y)
  expect_false(isTRUE(all.equal(x, y)))

dataset1=D1; dataset2=D2; MAF=NULL; p1=1e-4;p2=1e-4;p12=1e-5

test_that("adjust_prior doesn't mangle valid priors", {
  expect_equal(p1,adjust_prior(p1,50))
  expect_equal(rep(p1,50),adjust_prior(rep(p1,50),50))
  ## these are greater than 1/nsnp, should give warning and fix
  expect_warning(adjust_prior(0.1,50))
  expect_equal(1/11,suppressWarnings(adjust_prior(0.2,10)))
  expect_warning(adjust_prior(rep(0.1,50),50))
  expect_equal(rep(1/11,10),suppressWarnings(adjust_prior(rep(0.2,10),10)))
})

test_that("coloc.abf copes with variable priors", {
  result.const=coloc.abf(D1,D2)
  result.repped=coloc.abf(D1,D2,p1=rep(1e-4,length(D1$snp)),
                          p2=rep(1e-4,length(D1$snp)),
                          p12=rep(1e-5,length(D1$snp)))
  expect_equal(result.const$summary,result.repped$summary)
  result.varp=coloc.abf(D1,D2,p1=seq(1e-3,1e-4,length.out=length(D1$snp)))
  expect_not_equal(result.const$summary, result.varp$summary)
})

test_that("coloc.susie doesn't bail with variable priors", {
  sargs=susie.args=list(check_prior=FALSE)
  S3=do.call("runsusie",c(list(D3),sargs))
  S4=do.call("runsusie",c(list(D4),sargs))
  result.susie.data=coloc.susie(D3,D4,susie.args=sargs)
  result.susie=coloc.susie(S3,S4)
  expect_equal(result.susie.data$summary,result.susie$summary)
  result.susie.repped=coloc.susie(S3,S4,p2=rep(1e-4,length(D3$snp)))
  expect_equal(result.susie.repped$summary,result.susie$summary)
  result.susie.varp=coloc.susie(S3,S4,p2=seq(1e-4,1e-3,length.out=length(D3$snp)))
  expect_not_equal(result.susie.data$summary, result.susie.varp$summary)
})

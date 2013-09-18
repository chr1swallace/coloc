library(snpStats)
data(testdata)
X <- Autosomes[,101:201]
cs <- col.summary(X)
X <- X[, cs[,"MAF"]>0.05 & cs[,"Call.rate"]>0.9]
Y<-rnorm(nrow(X),mean=as(X[,8],"numeric")) + rnorm(nrow(X),sd=2)
eff<-snp.rhs.estimates(Y ~ 1, snp.data=X, family="gaussian")
vbeta <- sapply(eff@.Data, "[[", "Var.beta")
maf <- col.summary(X)[,"MAF"]
sd.est <- sdY.est(vbeta=vbeta, maf=maf, n=nrow(X))
test_that("sdY.est", {
  expect_that(abs(sd.est - sd(Y)) < 1, is_true())
})

test_that("process.dataset", {
  expect_that(process.dataset(list(), ""), throws_error())
  expect_that(process.dataset(list(beta=1,p=2,type="blah"), ""), throws_error())
})
          

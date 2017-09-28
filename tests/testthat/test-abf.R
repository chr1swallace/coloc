library(snpStats)
data(testdata)
X <- Autosomes[,101:201]
cs <- col.summary(X)
X <- X[, cs[,"MAF"]>0.05 & cs[,"Call.rate"]>0.9]
maf <- col.summary(X)[,"MAF"]
set.seed(42)

## quantitative trait
Y<-rnorm(nrow(X),mean=as(X[,8],"numeric"),sd=4) + rnorm(nrow(X),sd=2)
eff<-snp.rhs.estimates(Y ~ 1, snp.data=X, family="gaussian")
beta.q <- sapply(eff@.Data, "[[", "beta")
vbeta.q <- sapply(eff@.Data, "[[", "Var.beta")
p.q <- pchisq(beta.q^2/vbeta.q,df=1,lower.tail=FALSE)
sd.est <- coloc:::sdY.est(vbeta=vbeta.q, maf=maf, n=nrow(X))

## case-control trait
cc <- rbinom(nrow(X),1, p=(1+as(X[,8],"numeric"))/4)
eff<-snp.rhs.estimates(cc ~ 1, snp.data=X, family="binomial")
beta.cc <- sapply(eff@.Data, "[[", "beta")
vbeta.cc <- sapply(eff@.Data, "[[", "Var.beta")
p.cc <- pchisq(beta.cc^2/vbeta.cc,df=1,lower.tail=FALSE)

## general things
test_that("sdY.est", {
  expect_that(abs(sd.est - sd(Y)) < 0.1, is_true())
})
DQ <- list(beta=beta.q,
           varbeta=vbeta.q,
           type="quant",
           snp=colnames(X),
           sdY=sd.est,
           N=nrow(X))

DCC <- list(beta=beta.cc,
            varbeta=vbeta.cc,
            type="cc",
            snp=colnames(X),
            MAF=maf,
            s=mean(cc),
            N=nrow(X))
PQ <- list(pvalues=p.q,
           type="quant",
           MAF=maf,
           snp=colnames(X),
           sdY=sd.est,
           N=nrow(X))
PCC <- list(pvalues=p.cc,
            type="cc",
            MAF=maf,
            snp=colnames(X),
            s=mean(cc),
            N=nrow(X))
PCC.bad <- list(pvalues=p.cc,
            type="cc",
            MAF=maf,
            snp=colnames(X),
            s=sum(cc),
            N=nrow(X))


RESULTS <- list(dd = coloc.abf(dataset1=DQ,dataset2=DCC),
                dp = coloc.abf(dataset1=DQ,dataset2=PCC),
                pd = coloc.abf(dataset1=PQ,dataset2=DCC),
                pp = coloc.abf(dataset1=PQ,dataset2=PCC))
lapply(RESULTS,"[[","summary")

test_that("process.dataset", {
  expect_that(process.dataset(list(), ""), throws_error())
  expect_that(process.dataset(list(beta=1,p=2,type="blah"), ""), throws_error())
  expect_error(process.dataset(DQ,suffix=".q"), NA)
  expect_error(process.dataset(DCC,suffix=".q"), NA)
  expect_error(process.dataset(PQ,suffix=".q"), NA)
  expect_error(process.dataset(PCC,suffix=".q"), NA)
  expect_error(process.dataset(PCC.bad,suffix=".q"))
  pd.cc <- process.dataset(DCC,suffix=".cc")
  pd.q <- process.dataset(DQ,suffix=".q")
  expect_is(pd.q,"data.frame")
  expect_is(pd.cc,"data.frame")
  expect_equal(nrow(pd.q),ncol(X))
  expect_equal(nrow(pd.cc),ncol(X))
})

## coloc.abf with coefficients
test_that("coloc.abf", {
    expect_error(coloc.abf(dataset1=DQ,dataset2=DCC), NA)
    result <- coloc.abf(dataset1=DQ,dataset2=DCC)
    expect_true(which.max(result$summary[-1]) == 5)
    expect_true(result$summary[1] == ncol(X))
})


## alternative test data
## colocdata<- read.table("inst/tests/test.txt", sep="\t", header=T)
## N <- 18124
## result <- coloc.abf(dataset1=list(beta=colocdata$beta.dataset1,
##             varbeta=colocdata$varbeta.dataset1,
##             type="quant",
##             snp=colocdata$SNP,
##             N=N),
##           dataset2=list(beta=colocdata$beta.dataset1,
##             varbeta=colocdata$varbeta.dataset1,
##             type="quant",
##             snp=colocdata$SNP,
##             N=N),
##           MAF=colocdata$MAF)



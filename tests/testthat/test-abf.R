## library(snpStats)
## data(testdata)
set.seed(42)
X <- do.call("cbind", lapply(runif(100,0.05,0.5), function(maf) rbinom(400,2,maf)))
colnames(X)=paste0("S",1:ncol(X))
maf <- colMeans(X)/2

## quantitative trait
Y<-rnorm(nrow(X),mean=as(X[,1],"numeric"),sd=4) + rnorm(nrow(X),sd=2)
m=sapply(1:ncol(X), function(i) {
    mod=lm(Y ~ X[,i])
    c(coefficients(mod)[-1], vcov(mod)[-1,-1])
})
beta.q <- m[1,]
vbeta.q <- m[2,]
p.q <- pchisq(beta.q^2/vbeta.q,df=1,lower.tail=FALSE)
sd.est <- suppressWarnings(coloc:::sdY.est(vbeta=vbeta.q, maf=maf, n=nrow(X)))

## case-control trait
cc <- rbinom(nrow(X),1, p=(1+as(X[,1],"numeric"))/4)
m=sapply(1:ncol(X), function(i) {
    mod=glm(cc ~ X[,i], family="binomial")
    c(coefficients(mod)[-1], vcov(mod)[-1,-1])
})
beta.cc <- m[1,]
vbeta.cc <- m[2,]
p.cc <- pchisq(beta.cc^2/vbeta.cc,df=1,lower.tail=FALSE)

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

## general things
test_that("sdY.est", {
  expect_true(abs(sd.est - sd(Y)) < 0.1)
})

RESULTS <- suppressWarnings(list(dd = coloc.abf(dataset1=DQ,dataset2=DCC),
                                 dp = coloc.abf(dataset1=DQ,dataset2=PCC),
                                 pd = coloc.abf(dataset1=PQ,dataset2=DCC),
                                 pp = coloc.abf(dataset1=PQ,dataset2=PCC)))
lapply(RESULTS,"[[","summary")

test_that("process.dataset", {
  expect_error(process.dataset(DQ,suffix=".q"), NA)
  expect_error(process.dataset(DCC,suffix=".q"), NA)
  expect_error(process.dataset(PQ,suffix=".q"), NA)
  expect_error(process.dataset(PCC,suffix=".q"), NA)
  pd.cc <- process.dataset(DCC,suffix=".cc")
  pd.q <- process.dataset(DQ,suffix=".q")
  expect_is(pd.q,"data.frame")
  expect_is(pd.cc,"data.frame")
  expect_equal(nrow(pd.q),ncol(X))
  expect_equal(nrow(pd.cc),ncol(X))
})

## coloc.abf with coefficients
test_that("coloc.abf", {
    expect_error(suppressWarnings(coloc.abf(dataset1=DQ,dataset2=DCC), NA))
    expect_warning(result <- coloc.abf(dataset1=DQ,dataset2=DCC), "minimum p value")
    expect_true(which.max(result$summary[-1]) == 5)
    expect_true(result$summary[1] == ncol(X))
})



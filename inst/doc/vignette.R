### R code from vignette source 'vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: vignette.Rnw:82-119
###################################################
setClass("simdata",
         representation(df1="data.frame",df2="data.frame"))
setValidity("simdata", function(object) {
  n <- nrow(object@df1)
  if(nrow(object@df2)!=n)
    return("nrow of '@df1' should equal nrow of '@df2'")
})
setMethod("show", signature="simdata", function(object) {
  cat("pair of simulated datasets, with",ncol(object@df1)-1,"SNPs and",nrow(object@df1),"samples.\n")
})

sim.data <- function(nsnps=10,nsamples=200,causals=1:2,nsim=1) {
  cat("Generate",nsim,"small sets of data\n")
  ntotal <- nsnps * nsamples * nsim
  X1 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y1 <- rnorm(nsamples,rowSums(X1[,causals]),2)
  X2 <- matrix(rbinom(ntotal,1,0.4)+rbinom(ntotal,1,0.4),ncol=nsnps)
  Y2 <- rnorm(nsamples,rowSums(X2[,causals]),2)
  colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
  df1 <- cbind(Y=Y1,X1)
  df2 <- cbind(Y=Y2,X2)
  if(nsim==1) {
    return(new("simdata",
               df1=as.data.frame(df1),
               df2=as.data.frame(df2)))
  } else {
    index <- split(1:(nsamples * nsim), rep(1:nsim, nsamples))
    objects <- lapply(index, function(i) new("simdata", df1=as.data.frame(df1[i,]),
                                             df2=as.data.frame(df2[i,])))
    return(objects)
  }
}

## simulate some data and load the coloc library
set.seed(12345)
data <- sim.data(nsamples=1000,nsim=1)
library(coloc)


###################################################
### code chunk number 2: vignette.Rnw:133-138
###################################################
## run a coloc with pcs
pcs <- pcs.prepare(data@df1[,-1], data@df2[,-1])
pcs.1 <- pcs.model(pcs, group=1, Y=data@df1[,1], threshold=0.8)
pcs.2 <- pcs.model(pcs, group=2, Y=data@df2[,1], threshold=0.8)
ct.pcs <- coloc.test(pcs.1,pcs.2)


###################################################
### code chunk number 3: vignette.Rnw:148-150
###################################################
ct.pcs
str(summary(ct.pcs))


###################################################
### code chunk number 4: vignette.Rnw:164-166
###################################################
ct.pcs.bayes <- coloc.test(pcs.1,pcs.2, bayes=TRUE)
ci(ct.pcs.bayes)


###################################################
### code chunk number 5: vignette.Rnw:174-182
###################################################
ct.bma <- coloc.bma(data@df1, data@df2, 
                    family1="gaussian", family2="gaussian",
                    plot.coeff=TRUE)
ct.bma.bayes <- coloc.bma(data@df1, data@df2, 
                          family1="gaussian", family2="gaussian", 
                          bayes=TRUE)
ct.bma
ci(ct.bma.bayes)


###################################################
### code chunk number 6: vignette.Rnw:224-233
###################################################
## compare individual values of eta
ct.pcs <- coloc.test(pcs.1,pcs.2, bayes.factor=c(-1,0,1))
bf(ct.pcs)

## compare ranges of eta
ct.bma <- coloc.bma(data@df1, data@df2, 
                    family1="gaussian", family2="gaussian",
                    bayes.factor=list(c(-0.1,1), c(0.9,1.1)))
bf(ct.bma)


###################################################
### code chunk number 7: vignette.Rnw:283-294
###################################################
library(snpStats)

Y1 <- data@df1$Y
Y2 <- data@df2$Y

X1 <- new("SnpMatrix",as.matrix(data@df1[,-1]))
X2 <- new("SnpMatrix",as.matrix(data@df2[,-1]))

p1 <- snpStats::p.value(single.snp.tests(phenotype=Y1, snp.data=X1),df=1)
p2 <- snpStats::p.value(single.snp.tests(phenotype=Y2, snp.data=X2),df=1)
maf <- col.summary(X2)[,"MAF"]


###################################################
### code chunk number 8: vignette.Rnw:301-305
###################################################
my.res <- coloc.abf(dataset1=list(pvalues=p1,N=nrow(X1),type="quant"),
                    dataset2=list(pvalues=p2,N=nrow(X2),type="quant"),
                    MAF=maf)
print(my.res[[1]])


###################################################
### code chunk number 9: vignette.Rnw:313-315
###################################################
ct.abf <- coloc.abf.datasets(data@df1, data@df2, response1="Y", response2="Y",
                             type1="quant", type2="quant")



source("demo/common.R")

## simulate some data
set.seed(12345)
data <- sim.data(nsim=1)

## load code
library(devtools)
load_all()

## run a coloc with pcs
pcs <- pcs.prepare(data@df1[,-1], data@df2[,-1])
pcs.1 <- pcs.model(pcs, group=1, Y=data@df1[,1])
pcs.2 <- pcs.model(pcs, group=2, Y=data@df2[,1])
ct.pcs <- coloc.test(pcs.1,pcs.2)

## run a coloc with bma
ct.bma <- coloc.bma(data@df1, data@df2, family1="gaussian", family2="gaussian")


cb <- coloc.bma(data@df1, data@df2, family1="gaussian", family2="gaussian", snps=paste("s",1:10,sep=""),
                bayes.factor=c(0,1,Inf))
cbr <- coloc.bma(data@df1, data@df2, family1="gaussian", family2="gaussian", snps=paste("s",1:10,sep=""),
                bayes.factor=list(c(-0.1,0.1), c(0.9,1.1)))

bf(cb)
bf(cbr)

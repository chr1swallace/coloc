source("demo/common.R")

## simulate some data
data <- sim.data(nsim=1)

## load code
library(devtools)
load_all()

## run a coloc
cb <- coloc.bma(data@df1, data@df2, family1="gaussian", family2="gaussian", snps=paste("s",1:10,sep=""),
                bayes.factor=c(0,1,Inf))
cbr <- coloc.bma(data@df1, data@df2, family1="gaussian", family2="gaussian", snps=paste("s",1:10,sep=""),
                bayes.factor=list(c(-0.1,0.1), c(0.9,1.1)))

bf(cb)
bf(cbr)

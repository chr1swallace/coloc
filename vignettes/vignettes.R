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

## simulate some data
set.seed(12345)
data <- sim.data(nsim=1)

## load code
## library(devtools); load_all()

## run a coloc with pcs
pcs <- pcs.prepare(data@df1[,-1], data@df2[,-1])
pcs.1 <- pcs.model(pcs, group=1, Y=data@df1[,1])
pcs.2 <- pcs.model(pcs, group=2, Y=data@df2[,1])
ct.pcs <- coloc.test(pcs.1,pcs.2, bayes.factor=c(0,1,Inf))

## run a coloc with bma
ct.bma <- coloc.bma(data@df1, data@df2, family1="gaussian", family2="gaussian", bayes.factor=c(0,1,Inf))

## run a fully Bayesian coloc
ct.bayes <- 

bf(cb)
bf(cbr)

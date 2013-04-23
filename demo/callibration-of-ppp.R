### ppp are conservative when there is model uncertainty
### we attempt to callibrate them, and generate cppp, by simulating from the posterior
### this script checks that
### to do this in a reasonable time, it makes use of the parallel library
### set option mc.cores to something sensible if you have multiple cores available

require(parallel)
options(mc.cores=35)

require(ggplot2)
theme_set(theme_bw())

nsim <- max(getOption("mc.cores",1), 10)

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

cat("Generate",nsim,"small sets of data\n")
sim.data <- function(nsnps=10,nsamples=200,causals=1:2,nsim=1) {
  ntotal <- nsnps * nsamples * nsim
  X1 <- matrix(rbinom(ntotal,1,0.4),ncol=nsnps)
  Y1 <- rnorm(nsamples,rowSums(X1[,causals]),2)
  X2 <- matrix(rbinom(ntotal,1,0.4),ncol=nsnps)
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

data <- sim.data(nsim=nsim)
data[[1]]

cat("Do coloc testing on each dataset without callibration\n")
runone <- function(d) {
  coloc.bma(d@df1, d@df2, snps=colnames(d@df1)[-1],
            family1="gaussian", family2="gaussian",quiet=TRUE,thr=0.1,callibrate=0)
}
coloc.results <- lapply(data,runone,mc.cores=getOption("mc.cores",1))

cat("With small datasets, we expect the ppp values to be non-uniform, and peaked around 0.5:\n")
ppp <- unlist(lapply(coloc.results, ppp.value))
ggplot(data.frame(ppp=unlist(lapply(coloc.results, ppp.value))), aes(x=ppp)) +
  geom_histogram(binwidth=0.05, col="grey90",fill="grey70") +
  xlim(0,1) +
  geom_vline(xintercept=0.5,lty=2,col="red")

cat("Now do coloc testing on each dataset *with* callibration\n")
cat("Here, we use a moderate number of callibrations, but with real data and smaller significance thresholds you might prefer to use more.  Callibration is quite fast, anyway.\n")
runone <- function(d, calstyle=1) {
  coloc.bma(d@df1, d@df2, snps=colnames(d@df1)[-1],
            family1="gaussian", family2="gaussian",callibrate=100,callibrate.style=calstyle,thr=0.1,
            quiet=TRUE)
}
coloc.1 <- mclapply(data,runone,mc.cores=getOption("mc.cores",1), calstyle=1)
coloc.2 <- mclapply(data,runone,mc.cores=getOption("mc.cores",1), calstyle=2)
coloc.3 <- mclapply(data,runone,mc.cores=getOption("mc.cores",1), calstyle=3)

cat("With small datasets, we expect the ppp values to be non-uniform, and peaked around 0.5:\n")
cppp.1 <- do.call("rbind",lapply(coloc.1, function(x) x@cppp))
cppp.2 <- do.call("rbind",lapply(coloc.2, function(x) x@cppp))
cppp.3 <- do.call("rbind",lapply(coloc.3, function(x) x@cppp))

library(reshape)
myplot <- function(cppp) {
  df <- melt(cppp)

  ggplot(df, aes(x=value)) +
    geom_histogram(aes(y=..density..), binwidth=0.1, col="grey90",fill="grey70") +
      xlim(0,1) +
        geom_vline(xintercept=0.5,lty=2,col="red") +
          facet_wrap(~X2)
}

myplot(cppp.1)
myplot(cppp.2)
myplot(cppp.3)


                        

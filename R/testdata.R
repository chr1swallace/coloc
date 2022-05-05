#' Simulated data to use in testing and vignettes in the coloc package
#'
#' @docType data
#'
#' @usage data(coloc_test_data)
#'
#' @format A four of two coloc-style datasets. Elements D1 and D2 have a single
#'   shared causal variant, and 50 SNPs. Elements D3 and D4 have 100 SNPs, one
#'   shared causal variant, and one variant unique to D3. Use these as examples
#'   of what a coloc-style dataset for a quantitative trait should look like.
#'
#' @keywords datasets
#' @examples
#' data(coloc_test_data)
#' names(coloc_test_data)
#' str(coloc_test_data$D1)
#' check_dataset(coloc_test_data$D1) # should return NULL if data structure is ok
"coloc_test_data"

if(FALSE) {

  library(mvtnorm)
  library(bindata)
  simx <- function(nsamples,S,maf) {
    rmvbin(n=nsamples, margprob=maf, sigma=S)
    ## mu <- rep(0,nsnps)
    ## rawvars <- rmvnorm(n=nsamples, mean=mu, sigma=S)
    ## pvars <- pnorm(rawvars)
    ## x <- qbinom(1-pvars, 2, maf)
  }

  sim.data <- function(nsnps=50,nsamples=200,ncausals=c(1,1),sd.y=c(1,1.2)) {
    cat("Generate a small set of data\n")
    ntotal <- nsnps * nsamples
    ## S <- toeplitz(sample(1:10,nsnps,replace=TRUE)/10)
    S <- toeplitz((nsnps:1)/nsnps)
    R=rWishart(1,2*nsnps,S)[,,1]
    ## (1 - (abs(outer(1:nsnps,1:nsnps,`-`))/nsnps))
    maf=runif(nsnps,0.2,0.8)^2
    causals=sample(c(1:nsnps)[abs(maf-0.5) < 0.3], max(ncausals))
    X1 <- simx(nsamples,S,maf)
    X2 <- simx(nsamples,S,maf)
    X3 <- simx(nsamples,S,maf)
    Y1 <- rnorm(nsamples,rowSums(X1[,causals[1:ncausals[1]],drop=FALSE]/2),sd.y[1])

    Y2 <- rnorm(nsamples,rowSums(X2[,causals[1:ncausals[2]],drop=FALSE]/2),sd.y[2])
    colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
    df1 <- cbind(Y=Y1,X1)
    df2 <- cbind(Y=Y2,X2)
    return(list(
      df1=as.data.frame(df1),
      df2=as.data.frame(df2),
      X3=X3,
      causals=causals))
  }

  sim.null <- function(nsnps=50,nsamples=200,ncausals=0,sd.y=c(1,1.2)) {
    cat("Generate a small set of data\n")
    ntotal <- nsnps * nsamples
    ## S <- toeplitz(sample(1:10,nsnps,replace=TRUE)/10)
    S <- toeplitz((nsnps:1)/nsnps)
    R=rWishart(1,2*nsnps,S)[,,1]
    ## (1 - (abs(outer(1:nsnps,1:nsnps,`-`))/nsnps))
    maf=runif(nsnps,0.2,0.8)^2
    X1 <- simx(nsamples,S,maf)
    X2 <- simx(nsamples,S,maf)
    X3 <- simx(nsamples,S,maf)
    Y1 <- rnorm(nsamples,sd.y[1])
    Y2 <- rnorm(nsamples,sd.y[2])
    colnames(X1) <- colnames(X2) <- paste("s",1:nsnps,sep="")
    df1 <- cbind(Y=Y1,X1)
    df2 <- cbind(Y=Y2,X2)
    return(list(#new("simdata",
      df1=as.data.frame(df1),
      df2=as.data.frame(df2),
      X3=X3))
  }


  get.beta <- function(x,nm) {
    beta <- sapply(x,"[",1)
    varbeta <- sapply(x, "[", 2)^2
    names(beta) <- names(varbeta) <- nm
    return(list(beta=beta,varbeta=varbeta))
  }

make_data=function(data) {
  Y1 <- data$df1$Y
  Y2 <- data$df2$Y
  X1 <- as.matrix(data$df1[,-1])
  X2 <- as.matrix(data$df2[,-1])
  tests1 <- lapply(1:ncol(X1), function(i) summary(lm(Y1 ~ X1[,i]))$coefficients[2,])
  tests2 <- lapply(1:ncol(X2), function(i) summary(lm(Y2 ~ X2[,i]))$coefficients[2,])

  snpnames=make.unique(colnames(X2))
  maf <- colMeans(data$X3)/2
  names(maf) <- snpnames
  LD <- cor(data$X3)
  nsnp=length(maf)
  dimnames(LD)=list(snpnames,snpnames)
  b1 <- get.beta(tests1,colnames(LD))
  b2 <- get.beta(tests2,colnames(LD))

return(list(D1=list(beta=b1$beta,
             varbeta=b1$varbeta,
             N=nrow(X1),
             sdY=sd(Y1),
             type="quant",
             MAF=maf,
             LD=LD,
             snp=names(b1$beta),
             position=1:nsnp),
  D2=list(beta=b2$beta,
             varbeta=b2$varbeta,
             N=nrow(X2),
             sdY=sd(Y2),
             type="quant",
             MAF=maf,
             LD=LD,
             snp=names(b2$beta),
          position=1:nsnp)))
}


  ## set.seed(46411)
  set.seed(42)
  data=sim.data(nsamples=1000,nsnps=500,sd.y=c(1.1,1),ncausals=c(1,1))
  data$causals # 105
D1=make_data(data)[[1]]
  ## D1$causals=data$causals[1]
  D2=make_data(data)[[2]]
  causals=list(D1=data$causals, D2=data$causals)
  ## D2$causals=data$causals[1]
  par(mfrow=c(2,1)); plot_dataset(D1); abline(v=data$causals); plot_dataset(D2); abline(v=data$causals)
  set.seed(42)
  data=sim.data(nsamples=1000,nsnps=500,sd.y=c(1.2,1.3),ncausals=c(2,1))
  data$causals # 84, 379
  D3=make_data(data)[[1]]
  ## D3$causals=data$causals[1:2]
  D4=make_data(data)[[2]]
  causals=c(causals,list(D3=data$causals, D4=data$causals[1]))
  ## D4$causals=data$causals[1]
  par(mfrow=c(2,1)); plot_dataset(D3); abline(v=data$causals); plot_dataset(D4); abline(v=data$causals)
  null=sim.null(nsamples=1000,nsnps=500,sd.y=c(1.1,1))
N1=make_data(null)[[1]]
N2=make_data(null)[[2]]
  par(mfrow=c(2,1)); plot_dataset(N1) ; plot_dataset(N2);

  par(mfrow=c(2,2))
  plot_dataset(D1,main="D1")
  plot_dataset(D2,main="D2")
  plot_dataset(D3,main="D3")
  plot_dataset(D4,main="D4")

  S1=runsusie(D1)
  summary(S1)
  S2=runsusie(D2)
  summary(S2)
  S3=runsusie(D3)
  summary(S3)
  S4=runsusie(D4)
  summary(S4)

  D5=D3
  D5$varbeta=D5$varbeta * 2
  D5$N=D5$N / 2
  summary(runsusie(D5)) # default coverage 0.95
  summary(runsusie(D5,coverage=0.1))  # lower coverage
  summary(runsusie(D5,coverage=0.01)) # even lower

  coloc_test_data=list(D1=D1,D2=D2,D3=D3,D4=D4,N1=N1,N2=N2,causals=causals)
  save(coloc_test_data, file="data/coloc_test_data.rda", compress="xz")

}

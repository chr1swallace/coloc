library(devtools)
load_all("~/RP/coloc")
load('/scratch/wallace/twas/example_coloc_data.RData')
#raj-cd14, t1dgc+wtccc
#df1 = expr & genotype
#df2 = cc, stratum, gwas genotype
dim(df1)
#[1]  211 3652
dim(df2)
#[1] 13254  3653

response1 <- "expr"
response2 <- "cc"
stratum2 <- "stratum"
stratum1 <- NULL
family1 <- "gaussian"
family2 <- "binomial"
nsnps <- 2
thr <- 0.05
r2.trim <- 0.5
snps <- intersect(setdiff(colnames(df1),c(response1,stratum1)),
               setdiff(colnames(df2),c(response2,stratum2)))
 snps <- unique(snps)
  n.orig <- length(setdiff(snps,c(response1,response2,stratum1,stratum2)))
  df1 <- as.data.frame(df1[,c(response1,stratum1,snps)])
  df2 <- as.data.frame(df2[,c(response2,stratum2,snps)])

  ## remove missings and tag
  prep <- prepare.df(df1, df2, drop.cols=c(response1,response2,stratum1,stratum2),
                     r2.trim=r2.trim, dataset=1, quiet=quiet)  
dim(df1)
dim(prep$df1)
snps <- prep$snps
  df1 <- prep$df1[,c(response1,stratum1,snps)]
  df2 <- prep$df2[,c(response2,stratum2,snps)]
  
  if(!quiet) {
    message("Dropped ",n.orig - length(snps)," of ",n.orig," SNPs due to LD: r2 > ",r2.trim,".")
    message(length(snps)," SNPs remain.")
  }

  ## remove any completely predictive SNPs
  f1 <- as.formula(paste(response1, "~", paste(snps,collapse="+")))
f2 <- as.formula(paste(response2, "~", paste(snps,collapse="+")))
## use a subsample of rows if large
n1 <- nrow(df1)
if(n1 > 1000) {
    s1 <- df1[sample(1:nrow(df1),1000),]
} else {
    s1 <- df1
}
n2 <- nrow(df1)
if(n2 > 1000) {
    s2 <- df2[sample(1:nrow(df2),1000),]
} else {
    s2 <- df2
}

  lm1 <- lm(f1, data=s1) #, family=family1)
  lm2 <- lm(f2, data=s2) #, family=family2)
  while(any(is.na(coefficients(lm1))) || any(is.na(coefficients(lm2)))) {
    drop <- which((is.na(coefficients(lm1)) | is.na(coefficients(lm2)))[-1])
    if(!quiet)
      cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n",drop,"\n")
    snps <- snps[-drop]
    f1 <- as.formula(paste(response1, "~", paste(snps,collapse="+")))
    f2 <- as.formula(paste(response2, "~", paste(snps,collapse="+")))
    lm1 <- lm(f1, data=df1)#, family=family1)
    lm2 <- lm(f2, data=df2)#, family=family2)
  }

  x1 <- df1[,snps]
  x2 <- df2[,snps]
  n.clean <- length(snps)
  ## step2, evaluate all pairs of marginally interesting single SNPs
combs <- lapply(nsnps, function(n) {
    cmb <- combn(1:n.clean,n)
    idx <- apply(sapply(1:nrow(cmb), function(i) { cmb[i,] %in% which(use) }),1,any)
    t(cmb[,idx])
  }) 
  models <- lapply(combs, function(X) {
    tmp <- matrix(0,nrow(X),n.clean)
    for(j in 1:ncol(X)) {
      tmp[cbind(1:nrow(X),X[,j])] <- 1
    }
    tmp
  })
  if(length(nsnps)>1) {
    models <- do.call("rbind",models)
  } else {
    models <- models[[1]]
  }

system.time(ct.bma <- coloc.bma(df1, df2, response1 = "expr", response2 = "cc", stratum2 = 'stratum', family1 = "gaussian", family2 = "binomial", thr = 0.05, r2.trim = 0.8 ))

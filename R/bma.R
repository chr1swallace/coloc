#'Wrapper to use colocalization testing within a Bayesian model averaging
#'structure.
#'
#'Performs the colocalisation tests described in Plagnol et al (2009) and
#'Wallace et al (2012).
#'
#'This is a test for proportionality of regression coefficients from two
#'independent regressions.  Analysis can either be based on a profile
#'likelihood approach, where the proportionality coefficient, \code{eta}, is
#'replaced by its maximum likelihood value, and inference is based on a
#'chisquare test (\code{p.value}), or taking a hybrid-Bayesian approach and
#'integrating the p value over the posterior distribution of \code{eta}, which
#'gives a posterior predictive p value.  The Bayesian approach can also be used
#'to give a credible interval for \code{eta}.  See the references below for
#'further details.
#'
#'@inheritParams coloc.test
#'@param df1,df2 Each is a dataframe, containing response and potential explanatory variables for two independent datasets.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response1,response2 The names of the response variables in \code{df1} and \code{df2} respectively
#'@param family1,family2 the error family for use in \code{glm}
#'@param thr posterior probability threshold used to trim SNP list.  Only SNPs with a marginal posterior probability of inclusion greater than this with one or other trait will be included in the full BMA analysis
#'@param nsnps number of SNPs required to model both traits.  The BMA analysis will average over all possible \code{nsnp} SNP models, subject to \code{thr} above.
#'@param n.approx number of values at which to numerically approximate the posterior
#'@param bayes.factor Numeric vector, giving value(s) of \code{eta} which
#'should be compared to the null values 0 and Inf.
#'@param r2.trim for pairs SNPs with r2>\code{r2.trim}, only one SNP will be retained.  This avoids numerical instability problems caused by including two highly correlated SNPs in the model.
#'@param ... other parameters passed to \code{coloc.test}
#'@return a \code{colocBayes} object
#'@author Chris Wallace
#'@references Wallace et al (2012).  Statistical colocalisation of monocyte
#'gene expression and genetic risk variants for type 1 diabetes.  Hum Mol Genet
#'21:2815-2824.  \url{http://europepmc.org/abstract/MED/22403184}
#'
#'Plagnol et al (2009).  Statistical independence of the colocalized
#'association signals for type 1 diabetes and RPS26 gene expression on
#'chromosome 12q13. Biostatistics 10:327-34.
#'\url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
#'@examples
#'
#'  ## simulate covariate matrix (X) and continuous response vector (Y)
#'  ## for two populations/triats Y1 and Y2 depend equally on f1 and f2
#'  ## within each population, although their distributions differ between
#'  ## populations.  They are compatible with a null hypothesis that they
#'  ## share a common causal variant
#'set.seed(1)
#'  X1 <- matrix(rbinom(2000,1,0.4),ncol=4)
#'  Y1 <- rnorm(500,rowSums(X1[,1:2]),2)
#'  X2 <- matrix(rbinom(2000,1,0.6),ncol=4)
#'  Y2 <- rnorm(500,rowSums(X2[,1:2]),5)
#'  
#'  boxplot(list(Y1,Y2),names=c("Y1","Y2"))
#'  
#'  ## fit and store linear model objects
#'  colnames(X1) <- colnames(X2) <- sprintf("f%s",1:ncol(X1))
#'  summary(lm1 <- lm(Y1~f1+f2+f3+f4,data=as.data.frame(X1)))
#'  summary(lm2 <- lm(Y2~f1+f2+f3+f4,data=as.data.frame(X2)))
#'  
#'
#'## test colocalisation using bma
#'df1=cbind(Y1=Y1,X1)
#'df2=cbind(Y2=Y2,X2)
#'
#'coloc.bma( df1, df2, snps=colnames(X1), response1="Y1", response2="Y2",
#'family1="gaussian", family2="gaussian",
#'nsnps=2,bayes.factor=c(1,2,3),plot.coeff=TRUE)
#'
#'@export
coloc.bma <- function(df1,df2,snps,response1="Y", response2="Y", family1="binomial", family2="binomial",
                     thr=0.01,nsnps=2,n.approx=1001, bayes.factor=NULL,
                     plot.coeff=FALSE,r2.trim=0.95,...) {
  snps <- unique(snps)
  if(length(snps)<2)
    return(1)
  if(is.matrix(df1))
    df1 <- as.data.frame(df1)
  if(is.matrix(df2))
    df2 <- as.data.frame(df2)

  ## remove any rows with missing genotypes
  use1 <- apply(!is.na(df1),1,all)
  use2 <- apply(!is.na(df2),1,all)
  if(any(!use1)) {
    cat("dropping",sum(!use1),"observations from dataset 1 due to missingness.\n")
    df1 <- df1[use1,]
  }
  if(any(!use2)) {
    cat("dropping",sum(!use2),"observations from dataset 2 due to missingness.\n")
    df2 <- df2[use2,]
  }
  
  ## remove highly correlated SNPs - df1
  n.orig <- length(snps)
  r2 <- cor(df1[,snps])^2
  r2.high <- apply(r2, 1, function(x) which(x>r2.trim)[1])
  snps <- snps[ r2.high == 1:length(r2.high) ]

  ## remove additional highly correlated SNPs - df2
  r2 <- cor(df2[,snps])^2
  r2.high <- apply(r2, 1, function(x) which(x>r2.trim)[1])
  snps <- snps[ r2.high == 1:length(r2.high) ]

  cat("Dropped",n.orig - length(snps),"of",n.orig,"SNPs due to LD: r2 >",r2.trim,"\n",length(snps),"SNPs remain.\n")

  ## ## remove any completely predictive SNPs
  ## r2 <- cor(df1[,response1], df1[,snps])^2
  ## r2 <- cor(df2[,response2], df2[,snps])^2
 
  f1 <- as.formula(paste(response1, "~", paste(snps,collapse="+")))
  f2 <- as.formula(paste(response2, "~", paste(snps,collapse="+")))
  lm1 <- glm(f1, data=df1, family=family1)
  lm2 <- glm(f2, data=df2, family=family2)
  while(any(is.na(coefficients(lm1))) || any(is.na(coefficients(lm2)))) {
    drop <- which((is.na(coefficients(lm1)) | is.na(coefficients(lm2)))[-1])
    cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n")
    cat(drop,"\n")
    snps <- snps[-drop]
    f1 <- as.formula(paste(response1, "~", paste(snps,collapse="+")))
    f2 <- as.formula(paste(response2, "~", paste(snps,collapse="+")))
    lm1 <- glm(f1, data=df1, family=family1)
    lm2 <- glm(f2, data=df2, family=family2)
  }

  ## x <- rbind(df1[,snps],df2[,snps])
  ## y <- c(df1[,"Y"]+1,2*df2[,"Y"]+1)
  ## mm <- mlogit2logit(y ~ ., data=cbind(y=y,x), base.choice=1)

  x1 <- df1[,snps]
  x2 <- df2[,snps]
 
  ## step1, select marginally interesting single SNPs
  models <- matrix(0,ncol(x1),ncol(x1))
  diag(models) <- 1
  if(family1=="binomial") {
    mods1 <- glib(x1,df1$Y, error="binomial", link="logit",models=models)
  } else {
    mods1 <- glib(x1,df1$Y, error="gaussian", link="identity",models=models)
  }
  if(family2=="binomial") {
    mods2 <- glib(x2,df2$Y, error="binomial", link="logit",models=models)
  } else {
    mods2 <- glib(x2,df2$Y, error="gaussian", link="identity",models=models)
  }
  probs1 <- mods1$bf$postprob[,2]
  probs2 <- mods2$bf$postprob[,2]
  use <- probs1 > thr | probs2 > thr
  cat("Restricting model space.\n",length(snps) - sum(use), "SNPs have single SNP posterior probabilities <",thr,
      "\nModels containing only these SNPs will not be explored.\n")

  ## step2, evaluate all pairs of marginally interesting single SNPs
#  combs <- lapply(nsnps, function(n) t(combn(which(use),n)))
  combs <- lapply(nsnps, function(n) {
    cmb <- combn(1:ncol(x1),n)
    idx <- apply(sapply(1:nrow(cmb), function(i) { cmb[i,] %in% which(use) }),1,any)
    t(cmb[,idx])
  }) 
#  combs <- lapply(nsnps, function(n) t(combn(which(use),n)))
  models <- lapply(combs, function(X) {
    tmp <- matrix(0,nrow(X),ncol(x1))
    for(j in 1:ncol(X)) {
      tmp[cbind(1:nrow(X),X[,j])] <- 1
    }
    tmp
  })
  if(length(snps)>1) {
    models <- do.call("rbind",models)
  } else {
    models <- models[[1]]
  }
  cat("Fitting",nrow(models),"models each containing",paste(nsnps,collapse="/"),"SNPs to dataset 1.\n")
  if(family1=="binomial") {
    mods1 <- glib(x1,df1$Y, error="binomial", link="logit",models=models)
  } else {
    mods1 <- glib(x1,df1$Y, error="gaussian", link="identity",models=models)
  }
  cat("Fitting",nrow(models),"models each containing",paste(nsnps,collapse="/"),"SNPs to dataset 2.\n")
  if(family2=="binomial") {
    mods2 <- glib(x2,df2$Y, error="binomial", link="logit",models=models)
  } else {
    mods2 <- glib(x2,df2$Y, error="gaussian", link="identity",models=models)
  }

  probs <- mods1$bf$postprob[,2] * mods2$bf$postprob[,2]
  probs <- probs/sum(probs)

  models.l <- matrix(as.logical(models),nrow(models))
  post <- matrix(NA,length(probs),n.approx)
  if(length(bayes.factor)) {
    bf <- matrix(NA,length(probs),length(bayes.factor)+2)
  } else {
    bf <- numeric(0)
  }
  var.1 <- var.2 <- coef.1 <- coef.2 <- vector("list",nrow(models.l))
  p <- matrix(NA,length(probs),5)
  colnames(p) <- c("eta.hat","chisquare","n","p","ppp")
  wh <- which(probs>thr^2)
  if(length(wh) < min(3,length(snps)))
    wh <- order(probs,decreasing=TRUE)[1:3]
  cat("Evaluating",length(wh),"models with posterior probabilities >",signif(min(probs[wh]),digits=2),"\n")
  for(i in wh) {
    cat(".")
    lm1 <- glm(df1$Y ~ ., data=x1[,models.l[i,]], family=family1)
    lm2 <- glm(df2$Y ~ ., data=x2[,models.l[i,]], family=family2)
    coef.1[[i]] <- coefficients(lm1)[-1]
    coef.2[[i]] <- coefficients(lm2)[-1]
    var.1[[i]] <- diag(vcov(lm1))[-1]
    var.2[[i]] <- diag(vcov(lm1))[-1]
    tmp <-  coloc.test(lm1,lm2,plot.coeff=FALSE,bma=TRUE,n.approx=n.approx,bayes.factor=bayes.factor,...)
    post[i,] <- tmp@bma
    p[i,] <- c(tmp@result,p.value(tmp),ppp.value(tmp))
    if(length(bayes.factor)) {
      bf[i,] <- tmp@bayes.factor
    }
  }
  p <- colSums(p[wh,] * probs[wh] / sum(probs[wh]))
  post <- colSums(post[wh,] * probs[wh] / sum(probs[wh]))
  if(length(bayes.factor)) {
    bf <- colSums(bf[wh,] * probs[wh] / sum(probs[wh]))
    names(bf) <- c("0","Inf",bayes.factor)
  }
  ci <- credible.interval(post,interval=c(0,pi), n.approx=n.approx)
  
  if(plot.coeff) {
    coeff.plot(unlist(coef.1),unlist(coef.2),
               unlist(var.1),unlist(var.2),
               eta=ci$eta.mean,
               main="Coefficients",
                                        #         sub=paste("ppp =",format.pval(ppp$value,digits=2),"p =",format.pval(pchisq(X2,df=length(snps)-1,lower.tail=FALSE),digits=2)),
               xlab=expression(b[1]),ylab=expression(b[2]))
  }

  return(new("colocBayes",
             result=c(eta.hat=p[1],chisquare=NA,n=p[3]),
             ppp=p[5],
             credible.interval=ci,
             bayes.factor=bf))
}


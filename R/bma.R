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
#'@inheritParams coloc.test.summary
#'@param df1,df2 Each is a dataframe, containing response and potential explanatory variables for two independent datasets.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response1,response2 The names of the response variables in \code{df1} and \code{df2} respectively
#'@param family1,family2 the error family for use in \code{glm}
#'@param thr posterior probability threshold used to trim SNP list.  Only SNPs with a marginal posterior probability of inclusion greater than this with one or other trait will be included in the full BMA analysis
#'@param nsnps number of SNPs required to model both traits.  The BMA analysis will average over all possible \code{nsnp} SNP models, subject to \code{thr} above.
#'@param n.approx number of values at which to numerically approximate the posterior
#'@param r2.trim for pairs SNPs with r2>\code{r2.trim}, only one SNP will be retained.  This avoids numerical instability problems caused by including two highly correlated SNPs in the model.
#'@param quiet suppress messages about how the model spaced is trimmed for BMA
#'@param ... other parameters passed to \code{coloc.test}
#'@return a \code{coloc} or \code{colocBayes} object
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
coloc.bma <- function(df1,df2,snps=intersect(setdiff(colnames(df1),response1),
                                setdiff(colnames(df2),response2)),
                      response1="Y", response2="Y", family1="binomial", family2="binomial",
                     bayes=!is.null(bayes.factor),
                      thr=0.01,nsnps=2,n.approx=1001, bayes.factor=NULL,
                     plot.coeff=FALSE,r2.trim=0.95,quiet=FALSE,...) {
  snps <- unique(snps)
  n.orig <- length(snps)
  if(n.orig<2)
    return(1)
  prep <- prepare.df(df1, snps, r2.trim=r2.trim, dataset=1, quiet=quiet)
  df1 <- prep$df
  snps <- prep$snps
  prep <- prepare.df(df2, snps, r2.trim=r2.trim, dataset=2, quiet=quiet)
  df2 <- prep$df
  snps <- prep$snps
  
  if(!quiet)
    cat("Dropped",n.orig - length(snps),"of",n.orig,"SNPs due to LD: r2 >",r2.trim,"\n",length(snps),"SNPs remain.\n")

  ## remove any completely predictive SNPs
  f1 <- as.formula(paste(response1, "~", paste(snps,collapse="+")))
  f2 <- as.formula(paste(response2, "~", paste(snps,collapse="+")))
  lm1 <- glm(f1, data=df1, family=family1)
  lm2 <- glm(f2, data=df2, family=family2)
  while(any(is.na(coefficients(lm1))) || any(is.na(coefficients(lm2)))) {
    drop <- which((is.na(coefficients(lm1)) | is.na(coefficients(lm2)))[-1])
    if(!quiet)
      cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n",drop,"\n")
    snps <- snps[-drop]
    f1 <- as.formula(paste(response1, "~", paste(snps,collapse="+")))
    f2 <- as.formula(paste(response2, "~", paste(snps,collapse="+")))
    lm1 <- glm(f1, data=df1, family=family1)
    lm2 <- glm(f2, data=df2, family=family2)
  }

  x1 <- df1[,snps]
  x2 <- df2[,snps]
  n.clean <- length(snps)
 
  ## step1, select marginally interesting single SNPs
  models <- diag(1,n.clean,n.clean)
  use <- marg.bf(models, x=x1, y=df1$Y, family=family1)[,"pp"] > thr |
    marg.bf(models, x=x2, y=df2$Y, family=family2)[,"pp"] > thr
  if(!quiet)
    cat("Restricting model space.\n",n.clean - sum(use), "SNPs have single SNP posterior probabilities <",thr,
      "\nModels containing only these SNPs will not be explored.\n")

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

  ## fit the models to each dataset to get posterior probs
  probs <- multi.bf(models, x1, df1$Y, family1, dataset=1,quiet=quiet) * multi.bf(models, x2, df2$Y, family2, dataset=2,quiet=quiet)
  probs <- probs/sum(probs)

  ## run coloc.test on each model
  models.l <- matrix(as.logical(models),nrow(models))
  post <- matrix(NA,length(probs),n.approx)
  if(length(bayes.factor)) {
    bf <- matrix(NA,length(probs),length(bayes.factor))
  } else {
    bf <- numeric(0)
  }
  var.1 <- var.2 <- coef.1 <- coef.2 <- vector("list",nrow(models.l))
  p <- matrix(NA,length(probs),4 + as.integer(bayes))
  colnames(p) <- c("eta.hat","chisquare","n","p","ppp")[1:ncol(p)]
  wh <- which(probs>thr^2)
  min.models <- min(3,length(snps))
  if(length(wh) < min.models) 
    wh <- order(probs,decreasing=TRUE)[1:min.models]
  cat("Averaging coloc testing over",length(wh),
      "models with posterior probabilities >=",
      signif(min(probs[wh]),digits=2),"\n")
  for(i in wh) {
    if(!quiet)
    cat(".")
    lm1 <- glm(df1$Y ~ ., data=x1[,models.l[i,]], family=family1)
    lm2 <- glm(df2$Y ~ ., data=x2[,models.l[i,]], family=family2)
    coef.1[[i]] <- coefficients(lm1)[-1]
    coef.2[[i]] <- coefficients(lm2)[-1]
    var.1[[i]] <- vcov(lm1)[-1,-1]
    var.2[[i]] <- vcov(lm2)[-1,-1]
    this.coloc <-  coloc.test.summary(coef.1[[i]], coef.2[[i]], var.1[[i]], var.2[[i]], plot.coeff=FALSE,bma=TRUE,n.approx=n.approx,bayes.factor=bayes.factor,bayes=bayes,...)
    if(bayes) {
      post[i,] <- this.coloc@bma ## posterior distribution for eta
      p[i,] <- c(this.coloc@result,p.value(this.coloc),ppp.value(this.coloc))
    } else {
      p[i,] <- c(this.coloc@result,p.value(this.coloc))
    }      
    if(length(bayes.factor)) {
      bf[i,] <- this.coloc@bayes.factor
    }
  }

  stats <- colSums(p[wh,] * probs[wh] / sum(probs[wh]))
  if(bayes) {
    ## average posteriors over models to get a single posterior
    post <- colSums(post[wh,] * probs[wh] / sum(probs[wh]))
    ## then calculate the credible interval
    ci <- credible.interval(post,interval=c(0,pi), n.approx=n.approx)
  }
  if(length(bayes.factor)) {
    ## average bayes factors for particular values/ranges of eta over models
    bf <- colSums(bf[wh,] * probs[wh] / sum(probs[wh]))
    names(bf) <- c(bayes.factor)
  }
  var.1 <- lapply(var.1, diag)
  var.2 <- lapply(var.2, diag)
  
  if(plot.coeff) {
    coeff.plot(unlist(coef.1),unlist(coef.2),
               unlist(var.1),unlist(var.2),
               eta=stats["eta.hat"],
               main="Coefficients",
                                        #         sub=paste("ppp =",format.pval(ppp$value,digits=2),"p =",format.pval(pchisq(X2,df=length(snps)-1,lower.tail=FALSE),digits=2)),
               xlab=expression(b[1]),ylab=expression(b[2]))
  }

  if(!bayes) {
    return(new("coloc",
               result=c(stats["eta.hat"],chisquare=NA,stats["n"],stats["p"]),
               method="BMA"))
  } else {
    return(new("colocBayes",
               result=c(stats["eta.hat"],chisquare=NA,stats["n"],stats["p"]),
               method="BMA",
               ppp=stats["ppp"],
               credible.interval=ci,
               bayes.factor=bf))
  }
}

multi.bf <- function(models, x, y, family, dataset=1,quiet=FALSE) {
  if(!quiet)
  cat("Fitting",nrow(models),"multi SNP models to dataset",dataset,"\n")
  if(family=="binomial") {
    mods1 <- glib(x, y, error="binomial", link="logit",models=models)
  } else {
    mods1 <- glib(x, y, error="gaussian", link="identity",models=models)
  }
  return(mods1$bf$postprob[,2])
}

marg.bf <- function(models,x,y,family) {
  if(family=="binomial") {
    mods1 <- glib(x, y, error="binomial", link="logit",models=models)
  } else {
    mods1 <- glib(x, y, error="gaussian", link="identity",models=models)
  }
  cbind(pp=mods1$bf$postprob[,2],twologB10=mods1$bf$twologB10[,2])
}

prepare.df <- function(df1,snps,r2.trim,dataset=1,quiet=FALSE) {
  if(is.matrix(df1))
    df1 <- as.data.frame(df1)
  use1 <- apply(!is.na(df1),1,all)
  if(any(!use1)) {
    if(!quiet)
      cat("dropping",sum(!use1),"observations from dataset",dataset,"due to missingness.\n")
    df1 <- df1[use1,]
  }
  r2 <- cor(df1[,snps])^2
  r2.high <- apply(r2, 1, function(x) which(x>r2.trim)[1])
  snps <- snps[ r2.high == 1:length(r2.high) ]
  return(list(df=df1,snps=snps))
}


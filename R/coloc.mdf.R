#'Wrapper to use colocalization testing within a Bayesian model averaging
#'structure for datasets with common controls.
#'
#'This is a test for proportionality of regression coefficients from a multinomial regression.  
#'Analysis is based on taking a hybrid-Bayesian approach and
#'integrating the p value over the posterior distribution of \code{eta}, which
#'gives a posterior predictive p value.  The Bayesian approach can also be used
#'to give a credible interval for \code{eta}.
#'
#'@inheritParams coloc.test.summary
#'@export
#'@param df1 A dataframe, containing response and potential explanatory variables for the dataset.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response The name of the response variable in \code{df1}
#'@param thr posterior probability threshold used to trim SNP list.  Only SNPs with a marginal posterior probability of inclusion greater than this with one or other trait will be included in the full BMA analysis
#'@param nsnps number of SNPs required to model both traits.  The BMA analysis will average over all possible \code{nsnp} SNP models, subject to \code{thr} above.
#'@param n.approx number of values at which to numerically approximate the posterior
#'@param r2.trim for pairs SNPs with r2> \code{r2.trim}, only one SNP will be retained.  This avoids numerical instability problems caused by including two highly correlated SNPs in the model.
#'@param quiet suppress messages about how the model spaced is trimmed for BMA
#'@param ... other parameters passed to \code{coloc.test}
#'@return a \code{coloc} or \code{colocBayes} object
#'@author Mary Fortune
coloc.var.bma <- function(df1,snps=setdiff(colnames(df1),response),
                      response="Y", bayes=!is.null(bayes.factor),
                      thr=0.01,nsnps=2,n.approx=1001, bayes.factor=NULL,
                      plot.coeff=FALSE,r2.trim=0.95,quiet=FALSE,...) {
    snps <- unique(snps)
    n.orig <- length(snps)
    if(n.orig<2)
        return(1)
    prep <- prepare.df(df1, snps, r2.trim=r2.trim, dataset=1, quiet=quiet)
    df1 <- prep$df
    snps <- prep$snps
    
    if(!quiet)
        cat("Dropped",n.orig - length(snps),"of",n.orig,"SNPs due to LD: r2 >",r2.trim,"\n",length(snps),"SNPs remain.\n")
    
    ## remove any completely predictive SNPs
    f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    while(any(is.na(coefficients(lm1)))) {
        drop <- which(is.na(coefficients(lm1))[-1])
        if(!quiet)
            cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n",drop,"\n")
        snps <- snps[-drop]
        f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
        capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    }
    f1 <- as.formula(paste("Y ~ 1 | ", paste(snps,collapse="+")))
    x1 <- df1[,snps]
    n.clean <- length(snps)
    binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1
    binX<-binX[,-1]
    binY<-binmod[,"Y.star"]
    
    ## step1, select marginally interesting single SNPs
    models <- cbind(diag(1,n.clean,n.clean),diag(1,n.clean,n.clean),rep(1,n.clean))
    use <- marg.var.bf(models, x=binX, y=binY, family="binomial")[,"pp"] > thr
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
    npairs<-length(unlist(combs))/2
    if(length(nsnps)>1) {
        models <- do.call("rbind",models)
    } else {
        models <- models[[1]]
    }
    modelsbin<-cbind(models,models,rep(1,npairs))
    ## fit the models to each dataset to get posterior probs
    probs <- multi.var.bf(modelsbin,  x=binX, y=binY, family="binomial",quiet=quiet)
    probs <- probs/sum(probs)
    
    ## run coloc.test on each model
    models.l <- matrix(as.logical(models),nrow(models))
    post <- matrix(NA,length(probs),n.approx)
    if(length(bayes.factor)) {
        bf <- matrix(NA,length(probs),length(bayes.factor))
    } else {
        bf <- numeric(0)
    }
    var <- var.1 <-var.2 <- coef.1 <- coef.2 <- vector("list",nrow(models.l))
    p <- matrix(NA,length(probs),4 + as.integer(bayes))
    colnames(p) <- c("eta.hat","chisquare","n","p","ppp")[1:ncol(p)]
    wh <- which(probs>thr^2)
    min.models <- min(3,length(snps))
    if(length(wh) < min.models) 
        wh <- order(probs,decreasing=TRUE)[1:min.models]
    cat("Averaging coloc testing over",length(wh),"models with posterior probabilities >=",signif(min(probs[wh]),digits=2),"\n")
     
    for(i in wh) {
        modelsnps<-snps[which(models.l[i,])]
        snps1<-paste("1:",unlist(modelsnps),sep="")
        snps2<-paste("2:",unlist(modelsnps),sep="")
        if(!quiet)
            cat(".")
        capture.output(multiglm <- multinom(as.formula(paste("Y~",paste(modelsnps,collapse="+"))), data=df1,maxit=1000))
        coef.1[[i]] <- coefficients(multiglm)[1,modelsnps]
        coef.2[[i]] <- coefficients(multiglm)[2,modelsnps]
        var[[i]] <- vcov(multiglm)[c(snps1,snps2),c(snps1,snps2)]
        var.1[[i]]<- vcov(multiglm)[snps1,snps1]
        var.2[[i]]<- vcov(multiglm)[snps2,snps2]
        this.coloc <-coloc.var.test.summary(coef.1[[i]], coef.2[[i]], var[[i]], plot.coeff=FALSE,bma=TRUE,n.approx=n.approx,bayes.factor=bayes.factor,bayes=bayes,level.ci=0.95)
        if(bayes) {
            post[i,] <- this.coloc@bma
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
        post <- colSums(post[wh,] * probs[wh] / sum(probs[wh]))
        ci <- credible.interval(post,interval=c(0,pi), n.approx=n.approx)
    }
    if(length(bayes.factor)) {
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

marg.var.bf <- function(models,x,y,family) {
  if(family=="binomial") {
    mods1 <- glib(x, y, error="binomial", link="logit",models=models)
  } else {
    mods1 <- glib(x, y, error="gaussian", link="identity",models=models)
  }
  cbind(pp=mods1$bf$postprob[,2],twologB10=mods1$bf$twologB10[,2])
}

multi.var.bf <- function(models, x, y, family,quiet=FALSE) {
  if(!quiet)
    cat("Fitting",nrow(models),"multi SNP models to dataset \n")
  if(family=="binomial") {
    mods1 <- glib(x, y, error="binomial", link="logit",models=models)
  } else {
    mods1 <- glib(x, y, error="gaussian", link="identity",models=models)
  }
  return(mods1$bf$postprob[,2])
}


#'Prepares principal components of dataset with common control for
#'colocalisation testing.
#'
#'If \code{X1} is a \code{SnpMatrix} object, it are checked for
#'missing data, and any missing values imputed by repeated use of
#'\code{impute.snps} from the \code{snpStats} package.
#'
#'Principal components are calculated using \code{prcomp}.
#'
#'\code{pcs.var.model} can then be invoked to create \code{glm} objects.
#'
#'@aliases pcs.var.prepare
#'@param X Either a SnpMatrix or numeric matrix of genetic data.
#'Columns index SNPs, rows index samples.
#'@return a \code{colocPCs} object.
#'@export
#'@author Mary Fortune
pcs.var.prepare <- function(X) {
  #X gives the SNP matrix - all entries are 0,1,2
  snps <- colnames(X)
  if(length(snps)<2)
    stop("require at least 2 SNPs in object X")
  if(is(X,"SnpMatrix")) {
    if(any(X==as.raw("0"))) {
      X <- fillin(X)
    } else {
      X <- as(X,"numeric")
    }
  }
  rows.drop <- apply(is.na(X),1,any)
  if(sum(rows.drop))
    X <- X[!rows.drop,]
  pcs <- prcomp(X, centre=TRUE, scale.=TRUE)
  vars <- pcs$sdev^2
  cvars <- cumsum(vars)/sum(vars)
  
  return(new("colocPCs",
             pcs=pcs$x,
             use=!rows.drop,
             vars=cvars))    
}






#'Prepares models of response based on principal components of two datasets for
#'colocalisation testing.
#'
#'@aliases pcs.var.model
#'
#'@param object A colocPCs object, result of \code{pcs.var.prepare()}.
#'@param Y Numeric phenotype vector
#'@param threshold The minimum number of principal components which captures at
#'least threshold proportion of the variance will be selected.  Simulations
#'suggest \code{threshold=0.8} is a good default value.
#'@return \code{pcs.prepare} returns a \code{colocPCs} object.
#'@export
#'\code{pcs.model} returns a \code{glm} object.
#'@author Mary Fortune
pcs.var.model <- function(object, Y, threshold=0.8) {
  if(length(object@vars)<2)
    stop("require 2 or more principal components to test for proportionality")
  if (!(all(Y %in% c(0,1,2)))){
    stop("require the response values to be in 0, 1 or 2")
  }
  npc <- which(object@vars>threshold)[1]
  if(npc<2)
    npc <- 2
  cat("selecting",npc,"components out of",ncol(object@pcs),"to capture",object@vars[npc],"of total variance.\n")
  X <- object@pcs[, 1:npc ]
  if(nrow(X) != length(Y)){
    stop("length of Y not equal to the number of rows from X \n")
  }
  snps <- colnames(object@pcs)[1:npc]
  f <- as.formula(paste("Y ~", paste(snps,collapse="+")))
  data <-  as.data.frame(cbind(Y,X))
  m <- multinom(f,data=data,maxit=1000)
  while(any(is.na(coefficients(m)))) {
    drop <- which((is.na(coefficients(m)))[-1])
    cat("dropping",length(drop),"colinear PCs",snps[drop],"\n")
    snps <- snps[-drop]
    f <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    m <- multinom(f,data=data,maxit=1000)
  }
  return(m)  
}



#'Wrapper to use colocalization testing within a principle components
#'structure for datasets with common controls.
#'
#'This is a test for proportionality of regression coefficients from a multinomial regression.  
#'Analysis is based on taking a hybrid-Bayesian approach and
#'integrating the p value over the posterior distribution of \code{eta}, which
#'gives a posterior predictive p value.  The Bayesian approach can also be used
#'to give a credible interval for \code{eta}.
##'
##' @param X A glm object for the traits.  Any
##' Intercept term is dropped, but other covariates should have distinct names or
##' be listed in \code{vars.drop} to avoid them being included in the
##' colocalisation test.
##' @param vars.drop Character vector naming additional variables in either
##' regression which are not SNPs and should not be used in the colocalisation
##' test.  They should appear in
##' \code{c(names(coefficients(X)))}
##' @param ... other arguments passed to \code{\link{coloc.test.summary}()}.
##' @return a numeric vector with 3 named elements:
##' \item{eta.hat}{The estimated slope.}
##' \item{chisquare}{The chisquared test statistic}
##' \item{n}{The number of snps used in the test.  If eta were known, this
##' would be the degrees of freedom of the test. Because eta has been replaced by
##' its ML estimate, Plagnol et al suggest we expect the degrees of freedom to be
##' n-1, but this requires the likelihood to be well behaved which is not always
##' the case.  We prefer to consider the posterior predictive p value.}
##' \item{ppp}{The posterior predictive p value}
##' @export
##' @author Mary Fortune
coloc.var.test <- function(X,vars.drop=NULL, ...) {
  ## X is a multinom objects, fitted to the snps
  ## return values are
  ## return(c(eta.hat=eta.hat,chisq=X2,ppp=ppp$value))
  ## where
  ## eta.hat is the estimated slope
  ## chisq is the test statistic (degrees of freedom <= number of snps)
  ## ppp is the posterior predictive p value
  
  vars.drop <- c(vars.drop,"(Intercept)")
  snps <- setdiff(colnames(coefficients(X)),vars.drop)
  snps.dropX <- setdiff(colnames(coefficients(X)),c(snps,vars.drop))
  if(length(snps.dropX))
    cat("Warning:",length(snps.dropX),"variables dropped from regression:\n\t",snps.dropX,"\n")
  
  if(length(snps)<=1) { # 1 common coef => equal already
    cat("Only 1 factor,",snps," in regression.  Skipping\n")
    return(NULL)
  }
  b1 <- coefficients(X)[1,snps]
  b2 <- coefficients(X)[2,snps]
  
  if(any(is.na(b1)) | any(is.na(b2)))
    stop("Coefficients not all estimated.  Please refit the regressions.\n")
  snps1<-paste("1:",unlist(snps),sep="")
  snps2<-paste("2:",unlist(snps),sep="")
  #find the covariance matrix
  V <- vcov(X)[c(snps1,snps2),c(snps1,snps2)]
  coloc.var.test.summary(b1,b2,V,...)
}

##' Colocalisation testing supplying only regression coefficients and their variance-covariants matrices
##'
##' Typically this should be called from \code{\link{coloc.var.test}()}, but is left as a public function, to use at your own risk, if you have some other way to define the SNPs under test.
##' @title Colocalisation testing for shared controls using regression coefficients
##' @return an object of class coloc, colocBayes or colocBMA
##' @author Chris Wallace
##' @inheritParams coloc.test
##' @export
##' @param b1 regression coefficients for trait 1
##' @param b2 regression coefficients for trait 2
##' @param V variance-covariance matrix
##' @param k Theta has a Cauchy(0,k) prior.  The default, k=1, is equivalent to a
##' uniform (uninformative) prior.  We have found varying k to have little effect
##' on the results.
##' @param plot.coeff \code{TRUE} if you want to generate a plot showing the
##' coefficients from the two regressions together with confidence regions.
##' @param bma parameter set to \code{TRUE} when \code{coloc.test} is called by \code{coloc.bma}.  DO NOT SET THIS WHEN CALLING \code{coloc.test} DIRECTLY!
##' @param plots.extra list with 2 named elements, x and y, equal length
##' character vectors containing the names of the quantities to be plotted on the
##' x and y axes.
##' 
##' \code{x} is generally a sequence of \code{theta} and \code{eta}, with
##' \code{y} selected from \code{post.theta}, the posterior density of theta,
##' \code{chisq}, the chi-square values of the test, and \code{lhood}, the
##' likelihood function.
##' @param bayes Logical, indicating whether to perform Bayesian
##' inference for the coefficient of proportionality, eta.  If
##' \code{bayes.factor} is supplied, Bayes factors are additionally
##' computed for the specificed values.  This can add a little time as
##' it requires numerical integration, so can be set to FALSE to save
##' time in simulations, for example.
##' @param bayes.factor Either a numeric vector, giving single value(s) of \code{eta} or a list of numeric vectors, each of length two and specifying ranges of eta which should be compared to each other.  Thus, the vector or list needs to have length at least two.
##' @param level.ci,n.approx \code{level.ci} denotes the required level of the
##' credible interval for \code{eta}.  This is calculated numerically by
##' approximating the posterior distribution at \code{n.approx} distinct values.
##' @author Mary Fortune
coloc.var.test.summary <- function(b1,b2,V,k=1,plot.coeff=TRUE,plots.extra=NULL,bayes=!is.null(bayes.factor),
                               n.approx=1001, level.ci=0.95,
                               bayes.factor=NULL, bma=FALSE) {
  snpnum <- length(b1)
  b<-c(b1,b2)
  V11<-V[1:snpnum,1:snpnum]
  V12<-V[1:snpnum,(snpnum+1):(2*snpnum)]
  V21<-V[(snpnum+1):(2*snpnum),1:snpnum]
  V22<-V[(snpnum+1):(2*snpnum),(snpnum+1):(2*snpnum)]
  W<-solve(V)
  W11<-W[1:snpnum,1:snpnum]
  W12<-W[1:snpnum,(snpnum+1):(2*snpnum)]
  W21<-W[(snpnum+1):(2*snpnum),1:snpnum]
  W22<-W[(snpnum+1):(2*snpnum),(snpnum+1):(2*snpnum)]
  theta.min <- 0
  theta.max <- pi
  
  ## -2LL = Fieller's chisq
  d <- function(theta,b1,b2) { sin(theta) * b1 - cos(theta) * b2 }
  Vstar <- function(theta) { sin(theta)^2 * V11 -sin(theta)*cos(theta)*V12-sin(theta)*cos(theta)*V21+ cos(theta)^2 * V22 }
  chisq <- function(theta,b1,b2) { t(d(theta,b1,b2)) %*% solve(Vstar(theta)) %*% d(theta,b1,b2) }
  chisqV <- Vectorize(chisq, "theta")
  
  findmin <- function(b1,b2) {
    ## there are at most two minima, and never both on the same side of pi/2
    o.left <- optimize(chisq,interval=c(0,pi/2),b1=b1,b2=b2)
    o.right <- optimize(chisq,interval=c(pi/2,pi),b1=b1,b2=b2)
    if(o.left$objective < o.right$objective) {
      return(o.left)
    } else {
      return(o.right)
    }
  }
  fm <- findmin(b1,b2)
  theta.hat <- fm$minimum; eta.hat=tan(theta.hat)
  X2 <- fm$objective[1,1]
  
  ################################################################################
  
  ## Bayesian inference
  
  if(bayes) {
    ## cauchy prior for theta
    prior <- function(theta) { tt <- tan(theta);
                               k*(1+tt^2) / (2*pi*(1 + k^2 * tt^2)) }
    priorV <- Vectorize(prior,"theta")
    
    ## posterior for theta
    p <- length(b1)
    const <- ( sqrt(2*pi)^p * det(V))^(-0.5)
    M <- function(theta) { solve(cos(theta)^2 * W11 + cos(theta)*sin(theta)*W12 + cos(theta)*sin(theta)*W21+sin(theta)^2 * W22) }
    mu <- function(theta) { M(theta) %*% (  cos(theta)*(W11%*%b1) +cos(theta)*(W12%*%b2)+sin(theta)*(W21%*%b1)+sin(theta)*(W22%*%b2) )}
    L <- function(theta) {
      const * prior(theta) * det(M(theta))^(0.5) *
        exp( -0.5 * (t(b1) %*% W11 %*% b1 + t(b1) %*% W12 %*% b2 +t(b2) %*% W21 %*% b1 +t(b2) %*% W22 %*% b2 -
                       t(mu(theta)) %*% solve(M(theta)) %*% mu(theta)) )
    }
    LV <- Vectorize(L,"theta")
    LV.int <- integrate(LV,lower=0,upper=pi)
    post <- function(theta) { LV(theta) / LV.int$value }
    
    ##  posterior predictive p value
    pv <- function(theta) { pchisq(chisq(theta,b1,b2),df=p,lower.tail=FALSE) }
    pval <- Vectorize(pv,"theta")
    toint <- function(theta) { pval(theta) * post(theta) }
    ppp <- integrate(toint,lower=theta.min,upper=theta.max)
    
    if(bma) {
      ## numeric approx of lhood for BMA
      theta.bma <- seq(0,pi,length=n.approx)
      post.bma <- post(theta.bma)
    }
    
    ## bayes factors
    bf.calc <- function(eta) {
      if(length(eta)==1)
        return( LV(atan(eta))/priorV(atan(eta)) )      
      if(length(eta)==2) {
        theta.range <- sort(atan(eta))
        return(integrate(LV,lower=theta.range[1],upper=theta.range[2])$value/
                 integrate(priorV,lower=theta.range[1], upper=theta.range[2])$value)
      }
      warning(paste("Cannot calculate Bayes Factor for eta =",eta,"require length == 1 or 2."))
      return(NA)
    }
    
    if(!is.null(bayes.factor)) {
      if(length(bayes.factor)<2)
        warning("You are trying to prepare to calculate Bayes Factors for a single value or interval for eta.  bayes.factor should have length at least two.")
      post.bf <- sapply(bayes.factor, bf.calc)
      names(post.bf) <- c(bayes.factor)
    } else {
      post.bf <- numeric(0)
    }
    
    ## credible interval - only if not doing bma
    if(!bma) {
      cred.int <- credible.interval(post, interval=c(theta.min, theta.max), n.approx=n.approx,
                                    level.ci=level.ci)
    }
  }
  
  ################################################################################
  
  ## plots
  if(plot.coeff) {
    coeff.plot(b1,b2,diag(V11),diag(V22),eta=eta.hat,
               main="Coefficients",
               #         sub=paste("ppp =",format.pval(ppp$value,digits=2),"p =",format.pval(pchisq(X2,df=nsnps-1,lower.tail=FALSE),digits=2)),
               xlab=expression(b[1]),ylab=expression(b[2]))
  }
  
  x <- seq(theta.min,theta.max,length=1001)
  
  if(!is.null(plots.extra)) {
    plot.data <- list(theta=x,
                      eta=tan(x),
                      chisq=chisqV(x,b1,b2),
                      post.theta=if(bayes) { post(x) } else {rep(NA,length(x))},
                      lhood=chisqV(x,b1,b2))
    extra.plot(plot.data, plots.extra, theta.hat=theta.hat, eta.hat=eta.hat)   
  }
  
  ################################################################################
  
  ## return
  if(!bayes) {
    return(new("coloc",
               result=c(eta.hat=eta.hat,chisquare=X2,n=nsnps),
               method="single"))
  } else {
    if(!bma) {
      return(new("colocBayes",
                 result=c(eta.hat=eta.hat,chisquare=X2,n=nsnps),
                 method="single",
                 ppp=ppp$value,
                 credible.interval=cred.int,
                 bayes.factor=post.bf))
    } else {
      return(new("colocBayesBMA",
                 result=c(eta.hat=eta.hat,chisquare=X2,n=nsnps),
                 method="single",
                 ppp=ppp$value,
                 bma=post.bma,
                 bayes.factor=post.bf))
    }
  }  
}


##'Code to determine whether, for a specific variance-covariance matrix for shared controls
##'the form of the test statistic is the same under the two approaches
##'@title Checking equivalence of test statistics
##'@param V Variance-covariance matrix
##'@return a number giving the maximum difference between the two statistics
##'@author Mary Fortune
checklikechi<-function(V){
  #for a specific covariance matrix V, we check that the form of the chi squared statistic is that of the likelihood
    n=ncol(V)/2
    V11=V[1:n,1:n]
    V12=V[1:n,(n+1):(2*n)]
    V21=V[(n+1):(2*n),1:n]
    V22=V[(n+1):(2*n),(n+1):(2*n)]
    W=solve(V)
    W11=W[1:n,1:n]
    W12=W[1:n,(n+1):(2*n)]
    W21=W[(n+1):(2*n),1:n]
    W22=W[(n+1):(2*n),(n+1):(2*n)]
    chilike<-function(t){
      c=cos(t)
      s=sin(t)
    iX=solve(c*c*W11+c*s*W12+c*s*W21+s*s*W22)
    iY=solve(s*s*V11-c*s*V12-c*s*V21+c*c*V22)
    b1val=c*W11%*%iX%*%W12+s*W11%*%iX%*%W22-c*W12%*%iX%*%W11-s*W12%*%iX%*%W21-s*iY
    return(max(abs(b1val)))
  }
  minval<-optimize(chilike, interval=c(0, pi/2), maximum=TRUE)$objective
  return(minval)
}


#'@title finding the minimum single SNP p-value for each of the two traits
#'@param df1 A dataframe, containing response and potential explanatory variables for the dataset.
#'@param snps The SNPs to consider as potential explanatory variables
#'@param response The name of the response variable in \code{df1}
#'@param r2.trim for pairs SNPs with r2> \code{r2.trim}, only one SNP will be retained.  This avoids numerical instability problems caused by including two highly correlated SNPs in the model.
#'@param quiet suppress messages
#'@return the minimum single SNP p value for each trait
#'@author Mary Fortune
singlesnp.twotrait<-function(df1,response="Y",snps=setdiff(colnames(df1),response),r2.trim=0.95,quiet=TRUE){
    #returns the minimum single SNP p value for each trait
    snps <- unique(snps)
    n.orig <- length(snps)
    if(n.orig<2)
        return(1)
    prep <- prepare.df(df1, snps, r2.trim=r2.trim, dataset=1, quiet=quiet)
    df1 <- prep$df
    snps <- prep$snps
    
    if(!quiet)
        cat("Dropped",n.orig - length(snps),"of",n.orig,"SNPs due to LD: r2 >",r2.trim,"\n",length(snps),"SNPs remain.\n")
    
    ## remove any completely predictive SNPs
    f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    while(any(is.na(coefficients(lm1)))) {
        drop <- which(is.na(coefficients(lm1))[-1])
        if(!quiet)
            cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n",drop,"\n")
        snps <- snps[-drop]
        f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
        capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    }
    
    df.t1<-df1[which(df1[,response]!=2),]
    df.t2<-df1[which(df1[,response]!=1),]
    df.t2[,response]=df.t2[,response]/2
    capture.output(snp.data.t1<-new("SnpMatrix",as.matrix(df.t1[,snps])))
    capture.output(snp.data.t2<-new("SnpMatrix",as.matrix(df.t2[,snps])))
    
    results.t1<-single.snp.tests(df.t1[,response], snp.data=snp.data.t1)
    capture.output(show.results.t1<-show(results.t1))
    
    results.t2<-single.snp.tests(df.t2[,response], snp.data=snp.data.t2)
    capture.output(show.results.t2<-show(results.t2))
    
    cat("The minimum single SNP p value for RA is: ", min(show.results.t1[,4]),"\n")
    cat("The minimum single SNP p value for T1D is: ", min(show.results.t2[,4]),"\n")
    return(c(min(show.results.t1[,4]),min(show.results.t2[,4])))
}



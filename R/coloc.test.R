credible.interval <- function(post, interval, n.approx, level.ci=0.95) {
  if(is.function(post)) {
    ## minimum, so that we can reset theta.min, theta.max and centre on the mode
    minimum <- optimize(post, interval)
    t.min <- minimum$minimum
    t.max <- minimum$minimum + pi
    t.seq <- seq(t.min, t.max, length=n.approx)
    pval <- post(t.seq)
#    tpost <- function(theta) theta*post(theta)
#    t.mean <- integrate(tpost,lower=interval[1], upper=interval[2])$value
  #  epost <- function(theta) tan(theta)*post(theta)
#    e.mean <- integrate(epost,lower=interval[1], upper=interval[2])$value
#    e.mean <- tan(t.mean)
  } else {
    t.seq <- seq(interval[1],interval[2],length=n.approx)
    minimum <- which.min(post)
    t.min <- t.seq[minimum]
    t.max <- t.min + pi
    t.seq <- seq(t.min, t.max, length=n.approx)
    pval <- post[c(minimum:n.approx,1:(minimum - 1))]
#    e.mean <- tan(sum(pval * t.seq)/sum(pval))
  }    
  pval <- pval/sum(pval)
  
  ## search
  centre <- which.max(pval)
  t.mode <- t.seq[centre]
  v <- pval[centre-1]
  l <- which.min(abs(pval[1:(centre+1)] - v))
  u <- which.min(abs(pval[(centre+1):n.approx] - v)) + centre
  while(sum(pval[l:u])<level.ci) { # go down in twenties
    cat(".")
    v <- max(pval[l-20], pval[u+20])
    l <- which.min(abs(pval[1:(centre+1)] - v))
    u <- which.min(abs(pval[(centre+1):n.approx] - v)) + centre
  }
  while(sum(pval[l:u])>level.ci) { # go up in ones
    cat(".")
    v <- min(pval[l+1], pval[u-1])
    l <- which.min(abs(pval[1:(centre+1)] - v))
    u <- which.min(abs(pval[(centre+1):n.approx] - v)) + centre
  }

  return(list(#eta.mean=e.mean, #tan(t.mean %% pi),
              eta.mode=tan(t.mode), 
              lower=tan(t.seq[l]), upper=tan(t.seq[u]), level=level.ci,
              level.observed=if(is.function(post)) {
                integrate(post, lower=t.seq[l], upper=t.seq[u])$value
              } else { NA },
              interior=t.seq[l] %/% pi == t.seq[u] %/% pi))
}



##' Function to do colocalisation tests of two traits
##' 
##' Performs the colocalisation tests described in Plagnol et al (2009) and
##' Wallace et al (2012).
##' 
##' This is a test for proportionality of regression coefficients from two
##' independent regressions.  Analysis can either be based on a profile
##' likelihood approach, where the proportionality coefficient, \code{eta}, is
##' replaced by its maximum likelihood value, and inference is based on a
##' chisquare test (\code{p.value}), or taking a hybrid-Bayesian approach and
##' integrating the p value over the posterior distribution of \code{eta}, which
##' gives a posterior predictive p value.  The Bayesian approach can also be used
##' to give a credible interval for \code{eta}.  See the references below for
##' further details.
##'
##' @param X Either an lm or glm object for trait 1.  The intersection of
##' \code{names(coefficients(X))} and \code{names(coefficients(Y))} is used to
##' identify SNPs in common which will be tested for colocalisation.  Any
##' Intercept term is dropped, but other covariates should have distinct names or
##' be listed in \code{vars.drop} to avoid them being included in the
##' colocalisation test.
##' @param Y Either an lm or glm object for trait 2.
##' @param vars.drop Character vector naming additional variables in either
##' regression which are not SNPs and should not be used in the colocalisation
##' test.  They should appear in
##' \code{c(names(coefficients(X)),names(coefficients(Y)))}
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
##' @note Plagnol et al's original test was available in his R package
##' \code{QTLMatch v0.8} which now appears unavailable.  The numerically
##' identical test, extended to allow for more than two SNPs, can be found in
##' this package by looking at the chisquare statistic and the degrees of freedom
##' given by \code{chisquare()} and \code{df()} respectively.  %
##' \url{http://www-gene.cimr.cam.ac.uk/vplagnol/software.shtml}
##' @author Chris Wallace
##' @references Wallace et al (2012).  Statistical colocalisation of monocyte
##' gene expression and genetic risk variants for type 1 diabetes.  Hum Mol Genet
##' 21:2815-2824.  \url{http://europepmc.org/abstract/MED/22403184}
##' 
##' Plagnol et al (2009).  Statistical independence of the colocalized
##' association signals for type 1 diabetes and RPS26 gene expression on
##' chromosome 12q13. Biostatistics 10:327-34.
##' \url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
##' @examples
##' 
##'   ## simulate covariate matrix (X) and continuous response vector (Y)
##'   ## for two populations/triats Y1 and Y2 depend equally on f1 and f2
##'   ## within each population, although their distributions differ between
##'   ## populations.  They are compatible with a null hypothesis that they
##'   ## share a common causal variant
##' set.seed(1)
##'   X1 <- matrix(rbinom(1000,1,0.4),ncol=2)
##'   Y1 <- rnorm(500,apply(X1,1,sum),2)
##'   X2 <- matrix(rbinom(1000,1,0.6),ncol=2)
##'   Y2 <- rnorm(500,2*apply(X2,1,sum),5)
##'   
##'   boxplot(list(Y1,Y2),names=c("Y1","Y2"))
##'   
##'   ## fit and store linear model objects
##'   colnames(X1) <- colnames(X2) <- c("f1","f2")
##'   summary(lm1 <- lm(Y1~f1+f2,data=as.data.frame(X1)))
##'   summary(lm2 <- lm(Y2~f1+f2,data=as.data.frame(X2)))
##'   
##'   ## test whether the traits are compatible with colocalisation
##'   ### ppp should be large (>0.05, for example), indicating that they are.
##'   par(mfrow=c(2,2))
##'   coloc.test(lm1,lm2,plot.coeff=TRUE,
##'              plots.extra=list(x=c("eta","theta"),
##'                               y=c("lhood","lhood")))
##' 
##' @export
coloc.test <- function(X,Y,vars.drop=NULL, ...) {
  ## X and Y are glm objects, fitted to the same snps, with different outcome variables
  ## return values are
  ## return(c(eta.hat=eta.hat,chisq=X2,ppp=ppp$value))
  ## where
  ## eta.hat is the estimated slope
  ## chisq is the test statistic (degrees of freedom <= number of snps)
  ## ppp is the posterior predictive p value
  
  vars.drop <- c(vars.drop,"(Intercept)")
  snps <- setdiff(intersect(names(coefficients(X)),names(coefficients(Y))),
                  vars.drop)
  snps.dropX <- setdiff(names(coefficients(X)),c(snps,vars.drop))
  snps.dropY <- setdiff(names(coefficients(Y)),c(snps,vars.drop))
  if(length(snps.dropX))
    cat("Warning:",length(snps.dropX),"variables dropped from regression X:\n\t",snps.dropX,"\n")
  if(length(snps.dropY))
    cat("Warning:",length(snps.dropY),"variables dropped from regression Y:\n\t",snps.dropY,"\n")
  
  if(length(snps)<=1) { # 1 common coef => equal already
    cat("Only 1 factor,",snps," in common.  Skipping\n")
    return(NULL)
  }
  b1 <- coefficients(X)[snps]
  b2 <- coefficients(Y)[snps]
  
  if(any(is.na(b1)) | any(is.na(b2)))
    stop("Coefficients not all estimated.  Please refit the regressions.\n")
  
  V1 <- vcov(X)[snps,snps]
  V2 <- vcov(Y)[snps,snps]
  coloc.test.summary(b1,b2,V1,V2,...)
}

##' Colocalisation testing supplying only regression coefficients and their variance-covariants matrices
##'
##' Typically this should be called from \code{\link{coloc.test}()} or \code{\link{coloc.bma}()}, but is left as a public function, to use at your own risk, if you have some other way to define the SNPs under test.
##' @title Colocalisation testing using regression coefficients
##' @return an object of class coloc, colocBayes or colocBMA
##' @author Chris Wallace
##' @inheritParams coloc.test
##' @export
##' @param b1 regression coefficients for trait 1
##' @param b2 regression coefficients for trait 2
##' @param V1 variance-covariance matrix for trait 1
##' @param V2 variance-covariance matrix for trait 2
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
##' @param bayes.factor Calculate Bayes Factors to compare specific values of eta.  \code{bayes.factor} should either a numeric vector, giving single value(s) of \code{eta} or a list of numeric vectors, each of length two and specifying ranges of eta which should be compared to each other.  Thus, the vector or list needs to have length at least two.
##' @param level.ci,n.approx \code{level.ci} denotes the required level of the
##' credible interval for \code{eta}.  This is calculated numerically by
##' approximating the posterior distribution at \code{n.approx} distinct values.
coloc.test.summary <- function(b1,b2,V1,V2,k=1,plot.coeff=TRUE,plots.extra=NULL,bayes=!is.null(bayes.factor),
                               n.approx=1001, level.ci=0.95,
                               bayes.factor=NULL, bma=FALSE) {
  nsnps <- length(b1)
  S1 <- solve(V1)
  S2 <- solve(V2)
  theta.min <- 0
  theta.max <- pi
  
  ## -2LL = Fieller's chisq
  d <- function(theta,b1,b2) { sin(theta) * b1 - cos(theta) * b2 }
  Vstar <- function(theta) { sin(theta)^2 * V1 + cos(theta)^2 * V2 }
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
    const <- ( sqrt(2*pi)^p * det(V1) * det(V2) )^(-1)
    M <- function(theta) { solve(cos(theta)^2 * S1 + sin(theta)^2 * S2) }
    mu <- function(theta) { t( (cos(theta) * t(b1) %*% S1 +
                                sin(theta) * t(b2) %*% S2) %*% M(theta) ) }
    L <- function(theta) {
      const * prior(theta) * det(M(theta))^(-0.5) *
        exp( -0.5 * (t(b1) %*% S1 %*% b1 + t(b2) %*% S2 %*% b2 -
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
    coeff.plot(b1,b2,diag(V1),diag(V2),eta=eta.hat,
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



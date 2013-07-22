## given two snpMatrix objects or numeric matrices, find the common principal component structure, pick the top pcs and return regression objects ready for coloc testing

#' Impute missing genotypes
#'
#' Impute missing genotypes in a snpMatrix object in each SNP in turn, conditional on all the others.
#' @param X a snpMatrix object
#' @param bp optional vector giving basepair positions of the SNPs
#' @param strata optional vector giving stratificiation of the samples, one entry for each sample, and samples with the same value are assumed to come from a single strata
#' @return a numeric matrix of imputed genotypes, 0,2 = homs, 1 = het
fillin <- function(X,bp=1:ncol(X),strata=NULL) {
  N <- as(X,"numeric")
  if(!is.null(strata)) {
    strata <- as.factor(strata)
    if(length(levels(strata))>10)
      stop("too many levels in strata\n")
    for(i in levels(strata)) {
      cat("\nstrata",i,"\n")
      wh <- which(strata==i)
      N[wh,] <- fillin(X[wh,,drop=FALSE],bp)
    }
    return(N)
  }

  csumm <- col.summary(X)
  wh <- which(csumm[,"Call.rate"]<1)
  cat(length(wh),"to impute\n")
  for(i in wh) {
    cat(i,".")
    rule <- snp.imputation(X[,-i,drop=FALSE],X[,i,drop=FALSE],bp[-i],bp[i])
    if(is.null(rule@.Data[[1]]))
      next
    imp <- impute.snps(rule,X[,-i,drop=FALSE])
    wh.na <- which(is.na(N[,i]))
    N[wh.na,i] <- imp[wh.na]
  }
  cat("\n")
  N
}
   

#'Functions to prepare principle component models for colocalisation testing
#'
#'Prepares principal components of two datasets for
#'colocalisation testing.
#'
#'If \code{X1} and \code{X2} are \code{SnpMatrix} objects, they are checked for
#'missing data, and any missing values imputed by repeated use of
#'\code{impute.snps} from the \code{snpStats} package.
#'
#'Columns with common names are \code{rbind}ed together and principal
#'components calculated using \code{prcomp}.
#'
#' \code{pcs.model} can then be invoked to create \code{glm} objects.
#'
#'@aliases pcs.prepare
#'@param X1,X2 Each is either a SnpMatrix or numeric matrix of genetic data.
#'Columns index SNPs, rows index samples.
#'@return a \code{colocPCs} object.
#'
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
#'  ## share a common causal variant, with the effect twice as strong for
#'  ## Y2 as Y1
#'set.seed(1)
#'  X1 <- matrix(rbinom(5000,1,0.4),ncol=10)
#'  Y1 <- rnorm(500,apply(X1[,1:2],1,sum),2)
#'  X2 <- matrix(rbinom(5000,1,0.6),ncol=10)
#'  Y2 <- rnorm(500,2*apply(X2[,1:2],1,sum),5)
#'  
#'  ## generate principal components object
#'  colnames(X1) <- colnames(X2) <- make.names(1:ncol(X1))
#'  pcs <- pcs.prepare(X1,X2)
#'
#'  ## generate glm objects
#'  m1 <- pcs.model(pcs, group=1, Y=Y1)
#'  m2 <- pcs.model(pcs, group=2, Y=Y2)
#'
#'  ## test colocalisation using PCs
#'  coloc.test(m1,m2,plot.coeff=FALSE,bayes=FALSE)
#'
#'@export
pcs.prepare <- function(X1, X2) {
  snps.common <- intersect(colnames(X1),colnames(X2))
  if(length(snps.common)<2)
    stop("require at least 2 SNPs in common between objects X1, X2")
  if(!identical(class(X1),class(X2)))
    stop("require X1 and X2 to be of same class")
  if(is(X1,"SnpMatrix")) {
    if(any(X1==as.raw("0"))) {
      X1 <- fillin(X1)
    } else {
      X1 <- as(X1,"numeric")
    }
  }
  if(is(X2,"SnpMatrix")) {
    if(any(X2==as.raw("0"))) {
      X2 <- fillin(X2)
    } else {
      X2 <- as(X2,"numeric")
    }
  }
  X <- rbind(X1,X2)
  rows.drop <- apply(is.na(X),1,any)
  if(sum(rows.drop))
    X <- X[!rows.drop,]
  pcs <- prcomp(X, scale.=TRUE)
  vars <- pcs$sdev^2
  cvars <- cumsum(vars)/sum(vars)

  return(new("colocPCs",
             pcs=pcs$x,
             group=rep(c(1,2),times=c(nrow(X1),nrow(X2))),
             use=!rows.drop,
             vars=cvars))    
}

#'Functions to prepare principle component models for colocalisation testing
#'
#'Prepares models of response based on principal components of two datasets for
#'colocalisation testing.
#'
#'@title pcs.model
#'@param object A colocPCs object, result of \code{pcs.prepare()}.
#'@param group 1 or 2, indicating which group of samples to extract from
#'principal components matrix
#'@param Y Numeric phenotype vector, length equal to the number of samples from
#'the requested group
#'@param threshold The minimum number of principal components which captures at
#'least threshold proportion of the variance will be selected.  Simulations
#'suggest \code{threshold=0.8} is a good default value.
#'@param family Passed to \code{glm()} function.  \code{pcs.model} attempts to
#'guess, either "binomial" if \code{Y} contains only 0s and 1s, "gaussian"
#'otherwise.
##' @return \code{pcs.prepare} returns a \code{colocPCs} object, \code{pcs.model} returns a \code{glm} object.
##' @author Chris Wallace
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
#'  ## share a common causal variant, with the effect twice as strong for
#'  ## Y2 as Y1
#'set.seed(1)
#'  X1 <- matrix(rbinom(5000,1,0.4),ncol=10)
#'  Y1 <- rnorm(500,apply(X1[,1:2],1,sum),2)
#'  X2 <- matrix(rbinom(5000,1,0.6),ncol=10)
#'  Y2 <- rnorm(500,2*apply(X2[,1:2],1,sum),5)
#'  
#'  ## generate principal components object
#'  colnames(X1) <- colnames(X2) <- make.names(1:ncol(X1))
#'  pcs <- pcs.prepare(X1,X2)
#'
#'  ## generate glm objects
#'  m1 <- pcs.model(pcs, group=1, Y=Y1)
#'  m2 <- pcs.model(pcs, group=2, Y=Y2)
#'
#'  ## test colocalisation using PCs
#'  coloc.test(m1,m2,plot.coeff=FALSE,bayes=FALSE)
#'
#'@export
pcs.model <- function(object, group, Y, threshold=0.8, family=if(all(Y %in% c(0,1))) {"binomial"} else {"gaussian"}) {
  if(length(object@vars)<2)
    stop("require 2 or more principal components to test for proportionality")
  npc <- which(object@vars>threshold)[1]
  if(npc<2)
    npc <- 2
  cat("selecting",npc,"components out of",ncol(object@pcs),"to capture",object@vars[npc],"of total variance.\n")
  X <- object@pcs[ object@group[object@use]==group, 1:npc ]
  Y <- Y[ object@use[ object@group==group ] ]
  if(nrow(X) != length(Y))
    stop("length of Y not equal to the number of rows from group",group,"\n")
  snps <- colnames(object@pcs)[1:npc]
  f <- as.formula(paste("Y ~", paste(snps,collapse="+")))
  data <-  as.data.frame(cbind(Y,X))
  m <- glm(f,data=data,family=family)
  while(any(is.na(coefficients(m)))) {
    drop <- which((is.na(coefficients(m)))[-1])
    cat("dropping",length(drop),"colinear PCs",snps[drop],"\n")
    snps <- snps[-drop]
    f <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    m <- glm(f,data=data,family=family)
  }
  return(m)  
}

\name{coloc.bma}
\alias{coloc.bma}
\title{Wrapper to use colocalization testing within a Bayesian model averaging
structure.}
\usage{
  coloc.bma(df1, df2,
    snps = intersect(setdiff(colnames(df1), response1), setdiff(colnames(df2), response2)),
    response1 = "Y", response2 = "Y", family1 = "binomial",
    family2 = "binomial", bayes = !is.null(bayes.factor),
    thr = 0.01, nsnps = 2, n.approx = 1001,
    bayes.factor = NULL, plot.coeff = FALSE,
    r2.trim = 0.95, quiet = FALSE, ...)
}
\arguments{
  \item{df1,df2}{Each is a dataframe, containing response
  and potential explanatory variables for two independent
  datasets.}

  \item{snps}{The SNPs to consider as potential explanatory
  variables}

  \item{response1,response2}{The names of the response
  variables in \code{df1} and \code{df2} respectively}

  \item{family1,family2}{the error family for use in
  \code{glm}}

  \item{thr}{posterior probability threshold used to trim
  SNP list.  Only SNPs with a marginal posterior
  probability of inclusion greater than this with one or
  other trait will be included in the full BMA analysis}

  \item{nsnps}{number of SNPs required to model both
  traits.  The BMA analysis will average over all possible
  \code{nsnp} SNP models, subject to \code{thr} above.}

  \item{n.approx}{number of values at which to numerically
  approximate the posterior}

  \item{r2.trim}{for pairs SNPs with r2>\code{r2.trim},
  only one SNP will be retained.  This avoids numerical
  instability problems caused by including two highly
  correlated SNPs in the model.}

  \item{quiet}{suppress messages about how the model spaced
  is trimmed for BMA}

  \item{...}{other parameters passed to \code{coloc.test}}

  \item{bayes}{Logical, indicating whether to perform
  Bayesian inference for the coefficient of
  proportionality, eta.  If \code{bayes.factor} is
  supplied, Bayes factors are additionally computed for the
  specificed values.  This can add a little time as it
  requires numerical integration, so can be set to FALSE to
  save time in simulations, for example.}

  \item{bayes.factor}{Either a numeric vector, giving
  single value(s) of \code{eta} or a list of numeric
  vectors, each of length two and specifying ranges of eta
  which should be compared to each other.  Thus, the vector
  or list needs to have length at least two.}

  \item{plot.coeff}{\code{TRUE} if you want to generate a
  plot showing the coefficients from the two regressions
  together with confidence regions.}
}
\value{
  a \code{coloc} or \code{colocBayes} object
}
\description{
  Performs the colocalisation tests described in Plagnol et
  al (2009) and Wallace et al (2012).
}
\details{
  This is a test for proportionality of regression
  coefficients from two independent regressions.  Analysis
  can either be based on a profile likelihood approach,
  where the proportionality coefficient, \code{eta}, is
  replaced by its maximum likelihood value, and inference
  is based on a chisquare test (\code{p.value}), or taking
  a hybrid-Bayesian approach and integrating the p value
  over the posterior distribution of \code{eta}, which
  gives a posterior predictive p value.  The Bayesian
  approach can also be used to give a credible interval for
  \code{eta}.  See the references below for further
  details.
}
\examples{
## simulate covariate matrix (X) and continuous response vector (Y)
 ## for two populations/triats Y1 and Y2 depend equally on f1 and f2
 ## within each population, although their distributions differ between
 ## populations.  They are compatible with a null hypothesis that they
 ## share a common causal variant
set.seed(1)
 X1 <- matrix(rbinom(2000,1,0.4),ncol=4)
 Y1 <- rnorm(500,rowSums(X1[,1:2]),2)
 X2 <- matrix(rbinom(2000,1,0.6),ncol=4)
 Y2 <- rnorm(500,rowSums(X2[,1:2]),5)

 boxplot(list(Y1,Y2),names=c("Y1","Y2"))

 ## fit and store linear model objects
 colnames(X1) <- colnames(X2) <- sprintf("f\%s",1:ncol(X1))
 summary(lm1 <- lm(Y1~f1+f2+f3+f4,data=as.data.frame(X1)))
 summary(lm2 <- lm(Y2~f1+f2+f3+f4,data=as.data.frame(X2)))


## test colocalisation using bma
df1=cbind(Y1=Y1,X1)
df2=cbind(Y2=Y2,X2)

coloc.bma( df1, df2, snps=colnames(X1), response1="Y1", response2="Y2",
family1="gaussian", family2="gaussian",
nsnps=2,bayes.factor=c(1,2,3),plot.coeff=TRUE)
}
\author{
  Chris Wallace
}
\references{
  Wallace et al (2012).  Statistical colocalisation of
  monocyte gene expression and genetic risk variants for
  type 1 diabetes.  Hum Mol Genet 21:2815-2824.
  \url{http://europepmc.org/abstract/MED/22403184}

  Plagnol et al (2009).  Statistical independence of the
  colocalized association signals for type 1 diabetes and
  RPS26 gene expression on chromosome 12q13. Biostatistics
  10:327-34.
  \url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
}


\name{coloc.bayes.tag}
\alias{coloc.bayes.tag}
\title{bayesian colocalisation; two traits; with tagging}
\usage{
  coloc.bayes.tag(df1,
    snps = setdiff(colnames(df1), response),
    response = "Y", priors = list(c(1, 1, 1, 1, 1)),
    r2.trim = 0.99, pp.thr = 0.005, quiet = TRUE)
}
\arguments{
  \item{df1}{A dataframe, containing response and potential
  explanatory variables for the dataset.}

  \item{snps}{The SNPs to consider as potential explanatory
  variables}

  \item{response}{The name of the response variable in
  \code{df1}}

  \item{priors}{A list of priors over the hypotheses}

  \item{pp.thr}{posterior probability threshold used to
  trim SNP list.  Only SNPs with a marginal posterior
  probability of inclusion greater than this with one or
  other trait will be included in the full BMA analysis}

  \item{r2.trim}{If a pairs of SNPs has r2> \code{r2.trim},
  they are put in the same tag}

  \item{quiet}{suppress messages about how the model spaced
  is trimmed for BMA}
}
\value{
  a list of posterior probabilities that each tag is
  causitive to both traits, the tag names, and the
  corresponding SNPs
}
\description{
  Colocation for two traits with a common control Merges
  SNPs with high r2 into tags prior to analysis Generates
  bayes factors for each plausible one tag model
}
\author{
  Mary Fortune
}


\name{singlesnp.twotrait}
\alias{singlesnp.twotrait}
\title{finding the minimum single SNP p-value for each of the two traits}
\usage{
  singlesnp.twotrait(df1, response = "Y",
    snps = setdiff(colnames(df1), response),
    r2.trim = 0.95, quiet = TRUE)
}
\arguments{
  \item{df1}{A dataframe, containing response and potential
  explanatory variables for the dataset.}

  \item{snps}{The SNPs to consider as potential explanatory
  variables}

  \item{response}{The name of the response variable in
  \code{df1}}

  \item{r2.trim}{for pairs SNPs with r2> \code{r2.trim},
  only one SNP will be retained.  This avoids numerical
  instability problems caused by including two highly
  correlated SNPs in the model.}

  \item{quiet}{suppress messages}
}
\value{
  the minimum single SNP p value for each trait
}
\description{
  finding the minimum single SNP p-value for each of the
  two traits
}
\author{
  Mary Fortune
}


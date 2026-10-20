# Coloc data through Bayes factors

Colocalise two datasets represented by Bayes factors

## Usage

``` r
coloc.bf_bf(
  bf1,
  bf2,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 5e-06,
  overlap.min = 0.5,
  trim_by_posterior = TRUE,
  prior_weights1 = NULL,
  prior_weights2 = NULL
)
```

## Arguments

- bf1:

  named vector of log BF, or matrix of BF with colnames (cols=snps,
  rows=signals)

- bf2:

  named vector of log BF, or matrix of BF with colnames (cols=snps,
  rows=signals)

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

- overlap.min:

  see trim_by_posterior

- trim_by_posterior:

  it is important that the signals to be colocalised are covered by
  adequate numbers of snps in both datasets. If TRUE, signals for which
  snps in common do not capture least overlap.min proportion of their
  posteriors support are dropped and colocalisation not attempted.

- prior_weights1:

  Non-negative weights for the prior probability a SNP is associated
  with trait 1

- prior_weights2:

  Non-negative weights for the prior probability a SNP is asscoiated
  with trait 2

## Value

coloc.signals style result

## Details

This is the workhorse behind many coloc functions

## Author

Chris Wallace

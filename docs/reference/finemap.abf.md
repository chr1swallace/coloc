# Bayesian finemapping analysis

Bayesian finemapping analysis

## Usage

``` r
finemap.abf(dataset, p1 = 1e-04, prior_weights = NULL)
```

## Arguments

- dataset:

  a list with specifically named elements defining the dataset to be
  analysed. See
  [`check_dataset`](https://chr1swallace.github.io/coloc/reference/check_dataset.md)
  for details.

- p1:

  prior probability a SNP is associated with the trait 1, default 1e-4

- prior_weights:

  Non-negative weights for the prior probability a SNP is causal

## Value

a `data.frame`:

- an annotated version of the input data containing log Approximate
  Bayes Factors and intermediate calculations, and the posterior
  probability of the SNP being causal

## Details

This function calculates posterior probabilities of different causal
variant for a single trait.

If regression coefficients and variances are available, it calculates
Bayes factors for association at each SNP. If only p values are
available, it uses an approximation that depends on the SNP's MAF and
ignores any uncertainty in imputation. Regression coefficients should be
used if available.

## Author

Chris Wallace

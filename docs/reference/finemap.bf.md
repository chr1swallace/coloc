# Finemap data through Bayes factors

Finemap one dataset represented by Bayes factors

## Usage

``` r
finemap.bf(bf1, p1 = 1e-04)
```

## Arguments

- bf1:

  named vector of log BF, or matrix of log BF with colnames (cols=snps,
  rows=signals)

- p1:

  prior probability a SNP is associated with the trait 1, default 1e-4

## Value

finemap.signals style result

## Details

This is the workhorse behind many finemap functions

## Author

Chris Wallace

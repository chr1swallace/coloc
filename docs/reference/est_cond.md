# generate conditional summary stats

Internal helper function for est_all_cond

## Usage

``` r
est_cond(x, LD, YY, sigsnps, xtx = NULL)
```

## Arguments

- x:

  coloc dataset

- LD:

  named matrix of r

- YY:

  sum((Y-Ybar)^2)

- sigsnps:

  names of snps to jointly condition on

- xtx:

  optional, matrix X'X where X is the genotype matrix. If not available,
  will be estimated from LD, MAF, beta and sample size (the last three
  should be part of the coloc dataset)

## Value

data.table giving snp, beta and varbeta on remaining snps after
conditioning

## Author

Chris Wallace

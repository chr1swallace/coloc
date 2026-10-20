# logbf 2 pp

generic convenience function to convert logbf matrix to PP matrix

## Usage

``` r
logbf_to_pp(bf, pi, last_is_null)
```

## Arguments

- bf:

  an L by p or p+1 matrix of log Bayes factors

- pi:

  *either* a scalar representing the prior probability for any snp to be
  causal, *or* a full vector of per snp / null prior probabilities

- last_is_null:

  TRUE if last value of the bf vector or last column of a bf matrix
  relates to the null hypothesis of no association. This is standard for
  SuSiE results, but may not be for BF constructed in other ways.

## Value

matrix of posterior probabilities, same dimensions as bf

## Author

Chris Wallace

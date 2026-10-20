# Internal function, approx.bf.estimates

Internal function, approx.bf.estimates

## Usage

``` r
approx.bf.estimates(
  z,
  V,
  type,
  suffix = NULL,
  sdY = 1,
  effect_priors = c(quant = 0.15, cc = 0.2)
)
```

## Arguments

- z:

  normal deviate associated with regression coefficient and its variance

- V:

  its variance

- type:

  "quant" or "cc"

- suffix:

  suffix to append to column names of returned data.frame

- sdY:

  standard deviation of the trait. If not supplied, will be estimated.

- effect_priors:

  named vector with variance of prior effect sizes for quantitative and
  case-control traits. don't change unless you know what you are doing!

## Value

data.frame containing lABF and intermediate calculations

## Details

Calculate approximate Bayes Factors using supplied variance of the
regression coefficients

## Author

Vincent Plagnol, Chris Wallace

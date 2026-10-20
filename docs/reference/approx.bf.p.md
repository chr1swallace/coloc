# Internal function, approx.bf.p

Internal function, approx.bf.p

## Usage

``` r
approx.bf.p(p, f, type, N, s, suffix = NULL)
```

## Arguments

- p:

  p value

- f:

  MAF

- type:

  "quant" or "cc"

- N:

  sample size

- s:

  proportion of samples that are cases, ignored if type=="quant"

- suffix:

  suffix to append to column names of returned data.frame

## Value

data.frame containing lABF and intermediate calculations

## Details

Calculate approximate Bayes Factors

## Author

Claudia Giambartolomei, Chris Wallace

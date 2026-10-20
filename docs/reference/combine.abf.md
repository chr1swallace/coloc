# combine.abf

Internal function, calculate posterior probabilities for configurations,
given logABFs for each SNP and prior probs

## Usage

``` r
combine.abf(l1, l2, p1, p2, p12, quiet = FALSE)
```

## Arguments

- l1:

  merged.df\$lABF.df1

- l2:

  merged.df\$lABF.df2

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

- quiet:

  don't print posterior summary if TRUE. default=FALSE

## Value

named numeric vector of posterior probabilities

## Author

Claudia Giambartolomei, Chris Wallace

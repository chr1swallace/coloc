# Estimate trait variance, internal function

Estimate trait standard deviation given vectors of variance of
coefficients, MAF and sample size

## Usage

``` r
sdY.est(vbeta, maf, n)
```

## Arguments

- vbeta:

  vector of variance of coefficients

- maf:

  vector of MAF (same length as vbeta)

- n:

  sample size

## Value

estimated standard deviation of Y

## Details

Estimate is based on var(beta-hat) = var(Y) / (n \* var(X)) var(X) =
2*maf*(1-maf) so we can estimate var(Y) by regressing n\*var(X) against
1/var(beta)

## Author

Chris Wallace

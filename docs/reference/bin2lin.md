# binomial to linear regression conversion

Convert binomial to linear regression

## Usage

``` r
bin2lin(D, doplot = FALSE)
```

## Arguments

- D:

  standard format coloc dataset

- doplot:

  plot results if TRUE - useful for debugging

## Value

D, with original beta and varbeta in beta.bin, varbeta.bin, and beta and
varbeta updated to linear estimates

## Details

Estimate beta and varbeta if a linear regression had been run on a
binary outcome, given log OR and their variance + MAF in controls

sets beta = cov(x,y)/var(x) varbeta = (var(y)/var(x) -
cov(x,y)^2/var(x)^2)/N

## Author

Chris Wallace

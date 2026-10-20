# find the next most significant SNP, conditioning on a list of sigsnps

Internal helper function for finemap.signals

## Usage

``` r
map_cond(D, LD, YY, sigsnps = NULL)
```

## Arguments

- D:

  dataset in standard coloc format

- LD:

  named matrix of r

- YY:

  sum(y^2)

- sigsnps:

  names of snps to mask

## Value

named numeric - Z score named by snp

## Author

Chris Wallace

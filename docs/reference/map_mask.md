# find the next most significant SNP, masking a list of sigsnps

Internal helper function for finemap.signals

## Usage

``` r
map_mask(D, LD, r2thr = 0.01, sigsnps = NULL)
```

## Arguments

- D:

  dataset in standard coloc format

- LD:

  named matrix of r

- r2thr:

  mask all snps with r2 \> r2thr with any in sigsnps

- sigsnps:

  names of snps to mask

## Value

named numeric - Z score named by snp

## Author

Chris Wallace

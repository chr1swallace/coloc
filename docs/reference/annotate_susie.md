# annotate susie_rss output for use with coloc_susie

coloc functions need to be able to link summary stats from two different
datasets and they do this through snp identifiers. This function takes
the output of susie_rss() and adds snp identifiers. It is entirely the
user's responsibility to ensure snp identifiers are in the correct
order, coloc cannot make any sanity checks.

## Usage

``` r
annotate_susie(res, snp, LD)
```

## Arguments

- res:

  output of susie_rss()

- snp:

  vector of snp identifiers

- LD:

  matrix of LD (r) between snps in snp identifiers. Columns, rows should
  be named by a string that exists in the vector snp

## Value

res with column names added to some components

## Details

Note: this annotation step is not needed if you use runsusie() - this is
only required if you use the susieR functions directly

## Author

Chris Wallace

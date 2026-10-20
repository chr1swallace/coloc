# run coloc using susie to detect separate signals

coloc for susie output + a separate BF matrix

## Usage

``` r
coloc.susie_bf(
  dataset1,
  bf2,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 5e-06,
  susie.args = list(),
  ...
)
```

## Arguments

- dataset1:

  a list with specifically named elements defining the dataset to be
  analysed. See
  [`check_dataset`](https://chr1swallace.github.io/coloc/reference/check_dataset.md)
  for details.

- bf2:

  named vector of log BF, names are snp ids and will be matched to
  column names of susie object's alpha

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

- susie.args:

  named list of arguments to be passed to susieR::susie_rss()

- ...:

  other arguments passed to
  [coloc.bf_bf](https://chr1swallace.github.io/coloc/reference/coloc.bf_bf.md),
  in particular prior values for causal association with one trait (p1,
  p2) or both (p12) and prior weights.

## Value

coloc.signals style result

## Author

Chris Wallace

# Coloc with multiple signals per trait

New coloc function, builds on coloc.abf() by allowing for multiple
independent causal variants per trait through conditioning or masking.

## Usage

``` r
coloc.signals(
  dataset1,
  dataset2,
  MAF = NULL,
  LD = NULL,
  method = c("single", "cond", "mask"),
  mode = c("iterative", "allbutone"),
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = NULL,
  maxhits = 3,
  r2thr = 0.01,
  pthr = 1e-06
)
```

## Arguments

- dataset1:

  a list with specifically named elements defining the dataset to be
  analysed. See
  [`check_dataset`](https://chr1swallace.github.io/coloc/reference/check_dataset.md)
  for details.

- dataset2:

  as above, for dataset 2

- MAF:

  Common minor allele frequency vector to be used for both dataset1 and
  dataset2, a shorthand for supplying the same vector as parts of both
  datasets

- LD:

  required if method="cond". matrix of genotype *correlation* (ie r, not
  r^2) between SNPs. If dataset1 and dataset2 may have different LD, you
  can instead add LD=LD1 to the list of dataset1 and a different LD
  matrix for dataset2

- method:

  default "" means do no conditioning, should return similar to
  coloc.abf. if method="cond", then use conditioning to coloc multiple
  signals. if method="mask", use masking to coloc multiple signals. if
  different datasets need different methods (eg LD is only available for
  one of them) you can set method on a per-dataset basis by adding
  method="..." to the list for that dataset.

- mode:

  "iterative" or "allbutone". Easiest understood with an example.
  Suppose there are 3 signal SNPs detected for trait 1, A, B, C and only
  one for trait 2, D.

      Under "iterative" mode, 3 coloc will be performed:
      * trait 1 - trait 2
      * trait 1 conditioned on A - trait 2
      * trait 1 conditioned on A+B - trait 2

      Under "allbutone" mode, they would be
      * trait 1 conditioned on B+C - trait 2
      * trait 1 conditioned on A+C - trait 2
      * trait 1 conditioned on A+B - trait 2

      Only iterative mode is supported for method="mask".

      The allbutone mode is optimal if the signals are known with
      certainty (which they never are), because it allows each
      signal to be tested without influence of the others.  When
      there is uncertainty, it may make sense to use iterative mode,
      because the strongest signals aren't affected by conditioning
      incorrectly on weaker secondary and less certain signals.

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

- maxhits:

  maximum number of levels to condition/mask

- r2thr:

  if masking, the threshold on r2 should be used to call two signals
  independent. our experience is that this needs to be set low to avoid
  double calling the same strong signal.

- pthr:

  if masking or conditioning, what p value threshold to call a secondary
  hit "significant"

## Value

data.table of coloc results, one row per pair of lead snps detected in
each dataset

## Author

Chris Wallace

# Finemap multiple signals in a single dataset

This is an analogue to finemap.abf, adapted to find multiple signals
where they exist, via conditioning or masking - ie a stepwise procedure

## Usage

``` r
finemap.signals(
  D,
  LD = D$LD,
  method = c("single", "mask", "cond"),
  r2thr = 0.01,
  sigsnps = NULL,
  pthr = 1e-06,
  maxhits = 3,
  return.pp = FALSE
)
```

## Arguments

- D:

  list of summary stats for a single disease, see
  [check_dataset](https://chr1swallace.github.io/coloc/reference/check_dataset.md)

- LD:

  matrix of signed r values (not rsq!) giving correlation between SNPs

- method:

  if method="cond", then use conditioning to coloc multiple signals. The
  default is mask - this is less powerful, but safer because it does not
  assume that the LD matrix is properly allelically aligned to estimated
  effect

- r2thr:

  if mask==TRUE, all snps will be masked with r2 \> r2thr with any
  sigsnps. Otherwise ignored

- sigsnps:

  SNPs already deemed significant, to condition on or mask, expressed as
  a numeric vector, whose *names* are the snp names

- pthr:

  when p \> pthr, stop successive searching

- maxhits:

  maximum depth of conditioning. procedure will stop if p \> pthr OR
  abs(z)\<zthr OR maxhits hits have been found.

- return.pp:

  if FALSE (default), just return the hits. Otherwise return vectors of
  PP

- mask:

  use masking if TRUE, otherwise conditioning. defaults to TRUE

## Value

list of successively significant fine mapped SNPs, named by the SNPs

## Author

Chris Wallace

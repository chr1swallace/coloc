# check_dataset

Check coloc dataset inputs for errors

## Usage

``` r
check_dataset(d, suffix = "", req = c("type", "snp"), warn.minp = 1e-06)

check.dataset(...)
```

## Arguments

- d:

  dataset to check

- suffix:

  string to identify which dataset (1 or 2)

- req:

  names of elements that must be present

- warn.minp:

  print warning if no p value \< warn.minp

- ...:

  arguments passed to check_dataset()

## Value

NULL if no errors found

## Details

A coloc dataset is a list, containing a mixture of vectors capturing
quantities that vary between snps (these vectors must all have equal
length) and scalars capturing quantities that describe the dataset.

Coloc is flexible, requiring perhaps only p values, or z scores, or
effect estimates and standard errors, but with this flexibility, also
comes difficulties describing exactly the combinations of items
required.

Required vectors are some subset of

- beta:

  regression coefficient for each SNP from dataset 1

- varbeta:

  variance of beta

- pvalues:

  P-values for each SNP in dataset 1

- MAF:

  minor allele frequency of the variants

- snp:

  a character vector of snp ids, optional. It will be used to merge
  dataset1 and dataset2 and will be retained in the results.

Preferably, give `beta` and `varbeta`. But if these are not available,
sufficient statistics can be approximated from `pvalues` and `MAF`.

Required scalars are some subset of

- N:

  Number of samples in dataset 1

- type:

  the type of data in dataset 1 - either "quant" or "cc" to denote
  quantitative or case-control

- s:

  for a case control dataset, the proportion of samples in dataset 1
  that are cases

- sdY:

  for a quantitative trait, the population standard deviation of the
  trait. if not given, it can be estimated from the vectors of varbeta
  and MAF

You must always give `type`. Then,

- if `type`=="cc":

  `s`

- if `type`=="quant" and `sdY` known:

  `sdY`

- if beta, varbeta not known:

  `N`

If `sdY` is unknown, it will be approximated, and this will require

- summary data to estimate `sdY`:

  `beta`, `varbeta`, `N`, `MAF`

Optional vectors are

- position:

  a vector of snp positions, required for `plot_dataset`

`check_dataset` calls stop() unless a series of expectations on dataset
input format are met

This is a helper function for use by other coloc functions, but you can
use it directly to check the format of a dataset to be supplied to
coloc.abf(), coloc.signals(), finemap.abf(), or finemap.signals().

## Author

Chris Wallace

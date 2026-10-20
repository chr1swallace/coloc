# Bayesian colocalisation analysis with detailed output

Bayesian colocalisation analysis, detailed output

## Usage

``` r
coloc.detail(
  dataset1,
  dataset2,
  MAF = NULL,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 1e-05
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

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

## Value

a list of three `data.tables`s:

- summary is a vector giving the number of SNPs analysed, and the
  posterior probabilities of H0 (no causal variant), H1 (causal variant
  for trait 1 only), H2 (causal variant for trait 2 only), H3 (two
  distinct causal variants) and H4 (one common causal variant)

- df is an annotated version of the input data containing log
  Approximate Bayes Factors and intermediate calculations, and the
  posterior probability SNP.PP.H4 of the SNP being causal for the shared
  signal

- df3 is the same for all 2 SNP H3 models

## Details

This function replicates coloc.abf, but outputs more detail for further
processing using coloc.process

Intended to be called internally by coloc.signals

## See also

[`coloc.process`](https://chr1swallace.github.io/coloc/reference/coloc.process.md),
[`coloc.abf`](https://chr1swallace.github.io/coloc/reference/coloc.abf.md)

## Author

Chris Wallace

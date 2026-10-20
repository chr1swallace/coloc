# Post process a coloc.details result using masking

Internal helper function

## Usage

``` r
coloc.process(
  obj,
  hits1 = NULL,
  hits2 = NULL,
  LD = NULL,
  r2thr = 0.01,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 1e-06,
  LD1 = LD,
  LD2 = LD,
  mode = c("iterative", "allbutone")
)
```

## Arguments

- obj:

  object returned by coloc.detail()

- hits1:

  lead snps for trait 1. If length \> 1, will use masking

- hits2:

  lead snps for trait 2. If length \> 1, will use masking

- LD:

  named LD matrix (for masking)

- r2thr:

  r2 threshold at which to mask

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

- LD1:

  named LD matrix (for masking) for trait 1 only

- LD2:

  named LD matrix (for masking) for trait 2 only

- mode:

  either "iterative" (default) - successively condition on signals or
  "allbutone" - find all putative signals and condition on all but one
  of them in each analysis

## Value

data.table of coloc results

## Author

Chris Wallace

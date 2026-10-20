# convert CLPP to coloc posteriors

convert CLPP to coloc posteriors

## Usage

``` r
clPP_to_coloc(
  clPP,
  alpha1,
  alpha2,
  p1 = 1e-04,
  p2 = 1e-04,
  p12 = 5e-06,
  h4_only = FALSE
)
```

## Arguments

- clPP:

  CLPP value

- alpha1:

  size of credible set 1 (the sum of PP associated with variants in the
  credible set)

- alpha2:

  size of credible set 2

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

- h4_only:

  set TRUE if you only want PP.H4 returned

## Value

a named vector containing the posterior probabilities of H3 and H4

## Author

Jeffrey Pullin, Chris Wallace

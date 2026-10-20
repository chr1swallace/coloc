# process.dataset

Internal function, process each dataset list for coloc.abf.

## Usage

``` r
process.dataset(d, suffix, ...)
```

## Arguments

- d:

  list

- suffix:

  "df1" or "df2"

- ...:

  used to pass parameters to approx.bf.estimates, in particular the
  effect_priors parameter

## Value

data.frame with log(abf) or log(bf)

## Details

Made public for another package to use, but not intended for users to
use.

## Author

Chris Wallace

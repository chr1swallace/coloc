# trim a dataset to only peak(s)

tries to be smart about detecting the interesting subregion to
finemap/coloc.

## Usage

``` r
findpeaks(d, maxz = 4, maxr2 = 0.1, do.plot = FALSE)
```

## Arguments

- d:

  a coloc dataset

- maxz:

  keep all snps between the leftmost and rightmost snp with \|z\| \>
  maxz

- maxr2:

  expand window to keep all snps between snps with r2 \> maxr2 with the
  left/rightmost snps defined by the maxz threshold

- do.plot:

  if TRUE, plot dataset + boundaries

## Value

logical vector of length d\$position indicating which snps to keep

## Details

Differs from findends by finding multiple separate regions if there are
multiple peaks

## See also

findends

## Author

Chris Wallace

# Simulated data to use in testing and vignettes in the coloc package

Simulated data to use in testing and vignettes in the coloc package

## Usage

``` r
data(coloc_test_data)
```

## Format

A four of two coloc-style datasets. Elements D1 and D2 have a single
shared causal variant, and 50 SNPs. Elements D3 and D4 have 100 SNPs,
one shared causal variant, and one variant unique to D3. Use these as
examples of what a coloc-style dataset for a quantitative trait should
look like.

## Examples

``` r
data(coloc_test_data)
names(coloc_test_data)
#> [1] "D1"      "D2"      "D3"      "D4"      "causals"
str(coloc_test_data$D1)
#> List of 9
#>  $ beta    : num [1:500] 0.337 0.211 0.257 0.267 0.247 ...
#>  $ varbeta : num [1:500] 0.01634 0.00532 0.00748 0.01339 0.00664 ...
#>  $ N       : int 1000
#>  $ sdY     : num 1.1
#>  $ type    : chr "quant"
#>  $ MAF     : num [1:500] 0.031 0.166 0.0925 0.0405 0.118 ...
#>  $ LD      : num [1:500, 1:500] 1 0.365 0.54 0.851 0.463 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:500] "s1" "s2" "s3" "s4" ...
#>   .. ..$ : chr [1:500] "s1" "s2" "s3" "s4" ...
#>  $ snp     : chr [1:500] "s1" "s2" "s3" "s4" ...
#>  $ position: int [1:500] 1 2 3 4 5 6 7 8 9 10 ...
check_dataset(coloc_test_data$D1) # should return NULL if data structure is ok
#> NULL
```

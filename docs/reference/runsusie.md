# Run susie on a single coloc-structured dataset

run susie_rss storing some additional information for coloc

## Usage

``` r
runsusie(
  d,
  suffix = 1,
  maxit = 100,
  repeat_until_convergence = TRUE,
  s_init = NULL,
  ...
)
```

## Arguments

- d:

  coloc dataset, must include LD (signed correlation matrix) and N
  (sample size)

- suffix:

  suffix label that will be printed with any error messages

- maxit:

  maximum number of iterations for the first run of susie_rss(). If
  susie_rss() does not report convergence, runs will be extended
  assuming repeat_until_convergence=TRUE. Most users will not need to
  change this default.

- repeat_until_convergence:

  keep running until susie_rss() indicates convergence. Default TRUE. If
  FALSE, susie_rss() will run with maxit iterations, and if not
  converged, runsusie() will error. Most users will not need to change
  this default.

- s_init:

  used internally to extend runs that haven't converged. don't use.

- ...:

  arguments passed to susie_rss. In particular, if you want to match
  some coloc defaults, set

  - prior_variance=0.2^2 (if a case-control trait) or (0.15/sd(Y))^2 if
    a quantitative trait

  - estimate_prior_variance=FALSE

  otherwise susie_rss will estimate the prior variance itself

## Value

results of a susie_rss run, with some added dimnames

## Author

Chris Wallace

## Examples

``` r
library(coloc)
data(coloc_test_data)
result=runsusie(coloc_test_data$D1)
#> running max iterations: 100
#>  converged: TRUE
summary(result)
#> 
#> Variables in credible sets:
#> 
#>  variable variable_prob cs
#>       105     0.9593968  1
#> 
#> Credible sets summary:
#> 
#>  cs cs_log10bf cs_avg_r2 cs_min_r2 variable
#>   1    10.2121         1         1      105
```

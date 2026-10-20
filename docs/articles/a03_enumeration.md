# Coloc: under a single causal variant assumption

## Fine mapping under a single causal variant assumption

The Approximate Bayes Factor colocalisation analysis described in the
next section essentially works by fine mapping each trait under a single
causal variant assumption and then integrating over those two posterior
distributions to calculate probabilities that those variants are shared.
Of course, this means we can look at each trait on its own quite simply,
and we can do that with the function `finemap.abf`.

First we load some simulated data. See [the data
vignette](https://chr1swallace.github.io/coloc/articles/a02_data.md) to
understand how to format your datasets.

``` r

library(coloc)
```

    ## This is coloc version 7.0.2

``` r

data(coloc_test_data)
attach(coloc_test_data)
```

Then we analyse the statistics from a single study, asking about the
evidence that each SNP in turn is solely causal for any association
signal we see. As we might expect, that evidence is maximised at the SNP
with the smallest p value

``` r

plot_dataset(D1)
```

![](a03_enumeration_files/figure-html/unnamed-chunk-2-1.png)

``` r

my.res <- finemap.abf(dataset=D1)
my.res[21:30,]
```

    ##             V.       z.        r.       lABF. snp position prior       SNP.PP
    ## 21 0.006158578 3.105184 0.8155100  3.08656163 s21       21 1e-04 1.026731e-11
    ## 22 0.006047354 3.010816 0.8182363  2.85613791 s22       22 1e-04 8.154269e-12
    ## 23 0.017869410 1.247449 0.6037164  0.00691776 s23       23 1e-04 4.720461e-13
    ## 24 0.005427008 3.746694 0.8337827  4.95497248 s24       24 1e-04 6.651159e-11
    ## 25 0.005699008 3.168342 0.8268940  3.27341820 s25       25 1e-04 1.237678e-11
    ## 26 0.004712835 5.305229 0.8524283 11.03926980 s26       26 1e-04 2.919269e-08
    ## 27 0.006988649 2.570785 0.7957235  1.83530152 s27       27 1e-04 2.937930e-12
    ## 28 0.004730114 5.483924 0.8519673 11.85562507 s28       28 1e-04 6.604085e-08
    ## 29 0.004929736 4.771115 0.8466779  8.69908643 s29       29 1e-04 2.811547e-09
    ## 30 0.005566298 3.176653 0.8302407  3.30234473 s30       30 1e-04 1.274002e-11

The `SNP.PP` column shows the posterior probability that exactly that
SNP is causal. Note the last line in this data.frame does not correspond
to a SNP, but to the null model, that no SNP is causal.

``` r

tail(my.res,3)
```

Finally, if you do have full genotype data as here, while this is a fast
method for fine mapping, it can be sensible to consider multiple causal
variant models too. One package that allows you to do this is
[GUESSFM](https://github.com/chr1swallace/GUESSFM), described
in^([5](#fn.5))

## (Approximate) Bayes Factor colocalisation analyses

### Introduction

The idea behind the ABF analysis is that the association of each trait
with SNPs in a region may be summarised by a vector of 0s and at most a
single 1, with the 1 indicating the causal SNP (so, assuming a single
causal SNP for each trait). The posterior probability of each possible
configuration can be calculated and so, crucially, can the posterior
probabilities that the traits share their configurations. This allows us
to estimate the support for the following cases:

- $`H_0`$: neither trait has a genetic association in the region
- $`H_1`$: only trait 1 has a genetic association in the region
- $`H_2`$: only trait 2 has a genetic association in the region
- $`H_3`$: both traits are associated, but with different causal
  variants
- $`H_4`$: both traits are associated and share a single causal variant

### The basic `coloc.abf` function

The function `coloc.abf` is ideally suited to the case when only summary
data are available.

``` r

my.res <- coloc.abf(dataset1=D1,
                    dataset2=D2)
```

    ## PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
    ##  1.38e-18  2.94e-10  8.59e-12  8.34e-04  9.99e-01 
    ## [1] "PP abf for shared variant: 99.9%"

``` r

print(my.res) 
```

    ## Coloc analysis of trait 1, trait 2

    ## 
    ## SNP Priors

    ##    p1    p2   p12 
    ## 1e-04 1e-04 1e-05

    ## 
    ## Hypothesis Priors

    ##        H0   H1   H2       H3    H4
    ##  0.892505 0.05 0.05 0.002495 0.005

    ## 
    ## Posterior

    ##        nsnps           H0           H1           H2           H3           H4 
    ## 5.000000e+02 1.377000e-18 2.937336e-10 8.593226e-12 8.338917e-04 9.991661e-01

Note that if you do find strong evidence for H4, we can extract the
posterior probabilities for each SNP to be causal *conditional on H4
being true*. This is part of the calculation required by coloc, and
contained in the column SNP.PP.H4 in the “results” element of the
returned list. So we can extract the more likely causal variants by

``` r

subset(my.res$results,SNP.PP.H4>0.01)
```

    ##    snp position       V.df1    z.df1     r.df1 lABF.df1       V.df2    z.df2
    ## 6 s103      103 0.004846592 7.780900 0.8488730 24.75159 0.004102488 7.560002
    ## 8 s105      105 0.004808487 8.302463 0.8498828 28.34342 0.004100753 7.699676
    ##       r.df2 lABF.df2 internal.sum.lABF  SNP.PP.H4
    ## 6 0.8513388 23.37551           48.1271 0.01098066
    ## 8 0.8513923 24.28418           52.6276 0.98894507

or the 95% credible set by

``` r

o <- order(my.res$results$SNP.PP.H4,decreasing=TRUE)
cs <- cumsum(my.res$results$SNP.PP.H4[o])
w <- which(cs > 0.95)[1]
my.res$results[o,][1:w,]$snp
```

    ## [1] "s105"

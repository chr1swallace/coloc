coloc
=====



The coloc package can be used to perform genetic colocalisation
analysis of two potentially related phenotypes, to ask whether they
share common genetic causal variant(s) in a given region. 


# version 4

This is an updated version of coloc.  I have tested it, but there may be bugs. Please test it, and let me know whether it works or not (both kinds of feedback useful!).  

It is not yet on CRAN. To install the new version, do
```
if(!require("remotes"))
   install.packages("remotes") # if necessary
library(remotes)
install_github("chr1swallace/coloc")
```


# Background

The new ideas are described in 
> [Wallace C (2020) Eliciting priors and relaxing the single causal variant assumption in colocalisation analyses. PLOS Genetics 16(4): e1008720](https://doi.org/10.1371/journal.pgen.1008720)

For usage, please see the vignette at https://chr1swallace.github.io/coloc

Key previous references are:
- original propostion of proportional colocalisation [Plagnol et al (2009)](http://www.ncbi.nlm.nih.gov/pubmed/19039033)
- proportional colocalisation with type 1 error rate control [Wallace et al (2013)](http://onlinelibrary.wiley.com/doi/10.1002/gepi.21765/abstract)
- colocalisation by enumerating all the possible causal SNP configurations between two traits, assuming at most one causal variant per trait [Giambartolomei et al (2013)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)

[![Build Status](https://travis-ci.org/chr1swallace/coloc.svg?branch=master)](https://travis-ci.org/chr1swallace/coloc)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/coloc)](https://cran.r-project.org/package=coloc)


# Frequently Asked Questions

- [If I understand correctly, coloc.abf() can be run with correlated variants, that is, no prerequisite for taking through LD pruning/clumping is required. Am I correct in my understanding ?](#user-content-if-i-understand-correctly-coloc.abf-can-be-run-with-correlated-variants-that-is-no-prerequisite-for-taking-through-ld-pruning-clumping-is-required-am-i-correct-in-my-understanding)
- [Assume I identify a sentinel variant for a block of genome, can I do a comparison with just one variant using coloc.abf()?](#user-content-assume-i-identify-a-sentinel-variant-for-a-block-of-genome-can-i-do-a-comparison-with-just-one-variant-using-coloc.abf)
- [Can the process of identifying colocalized variants be carried out genome wide or is it meant to be done in defined small regions?](#user-content-can-the-process-of-identifying-colocalized-variants-be-carried-out-genome-wide-or-is-it-meant-to-be-done-in-defined-small-regions)

## If I understand correctly, coloc.abf() can be run with correlated variants, that is, no prerequisite for taking through LD pruning/clumping is required. Am I correct in my understanding ?

Yes, coloc.abf() and coloc.signals() assume they are given a dense map of all SNPs in a region that could be causal.   Do not prune and clump.

## Assume I identify a sentinel variant for a block of genome, can I do a comparison with just one variant using coloc.abf()?

No, coloc.abf() and coloc.signals() assume they are given a dense map of all SNPs in a region that could be causal. This means you need to give all SNPs in a region. You can imagine they ask whether the patterns "match" across this region of SNPs, and a single variant does not represent a pattern. 

## Can the process of identifying colocalized variants be carried out genome wide or is it meant to be done in defined small regions?

You need to break the genome into smaller regions, within which it is reasonable to assume there is at most one (coloc.abf) or a small number (coloc.signals) of causal variants per trait.  One way to do this is to use the boundaries defined by recombination hotspots, proxied by [this map](https://bitbucket.org/nygcresearch/ldetect-data/src/master/) created by [lddetect](https://academic.oup.com/bioinformatics/article/32/2/283/1743626).

# Notes to self

### Note to self: to generate vignettes:
```
cp vignettes/colocqq-tests.R.tospin vignettes/colocqq-tests.R && Rscript -e 'knitr::spin("vignettes/colocqq-tests.R",knit=FALSE); devtools::build_vignettes()'
```

### Note to self: to generate website:
https://chr1swallace.github.io/coloc/
```
Rscript -e "pkgdown::build_site()"
```

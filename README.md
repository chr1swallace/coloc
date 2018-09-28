coloc
=====

Repo for the R package coloc

For usage, see the vignette at https://chr1swallace.github.io/coloc

[![Build Status](https://travis-ci.org/chr1swallace/coloc.svg?branch=master)](https://travis-ci.org/chr1swallace/coloc)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/coloc)](https://cran.r-project.org/package=coloc)

### To generate vignettes:
```
cp vignettes/colocqq-tests-tospin.R vignettes/colocqq-tests.R && Rscript -e 'knitr::spin("vignettes/colocqq-tests.R",knit=FALSE); devtools::build_vignettes()'
```

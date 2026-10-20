# Coloc: a package for colocalisation analyses

## A brief outline of colocalisation analysis

The coloc package can be used to perform genetic colocalisation analysis
of two potentially related phenotypes, to ask whether they share common
genetic causal variant(s) in a given region. There are a few key
references which this vignette will not duplicate (see below).

In brief, two approaches can be implemented. The proportional testing
approach ^([1](#fn.1)[2](#fn.2)) has now been moved to its own package,
coloc.prop.

This package implements the more population enumeration approach.

You can read about how to prepare your data in [this
vignette](https://chr1swallace.github.io/coloc/articles/a02_data.md) or
read the vignettes listed in one of the sections below to understand how
coloc works.

### A single causal variant assumption

Claudia Giambartolomei and Vincent Plagnol proposed the enumeration
method, which makes use of Jon Wakefield’s work on determining
approximate Bayes Factors from p values ^([3](#fn.3)) to generate a
colocalisation analysis ^([4](#fn.4)), implemented in the function
[`coloc.abf()`](https://chr1swallace.github.io/coloc/reference/coloc.abf.md).
By assuming there is at most one causal variant per trait, every
possible configuration can be individually enumerated and evaluated, and
aggregating over these allows us to gauge the relative support for
models which support colocalisation to those that don’t.

You can see more about the enumeration approach on [this
blogpost](https://haldanessieve.org/2013/05/21/our-paper-bayesian-test-for-co-localisation-between-pairs-of-genetic-association-studies-using-summary-statistics/).

See vignette:
[enumeration](https://chr1swallace.github.io/coloc/articles/a03_enumeration.md)

### Sensitivity analysis

As a Bayesian method, coloc.abf() requires the user to specify prior
probabilities of SNP causality and colocalisation. Post-hoc sensitivity
analysis can be used to assess whether results are robust across a range
of plausible priors.

See vignette:
[sensitivity](https://chr1swallace.github.io/coloc/articles/a04_sensitivity.md)

### Deprecated: relaxing the single causal variant assumption through conditioning

The single variant assumption can be relaxed through conditioning. We
have implemented a conditioning step within coloc, which we hope will
increase use of conditioning, and proposed an alternative, masking.

See vignette:
[conditioning/masking](https://chr1swallace.github.io/coloc/articles/a05_conditioning.md)

### An improved approach to relaxing the single causal variant assumption: SuSiE

The sum of single effects regression method proposed by Wang et al
^([6](#fn.6)) can instead simultaneously decompose multiple signals from
the marginal summary stats, and appears to work better than
conditioning.

See vignette:
[SuSiE](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.md)

## Footnotes

^([1](#fnr.1))<https://pubmed.ncbi.nlm.nih.gov/19039033/>

^([2](#fnr.2))<https://doi.org/10.1002/gepi.21765>

^([3](#fnr.3))<https://pubmed.ncbi.nlm.nih.gov/18642345/>

^([4](#fnr.4))<https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383>

^([5](#fnr.5))<https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005272>

^([6](#fnr.6))<https://rss.onlinelibrary.wiley.com/doi/10.1111/rssb.12388>

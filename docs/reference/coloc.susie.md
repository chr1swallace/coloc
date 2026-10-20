# run coloc using susie to detect separate signals

colocalisation with multiple causal variants via SuSiE

## Usage

``` r
coloc.susie(
  dataset1,
  dataset2,
  back_calculate_lbf = FALSE,
  susie.args = list(),
  ...
)
```

## Arguments

- dataset1:

  *either* a coloc-style input dataset (see
  [check_dataset](https://chr1swallace.github.io/coloc/reference/check_dataset.md)),
  or the result of running
  [runsusie](https://chr1swallace.github.io/coloc/reference/runsusie.md)
  on such a dataset

- dataset2:

  *either* a coloc-style input dataset (see
  [check_dataset](https://chr1swallace.github.io/coloc/reference/check_dataset.md)),
  or the result of running
  [runsusie](https://chr1swallace.github.io/coloc/reference/runsusie.md)
  on such a dataset

- back_calculate_lbf:

  by default, use the log Bayes factors returned by susie_rss. It is
  also possible to back-calculate these from the posterior
  probabilities. It is not advised to set this to TRUE, the option
  exists really for testing purposes only.

- susie.args:

  a named list of additional arguments to be passed to
  [runsusie](https://chr1swallace.github.io/coloc/reference/runsusie.md)

- ...:

  other arguments passed to
  [coloc.bf_bf](https://chr1swallace.github.io/coloc/reference/coloc.bf_bf.md),
  in particular prior values for causal association with one trait (p1,
  p2) or both (p12) and and prior weights.

## Value

a list, containing elements \* summary a data.table of posterior
probabilities of each global hypothesis, one row per pairwise comparison
of signals from the two traits \* results a data.table of detailed
results giving the posterior probability for each snp to be jointly
causal for both traits *assuming H4 is true*. Please ignore this column
if the corresponding posterior support for H4 is not high. \* priors a
vector of the priors used for the analysis

## Author

Chris Wallace

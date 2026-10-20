# Query Open Targets or eQTL Catalogue for many credible sets

Runs
[`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md)
or
[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md)
for each credible set in a list. The results can then be colocalised
with the same credible sets using
[`coloc_query_many()`](https://chr1swallace.github.io/coloc/reference/coloc_query_many.md).
Keeping the two steps separate means the queries, which need network
access, can be run once (eg on an HPC login node) and saved with
saveRDS(), and the colocalisation run later or repeated without
downloading again.

## Usage

``` r
query_ot_many(cs_list, ..., max_fail = 5)

query_eqtl_cat_many(cs_list, ..., max_fail = 5)
```

## Arguments

- cs_list:

  a list of index credible sets, each as described for
  [`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md).
  Names are used to label the results; if unnamed, the credible sets are
  numbered.

- ...:

  further arguments passed to
  [`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md)
  or
  [`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md),
  eg type, quant_method, version, window

- max_fail:

  stop after this many credible sets in a row fail

## Value

a list, named as cs_list, of the query results for each credible set
(NULL for those that failed). Attribute "failed" gives the error message
for each failed credible set.

## Details

A credible set whose query fails is skipped, with its error recorded,
and the others continue. If `max_fail` credible sets in a row fail, or
all of them do, the run stops, as that suggests a systematic problem
such as no network access. Warnings are collected and reported once each
at the end, with the names of the credible sets they came from.

## See also

[`coloc_query_many()`](https://chr1swallace.github.io/coloc/reference/coloc_query_many.md)

## Author

Chris Wallace

## Examples

``` r
if (FALSE) { # \dontrun{
cs_list <- list(BACH2 = q, other = q2) # index credible sets, as for query_ot()
ot <- query_ot_many(cs_list, type = "qtl")
ecat <- query_eqtl_cat_many(cs_list, quant_method = c("ge", "leafcutter"))
res_ot <- coloc_query_many(cs_list, ot)
res_ecat <- coloc_query_many(cs_list, ecat)
} # }
```

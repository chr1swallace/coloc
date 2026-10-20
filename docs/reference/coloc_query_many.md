# Colocalise many credible sets with their query results

Runs
[`coloc_query()`](https://chr1swallace.github.io/coloc/reference/coloc_query.md)
for each credible set in cs_list against its own query results from
[`query_ot_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md)
or
[`query_eqtl_cat_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md),
matched by name, and combines the results. Credible sets whose query
failed, or found nothing, are skipped. Failures and warnings are handled
as for
[`query_ot_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md).

## Usage

``` r
coloc_query_many(cs_list, query_list, max_fail = 5)
```

## Arguments

- cs_list:

  a list of index credible sets, each as described for
  [`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md).
  Names are used to label the results; if unnamed, the credible sets are
  numbered.

- query_list:

  output of
  [`query_ot_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md)
  or
  [`query_eqtl_cat_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md)
  for the same cs_list

- max_fail:

  stop after this many credible sets in a row fail

## Value

a data.frame of coloc results, as from
[`coloc_query()`](https://chr1swallace.github.io/coloc/reference/coloc_query.md),
with an extra first column query_cs giving the name of the index
credible set each row refers to. Attribute "failed" gives the error
message for each credible set that failed, in the query or in the
colocalisation.

## See also

[`query_ot_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md),
[`coloc_query()`](https://chr1swallace.github.io/coloc/reference/coloc_query.md)

## Author

Chris Wallace

# query Open Targets credible sets

This function takes an input credible set, uses that to define a region
of interest, and searches for Open Targets credible sets that may
overlap - ie contain at least one snp within the region +/- window. The
resulting list of Open Targets credible sets may then be colocalised
with the index credible set using the coloc_query() function.

## Usage

``` r
query_ot(
  cs,
  window = 1e+05,
  type = c("all", "gwas", "qtl"),
  extra = c("traitFromSource", "pubmedId")
)
```

## Arguments

- cs:

  an index data.frame containing columns "chr" and "pos", for a single
  chromosome. chr may be numeric or character, eg 6, "6" or "X".
  Credible sets returned by query_ot() can be used directly.

- window:

  numeric value saying how far outside the index credible set we should
  search. Given that only overlapping credible sets can show
  colocalisation there is no reason to set this particularly large. The
  default, 1e+5, has been adequate in all our tests.

- type:

  all, gwas or qtl. Depending on your use call, you may wish to
  colocalise against all study types in the Open Targets dataset, just
  gwas datasets (you want to colocalise with human level traits eg to
  look for pleiotropy) or qtl datasets (you want to colocalise with
  molecular datasets eg to nominate causal genes).

- extra:

  extra columns from the Open Targets manifest

## Value

a list of credible sets with an attribute "annot" containing a
data.frame where one row per element in the list describes the
annotation of the list (most importantly, studyId).

## Details

The Open Targets credible sets are stored in 200 files, each covering a
contiguous block of the genome. By default, query_ot() downloads only
the file(s) overlapping the query region (typically one or two, ~15 MB
each) to the cache dir, and reuses them in later queries. To download
everything up front (eg before working on HPC compute nodes without
internet access), use
[`ot_fetch()`](https://chr1swallace.github.io/coloc/reference/ot_fetch.md).

Files are kept in the directory given by
`options(coloc.cache_dir = "path/to/dir")`, so they persist between R
sessions. Setting this is strongly recommended: otherwise a temporary
directory is used, and files are downloaded again in every session. See
[coloc-options](https://chr1swallace.github.io/coloc/reference/coloc-options.md).

A few Open Targets credible sets include variants on alt or patch
contigs (eg 17_KI270857v1_alt_811563_C_G). These are returned as Open
Targets provides them, with pos set to NA. Because they cannot match
variants on the primary chromosomes, their posterior probability counts
against colocalisation, so query_ot() warns when it finds them: assess
any lack of colocalisation with those credible sets carefully.

## See also

[`coloc_query()`](https://chr1swallace.github.io/coloc/reference/coloc_query.md),
[`ot_fetch()`](https://chr1swallace.github.io/coloc/reference/ot_fetch.md)

## Author

Chris Wallace, Jeffrey Pullin

## Examples

``` r
if (FALSE) { # \dontrun{
q <- data.table(chr = 6,
   pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
   snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T",
           "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
   PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))

ot <- query_ot(q)
out <- coloc_query(q, ot)

out[out$PP.H4 > 0 & out$studyType != "gwas",]
} # }
```

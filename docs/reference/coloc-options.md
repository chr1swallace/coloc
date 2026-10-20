# Package options

coloc reads the following options, set with
[`options()`](https://rdrr.io/r/base/options.html), for example in your
`.Rprofile`. They mainly control the functions that download public
credible sets
([`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md),
[`ot_fetch()`](https://chr1swallace.github.io/coloc/reference/ot_fetch.md),
[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md),
[`ensembl_to_symbol()`](https://chr1swallace.github.io/coloc/reference/ensembl_to_symbol.md)).

## Details

- `coloc.cache_dir`:

  Directory where downloaded files are kept, so they are reused between
  R sessions. If unset, files go to a temporary directory that is
  deleted when R exits, with a warning. Default: unset.

- `coloc.max_temp_download`:

  Largest download, in bytes, allowed into that temporary directory.
  Larger downloads (eg
  [`ot_fetch()`](https://chr1swallace.github.io/coloc/reference/ot_fetch.md),
  or several eQTL Catalogue quantification methods at once) stop with a
  request to set a cache directory. Default: `1e9` (1 GB).

## Examples

``` r
if (FALSE) { # \dontrun{
## in ~/.Rprofile, to keep downloads in one place
options(coloc.cache_dir = "~/data/coloc-cache")
} # }
```

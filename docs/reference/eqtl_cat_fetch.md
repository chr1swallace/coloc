# Download eQTL Catalogue credible sets

Downloads all the eQTL Catalogue credible set files for the chosen
quantification methods and release, plus the dataset metadata, to the
cache directory. query_eqtl_cat() downloads only the files it needs for
each query, so you only need this to download everything in advance, for
example on an HPC login node before running query_eqtl_cat() on compute
nodes without internet access, using the same cache directory. Files
already present are not downloaded again.

## Usage

``` r
eqtl_cat_fetch(quant_method = "ge", version = .eqtl_cat_default_version)
```

## Arguments

- quant_method:

  the eQTL Catalogue contains results for up to six different ways of
  quantifying gene expression: exon, ge, leafcutter, majiq, tx, txrev.
  See <https://www.ebi.ac.uk/eqtl/Methods/>.
  [majiq](https://www.biorxiv.org/content/10.1101/2024.07.02.601792v1)
  is a recent splicing quantification that was introduced to eQTL
  Catalogue in r8, so is not available in r7. The default is ge, but you
  can choose more than one if you wish.

- version:

  eQTL Catalogue release: "r8-beta" (default,
  <https://zenodo.org/records/23398345>) or "r7"
  (<https://zenodo.org/records/23397994>). Each release is kept in its
  own subdirectory of the cache directory.

## Value

named character vector of local file paths, named by quant_method, with
the dataset metadata file named "metadata"

## See also

[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md)

## Author

Chris Wallace

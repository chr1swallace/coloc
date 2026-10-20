# Translate Ensembl gene ids to HGNC gene symbols

Uses the HGNC complete gene set, downloaded on first use to the same
directory as the eQTL Catalogue files. Version suffixes (eg
ENSG00000160712.13) are ignored.

## Usage

``` r
ensembl_to_symbol(ids)
```

## Arguments

- ids:

  character vector of Ensembl gene ids

## Value

character vector of gene symbols, the same length as ids, with NA where
no symbol is found

# Draw extended plot of summary statistics for two coloc datasets

Draw Manhattan-style locus plots for p-values in each dataset, gene
annotations, a scatter plot of z-scores, and a table of coloc summary
statistics.

## Usage

``` r
plot_extended_datasets(
  d1,
  d2,
  x,
  first_highlight_list = NULL,
  second_highlight_list = NULL,
  first_index_snp = NULL,
  second_index_snp = NULL,
  first_trait = "Trait 1",
  second_trait = "Trait 2",
  snp_label = "snp",
  ld_label = NULL,
  show_ld = FALSE,
  locus_plot_ylim = NULL,
  ens_db = "EnsDb.Hsapiens.v86"
)
```

## Arguments

- d1:

  a coloc dataset

- d2:

  another coloc dataset

- x:

  object of class `coloc_abf` returned by
  [`coloc.abf()`](https://chr1swallace.github.io/coloc/reference/coloc.abf.md)

- first_highlight_list:

  character vector of snps to highlight in the first dataset; `'index'`
  can be passed to highlight the most significant SNP

- second_highlight_list:

  character vector of snps to highlight in the second dataset; `'index'`
  can be passed to highlight the most significant SNP

- first_index_snp:

  snp to designate as index snp in the first dataset for the purpose of
  LD visualisation

- second_index_snp:

  snp to designate as index snp in the second dataset for the purpose of
  LD visualisation

- first_trait:

  label for the first trait

- second_trait:

  label for the second trait

- snp_label:

  label for the snp column

- ld_label:

  label for the LD column

- show_ld:

  logical, whether to show LD in the locus plots

- locus_plot_ylim:

  numeric vector of length 2 specifying the y-axis limits for the locus
  plots on the -log10 scale

- ens_db:

  character string specifying Ensembl database package from which to get
  gene positions. Current (at time of writing this documentation!)
  sensible values are "EnsDb.Hsapiens.v86" for build 38 and
  "EnsDb.Hsapiens.v75" for build 37.

## Value

a gtable object

## Details

Plot a pair of coloc datasets in an extended format

This function expects that the first two elements of the coloc dataset
list d contain summary statistics. Columns labelled 'chromosome',
'position', and 'p' (p-values) are expected in each dataset. The
packages `locuszoomr` and one of `EnsDb.Hsapiens.v86` or
`EnsDb.Hsapiens.v75` are required for this function. One of the
`EnsDb.Hsapiens` libraries containing the Ensembl database specified by
ens_db must be loaded prior to use (see the example). EnsDb.Hsapiens.v86
is the default database (GRCh38/hg19); for GRCh37/hg19, use
EnsDb.Hsapiens.v75.

## Author

Tom Willis

## Examples

``` r
if (FALSE) { # \dontrun{
library(coloc)
library(EnsDb.Hsapiens.v86)
plot(plot_extended_dataset(list(first_dataset, second_dataset), coloc,
first_highlight_list = c("rs123", "rs456"), 
second_highlight_list = c("rs789", "rs1011"), 
ens_db = "EnsDb.Hsapiens.v86"))
} # }
```

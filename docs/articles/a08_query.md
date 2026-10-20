# coloc_pip and using coloc_query to query against eQTL Catalogue and Open Targets

## Using credible sets to approximate full coloc

As GWAS and QTL datasets continue to grow in number, it is becoming
harder to store complete summary statistics. Organisations like eQTL
Catalogue and Open Targets are instead moving to store/make publicly
available just the credible sets at significant loci.

A credible set is defined for coloc as a data.frame() containing at
minimum character column snp and a numeric PP column containing the
posterior probabilities under a fine mapping. Here we will generate a
couple of credible sets from our example data:

``` r

library(coloc)
```

    ## This is coloc version 7.0.2

``` r

data(coloc_test_data)
attach(coloc_test_data) ## datasets D1, D2, D3 and D4
cs1 <- finemap.abf(D1) |> credible.sets()
cs2 <- finemap.abf(D2) |> credible.sets()
cs1
```

    ##      snp        PP
    ## 105 s105 0.9558118

``` r

cs2
```

    ##      snp         PP
    ## 105 s105 0.56397680
    ## 103 s103 0.22731698
    ## 109 s109 0.11179820
    ## 116 s116 0.05601633

Now we can compare performing colocalisation using the full data using
the usual
[`coloc.abf()`](https://chr1swallace.github.io/coloc/reference/coloc.abf.md)
function or just the credible sets using the new
[`coloc_pip()`](https://chr1swallace.github.io/coloc/reference/coloc_pip.md)
function.

``` r

coloc.abf(D1,D2)
```

    ## PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf 
    ##  1.38e-18  2.94e-10  8.59e-12  8.34e-04  9.99e-01 
    ## [1] "PP abf for shared variant: 99.9%"

    ## Coloc analysis of trait 1, trait 2

    ## 
    ## SNP Priors

    ##    p1    p2   p12 
    ## 1e-04 1e-04 1e-05

    ## 
    ## Hypothesis Priors

    ##        H0   H1   H2       H3    H4
    ##  0.892505 0.05 0.05 0.002495 0.005

    ## 
    ## Posterior

    ##        nsnps           H0           H1           H2           H3           H4 
    ## 5.000000e+02 1.377000e-18 2.937336e-10 8.593226e-12 8.338917e-04 9.991661e-01

``` r

coloc_pip(cs1,cs2)
```

    ##          nsnps1          nsnps2       nsnps_int           PP.H3           PP.H4 
    ##     1.000000000     4.000000000     1.000000000     0.001399272     0.998600728 
    ## freq_sign_agree 
    ##              NA

Being able to work with credible sets opens up the chance to rapidly run
colocalisation with the large public resources like eQTL Catalogue and
Open Targets.

We first generate our own credible set. Here I’m using one created
earlier:

``` r

q <- data.frame(chr = 6,
   pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
   snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T",
           "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
   PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))
```

You will see this credible set has additional chr and pos columns. We
will need these to define the lookup region in the public data
resources. The SNP identifiers also have a particular format:
chr_pos_ref_alt. This is how snp identifiers are formatted in the public
resources, so yours will need to be too in order to match identifiers.

### Set up a cache directory first

Both query functions download files: Open Targets credible sets and a
study manifest for
[`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md),
and eQTL Catalogue files for
[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md).
coloc keeps these in a *cache directory* and reuses them, so each file
is only downloaded once.

Set a cache directory before you start with the `coloc.cache_dir`. You
can put this at the top of a script (more reproducible) or put it in
your `~/.Rprofile`, so that every coloc session uses it:

``` r

options(coloc.cache_dir = "~/data/coloc-cache")
```

If you do not, coloc falls back to a temporary directory with a warning.
This works, but the temporary directory is deleted when R exits, so
everything is downloaded again in every new session. For this reason,
downloads larger than 1 GB (such as downloading all of Open Targets with
[`ot_fetch()`](https://chr1swallace.github.io/coloc/reference/ot_fetch.md),
or all of several eQTL Catalogue quantification methods with
[`eqtl_cat_fetch()`](https://chr1swallace.github.io/coloc/reference/eqtl_cat_fetch.md))
will stop and ask you to set a cache directory.

The options that control downloading are listed in `?coloc-options`:

- `coloc.cache_dir` - where downloaded files are kept (default: unset,
  so a temporary directory is used, with a warning)
- `coloc.max_temp_download` - the largest download, in bytes, allowed
  into that temporary directory (default: `1e9`, 1 GB). Set
  `coloc.cache_dir` to download larger files.
- `coloc.ot_remote` - whether
  [`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md)
  reads Open Targets files directly from S3 instead of downloading them
  (default: `FALSE`; see below)

In this vignette we use a cache directory inside R’s temporary
directory, so that building the vignette leaves no files behind. For
your own work, use a directory that persists, as above.

``` r

options(coloc.cache_dir = file.path(tempdir(), "coloc-cache"))
```

## Query public resources

To query a public resource, you need two steps:

1.  `x <- query_*(cs)` will download the data overlapping your credible
    set

2.  `coloc_query(cs, x)` will colocalise your credible set against the
    downloaded data.

Here are the steps in detail for Open Targets and eQTL Catalogue, in
turn.

### Query Open Targets

Open Targets stores its fine-mapped credible sets as 200 files in a
public S3 bucket, about 3.2 GB in total for release 26.06. Each file
covers a contiguous block of the genome.
[`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md)
works out which file or files overlap your region (usually just one,
around 15 MB), downloads only those to the cache directory, and reads
them locally. Later queries in the same region reuse the downloaded
files, so they are fast and need no network access. It also downloads a
study manifest (about 100 MB, once) so it can describe what each Open
Targets study relates to.

``` r

ot <- query_ot(q)
```

    ## Downloading 1 Open Targets credible set file(s) (24 MB) to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/credible_set_26.06

    ## Downloading part-00166-8ee8599e-d3fd-478e-bb9b-6a3c3243641d-c000.snappy.parquet to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/credible_set_26.06

    ## Downloading part-00000-1c0a4198-4d65-425e-82a1-d9d8b54d8549-c000.snappy.parquet to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache

``` r

out <- coloc_query(q, ot)
head(out) 
```

    ##                                     leadVariant
    ## 6743c2da9f466cf07a57e5c4c37de167 6_90177105_C_T
    ## bfc7ad872426a5b36b0c8e3b6a4c4450 6_90177105_C_T
    ## b0426b8b92d8e1459f9decdfa27dec2b 6_90177105_C_T
    ## fc3d89e923b8804799a4353dac6d5065 6_90177105_C_T
    ## 8c45c9fef65991d8ad549676dc385768 6_90178232_C_T
    ## 5c515b6d48a02f4a435592567ff98c34 6_90182225_T_C
    ##                                                      studyLocusId
    ## 6743c2da9f466cf07a57e5c4c37de167 6743c2da9f466cf07a57e5c4c37de167
    ## bfc7ad872426a5b36b0c8e3b6a4c4450 bfc7ad872426a5b36b0c8e3b6a4c4450
    ## b0426b8b92d8e1459f9decdfa27dec2b b0426b8b92d8e1459f9decdfa27dec2b
    ## fc3d89e923b8804799a4353dac6d5065 fc3d89e923b8804799a4353dac6d5065
    ## 8c45c9fef65991d8ad549676dc385768 8c45c9fef65991d8ad549676dc385768
    ## 5c515b6d48a02f4a435592567ff98c34 5c515b6d48a02f4a435592567ff98c34
    ##                                                                                   studyId
    ## 6743c2da9f466cf07a57e5c4c37de167                                             GCST90479706
    ## bfc7ad872426a5b36b0c8e3b6a4c4450                                             GCST90479705
    ## b0426b8b92d8e1459f9decdfa27dec2b                                             GCST90475514
    ## fc3d89e923b8804799a4353dac6d5065                                             GCST90475511
    ## 8c45c9fef65991d8ad549676dc385768 ipscore_exon_ipsc_ensg00000112159_12_6_89684876_89684985
    ## 5c515b6d48a02f4a435592567ff98c34                                             GCST90479855
    ##                                  chromosome position studyType
    ## 6743c2da9f466cf07a57e5c4c37de167          6 90177105      gwas
    ## bfc7ad872426a5b36b0c8e3b6a4c4450          6 90177105      gwas
    ## b0426b8b92d8e1459f9decdfa27dec2b          6 90177105      gwas
    ## fc3d89e923b8804799a4353dac6d5065          6 90177105      gwas
    ## 8c45c9fef65991d8ad549676dc385768          6 90178232      eqtl
    ## 5c515b6d48a02f4a435592567ff98c34          6 90182225      gwas
    ##                                                                                                       traitFromSource
    ## 6743c2da9f466cf07a57e5c4c37de167                                   monocyte (fraction, minimum, inv-norm transformed)
    ## bfc7ad872426a5b36b0c8e3b6a4c4450                                      monocyte (fraction, mean, inv-norm transformed)
    ## b0426b8b92d8e1459f9decdfa27dec2b                                   monocyte (fraction, minimum, inv-norm transformed)
    ## fc3d89e923b8804799a4353dac6d5065                                      monocyte (fraction, mean, inv-norm transformed)
    ## 8c45c9fef65991d8ad549676dc385768                                               ENSG00000112159.12_6_89684876_89684985
    ## 5c515b6d48a02f4a435592567ff98c34 Benign neoplasm of pituitary gland and craniopharyngeal duct (pouch) (PheCode 227.3)
    ##                                  pubmedId nsnps1 nsnps2 nsnps_int PP.H3 PP.H4
    ## 6743c2da9f466cf07a57e5c4c37de167 39024449      6     43         0     1     0
    ## bfc7ad872426a5b36b0c8e3b6a4c4450 39024449      6     20         0     1     0
    ## b0426b8b92d8e1459f9decdfa27dec2b 39024449      6     12         0     1     0
    ## fc3d89e923b8804799a4353dac6d5065 39024449      6      6         0     1     0
    ## 8c45c9fef65991d8ad549676dc385768 28410642      6     32         0     1     0
    ## 5c515b6d48a02f4a435592567ff98c34 39024449      6      1         0     1     0
    ##                                  freq_sign_agree
    ## 6743c2da9f466cf07a57e5c4c37de167              NA
    ## bfc7ad872426a5b36b0c8e3b6a4c4450              NA
    ## b0426b8b92d8e1459f9decdfa27dec2b              NA
    ## fc3d89e923b8804799a4353dac6d5065              NA
    ## 8c45c9fef65991d8ad549676dc385768              NA
    ## 5c515b6d48a02f4a435592567ff98c34              NA

``` r

out[out$PP.H4 > 0 & out$studyType != "gwas",]
```

    ##                                     leadVariant
    ## aac59cadd0f59b18622471f63455c113 6_90227175_C_A
    ## d99924866bffc86e4d33b655e076c554 6_90241017_T_C
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 6_90279406_G_C
    ## 056181f6aca3d589248bcc5bf3a71ef8 6_90279406_G_C
    ## 906056ee179f3e91bf0e27e6b8765629 6_90281769_T_C
    ##                                                      studyLocusId
    ## aac59cadd0f59b18622471f63455c113 aac59cadd0f59b18622471f63455c113
    ## d99924866bffc86e4d33b655e076c554 d99924866bffc86e4d33b655e076c554
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 3402ca2aaa4e1c7cd61fd49e0ecb7e48
    ## 056181f6aca3d589248bcc5bf3a71ef8 056181f6aca3d589248bcc5bf3a71ef8
    ## 906056ee179f3e91bf0e27e6b8765629 906056ee179f3e91bf0e27e6b8765629
    ##                                                                                                       studyId
    ## aac59cadd0f59b18622471f63455c113                                            onek1k_ge_cd4_tcm_ensg00000112182
    ## d99924866bffc86e4d33b655e076c554                  blueprint_leafcutter_t-cell_6_90206681_90296480_clu_39598__
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 bossini-castillo_2019_exon_treg_naive_ensg00000112182_15_6_89926529_89931484
    ## 056181f6aca3d589248bcc5bf3a71ef8                          bossini-castillo_2019_ge_treg_naive_ensg00000112182
    ## 906056ee179f3e91bf0e27e6b8765629                                          onek1k_ge_cd4_naive_ensg00000112182
    ##                                  chromosome position studyType
    ## aac59cadd0f59b18622471f63455c113          6 90227175    sceqtl
    ## d99924866bffc86e4d33b655e076c554          6 90241017      sqtl
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48          6 90279406      eqtl
    ## 056181f6aca3d589248bcc5bf3a71ef8          6 90279406      eqtl
    ## 906056ee179f3e91bf0e27e6b8765629          6 90281769    sceqtl
    ##                                                         traitFromSource
    ## aac59cadd0f59b18622471f63455c113                        ENSG00000112182
    ## d99924866bffc86e4d33b655e076c554        6:90206681:90296480:clu_39598_+
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 ENSG00000112182.15_6_89926529_89931484
    ## 056181f6aca3d589248bcc5bf3a71ef8                        ENSG00000112182
    ## 906056ee179f3e91bf0e27e6b8765629                        ENSG00000112182
    ##                                  pubmedId nsnps1 nsnps2 nsnps_int      PP.H3
    ## aac59cadd0f59b18622471f63455c113 35389779      6     12         5 0.34481976
    ## d99924866bffc86e4d33b655e076c554 27863251      6     14         5 0.12097880
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 35591976      6      4         1 0.00104726
    ## 056181f6aca3d589248bcc5bf3a71ef8 35591976      6     53         1 0.05568514
    ## 906056ee179f3e91bf0e27e6b8765629 35389779      6     14         2 0.70183927
    ##                                      PP.H4 freq_sign_agree
    ## aac59cadd0f59b18622471f63455c113 0.6551802              NA
    ## d99924866bffc86e4d33b655e076c554 0.8790212              NA
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 0.9989527              NA
    ## 056181f6aca3d589248bcc5bf3a71ef8 0.9443149              NA
    ## 906056ee179f3e91bf0e27e6b8765629 0.2981607              NA

In this case, we are only interested in seeing what evidence there is
for a causal gene, so we exclude studyType=“gwas”. We could have done
this directly using the type argument:

``` r

ot <- query_ot(q, type="qtl")
out <- coloc_query(q, ot)
out[out$PP.H4 > 0,]
```

    ##                                     leadVariant
    ## aac59cadd0f59b18622471f63455c113 6_90227175_C_A
    ## d99924866bffc86e4d33b655e076c554 6_90241017_T_C
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 6_90279406_G_C
    ## 056181f6aca3d589248bcc5bf3a71ef8 6_90279406_G_C
    ## 906056ee179f3e91bf0e27e6b8765629 6_90281769_T_C
    ##                                                      studyLocusId
    ## aac59cadd0f59b18622471f63455c113 aac59cadd0f59b18622471f63455c113
    ## d99924866bffc86e4d33b655e076c554 d99924866bffc86e4d33b655e076c554
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 3402ca2aaa4e1c7cd61fd49e0ecb7e48
    ## 056181f6aca3d589248bcc5bf3a71ef8 056181f6aca3d589248bcc5bf3a71ef8
    ## 906056ee179f3e91bf0e27e6b8765629 906056ee179f3e91bf0e27e6b8765629
    ##                                                                                                       studyId
    ## aac59cadd0f59b18622471f63455c113                                            onek1k_ge_cd4_tcm_ensg00000112182
    ## d99924866bffc86e4d33b655e076c554                  blueprint_leafcutter_t-cell_6_90206681_90296480_clu_39598__
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 bossini-castillo_2019_exon_treg_naive_ensg00000112182_15_6_89926529_89931484
    ## 056181f6aca3d589248bcc5bf3a71ef8                          bossini-castillo_2019_ge_treg_naive_ensg00000112182
    ## 906056ee179f3e91bf0e27e6b8765629                                          onek1k_ge_cd4_naive_ensg00000112182
    ##                                  chromosome position studyType
    ## aac59cadd0f59b18622471f63455c113          6 90227175    sceqtl
    ## d99924866bffc86e4d33b655e076c554          6 90241017      sqtl
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48          6 90279406      eqtl
    ## 056181f6aca3d589248bcc5bf3a71ef8          6 90279406      eqtl
    ## 906056ee179f3e91bf0e27e6b8765629          6 90281769    sceqtl
    ##                                                         traitFromSource
    ## aac59cadd0f59b18622471f63455c113                        ENSG00000112182
    ## d99924866bffc86e4d33b655e076c554        6:90206681:90296480:clu_39598_+
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 ENSG00000112182.15_6_89926529_89931484
    ## 056181f6aca3d589248bcc5bf3a71ef8                        ENSG00000112182
    ## 906056ee179f3e91bf0e27e6b8765629                        ENSG00000112182
    ##                                  pubmedId nsnps1 nsnps2 nsnps_int      PP.H3
    ## aac59cadd0f59b18622471f63455c113 35389779      6     12         5 0.34481976
    ## d99924866bffc86e4d33b655e076c554 27863251      6     14         5 0.12097880
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 35591976      6      4         1 0.00104726
    ## 056181f6aca3d589248bcc5bf3a71ef8 35591976      6     53         1 0.05568514
    ## 906056ee179f3e91bf0e27e6b8765629 35389779      6     14         2 0.70183927
    ##                                      PP.H4 freq_sign_agree
    ## aac59cadd0f59b18622471f63455c113 0.6551802              NA
    ## d99924866bffc86e4d33b655e076c554 0.8790212              NA
    ## 3402ca2aaa4e1c7cd61fd49e0ecb7e48 0.9989527              NA
    ## 056181f6aca3d589248bcc5bf3a71ef8 0.9443149              NA
    ## 906056ee179f3e91bf0e27e6b8765629 0.2981607              NA

### Query eQTL Catalogue

Querying the eQTL Catalogue works the same way as Open Targets. The eQTL
Catalogue publishes its fine-mapped credible sets as one large file per
quantification method, each covering the whole genome. For coloc, these
have been split into smaller files (about 25 MB each) sorted by genomic
position, with an index recording which part of the genome each file
covers. No credible set is split across files, and the data are
otherwise unchanged.
[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md)
uses the index to download only the file(s) overlapping your region -
usually one per quantification method - to the cache directory, and
reuses them in later queries.

#### Releases and quantification methods

Two eQTL Catalogue releases are available, chosen with the `version`
argument: “r8-beta” and “r7” (the default). The split files for each are
on Zenodo, for [r7](https://zenodo.org/records/23397994) and
[r8-beta](https://zenodo.org/records/23398345). Each release is kept in
its own subdirectory of the cache directory, so they never overwrite
each other and you can use both side by side.

Within a release, the data are divided by the method used to quantify
gene expression, chosen with `quant_method` (default “ge”). You can ask
for more than one at once, eg `quant_method = c("ge", "leafcutter")`.
See <https://www.ebi.ac.uk/eqtl/Methods/> for details of each method.

| quant_method | what it measures | r7 credible sets | r8-beta credible sets |
|----|----|----|----|
| ge | total gene expression | 503,000 | 721,000 |
| exon | exon expression | 626,000 | 416,000 |
| tx | transcript usage | 368,000 | 353,000 |
| txrev | transcriptional event usage | 411,000 | 403,000 |
| leafcutter | splicing (LeafCutter) | 234,000 | 251,000 |
| majiq | splicing (MAJIQ), new in r8 | \- | 713,000 |

The two releases differ slightly in their columns, and this shows in the
annotation of the results: r8-beta credible sets carry a
`molecular_trait_object_id` column, while r7 has a `gene_id` column
instead.

#### Example

We will use “exon” for this vignette. In this case, there is only one
exon credible set in the region but it doesn’t overlap our credible set
and so doesn’t colocalise.

``` r

ecat <- query_eqtl_cat(q, quant_method="exon")
```

    ## Downloading 1 eQTL Catalogue r8-beta file(s) (23 MB) to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/eQTL_Catalogue_r8-beta

    ## Downloading eQTL_Catalogue_r8-beta_exon_004.parquet to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/eQTL_Catalogue_r8-beta

    ## Downloading dataset_metadata.tsv to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/eQTL_Catalogue_r8-beta

``` r

out <- coloc_query(q, ecat)
out
```

    ##                                                     dataset_id
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2  QTD000357
    ##                                                                         molecular_trait_id
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 ENSG00000083099.11_6_89636857_89636912
    ##                                                     molecular_trait_object_id
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2           ENSG00000083099
    ##                                                             gene_id
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 ENSG00000083099
    ##                                                                                                   cs_id
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 QTD000357_ENSG00000083099.11_6_89636857_89636912_L2
    ##                                                                     region
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 chr6:88636884-90636884
    ##                                                      study_id study_label
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 QTS000015    GTEx_v10
    ##                                                     sample_group      tissue_id
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2        blood UBERON_0000178
    ##                                                     tissue_label
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2        blood
    ##                                                     condition_label sample_size
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2           naive         853
    ##                                                     quant_method     pmid
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2         exon 32913098
    ##                                                     study_type nsnps1 nsnps2
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2       bulk      6      1
    ##                                                     nsnps_int PP.H3 PP.H4
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2         0     1     0
    ##                                                     freq_sign_agree
    ## QTD000357_ENSG00000083099.11_6_89636857_89636912_L2              NA

The same query against the earlier release:

``` r

ecat7 <- query_eqtl_cat(q, quant_method="exon", version = "r7")
```

    ## Downloading 1 eQTL Catalogue r7 file(s) (21 MB) to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/eQTL_Catalogue_r7

    ## Downloading eQTL_Catalogue_r7_exon_009.parquet to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/eQTL_Catalogue_r7

    ## Downloading dataset_metadata_r7.tsv to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache/eQTL_Catalogue_r7

    ## Warning in query_eqtl_cat(q, quant_method = "exon", version = "r7"): 1 / 10
    ## credsets have sum(PP) not in [0,1]. Removing:
    ## QTD000032_ENSG00000112182.15_6_90089066_90089109_L1

``` r

coloc_query(q, ecat7)
```

    ##                                                     dataset_id
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1  QTD000037
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1  QTD000057
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1  QTD000122
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1  QTD000217
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1  QTD000282
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1  QTD000302
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1  QTD000367
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1  QTD000445
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1  QTD000495
    ##                                                                         molecular_trait_id
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 ENSG00000112182.15_6_89926529_89931484
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 ENSG00000112159.12_6_89745492_89745626
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 ENSG00000112159.12_6_89751431_89751582
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 ENSG00000112159.12_6_89693149_89693153
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 ENSG00000112159.12_6_89686902_89687023
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 ENSG00000112159.12_6_89693149_89693153
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 ENSG00000112159.12_6_89684876_89684985
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 ENSG00000112159.12_6_89658610_89658917
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 ENSG00000112159.12_6_89716650_89716809
    ##                                                             gene_id
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 ENSG00000112182
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 ENSG00000112159
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 ENSG00000112159
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 ENSG00000112159
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 ENSG00000112159
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 ENSG00000112159
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 ENSG00000112159
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 ENSG00000112159
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 ENSG00000112159
    ##                                                                                                   cs_id
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 QTD000037_ENSG00000112182.15_6_89926529_89931484_L1
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 QTD000057_ENSG00000112159.12_6_89745492_89745626_L1
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 QTD000122_ENSG00000112159.12_6_89751431_89751582_L1
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 QTD000217_ENSG00000112159.12_6_89693149_89693153_L1
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 QTD000282_ENSG00000112159.12_6_89686902_89687023_L1
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 QTD000302_ENSG00000112159.12_6_89693149_89693153_L1
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 QTD000367_ENSG00000112159.12_6_89684876_89684985_L1
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 QTD000445_ENSG00000112159.12_6_89658610_89658917_L1
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 QTD000495_ENSG00000112159.12_6_89716650_89716809_L1
    ##                                                                     region
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 chr6:88929006-90929006
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 chr6:88745559-90745559
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 chr6:88751506-90751506
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 chr6:88693151-90693151
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 chr6:88686962-90686962
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 chr6:88693151-90693151
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 chr6:88684930-90684930
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 chr6:88658763-90658763
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 chr6:88716729-90716729
    ##                                                      study_id
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 QTS000003
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 QTS000006
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 QTS000015
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 QTS000015
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 QTS000015
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 QTS000015
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 QTS000017
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 QTS000026
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 QTS000026
    ##                                                               study_label
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 Bossini-Castillo_2019
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1                   CAP
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1                  GTEx
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1                  GTEx
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1                  GTEx
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1                  GTEx
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1               iPSCORE
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1        Schmiedel_2018
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1        Schmiedel_2018
    ##                                                                 sample_group
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1               Treg_naive
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1                LCL_naive
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1         adipose_visceral
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1               fibroblast
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1                   muscle
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1                pituitary
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1                     iPSC
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1              Th17_memory
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 CD8_T-cell_anti-CD3-CD28
    ##                                                          tissue_id
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1     CL_0002677
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1    EFO_0005292
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 UBERON_0010414
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1     CL_0000057
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 UBERON_0001134
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 UBERON_0000007
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1    EFO_0004905
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1     CL_0000899
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1     CL_0000625
    ##                                                           tissue_label
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1        Treg memory
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1                LCL
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 adipose (visceral)
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1         fibroblast
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1             muscle
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1          pituitary
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1               iPSC
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1          Th17 cell
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1        CD8+ T cell
    ##                                                      condition_label
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1            naive
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1            naive
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1            naive
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1            naive
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1            naive
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1            naive
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1            naive
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1            naive
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 anti-CD3-CD28_4h
    ##                                                     sample_size quant_method
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1         119         exon
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1         148         exon
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1         469         exon
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1         483         exon
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1         702         exon
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1         237         exon
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1         101         exon
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1          89         exon
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1          88         exon
    ##                                                         pmid study_type nsnps1
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 35591976       bulk      6
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 32787775       bulk      6
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 32913098       bulk      6
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 32913098       bulk      6
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 32913098       bulk      6
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 32913098       bulk      6
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 28410642       bulk      6
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 30449622       bulk      6
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 30449622       bulk      6
    ##                                                     nsnps2 nsnps_int      PP.H3
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1      4         1 0.00104726
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1      3         0 1.00000000
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1     23         0 1.00000000
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1      4         0 1.00000000
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1     13         0 1.00000000
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1      5         0 1.00000000
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1     28         0 1.00000000
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1     35         0 1.00000000
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1      7         0 1.00000000
    ##                                                         PP.H4 freq_sign_agree
    ## QTD000037_ENSG00000112182.15_6_89926529_89931484_L1 0.9989527              NA
    ## QTD000057_ENSG00000112159.12_6_89745492_89745626_L1 0.0000000              NA
    ## QTD000122_ENSG00000112159.12_6_89751431_89751582_L1 0.0000000              NA
    ## QTD000217_ENSG00000112159.12_6_89693149_89693153_L1 0.0000000              NA
    ## QTD000282_ENSG00000112159.12_6_89686902_89687023_L1 0.0000000              NA
    ## QTD000302_ENSG00000112159.12_6_89693149_89693153_L1 0.0000000              NA
    ## QTD000367_ENSG00000112159.12_6_89684876_89684985_L1 0.0000000              NA
    ## QTD000445_ENSG00000112159.12_6_89658610_89658917_L1 0.0000000              NA
    ## QTD000495_ENSG00000112159.12_6_89716650_89716809_L1 0.0000000              NA

The credible sets returned contain those of their variants that lie
within the query window (the region spanned by your credible set, plus
`window` either side, 100 kb by default). Because only overlapping
credible sets can colocalise, this rarely matters, but you can increase
`window` if you want more of each credible set.

Because ecat is a small object, we can examine its structure. It is a
list of credible sets, with an attribute “annot”. annot is a data.frame
containing one row per credible set giving information about the
analysis it refers to.

``` r

ecat
```

    ## $QTD000357_ENSG00000083099.11_6_89636857_89636912_L2
    ##              snp        PP      beta
    ## 1 6_90383144_A_T 0.5050769 -0.714343
    ## 
    ## attr(,"annot")
    ##   dataset_id                     molecular_trait_id molecular_trait_object_id
    ## 1  QTD000357 ENSG00000083099.11_6_89636857_89636912           ENSG00000083099
    ##           gene_id                                               cs_id
    ## 1 ENSG00000083099 QTD000357_ENSG00000083099.11_6_89636857_89636912_L2
    ##                   region  study_id study_label sample_group      tissue_id
    ## 1 chr6:88636884-90636884 QTS000015    GTEx_v10        blood UBERON_0000178
    ##   tissue_label condition_label sample_size quant_method     pmid study_type
    ## 1        blood           naive         853         exon 32913098       bulk

#### Gene symbols

The eQTL Catalogue identifies genes by Ensembl gene id.
[`ensembl_to_symbol()`](https://chr1swallace.github.io/coloc/reference/ensembl_to_symbol.md)
translates these to HGNC gene symbols, using the HGNC gene list
(downloaded once to the cache directory, about 15 MB):

``` r

annot <- attr(ecat, "annot")
gene_col <- intersect(c("gene_id", "molecular_trait_object_id"), names(annot))[1] # depends on the release
annot$symbol <- ensembl_to_symbol(annot[[gene_col]])
```

    ## Downloading hgnc_complete_set.txt to /private/var/folders/p7/d3q7qs6d7fl0dpb73g4rr7jh0000gn/T/Rtmp12HWGR/coloc-cache

``` r

annot[, c("cs_id", gene_col, "symbol")]
```

    ##                                                 cs_id         gene_id symbol
    ## 1 QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 ENSG00000083099  LYRM2

Ids without an approved HGNC symbol (eg some lncRNAs) are returned as
`NA`.

### Downloading data in advance

Normally you don’t need to think about downloading: each query downloads
the files it needs to the cache directory, and later queries reuse them.
Two functions download data in advance and populate the cache directory:

- [`ot_fetch()`](https://chr1swallace.github.io/coloc/reference/ot_fetch.md)
  downloads all the Open Targets credible sets (200 files, about 3.2 GB
  for release 26.06).
- [`eqtl_cat_fetch()`](https://chr1swallace.github.io/coloc/reference/eqtl_cat_fetch.md)
  downloads all the eQTL Catalogue files for the quantification methods
  and release you choose, plus the dataset metadata. A whole release is
  about 2.2 GB; “ge” alone is about 0.7 GB.

This is useful when:

- **working on an HPC whose compute nodes have no internet access.** Run
  the fetch functions once on a login node (or any node with internet
  access), then run your queries on compute nodes using the same cache
  directory. The queries find the files already there and never touch
  the network.
- **running many queries**, so that every region’s files are downloaded
  in one go, rather than piecemeal during the run. See [Many credible
  sets](#many-credible-sets) below.

``` r

options(coloc.cache_dir = "/shared/project/coloc-cache")   # the same directory you will query from
ot_fetch()                                                # all Open Targets credible sets
eqtl_cat_fetch("ge")                                      # eQTL Catalogue r8-beta, total gene expression
eqtl_cat_fetch(c("ge", "leafcutter"), version = "r7")     # or several methods, from another release
```

Because these downloads are large, they need a persistent cache
directory: set `coloc.cache_dir`.

### Multiple credible sets

Often you will have several index credible sets to test.
[`query_ot_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md),
[`query_eqtl_cat_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md)
and
[`coloc_query_many()`](https://chr1swallace.github.io/coloc/reference/coloc_query_many.md)
do the same as
[`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md),
[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md)
and
[`coloc_query()`](https://chr1swallace.github.io/coloc/reference/coloc_query.md),
for a named list of credible sets. They work in two steps: first query,
then colocalise.

``` r

cs_list <- list(locus1 = q, locus2 = q)   # add more credible sets to the list, each named
ot_list <- query_ot_many(cs_list, type = "qtl")
ecat_list <- query_eqtl_cat_many(cs_list, quant_method = "exon")
res <- coloc_query_many(cs_list, ot_list)
head(res)
```

    ##   query_cs    leadVariant                     studyLocusId
    ## 1   locus1 6_90178232_C_T 8c45c9fef65991d8ad549676dc385768
    ## 2   locus1 6_90189436_G_A 6834396b03b49531eb226a3f68360557
    ## 3   locus1 6_90216893_C_A 416f17664700bc789e3317cbc8d760d3
    ## 4   locus1 6_90222139_C_T 46507521195da882087038dd500d4a7f
    ## 5   locus1 6_90227175_C_A aac59cadd0f59b18622471f63455c113
    ## 6   locus1 6_90231570_C_T 5631830e7927e93c292884b5ac17a8ee
    ##                                                      studyId chromosome
    ## 1   ipscore_exon_ipsc_ensg00000112159_12_6_89684876_89684985          6
    ## 2 gtex_exon_pituitary_ensg00000112159_12_6_89693149_89693153          6
    ## 3                   UKB_PPP_EUR_TNFRSF13B_O14836_OID20702_v1          6
    ## 4                       UKB_PPP_EUR_CD207_Q9UJ71_OID21298_v1          6
    ## 5                          onek1k_ge_cd4_tcm_ensg00000112182          6
    ## 6                  nathan_2022_ge_cd4__cd27__ensg00000112182          6
    ##   position studyType                        traitFromSource pubmedId nsnps1
    ## 1 90178232      eqtl ENSG00000112159.12_6_89684876_89684985 28410642      6
    ## 2 90189436      eqtl ENSG00000112159.12_6_89693149_89693153 32913098      6
    ## 3 90216893      pqtl           TNFRSF13B:O14836:OID20702:v1     <NA>      6
    ## 4 90222139      pqtl               CD207:Q9UJ71:OID21298:v1     <NA>      6
    ## 5 90227175    sceqtl                        ENSG00000112182 35389779      6
    ## 6 90231570    sceqtl                        ENSG00000112182 35545678      6
    ##   nsnps2 nsnps_int     PP.H3     PP.H4 freq_sign_agree
    ## 1     32         0 1.0000000 0.0000000              NA
    ## 2      7         0 1.0000000 0.0000000              NA
    ## 3     29         0 1.0000000 0.0000000              NA
    ## 4     22         0 1.0000000 0.0000000              NA
    ## 5     12         5 0.3448198 0.6551802              NA
    ## 6     50         0 1.0000000 0.0000000              NA

``` r

res_ecat <- coloc_query_many(cs_list, ecat_list)
res_ecat
```

    ##   query_cs dataset_id                     molecular_trait_id
    ## 1   locus1  QTD000357 ENSG00000083099.11_6_89636857_89636912
    ## 2   locus2  QTD000357 ENSG00000083099.11_6_89636857_89636912
    ##   molecular_trait_object_id         gene_id
    ## 1           ENSG00000083099 ENSG00000083099
    ## 2           ENSG00000083099 ENSG00000083099
    ##                                                 cs_id                 region
    ## 1 QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 chr6:88636884-90636884
    ## 2 QTD000357_ENSG00000083099.11_6_89636857_89636912_L2 chr6:88636884-90636884
    ##    study_id study_label sample_group      tissue_id tissue_label
    ## 1 QTS000015    GTEx_v10        blood UBERON_0000178        blood
    ## 2 QTS000015    GTEx_v10        blood UBERON_0000178        blood
    ##   condition_label sample_size quant_method     pmid study_type nsnps1 nsnps2
    ## 1           naive         853         exon 32913098       bulk      6      1
    ## 2           naive         853         exon 32913098       bulk      6      1
    ##   nsnps_int PP.H3 PP.H4 freq_sign_agree
    ## 1         0     1     0              NA
    ## 2         0     1     0              NA

[`query_ot_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md)
and
[`query_eqtl_cat_many()`](https://chr1swallace.github.io/coloc/reference/query_ot_many.md)
take the same arguments as
[`query_ot()`](https://chr1swallace.github.io/coloc/reference/query_ot.md)
and
[`query_eqtl_cat()`](https://chr1swallace.github.io/coloc/reference/query_eqtl_cat.md),
and return a list of query results, one per credible set, named as
`cs_list`.
[`coloc_query_many()`](https://chr1swallace.github.io/coloc/reference/coloc_query_many.md)
matches each credible set to its own query results by name, and returns
one table of results, with a `query_cs` column first saying which of
your credible sets each row refers to.

Keeping the two steps separate means the slow step, querying (which may
need to download), can be done once and saved, and the fast step,
colocalisation, run or repeated later without the network. For example,
on an HPC:

``` r

## on a node with internet access
ot_list <- query_ot_many(cs_list)
saveRDS(ot_list, "ot_list.rds")

## later, anywhere
res <- coloc_query_many(cs_list, readRDS("ot_list.rds"))
```

With many credible sets, a few may fail - for example a credible set
with variants in a format the query can’t use. Rather than stop, the
`_many` functions skip a credible set that fails, record why in the
“failed” attribute of their result, and carry on:

``` r

attr(ot_list, "failed")   # error message for each credible set that failed
attr(res, "failed")       # failures in the query or the colocalisation
```

If several credible sets in a row fail (5 by default, set by
`max_fail`), or all of them do, the run stops, as that suggests a
problem affecting everything, such as no network access. Warnings - for
example about a temporary cache directory, or Open Targets variants on
alt/patch contigs - are collected and reported once each at the end,
with the names of the credible sets they came from, rather than repeated
for every credible set.

The `_many` functions run one credible set after another. To speed up a
large run, download the data first with [the fetch
functions](#downloading-data-in-advance), so the queries read local
files. If you do run queries in parallel yourself, for example with
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
download first in the same way, so that workers don’t each download the
same files at once. Forked workers (`mclapply()`) inherit your options,
but socket-based workers (eg
[`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html),
or `future`’s multisession plan) start fresh R sessions without them, so
set the cache directory in each worker too, eg
`parallel::clusterEvalQ(cl, options(coloc.cache_dir = "/shared/project/coloc-cache"))`.

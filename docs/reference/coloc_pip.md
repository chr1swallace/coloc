# An approximation to coloc using only information from fine-mapped credible sets.

A credible set should be passed in the form of a list, data.frame, or
equivalent tabular object containing at least two equal-length vectors.
The object should contain at a minimim

## Usage

``` r
coloc_pip(cs1, cs2, p1 = 1e-04, p2 = 1e-04, p12 = 5e-06)
```

## Arguments

- cs1:

  credible set 1

- cs2:

  credible set 2

- p1:

  prior probability a SNP is associated with trait 1, default 1e-4

- p2:

  prior probability a SNP is associated with trait 2, default 1e-4

- p12:

  prior probability a SNP is associated with both traits, default 1e-5

## Value

a named vector containing the number of SNPs found in each credible set
(nsnps1, nsnps2), the number in both credible sets (nsnps_int) and the
posterior probabilities of H3 and H4 (PP.H3 and PP.H4)

## Details

- snp: a character vector of snp identifiers

- PP: a numeric vector of posterior probabilities

Each cs may additionally contain a column

- beta: the estimated per all snp effect on the trait

If they both do, then the returned object will also include
"freq_sign_agree", which is the proportion of snps in both credible sets
for which the difrection of effects (sign of beta) agrees. This is a
useful quantity to be able to say "higher gene expression associates
with higher (lower) disease risk" when freq_sign_agree=1 (0). If 0 \<
freq_sign_agree \< 1 then there is some disagreement between SNPs and
users are invited to investigate potential causes of this further.

## Author

Jeffrey Pullin, Chris Wallace

## Examples

``` r
cs1 <- data.frame(snp=c("s1","s2"), PP=c(0.05,0.9))
cs2 <- data.frame(snp=c("s1","s2","s3"), PP=c(0.1,0.8,0.05))
coloc_pip(cs1,cs2)
#>          nsnps1          nsnps2       nsnps_int           PP.H3           PP.H4 
#>    2.0000000000    3.0000000000    2.0000000000    0.0004894155    0.9995105845 
#> freq_sign_agree 
#>              NA 
```

<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [If I understand correctly, coloc.abf() can be run with correlated variants, that is, no prerequisite for taking through LD pruning/clumping is required. Am I correct in my understanding ?](#if-i-understand-correctly-colocabf-can-be-run-with-correlated-variants-that-is-no-prerequisite-for-taking-through-ld-pruningclumping-is-required-am-i-correct-in-my-understanding-)
- [Assume I identify a sentinel variant for a block of genome, can I do a comparison with just one variant using coloc.abf()?](#assume-i-identify-a-sentinel-variant-for-a-block-of-genome-can-i-do-a-comparison-with-just-one-variant-using-colocabf)
- [Can the process of identifying colocalized variants be carried out genome wide?](#can-the-process-of-identifying-colocalized-variants-be-carried-out-genome-wide)
- [How is coloc.abf accounting for the correlated variants?](#how-is-colocabf-accounting-for-the-correlated-variants)
- [How to define priors, is it dependent on sample size or any other parameters?](#how-to-define-priors-is-it-dependent-on-sample-size-or-any-other-parameters)
- [What does high PP4 mean?](#what-does-high-pp4-mean)

<!-- markdown-toc end -->

# If I understand correctly, coloc.abf() can be run with correlated variants, that is, no prerequisite for taking through LD pruning/clumping is required. Am I correct in my understanding ?

Yes, coloc.abf() and coloc.signals() assume they are given a dense map of all SNPs in a region that could be causal.   Do not prune and clump.

# Assume I identify a sentinel variant for a block of genome, can I do a comparison with just one variant using coloc.abf()?

No, coloc.abf() and coloc.signals() assume they are given a dense map of all SNPs in a region that could be causal. This means you need to give all SNPs in a region. You can imagine they ask whether the patterns "match" across this region of SNPs, and a single variant does not represent a pattern. 

# Can the process of identifying colocalized variants be carried out genome wide?

Coloc is designed to address whether two traits share causal variant(s) in a genomic region.  It leaves the definition of "region" up to the user.  You need to break the genome into smaller regions, within which it is reasonable to assume there is at most one (coloc.abf) or a small number (coloc.signals) of causal variants per trait.  One way to do this is to use the boundaries defined by recombination hotspots, proxied by [this map](https://bitbucket.org/nygcresearch/ldetect-data/src/master/) created by [lddetect](https://academic.oup.com/bioinformatics/article/32/2/283/1743626).

How big should a region be?  Big enough that all variants in LD with a lead SNP are included; small enough that only one or a small number of causal variants might exist in it.  I have found using densely genotyped studies and the lddetect boundaries above that regions typically contain 1,000-10,000 SNPs.

# How is coloc.abf accounting for the correlated variants?

That's really how coloc works - by exploiting a dense SNP map - please see the [original paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383)

# How to define priors, is it dependent on sample size or any other parameters?

This is described in detail in the [latest paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008720)

# What does high PP4 mean?

The summary printed on the screen by coloc.abf() and coloc.signals() shows the posterior probability of whether a shared causal variant exists in the region. High PP4 does not mean all variants are causal and shared - to check which variants are most likely to be causal look at the SNP.PP column in the returned detailed results data.frame.

## ---------------------------------------------------------------------------
## The allowlist: every way this branch is permitted to differ from 6.0.3.
##
## One registry of records, not scattered expect_known_difference() calls.
## Scattering spreads the knowledge across every call site, cannot answer
## "list every way we differ", and cannot detect staleness. A registry can
## generate the NEWS entries and the migration guide, and can be tested
## against the migration notes later.
##
## Each record carries:
##   id          AL01..AL13, matching rewrite-design section 6
##   why         the behaviour change in user terms, as long as it needs to be
##   ref         where the 6.0.3 behaviour lives, cited by line at 8f20f0b
##   step        the build step at which this entry starts firing
##   applies_to  function(case_id, layer) -> logical
##   fields      the field names it is allowed to move (see diff_fields())
##   news        the NEWS.md bullet this becomes
##
## Entries are added by the step that CAUSES them, never in advance, so that
## this list is a deliberate act rather than a rubber stamp. Steps 0 and 1 ran
## with the registry empty, which is what proved the harness before it had
## anything to catch. AL10, AL11 and AL12 all fire at step 2, when
## coloc.bf_bf() became a wrapper over the engine and the grid could finally
## see the two implementations disagree.
##
## AL08 is the first entry of the other kind. Steps 2 to 5 only MOVED 6.0.3's
## arithmetic, so the differences they caused were ones the grid found; AL08 is
## a difference we chose to introduce, and it exists from the moment the code
## is written rather than from the moment a case covers it. AL13 is the second.
##
## AL08 and AL13 are also the first entries to OVERLAP: both name
## summary$PP.H0.abf, and AL13 covers every case AL08 does. allowlist_fields()
## unions, so neither can go stale on its own any more. Both directions were
## measured rather than assumed:
##
##   drop AL13 from this list  -> only t0-all-trimmed-null reports an
##                                unexpected difference. The other two stay
##                                quiet, because AL08 still permits the column
##   revert AL08 in the CODE,  -> the grid stays quiet on all three. PP.H0
##   keeping both records         becomes 2.5e-32, noise squared, and nothing
##                                in test-regression.R objects
##
## What catches that second one is test-engine.R, which fails in three tests
## because it names the expected value from an independent closed form rather
## than from the allowlist. Those assertions are load-bearing here, not belt
## and braces. Staleness detection is this registry's whole advantage over
## scattered expect_known_difference() calls, and this one field has lost it.
## ---------------------------------------------------------------------------

coloc_allowlist <- list(

    list(id = "AL08",
         step = 6,
         why = paste("PP.H0 in the all-trimmed branch was an unearned",
                     "certainty. ph0 came from 1 - rowSums(pp), but",
                     "logbf_to_pp() normalises each row over the columns it",
                     "was given, so without a null column that is 1 - 1: not a",
                     "small probability but noise, of either sign. A signal",
                     "whose log BFs are all 0, carrying no evidence for any",
                     "snp, was reported at PP.H0 = 2.2e-16 when its actual",
                     "probability of no causal variant is the prior, 0.9996.",
                     "The engine appends the null at log BF 0, which is what",
                     "these BFs are already measured against. Only the",
                     "no-null-column shape moves: where there is a null column",
                     "the expression is 6.0.3's, unchanged, which is what",
                     "t0-all-trimmed-null is in the grid to demonstrate"),
         ref = "R/susie.R:277-278, reported at :292",
         ## one column, named: nsnps, hit1, hit2, idx1 and idx2 must NOT move,
         ## and PP.H1..H4 are NA on both sides. Naming the column is what makes
         ## the harness say so if any of that changes
         fields = "summary$PP.H0.abf",
         ## every case that reaches the all-trimmed branch, which is these two
         ## and no others: t0-overlap-min gets there by lowering overlap.min
         ## to 0.1 against a trait 2 whose posterior sits almost entirely on
         ## an snp outside the intersection
         applies_to = function(case_id, layer)
             case_id %in% c("t0-all-trimmed", "t0-overlap-min"),
         news = paste("the all-trimmed result now reports a real PP.H0.",
                      "Previously it was ~0 for any input without a null",
                      "column, which is every SuSiE run without null_weight",
                      "and every coloc.abf() route.")),

    list(id = "AL10",
         step = 2,
         why = paste("per-snp p1/p2 with no prior weights: combine.abf() is",
                     "scalar only, so 6.0.3 recycles them and reports 2*Q+3",
                     "hypotheses (PP.H0.abf..PP.H32.abf for Q=8) instead of",
                     "five. The engine routes per-snp priors through",
                     "combine_abf_weighted(), which collapses them with",
                     "logsum(). $priors also differs, but note which way:",
                     "6.0.3 returns list(p1, p2, p12) with p12 reconciled to a",
                     "vector, while the engine reports the scalar p12 it was",
                     "given. Both sides are lists; they are not the same list"),
         ref = "R/claudia.R:123-131 combine.abf(), reached from R/susie.R:346",
         fields = c("summary", "priors"),
         applies_to = function(case_id, layer) case_id == "t0-vector-priors",
         news = paste("coloc.bf_bf() with per-snp p1/p2 and no prior weights now",
                      "returns five posteriors. Previously the per-snp priors",
                      "were recycled and the result carried one PP.Hn.abf",
                      "column per snp per hypothesis.")),

    list(id = "AL11",
         step = 2,
         why = paste("an empty snp intersection now warns. 6.0.3 returned",
                     "data.table(nsnps=NA) silently, which is hard to notice in",
                     "a loop over thousands of gene-trait pairs"),
         ref = "R/susie.R:263-264 in 6.0.3, the bare return with no condition",
         fields = "conditions",
         applies_to = function(case_id, layer) case_id == "t0-no-overlap",
         news = paste("coloc.bf_bf() now warns when the two traits share no",
                      "snps, instead of returning nsnps=NA silently.")),

    list(id = "AL12",
         step = 2,
         why = paste("a signal row that is entirely NA used to come back as a",
                     "CERTAINTY. 6.0.3 sets pp.abf[1:5] <- c(1,0,0,0,0), then",
                     "overwrites it with pp.abf[c(1,3)] <- c(0,1), because",
                     "which.max() on an all-NA row returns nothing. The pair",
                     "therefore reports PP.H2 = 1 with hit1 = '-': certainty",
                     "that only trait 2 is associated, from a trait 1 row that",
                     "carried no information at all. The engine reports NA and",
                     "sets status = degenerate. Only the degenerate PAIR moves;",
                     "the ordinary row beside it and the whole results table",
                     "stay bit-identical. Reachable only with",
                     "trim_by_posterior = FALSE, which is why step 2 passed",
                     "without it: under the default the NA trips 6.0.3's own",
                     "drop guards and it errors before reaching this branch"),
         ref = "R/susie.R:350-359",
         ## named per column, not just "summary": hit2, nsnps, idx1, idx2 and
         ## the whole $results table must NOT move, and listing the columns
         ## individually is what makes the harness say so if they ever do
         fields = c("summary$hit1", paste0("summary$PP.H", 0:4, ".abf")),
         applies_to = function(case_id, layer) case_id == "t0-degenerate",
         news = paste("coloc.bf_bf() with trim_by_posterior = FALSE now returns",
                      "NA for a signal pair whose log Bayes factors are all NA.",
                      "Previously such a pair was reported as PP.H2 = 1, an",
                      "unearned certainty.")),

    list(id = "AL13",
         step = 6,
         why = paste("PP.H0 in the all-trimmed branch was pmin(ph0.1, ph0.2).",
                     "H0 is no causal variant in EITHER trait, so under",
                     "independence it is the product; pmin is an upper bound",
                     "that reports whichever trait is more confidently",
                     "associated and discards the other. On every case in the",
                     "grid the two signal rows collapse to one number, hiding",
                     "that the two signals differ in how likely they are to be",
                     "nothing: t0-all-trimmed reports 0.898988 twice where the",
                     "joint is 0.898629 and 0.809768. Unlike AL08 this moves",
                     "every all-trimmed result, null column or not, which is",
                     "why t0-all-trimmed-null is listed here and why its role",
                     "as AL08's control now rests on an assertion about ph0.1",
                     "and ph0.2 rather than about the output column"),
         ref = "R/susie.R:292",
         fields = "summary$PP.H0.abf",
         applies_to = function(case_id, layer)
             case_id %in% c("t0-all-trimmed", "t0-overlap-min",
                            "t0-all-trimmed-null"),
         news = paste("the all-trimmed result now reports PP.H0 as the joint",
                      "probability that neither trait has a causal variant,",
                      "ph0.1 * ph0.2. Previously it was pmin(ph0.1, ph0.2),",
                      "which reports one trait's marginal for every pair."))

)

names(coloc_allowlist) <- vapply(coloc_allowlist, function(r) r$id, character(1))

## Fields this case is permitted to differ in, for this layer.
allowlist_fields <- function(case_id, layer = "B") {
    hits <- Filter(function(r) isTRUE(r$applies_to(case_id, layer)), coloc_allowlist)
    ## always character(0), never NULL: unique(unlist(list())) is NULL, which
    ## setdiff() happens to tolerate but anything comparing against this does
    ## not
    as.character(unique(unlist(lapply(hits, `[[`, "fields"))))
}

## Used by the "allowlist is documented" test, and later to generate NEWS.
allowlist_table <- function() {
    if (!length(coloc_allowlist))
        return(data.frame(id = character(0), step = integer(0), why = character(0),
                          stringsAsFactors = FALSE))
    data.frame(id   = vapply(coloc_allowlist, `[[`, character(1), "id"),
               step = vapply(coloc_allowlist, `[[`, numeric(1), "step"),
               why  = vapply(coloc_allowlist, `[[`, character(1), "why"),
               stringsAsFactors = FALSE, row.names = NULL)
}

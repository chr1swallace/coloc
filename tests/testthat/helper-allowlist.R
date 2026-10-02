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
##   id          AL01..AL14, matching rewrite-design section 6
##   why         the behaviour change in user terms, as long as it needs to be.
##               Its first sentence labels the check that it still differs
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
## unions them. Both directions were measured rather than assumed, and the
## second changed when the golden comparisons moved from exact to testthat's
## default tolerance:
##
##   drop AL13 from this list  -> only t0-all-trimmed-null reports an
##                                unexpected difference. The other two stay
##                                quiet, because AL08 still permits the column
##   revert AL08 in the CODE,  -> t0-all-trimmed and t0-overlap-min fail their
##   keeping both records         "still differs" check. PP.H0 becomes 2.5e-32,
##                                noise squared, against 6.0.3's noise of about
##                                1e-16. Both are below the tolerance, where the
##                                comparison is absolute, so they count as
##                                equal. Compared exactly, the grid stayed
##                                quiet on all three
##
## That second catch rests on the tolerance calling two kinds of noise equal,
## which is not something to lean on. test-engine.R catches it on purpose: it
## fails in three tests, because it names the expected value from an
## independent closed form rather than from the allowlist. Those assertions
## are load-bearing here, not belt and braces.
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
         why = paste("a vector p1/p2 is rejected. p1, p2 and p12 are the prior",
                     "for any one snp, so each is a single number, and prior",
                     "weights are what vary it across snps. 6.0.3's",
                     "coloc.bf_bf() accepted a vector p1/p2 with no weights and",
                     "passed it to combine.abf(), which is scalar only, so it",
                     "was recycled and the result reported 2*Q+3 hypotheses",
                     "(PP.H0.abf..PP.H32.abf for Q=8) instead of five. coloc()",
                     "now stops in coloc_prior(), and coloc.bf_bf() with the",
                     "same message, so there is no result to compare: every",
                     "field differs"),
         ref = "R/claudia.R:123-131 combine.abf(), reached from R/susie.R:346",
         fields = c("shape", "class", "summary", "results", "priors", "conditions"),
         applies_to = function(case_id, layer) case_id == "t0-vector-priors",
         news = paste("coloc.bf_bf(), coloc.susie() and coloc.susie_bf() now",
                      "reject a vector p1 or p2, as coloc_prior() does; use",
                      "prior weights to vary the prior across snps. Previously",
                      "coloc.bf_bf() recycled the vector and the result carried",
                      "one PP.Hn.abf column per snp per hypothesis.")),

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
         why = paste("log BFs with missing values are rejected. A signal row",
                     "that is entirely NA used to come back as a CERTAINTY:",
                     "6.0.3 sets pp.abf[1:5] <- c(1,0,0,0,0), then overwrites",
                     "it with pp.abf[c(1,3)] <- c(0,1), because which.max() on",
                     "an all-NA row returns nothing, so the pair reports",
                     "PP.H2 = 1 with hit1 = '-': certainty that only trait 2 is",
                     "associated, from a trait 1 row that carried no",
                     "information at all. coloc() now stops, as check_dataset()",
                     "does on a dataset with missing values, so there is no",
                     "result to compare: every field differs. The engine",
                     "itself reports such a pair as NA with status degenerate,",
                     "confined to that pair, which test-engine.R checks and",
                     "coloc.bf_bf() reaches until step 7. Reachable in 6.0.3",
                     "only with trim_by_posterior = FALSE: under the default",
                     "the NA trips its own drop guards and it errors first"),
         ref = "R/susie.R:350-359",
         fields = c("shape", "class", "summary", "results", "priors", "conditions"),
         applies_to = function(case_id, layer) case_id == "t0-degenerate",
         news = paste("coloc() rejects log Bayes factors with missing or",
                      "infinite values. coloc.bf_bf() with",
                      "trim_by_posterior = FALSE now returns NA for a signal",
                      "pair whose log Bayes factors are all NA. Previously such",
                      "a pair was reported as PP.H2 = 1, an unearned certainty.")),

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
                      "which reports one trait's marginal for every pair.")),

    list(id = "AL14",
         step = 6,
         why = paste("coloc() warns that method cannot apply to a trait given as",
                     "precomputed log BFs. 6.0.3's coloc.susie_bf() took a",
                     "dataset and a BF matrix and had no method argument, so",
                     "there was nothing to warn about. The coloc() call that",
                     "replaces it has to name \"susie\" for the dataset, and a",
                     "length-1 method applies to both traits, so it warns",
                     "'trait 2: method is inapplicable to precomputed log BFs'.",
                     "Only the conditions differ; the numbers are 6.0.3's"),
         ref = "R/susie.R:132, coloc.susie_bf(), which takes no method",
         fields = "conditions",
         applies_to = function(case_id, layer)
             case_id == "t2-susie-bf" && layer == "C",
         news = paste("coloc() warns when method is supplied for a trait given",
                      "as precomputed log Bayes factors, where it cannot apply.",
                      "That includes coloc(dataset, bf, method = \"susie\"), the",
                      "replacement for coloc.susie_bf()."))

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

## The entries that permit this field to differ for this case and layer, each
## as its id and the first sentence of its reason: the label of the check that
## the field still differs.
allowlist_reason <- function(case_id, layer, field) {
    hits <- Filter(function(r) isTRUE(r$applies_to(case_id, layer)) && field %in% r$fields,
                   coloc_allowlist)
    paste(vapply(hits, function(r) paste0(r$id, ", ", sub("\\. .*$", ".", r$why)),
                 character(1)), collapse = " ")
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

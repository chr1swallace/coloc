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
##   id          AL01..AL12, matching rewrite-design section 6
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
## anything to catch. Everything here so far fires at step 2, when coloc.bf_bf()
## became a wrapper over the engine and the grid could finally see the two
## implementations disagree.
## ---------------------------------------------------------------------------

coloc_allowlist <- list(

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
                      "unearned certainty."))

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

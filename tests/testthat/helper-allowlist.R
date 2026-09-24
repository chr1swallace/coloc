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
##   id          AL01..AL09, matching rewrite-design section 6
##   why         one line, the behaviour change in user terms
##   ref         where the 6.0.3 behaviour lives
##   step        the build step at which this entry starts firing
##   applies_to  function(case_id, layer) -> logical
##   fields      the field names it is allowed to move (see diff_fields())
##   news        the NEWS.md bullet this becomes
##
## Steps 0 to 3 must run with this registry EMPTY. The first entry fires at
## step 4, when coloc() introduces a second default for p12.
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
                      "snps, instead of returning nsnps=NA silently."))


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

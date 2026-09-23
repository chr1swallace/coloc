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

    ## Nothing yet. Entries are added by the build step that causes them, never
    ## in advance, so that an empty allowlist is a meaningful gate.

)

names(coloc_allowlist) <- vapply(coloc_allowlist, function(r) r$id, character(1))

## Fields this case is permitted to differ in, for this layer.
allowlist_fields <- function(case_id, layer = "B") {
    if (!length(coloc_allowlist)) return(character(0))
    hits <- Filter(function(r) isTRUE(r$applies_to(case_id, layer)), coloc_allowlist)
    unique(unlist(lapply(hits, `[[`, "fields")))
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

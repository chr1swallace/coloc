## ---------------------------------------------------------------------------
## The coloc_result contract.
##
## 6.0.3's problem is that THE SHAPE ENCODES THE OUTCOME: coloc.bf_bf() alone
## has five return shapes, and a consumer cannot know what it is holding until
## it inspects it. Two of those shapes are already broken, failing
## sensitivity()'s stopifnot("priors" %in% names(obj)) with an opaque message.
##
## One class, one shape, outcome in a field. The governing rule:
##
##   Rows exist iff a signal pair was considered.
##   Posteriors are NA iff they could not be computed.
##
## THE SHAPE IS FROZEN HERE. Any later change to it is stop-the-line, which is
## why the contract is a function rather than a roxygen block: a contract in a
## comment drifts, a contract that is asserted on every return cannot.
## ---------------------------------------------------------------------------

## The one summary shape. Zero rows is legal; a different column set is not.
.empty_summary <- function() {
    data.table(nsnps = integer(0), hit1 = character(0), hit2 = character(0),
               PP.H0.abf = numeric(0), PP.H1.abf = numeric(0),
               PP.H2.abf = numeric(0), PP.H3.abf = numeric(0),
               PP.H4.abf = numeric(0),
               idx1 = integer(0), idx2 = integer(0))
}

.coloc_result_statuses <- c("ok", "no_overlap", "all_trimmed", "no_signal", "degenerate")

##' Validate a coloc_result against its frozen contract
##'
##' Called on every return from the engine and asserted across the whole
##' regression grid. Not exported: this is an internal invariant, not an API.
##'
##' @param x a coloc_result
##' @return x invisibly, or an error naming the breach
validate_coloc_result <- function(x) {
    stopifnot(inherits(x, "coloc_result"))
    required <- c("summary", "results", "priors", "status", "messages", "meta")
    missing_fields <- setdiff(required, names(x))
    if (length(missing_fields))
        stop("coloc_result: missing required field(s) ",
             paste(missing_fields, collapse = ", "))
    if (!identical(names(x$summary), names(.empty_summary())))
        stop("coloc_result: $summary columns have drifted from the frozen contract")
    if (!x$status %in% .coloc_result_statuses)
        stop("coloc_result: unknown status '", x$status, "'")
    if (!is.numeric(x$meta$schema))
        stop("coloc_result: $meta$schema must be present, for migrating saved objects")

    ## the governing rule, as assertions
    if (x$status %in% c("no_overlap", "no_signal") && nrow(x$summary))
        stop("coloc_result: status '", x$status, "' must carry zero rows")
    if (x$status == "all_trimmed" && !nrow(x$summary))
        stop("coloc_result: status 'all_trimmed' keeps its rows")
    if (x$status == "ok" && anyNA(x$summary$PP.H4.abf))
        stop("coloc_result: status 'ok' with NA posteriors")

    ## optional, because the results the engine builds for the old entry
    ## points do not carry it; coloc() always sets it
    if (!is.null(x$bf)) {
        is_bf <- function(m) is.matrix(m) && is.numeric(m) && !is.null(colnames(m))
        if (!is.list(x$bf) || !identical(names(x$bf), c("bf1", "bf2")) ||
            !all(vapply(x$bf, is_bf, logical(1))))
            stop("coloc_result: $bf must be list(bf1, bf2) of numeric matrices ",
                 "with snp ids as colnames")
    }
    invisible(x)
}

## Every exit from the engine goes through here, and validates on the way out.
.new_coloc_result <- function(summary = .empty_summary(),
                              results = data.table(snp = character(0)),
                              prior, weights = NULL, status = "ok",
                              messages = character(0), meta = list()) {
    ## 2 from when coloc() kept each trait's log BFs in $bf; a result saved
    ## at schema 1 has none
    meta$schema <- 2L
    out <- list(summary = summary,
                results = results,
                priors  = c(p1 = prior$p1, p2 = prior$p2, p12 = prior$p12),
                status  = status,
                messages = messages,
                meta    = meta)
    ## present iff prior weights were actually used
    if (!is.null(weights)) out$weights <- weights
    validate_coloc_result(structure(out, class = "coloc_result"))
}

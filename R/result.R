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

##' Check a coloc_result has the shape every result must have
##'
##' Called on every return from the engine and asserted across the whole
##' regression grid. A failure is a bug in coloc, not in the caller's input,
##' and the message says so.
##'
##' @param x a coloc_result
##' @return x invisibly, or an error saying what is wrong
##' @noRd
validate_coloc_result <- function(x) {
    stopifnot(inherits(x, "coloc_result"))
    ## every message says what is wrong, what it should be, and whose fault it is
    bug <- function(...)
        stop("internal error in coloc: ", ..., ". This is a bug in coloc; please report it",
             call. = FALSE)
    shown <- function(v) paste0("[", paste(v, collapse = ", "), "]")

    required <- c("summary", "results", "priors", "status", "messages", "meta")
    if (length(setdiff(required, names(x))))
        bug("the result has fields ", shown(names(x)), " but should have at least ",
            shown(required))
    if (!identical(names(x$summary), names(.empty_summary())))
        bug("result $summary has columns ", shown(names(x$summary)),
            " but should have ", shown(names(.empty_summary())))
    if (!x$status %in% .coloc_result_statuses)
        bug("result $status is '", x$status, "' but should be one of ",
            shown(.coloc_result_statuses))
    if (!is.numeric(x$meta$schema))
        bug("result $meta$schema is ",
            if (is.null(x$meta$schema)) "missing" else class(x$meta$schema)[1],
            " but should be a number, the version of the result's shape")

    ## rows exist iff a signal pair was considered; posteriors are NA iff they
    ## could not be computed
    if (x$status %in% c("no_overlap", "no_signal") && nrow(x$summary))
        bug("a result with status '", x$status, "' has ", nrow(x$summary),
            " row(s) in $summary but should have none, as no signal pair was considered")
    if (x$status == "all_trimmed" && !nrow(x$summary))
        bug("a result with status 'all_trimmed' has no rows in $summary but should ",
            "keep one for each signal pair it trimmed")
    if (x$status == "ok" && anyNA(x$summary$PP.H4.abf))
        bug("a result with status 'ok' has NA in $summary$PP.H4.abf, row(s) ",
            shown(which(is.na(x$summary$PP.H4.abf))), ", but should have none")

    ## optional, because the results the engine builds for the old entry
    ## points do not carry it; coloc() always sets it
    if (!is.null(x$bf)) {
        what <- function(m) if (!is.matrix(m)) class(m)[1]
                            else if (!is.numeric(m)) paste(typeof(m), "matrix")
                            else if (is.null(colnames(m))) "matrix without snp ids"
                            else "numeric matrix"
        is_bf <- function(m) identical(what(m), "numeric matrix")
        if (!is.list(x$bf) || !identical(names(x$bf), c("bf1", "bf2")) ||
            !all(vapply(x$bf, is_bf, logical(1))))
            bug("result $bf is ",
                if (is.list(x$bf)) paste0("list(", paste0(names(x$bf), " = ",
                                          vapply(x$bf, what, ""), collapse = ", "), ")")
                else class(x$bf)[1],
                " but should be list(bf1, bf2) of numeric matrices with snp ids as ",
                "column names")
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

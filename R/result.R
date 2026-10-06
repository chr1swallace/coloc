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

##' The result of coloc()
##'
##' [coloc()] returns an object of class `coloc_result`: a list with the
##' fields below. Every result has them, apart from `weights`, and the same
##' `$summary` columns, whatever the inputs were and whatever happened;
##' `$status` says what happened.
##'
##' The posterior probabilities are for five hypotheses about the region:
##' H0, neither trait has a causal variant; H1, only trait 1 does; H2, only
##' trait 2 does; H3, both do, but different variants; H4, both do, and share
##' one causal variant.
##'
##' @section Fields:
##' \describe{
##'   \item{`summary`}{a data.table with one row per pair of signals
##'     compared, one signal from each trait. Its columns are `nsnps`, the
##'     number of SNPs the two traits share; `hit1` and `hit2`, each signal's
##'     top SNP; `PP.H0.abf` to `PP.H4.abf`, the posterior probabilities of
##'     the five hypotheses; and `idx1` and `idx2`, which signal of each trait
##'     the row compares, as listed in `meta$idx1` and `meta$idx2`.}
##'   \item{`results`}{a data.table with one row per shared SNP: `snp`, and
##'     for each row of `summary` the posterior probability that this SNP is
##'     the shared causal variant, if H4 is true. The column is
##'     `SNP.PP.H4.abf` when there is one row, and `SNP.PP.H4.row1`,
##'     `SNP.PP.H4.row2` and so on when there are several.}
##'   \item{`priors`}{the prior probabilities used, `c(p1 =, p2 =, p12 =)`.
##'     For a large region these may be smaller than those given to
##'     [coloc_prior()], and [coloc()] warns when they are.}
##'   \item{`weights`}{present only if prior weights were given:
##'     `list(w1 =, w2 =)`, restricted to the shared SNPs.}
##'   \item{`status`}{what happened; see Status.}
##'   \item{`messages`}{the warnings about the comparison itself (a trait
##'     with no signal, no shared SNPs, pairs left out), as text, so that a
##'     loop over many regions can check them afterwards.}
##'   \item{`meta`}{details of how the result was made: `method`, how each
##'     trait's log Bayes factors were made (`"single"`, `"susie"`, or `NA`
##'     for log Bayes factors given directly); `idx1` and `idx2`, the ids of
##'     each trait's signals (for a SuSiE fit, the SuSiE component behind
##'     each credible set);
##'     `n_pairs_in` and `n_pairs_kept`, the pairs of signals before and after
##'     trimming, and `dropped`, why pairs were left out, when there were any
##'     pairs; and `schema`, the version of this structure.}
##'   \item{`bf`}{each trait's log Bayes factors, `list(bf1 =, bf2 =)`, as
##'     they went into the analysis; see [log_bf()].}
##' }
##'
##' @section Status:
##' \describe{
##'   \item{`"ok"`}{the signal pairs were compared. Pairs left out because
##'     too little of a signal's posterior falls on SNPs the other trait has
##'     (`overlap.min` in [coloc_control()]) have no row; `meta` counts
##'     them.}
##'   \item{`"no_signal"`}{a trait had no signal to compare, for example
##'     because SuSiE found no credible set. `summary` has no rows.}
##'   \item{`"no_overlap"`}{the two traits share no SNPs. `summary` has no
##'     rows.}
##'   \item{`"all_trimmed"`}{every pair was left out under `overlap.min`.
##'     `summary` keeps a row for each, with `PP.H1.abf` to `PP.H4.abf` `NA`,
##'     because telling those apart needs the shared SNPs; `PP.H0.abf` is the
##'     probability that neither trait has a causal variant, from each
##'     trait's own SNPs.}
##'   \item{`"degenerate"`}{the posteriors of some pair could not be
##'     computed, for example from a SuSiE fit whose log Bayes factors are
##'     `NA`. That pair's posteriors are `NA`.}
##' }
##'
##' Two rules hold for every result: `summary` has no rows only when no pair
##' of signals could be formed (`"no_signal"` and `"no_overlap"`), and a
##' posterior probability is `NA` only when it could not be computed.
##'
##' @section What you can do with it:
##' [pp()] gives one row's posterior probabilities as a named vector,
##' [is_ok()] says whether the posteriors were computed, [n_signals()]
##' counts the pairs compared, and [log_bf()] returns a trait's log Bayes
##' factors, to fine-map from the same numbers. `print()` shows a summary.
##' The fields can also be read directly: `summary` and `results` are
##' data.tables.
##'
##' @seealso [coloc()], [pp()], [is_ok()], [n_signals()], [log_bf()]
##' @examples
##' data(coloc_test_data)
##' res <- coloc(coloc_test_data$D1, coloc_test_data$D2)
##' names(res)
##' res$status
##' res$summary
##' pp(res)
##' head(res$results)
##' @name coloc_result
NULL

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
##' Called on every return from the engine, so every test that makes a result
##' asserts it. A failure is a bug in coloc, not in the caller's input, and the
##' message says so.
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

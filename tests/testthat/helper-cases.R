## ---------------------------------------------------------------------------
## The regression case table.
##
## Every case is a deterministic transform of coloc_test_data, or a literal
## matrix, declared as DATA. Nothing here is generated at test time, so the
## same table drives data-raw/make_golden.R and the tests that read its output.
##
## Each case is a list with:
##   id     unique key; names the golden file and keys the allowlist
##   tier   0 literal BF matrices | 1 single-variant route | 2 susie route
##   args   extra arguments passed to coloc.bf_bf() for layer B
##   build  function() returning list(bf1=, bf2=) for tier 0,
##          or list(d1=, d2=) for tiers 1 and 2
##   entry  NULL, or a function(d1, d2) snapshotted end to end as layer C
## ---------------------------------------------------------------------------

## ---- dataset transforms, all deterministic --------------------------------

as_cc <- function(d, s = 0.3) {
    d$type <- "cc"
    d$s <- s
    d$sdY <- NULL
    d
}

to_pvalues <- function(d) {
    d$pvalues <- 2 * pnorm(-abs(d$beta / sqrt(d$varbeta)))
    d$beta <- NULL
    d$varbeta <- NULL
    d$sdY <- NULL
    d
}

drop_sdY <- function(d) {
    d$sdY <- NULL
    d
}

reverse <- function(d) {
    n <- length(d$snp)
    ix <- rev(seq_len(n))
    for (v in c("beta", "varbeta", "pvalues", "MAF", "snp", "position"))
        if (!is.null(d[[v]]) && length(d[[v]]) == n) d[[v]] <- d[[v]][ix]
    if (!is.null(d$LD)) d$LD <- d$LD[d$snp, d$snp, drop = FALSE]
    d
}

## ---- small literal BF matrices for tier 0 ---------------------------------

.snps <- function(n) paste0("s", seq_len(n))

.mat <- function(..., rows = NULL) {
    m <- do.call(rbind, list(...))
    if (is.null(rows)) rows <- paste0("sig", seq_len(nrow(m)))
    dimnames(m) <- list(rows, .snps(ncol(m)))
    m
}

.bf_a <- .mat(c(0, 0, 8, 2, 0, 0, 0, 0),                ## peak at s3
              c(0, 0, 0, 0, 0, 7, 1, 0))                ## peak at s6
.bf_b <- .mat(c(0, 0, 7, 3, 0, 0, 0, 0))                ## peak at s3
.bf_v <- setNames(c(0, 0, 8, 2, 0, 0, 0, 0), .snps(8))  ## named vector

## row 1 carries no information at all. Pairing it against .bf_b is the only
## way to reach 6.0.3 R/susie.R:350-359, and it needs trim_by_posterior = FALSE
## to get there: under the default, the NA reaches 6.0.3's own drop guards and
## errors instead, so there is nothing to freeze.
.bf_na <- .bf_a; .bf_na[1, ] <- NA_real_

.with_null <- function(m, value = 0) {
    out <- cbind(m, null = value)
    rownames(out) <- rownames(m)
    out
}

.rename <- function(m, prefix) {
    colnames(m) <- paste0(prefix, colnames(m))
    m
}

## ---- the table ------------------------------------------------------------

coloc_cases <- list(

    ## ---- tier 0: literal BF matrices, no dataset machinery at all ----------

    list(id = "t0-basic", tier = 0, args = list(),
         build = function() list(bf1 = .bf_a, bf2 = .bf_b)),

    list(id = "t0-vector", tier = 0, args = list(),
         build = function() list(bf1 = .bf_v, bf2 = .bf_v)),

    list(id = "t0-null-both", tier = 0, args = list(),
         build = function() list(bf1 = .with_null(.bf_a), bf2 = .with_null(.bf_b))),

    list(id = "t0-null-one", tier = 0, args = list(),
         build = function() list(bf1 = .with_null(.bf_a), bf2 = .bf_b)),

    ## one row carrying a "null" column: the shape that broke in #210
    list(id = "t0-onerow-null", tier = 0, args = list(),
         build = function() list(bf1 = .with_null(.bf_b), bf2 = .with_null(.bf_b))),

    ## a "null" column that is NOT last, so nulllast() has to move it
    list(id = "t0-null-notlast", tier = 0, args = list(),
         build = function() {
             m <- .with_null(.bf_a)
             list(bf1 = m[, c("null", .snps(8)), drop = FALSE], bf2 = .bf_b)
         }),

    list(id = "t0-no-overlap", tier = 0, args = list(),
         build = function() list(bf1 = .bf_a, bf2 = .rename(.bf_b, "x"))),

    ## bf1 offers only SNPs far from bf2's peak: every pair trimmed
    list(id = "t0-all-trimmed", tier = 0, args = list(),
         build = function() list(bf1 = .bf_a[, 5:8, drop = FALSE], bf2 = .bf_b)),

    list(id = "t0-partial-overlap", tier = 0, args = list(),
         build = function() list(bf1 = .bf_a[, 1:6, drop = FALSE], bf2 = .bf_b)),

    list(id = "t0-trim-off", tier = 0, args = list(trim_by_posterior = FALSE),
         build = function() list(bf1 = .bf_a[, 5:8, drop = FALSE], bf2 = .bf_b)),

    ## An all-NA signal row, which 6.0.3 rewrites to a certainty and the engine
    ## reports as NA. Two rows on purpose: pair (1,1) is the degenerate one and
    ## pair (2,1) is ordinary, so the case pins that the change is per pair and
    ## does not leak into the row beside it. Allowlist AL12.
    list(id = "t0-degenerate", tier = 0, args = list(trim_by_posterior = FALSE),
         build = function() list(bf1 = .bf_na, bf2 = .bf_b)),

    list(id = "t0-overlap-min", tier = 0, args = list(overlap.min = 0.1),
         build = function() list(bf1 = .bf_a[, 4:8, drop = FALSE], bf2 = .bf_b)),

    list(id = "t0-vector-priors", tier = 0,
         args = list(p1 = rep(1e-4, 8), p2 = rep(1e-4, 8)),
         build = function() list(bf1 = .bf_a, bf2 = .bf_b)),

    list(id = "t0-weights", tier = 0,
         args = list(prior_weights1 = c(1, 1, 1, 1, 1, 1, 1, 100),
                     prior_weights2 = c(1, 1, 1, 1, 1, 1, 1, 100)),
         build = function() list(bf1 = .bf_a, bf2 = .bf_b)),

    list(id = "t0-p12", tier = 0, args = list(p12 = 1e-5),
         build = function() list(bf1 = .bf_a, bf2 = .bf_b)),

    ## ---- tier 1: the single-variant route ----------------------------------

    list(id = "t1-quant-beta", tier = 1, args = list(),
         build = function() list(d1 = coloc_test_data$D1, d2 = coloc_test_data$D2),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    list(id = "t1-cc", tier = 1, args = list(),
         build = function() list(d1 = as_cc(coloc_test_data$D1),
                                 d2 = as_cc(coloc_test_data$D2)),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    list(id = "t1-pvalues", tier = 1, args = list(),
         build = function() list(d1 = to_pvalues(coloc_test_data$D1),
                                 d2 = to_pvalues(coloc_test_data$D2)),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    list(id = "t1-no-sdY", tier = 1, args = list(),
         build = function() list(d1 = drop_sdY(coloc_test_data$D1),
                                 d2 = drop_sdY(coloc_test_data$D2)),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    ## partial SNP overlap between the two traits
    list(id = "t1-subset", tier = 1, args = list(),
         build = function() list(d1 = subset_dataset(coloc_test_data$D1, 1:400),
                                 d2 = subset_dataset(coloc_test_data$D2, 101:500)),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    ## the order-invariance case that b63917a fixed
    list(id = "t1-reversed", tier = 1, args = list(),
         build = function() list(d1 = coloc_test_data$D1,
                                 d2 = reverse(coloc_test_data$D2)),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    list(id = "t1-mixed-type", tier = 1, args = list(),
         build = function() list(d1 = coloc_test_data$D1,
                                 d2 = as_cc(coloc_test_data$D2)),
         entry = function(d1, d2) coloc.abf(d1, d2)),

    ## ---- tier 2: the susie route -------------------------------------------

    list(id = "t2-d3d4", tier = 2, args = list(),
         build = function() list(d1 = coloc_test_data$D3, d2 = coloc_test_data$D4),
         entry = function(d1, d2) coloc.susie(d1, d2)),

    ## D4 has exactly one credible set: the dimension-collapse shape from #210
    list(id = "t2-single-cs-first", tier = 2, args = list(),
         build = function() list(d1 = coloc_test_data$D4, d2 = coloc_test_data$D3),
         entry = function(d1, d2) coloc.susie(d1, d2)),

    list(id = "t2-d3d3", tier = 2, args = list(),
         build = function() list(d1 = coloc_test_data$D3, d2 = coloc_test_data$D3),
         entry = function(d1, d2) coloc.susie(d1, d2)),

    ## susie against a raw BF matrix, the existing asymmetric entry point
    list(id = "t2-susie-bf", tier = 2, args = list(),
         build = function() list(d1 = coloc_test_data$D3, d2 = coloc_test_data$D4),
         entry = function(d1, d2) {
             s2 <- runsusie(d2)
             coloc.susie_bf(d1, s2$lbf_variable[s2$sets$cs_index, , drop = FALSE])
         })
)

names(coloc_cases) <- vapply(coloc_cases, `[[`, character(1), "id")

## ---- layer A: the matrices 6.0.3 would have fed the engine ----------------

case_layerA <- function(case) {
    b <- case$build()
    if (case$tier == 0)
        return(list(bf1 = b$bf1, bf2 = b$bf2))
    if (case$tier == 1) {
        ## what coloc.abf() computes internally, as a one-row matrix per trait
        one_row <- function(d, suffix) {
            pd <- process.dataset(d = d, suffix = suffix)
            matrix(pd[[paste0("lABF.", suffix)]], nrow = 1L,
                   dimnames = list("sig1", pd$snp))
        }
        return(list(bf1 = one_row(b$d1, "df1"), bf2 = one_row(b$d2, "df2")))
    }
    ## tier 2: the credible-set rows of each susie fit
    cs_rows <- function(d) {
        s <- runsusie(d)
        s$lbf_variable[s$sets$cs_index, , drop = FALSE]
    }
    list(bf1 = cs_rows(b$d1), bf2 = cs_rows(b$d2))
}

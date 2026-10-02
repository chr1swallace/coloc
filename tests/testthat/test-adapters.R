## ---------------------------------------------------------------------------
## Step 3: the adapters, tested one at a time against the frozen layer A
## matrices. Layer A is what 6.0.3 would have fed the engine, so an adapter
## that reproduces it is a faithful replacement for the code path it lifts.
## ---------------------------------------------------------------------------

test_that("stage 1 classifies from class alone, and the order is load-bearing", {
    ## a class-"susie" stub exercises stage 1 without fitting anything
    fit <- structure(list(), class = "susie")
    m <- matrix(0, 2, 3, dimnames = list(NULL, paste0("s", 1:3)))
    v <- setNames(c(0, 0, 0), paste0("s", 1:3))

    expect_equal(.classify_input(coloc_test_data$D1), "dataset")
    expect_equal(.classify_input(fit), "susie_fit")
    expect_equal(.classify_input(m), "bf")
    expect_equal(.classify_input(v), "bf")

    ## Stage 1 answers KIND only. An unnamed vector is still a vector, so it
    ## classifies and the adapter is what rejects it, with a message that says
    ## what is actually wrong. Asking for names here would send an obviously
    ## vector-shaped input to "not any of these".
    expect_equal(.classify_input(1:3), "bf")
    expect_equal(.classify_input(matrix(0, 2, 2)), "bf")
    expect_error(.adapt_bf(1:3), "must be named by snp")
    expect_error(.adapt_bf(matrix(0, 2, 2)), "colnames")

    ## and the generic message still fires for something genuinely unclassifiable
    expect_error(.classify_input(list(a = 1)), "not a susie fit")

    ## a susie fit is also a list: the class test has to come first
    fit_like_dataset <- structure(list(snp = "s1"), class = "susie")
    expect_equal(.classify_input(fit_like_dataset), "susie_fit")
})

test_that("stage 2 resolves each trait against its own input class", {
    expect_equal(.resolve_method(c("single", "susie"), c("dataset", "dataset")),
                 c("single", "susie"))
    expect_equal(.resolve_method("susie", c("dataset", "dataset")),
                 c("susie", "susie"))

    ## a fit carries no beta/varbeta, so "single" is an error not a fallback
    expect_error(.resolve_method("single", c("susie_fit", "dataset")), "susie fit")

    ## method is inapplicable to precomputed BFs
    expect_warning(m <- .resolve_method("single", c("bf", "dataset")), "inapplicable")
    expect_true(is.na(m[1]))
    expect_equal(m[2], "single")
    expect_silent(.resolve_method("single", c("bf", "dataset"), supplied = FALSE))

    ## the USE.NAMES trap: a named result makes identical(m, "susie") FALSE
    expect_null(names(.resolve_method(c("single", "susie"), c("dataset", "dataset"))))
})

test_that(".adapt_single reproduces the matrices coloc.abf() feeds the engine", {
    for (case in coloc_cases) {
        if (case$tier != 1) next
        g <- read_golden(case$id)
        b <- case$build()
        a1 <- suppressWarnings(.adapt_single(b$d1))
        a2 <- suppressWarnings(.adapt_single(b$d2))
        expect_equal(a1$bf, g$layerA$bf1, label = paste0(case$id, " bf1"))
        expect_equal(a2$bf, g$layerA$bf2, label = paste0(case$id, " bf2"))
        expect_equal(a1$meta$method, "single")
        expect_equal(nrow(a1$bf), 1L)
        expect_false(a1$meta$has_null)
    }
})

test_that(".adapt_susie reproduces the credible-set rows", {
    skip_if_susieR_differs()
    for (case in coloc_cases) {
        if (case$tier != 2) next
        g <- read_golden(case$id)
        b <- case$build()
        a1 <- suppressMessages(.adapt_susie(b$d1))
        ## values and snps must match 6.0.3; the row NAMES are a deliberate
        ## addition, since the adapter contract is
        ## dimnames = list(signal_id, snp) and lbf_variable carries none
        got <- a1$bf; rownames(got) <- NULL
        expect_equal(got, g$layerA$bf1, label = paste0(case$id, " bf1"))
        expect_equal(rownames(a1$bf), paste0("sig", a1$meta$idx))
        expect_equal(a1$meta$method, "susie")
    }
})

test_that(".adapt_susie selects the credible-set rows", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    ## the row-selection policy, made visible: ten components, two sets
    fit <- suppressMessages(runsusie(coloc_test_data$D3))
    sel <- .adapt_susie(fit)
    all_rows <- .adapt_bf(fit$lbf_variable)
    expect_equal(nrow(sel$bf), length(fit$sets$cs_index))
    expect_equal(nrow(all_rows$bf), nrow(fit$lbf_variable))
    expect_lt(nrow(sel$bf), nrow(all_rows$bf))
    expect_equal(sel$meta$idx, fit$sets$cs_index)
})

test_that(".adapt_susie survives a fit with no credible sets", {
    fit <- structure(list(lbf_variable = matrix(0, 3, 4,
                                                dimnames = list(NULL, paste0("s", 1:4))),
                          sets = list(cs_index = integer(0))),
                     class = "susie")
    a <- .adapt_susie(fit)
    expect_equal(nrow(a$bf), 0L)
    expect_equal(colnames(a$bf), paste0("s", 1:4))
    ## and the engine turns that into a status rather than an error
    res <- suppressWarnings(.coloc_engine(a$bf, .adapt_bf(fit$lbf_variable)$bf))
    expect_equal(res$status, "no_signal")
})

test_that(".adapt_bf validates and takes every row", {
    v <- setNames(c(1, 2, 3), paste0("s", 1:3))
    a <- .adapt_bf(v)
    expect_equal(dim(a$bf), c(1L, 3L))
    expect_equal(rownames(a$bf), "sig1")
    expect_true(is.na(a$meta$method))

    m <- matrix(0, 3, 2, dimnames = list(NULL, c("s1", "null")))
    expect_true(.adapt_bf(m)$meta$has_null)
    expect_equal(.adapt_bf(m)$meta$idx, 1:3)

    expect_error(.adapt_bf(c(1, 2, 3)), "named by snp")
    expect_error(.adapt_bf(matrix(0, 2, 2)), "colnames")

    ## is.matrix() is true of a character matrix, so without a type check it
    ## reaches the engine and dies on arithmetic instead
    chr <- matrix("a", 2, 3, dimnames = list(NULL, paste0("s", 1:3)))
    expect_error(.adapt_bf(chr), "must be numeric")
    expect_error(.adapt_bf(setNames(letters[1:3], paste0("s", 1:3))), "must be numeric")
    ## and the same via coloc(), so a character VECTOR gets the type message
    ## rather than the generic one a character matrix already avoided
    expect_equal(.classify_input(setNames(letters[1:3], paste0("s", 1:3))), "bf")
})

test_that("coloc() reports what is wrong with a bad BF input, not just that it is bad", {
    bf <- matrix(0, 1, 3, dimnames = list("sig1", paste0("s", 1:3)))
    expect_error(coloc(c(1, 2, 3), bf), "must be named by snp")
    expect_error(coloc(matrix(0, 2, 3), bf), "colnames")
    expect_error(coloc(matrix("a", 1, 3, dimnames = list(NULL, paste0("s", 1:3))), bf),
                 "must be numeric")
})

test_that("the old doors keep 6.0.3's behaviour on the same bad inputs", {
    ## D6 freezes coloc.bf_bf(), which never calls .adapt_bf(): it does its own
    ## is.vector promotion at R/susie.R:247-250. A character matrix must still
    ## fail the way it always has, deep in the arithmetic. If this test starts
    ## failing, the tightening above has leaked into a frozen door.
    chr <- matrix("a", 2, 3, dimnames = list(NULL, paste0("s", 1:3)))
    num <- matrix(0,   2, 3, dimnames = list(NULL, paste0("s", 1:3)))
    expect_error(coloc.bf_bf(chr, num), "non-numeric argument")
})

test_that("meta$has_null never contradicts the matrix it describes", {
    ## The zero-credible-set branch keeps the fit's colnames, so it can carry a
    ## "null" column just like the normal branch. It used to hardcode FALSE,
    ## which nothing noticed because the engine detects the null column itself
    ## rather than trusting meta. That makes this a latent wrongness: it only
    ## bites once something starts reading the field.
    mk <- function(snps, empty) structure(list(
        lbf_variable = matrix(0, 3, length(snps), dimnames = list(NULL, snps)),
        sets = list(cs = if (empty) list() else list(1),
                    cs_index = if (empty) integer(0) else 1L)), class = "susie")

    for (snps in list(c("s1", "s2", "null"), c("s1", "s2", "s3")))
        for (empty in c(TRUE, FALSE)) {
            a <- .adapt_susie(mk(snps, empty))
            expect_equal(a$meta$has_null, "null" %in% colnames(a$bf),
                         label = paste0("snps=", paste(snps, collapse = ","),
                                        " empty=", empty))
        }

    ## .adapt_bf() likewise describes what it was handed
    for (cn in list(c("s1", "null"), c("s1", "s2"))) {
        m <- matrix(0, 2, 2, dimnames = list(NULL, cn))
        expect_equal(.adapt_bf(m)$meta$has_null, "null" %in% cn)
    }
})

test_that("has_null is detected where it can be, and FALSE by construction where it cannot", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    ## runsusie() only produces a "null" column when susie_rss() was given a
    ## null_weight (annotate_susie(), R/susie.R:502); on the default route it
    ## does not, so the engine's null handling stays inert
    fit <- suppressMessages(runsusie(coloc_test_data$D4))
    expect_false(.adapt_susie(fit)$meta$has_null)

    ## .adapt_single() asserts FALSE rather than detecting, and that is correct:
    ## process.dataset() never appends a null pseudo-snp. finemap.abf() does
    ## (R/claudia.R:329-333), coloc.abf() does not, and this reproduces
    ## coloc.abf(). So the claim is about construction, not inspection.
    a <- suppressWarnings(.adapt_single(coloc_test_data$D1))
    expect_false(a$meta$has_null)
    expect_false("null" %in% colnames(a$bf))
})

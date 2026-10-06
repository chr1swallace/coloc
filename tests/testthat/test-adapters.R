## The adapters, one at a time. Each hands the engine the log BFs that 6.0.3
## computed on the code path it lifts: process.dataset() for .adapt_single(),
## the credible-set rows of a SuSiE fit for .adapt_susie().

snps <- paste0("s", 1:3)

## a susie fit made by hand, with or without credible sets
stub_fit <- function(cols, empty = FALSE)
    structure(list(lbf_variable = matrix(0, 3, length(cols), dimnames = list(NULL, cols)),
                   sets = list(cs = if (empty) list() else list(1),
                               cs_index = if (empty) integer(0) else 1L)),
              class = "susie")

test_that(".classify_input() tells the kind of input from its class", {
    expect_equal(.classify_input(coloc_test_data$D1), "dataset")
    expect_equal(.classify_input(stub_fit(snps)), "susie_fit")
    expect_equal(.classify_input(matrix(0, 2, 3, dimnames = list(NULL, snps))), "bf")
    expect_equal(.classify_input(setNames(c(0, 0, 0), snps)), "bf")
    ## a fit is also a list with a snp element: the class test comes first
    expect_equal(.classify_input(structure(list(snp = "s1"), class = "susie")), "susie_fit")
    ## kind only: an unnamed vector is still log BFs, for .adapt_bf() to reject
    expect_equal(.classify_input(1:3), "bf")
    expect_equal(.classify_input(setNames(letters[1:3], snps)), "bf")
    expect_error(.classify_input(list(a = 1)), "not a susie fit")
})

test_that(".resolve_method() resolves each trait against its own input", {
    expect_equal(.resolve_method("susie", c("dataset", "dataset")), c("susie", "susie"))
    expect_equal(.resolve_method(c("single", "susie"), c("dataset", "dataset")),
                 c("single", "susie"))
    ## vapply would name the result after the values, and break identical()
    expect_null(names(.resolve_method(c("single", "susie"), c("dataset", "dataset"))))
    ## a fit implies "susie"; only a "single" the caller asked for is an error
    expect_error(.resolve_method("single", c("susie_fit", "dataset")), "susie fit")
    expect_equal(.resolve_method("single", c("susie_fit", "dataset"), supplied = FALSE),
                 c("susie", "single"))
    ## log BFs take no method, and naming one warns
    expect_warning(m <- .resolve_method("single", c("bf", "dataset")), "inapplicable")
    expect_equal(m, c(NA, "single"))
    expect_silent(.resolve_method("single", c("bf", "dataset"), supplied = FALSE))
})

test_that(".adapt_single() hands the engine process.dataset()'s log ABFs", {
    for (trait in 1:2) {
        d <- coloc_test_data[[trait]]                 # D1, D2
        pd <- process.dataset(d, suffix = paste0("df", trait))
        a <- .adapt_single(d, trait)
        expect_equal(a$bf, matrix(pd[[paste0("lABF.df", trait)]], nrow = 1,
                                  dimnames = list("sig1", pd$snp)))
        expect_equal(a$meta, list(method = "single", idx = 1L, has_null = FALSE))
    }
})

test_that(".adapt_single() checks the dataset first, naming the trait", {
    bad <- coloc_test_data$D1
    bad$beta[1] <- NA
    expect_error(.adapt_single(bad, 2), "dataset 2: beta contains missing values")
})

test_that(".adapt_susie() hands the engine the credible-set rows", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    for (d in coloc_test_data[c("D3", "D4")]) {
        fit <- suppressMessages(runsusie(d))
        a <- suppressMessages(.adapt_susie(d))
        expect_equal(unname(a$bf), unname(fit$lbf_variable[fit$sets$cs_index, , drop = FALSE]))
        expect_equal(colnames(a$bf), colnames(fit$lbf_variable))
        ## lbf_variable has no row names; the adapter names rows by component
        expect_equal(rownames(a$bf), paste0("sig", fit$sets$cs_index))
        expect_equal(a$meta$idx, fit$sets$cs_index)
        expect_equal(a$meta$method, "susie")
        ## the same from the fit itself, and fewer rows than the fit has
        expect_equal(.adapt_susie(fit), a)
        expect_lt(nrow(a$bf), nrow(fit$lbf_variable))
    }
})

test_that(".adapt_susie() takes a fit with no credible sets", {
    a <- .adapt_susie(stub_fit(snps, empty = TRUE))
    expect_equal(dim(a$bf), c(0L, 3L))
    expect_equal(colnames(a$bf), snps)
    expect_equal(a$meta$idx, integer(0))
})

test_that(".adapt_susie() rejects a broken fit, and a suffix in susie.args", {
    no_lbf <- structure(list(sets = list(cs_index = 1L)), class = "susie")
    no_names <- structure(list(lbf_variable = matrix(0, 2, 4), sets = list(cs_index = 1L)),
                          class = "susie")
    expect_error(.adapt_susie(no_lbf, 2), "trait 2: the susie fit has no lbf_variable")
    expect_error(.adapt_susie(no_names, 1), "trait 1: the susie fit's lbf_variable needs snp ids")
    ## coloc() sets suffix to the trait number
    expect_error(.adapt_susie(coloc_test_data$D3, susie.args = list(suffix = 3)),
                 "susie.args cannot set suffix")
})

test_that(".adapt_bf() takes a named vector or a matrix, every row", {
    v <- .adapt_bf(setNames(c(1, 2, 3), snps))
    expect_equal(v$bf, matrix(c(1, 2, 3), nrow = 1, dimnames = list("sig1", snps)))
    expect_true(is.na(v$meta$method))
    m <- .adapt_bf(matrix(0, 3, 2, dimnames = list(NULL, c("s1", "s2"))))
    expect_equal(rownames(m$bf), c("sig1", "sig2", "sig3"))
    expect_equal(m$meta$idx, 1:3)
    ## a single snp is still a matrix
    expect_equal(dim(.adapt_bf(c(s1 = 2))$bf), c(1L, 1L))
})

test_that(".adapt_bf() takes a matrix with no rows, which coloc() reports as no signal", {
    none <- matrix(0, 0, 3, dimnames = list(NULL, snps))
    expect_equal(dim(.adapt_bf(none)$bf), c(0L, 3L))
    bf <- matrix(c(0, 8, 0), nrow = 1, dimnames = list("sig1", snps))
    expect_warning(res <- coloc(none, bf), "no signals")
    expect_equal(res$status, "no_signal")
})

test_that(".adapt_bf() rejects what are not log BFs, naming the trait", {
    good <- matrix(c(0, 8, 2), nrow = 1, dimnames = list("sig1", snps))
    na <- good; na[1, 2] <- NA
    inf <- good; inf[1, 2] <- Inf
    dup <- good; colnames(dup)[2] <- "s1"
    expect_error(.adapt_bf(c(1, 2, 3), 2), "trait 2: a log-BF vector must be named by snp")
    expect_error(.adapt_bf(matrix(0, 2, 2), 2), "trait 2: a log-BF matrix needs snp ids")
    ## is.matrix() is TRUE of a character matrix too
    expect_error(.adapt_bf(matrix("a", 1, 3, dimnames = list(NULL, snps)), 2),
                 "trait 2: log Bayes factors must be numeric")
    expect_error(.adapt_bf(na, 2), "trait 2: log Bayes factors contain missing or infinite")
    expect_error(.adapt_bf(inf, 2), "trait 2: log Bayes factors contain missing or infinite")
    expect_error(.adapt_bf(dup, 2), "trait 2: duplicated snps found")
})

test_that("meta$has_null says whether the matrix has a null column", {
    for (cols in list(c("s1", "s2", "null"), snps)) {
        expect_equal(.adapt_bf(matrix(0, 1, 3, dimnames = list(NULL, cols)))$meta$has_null,
                     "null" %in% cols)
        for (empty in c(TRUE, FALSE))
            expect_equal(.adapt_susie(stub_fit(cols, empty))$meta$has_null, "null" %in% cols)
    }
    ## runsusie() adds a null column only when given null_weight
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    expect_false(.adapt_susie(suppressMessages(runsusie(coloc_test_data$D4)))$meta$has_null)
})

## ---------------------------------------------------------------------------
## Step 3: the adapters, tested one at a time against the frozen layer A
## matrices. Layer A is what 6.0.3 would have fed the engine, so an adapter
## that reproduces it is a faithful replacement for the code path it lifts.
## The frozen matrices are in inst/testdata/golden/; get_match_6_0_3() in
## helper-golden.R reads one.
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
        a1 <- suppressWarnings(.adapt_single(b$d1, trait = 1))
        a2 <- suppressWarnings(.adapt_single(b$d2, trait = 2))
        expect_bf_equal(a1$bf, g$layerA$bf1, paste0(case$id, " bf1"))
        expect_bf_equal(a2$bf, g$layerA$bf2, paste0(case$id, " bf2"))
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
        expect_bf_equal(got, g$layerA$bf1, paste0(case$id, " bf1"))
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

## ---- coloc() checks its inputs ---------------------------------------------

## One mistake each, made to a copy of a good dataset, so the message is
## about this mistake and no other.
.mistakes <- list(
    "no type"            = function(d) { d$type <- NULL; d },
    "a wrong type"       = function(d) { d$type <- "quantitative"; d },
    "an NA beta"         = function(d) { d$beta[5] <- NA; d },
    "a duplicated snp"   = function(d) { d$snp[2] <- d$snp[1]; d },
    "mismatched lengths" = function(d) { d$beta <- d$beta[-1]; d },
    "a zero varbeta"     = function(d) { d$varbeta[5] <- 0; d })

test_that("coloc() stops on a bad dataset with coloc.abf()'s message, for either trait", {
    ## .adapt_single() runs check_dataset() before process.dataset(), as
    ## coloc.abf() does. Without it, one NA in beta came back as NA posteriors
    ## with a named hit, where coloc.abf() stops.
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    for (what in names(.mistakes)) for (trait in 1:2) {
        args <- list(D1, D2)
        args[[trait]] <- .mistakes[[what]](args[[trait]])
        expected <- tryCatch(do.call(coloc.abf, args), error = conditionMessage)
        expect_match(expected, paste0("dataset ", trait), fixed = TRUE,
                     label = paste0("coloc.abf()'s message for ", what, " in trait ", trait))
        expect_error(do.call(coloc, c(args, method = "single")), expected, fixed = TRUE,
                     label = paste0("coloc() on ", what, " in trait ", trait))
    }
})

test_that("on the SuSiE route a bad dataset is reported with its own trait number", {
    ## runsusie() checks the dataset, and coloc() now tells it which trait it
    ## is: it used to report a bad trait 2 as dataset 1. A dataset with no type
    ## used to fail on a zero-length condition, because runsusie() leaves type
    ## out of check_dataset()'s req.
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    fit3 <- suppressMessages(runsusie(D3)); fit4 <- suppressMessages(runsusie(D4))
    problem <- c("no type"            = "variable type not set",
                 "a wrong type"       = "type must be quant or cc",
                 "an NA beta"         = "beta contains missing values",
                 "a duplicated snp"   = "duplicated snps found",
                 "mismatched lengths" = "lengths of inputs don't match",
                 "a zero varbeta"     = "Zero values in varbeta")
    for (what in names(.mistakes)) {
        expect_error(suppressMessages(coloc(.mistakes[[what]](D3), fit4, method = "susie")),
                     paste0("dataset 1: ", problem[[what]]), fixed = TRUE,
                     label = paste0(what, " in trait 1"))
        expect_error(suppressMessages(coloc(fit3, .mistakes[[what]](D4), method = "susie")),
                     paste0("dataset 2: ", problem[[what]]), fixed = TRUE,
                     label = paste0(what, " in trait 2"))
    }
})

test_that("susie.args cannot set suffix, which coloc() sets per trait", {
    D3 <- coloc_test_data$D3; D4 <- coloc_test_data$D4
    expect_error(coloc(D3, D4, method = "susie", susie.args = list(suffix = 3)),
                 "susie.args cannot set suffix")
})

test_that("coloc() stops on bad log BFs and a bad susie fit, naming the trait", {
    good <- matrix(c(0, 0, 8, 2), nrow = 1, dimnames = list("sig1", paste0("s", 1:4)))
    na  <- good; na[1, 2]  <- NA
    inf <- good; inf[1, 2] <- Inf
    dup <- good; colnames(dup)[2] <- "s1"
    for (trait in 1:2) {
        pair <- function(bad) if (trait == 1) list(bad, good) else list(good, bad)
        lead <- paste0("trait ", trait, ": ")
        expect_error(do.call(coloc, pair(na)),
                     paste0(lead, "log Bayes factors contain missing"), fixed = TRUE)
        expect_error(do.call(coloc, pair(inf)),
                     paste0(lead, "log Bayes factors contain missing or infinite"), fixed = TRUE)
        expect_error(do.call(coloc, pair(dup)),
                     paste0(lead, "duplicated snps found"), fixed = TRUE)
    }

    no_lbf   <- structure(list(sets = list(cs_index = 1L)), class = "susie")
    no_names <- structure(list(lbf_variable = matrix(0, 2, 4), sets = list(cs_index = 1L)),
                          class = "susie")
    ## a fit needs method = "susie", which then warns that it cannot apply to
    ## the log BFs beside it (AL14); that warning is not what is tested here
    expect_error(.without_method_warning(coloc(no_lbf, good, method = "susie")),
                 "trait 1: the susie fit has no lbf_variable", fixed = TRUE)
    expect_error(.without_method_warning(coloc(good, no_names, method = "susie")),
                 "trait 2: the susie fit's lbf_variable needs snp ids", fixed = TRUE)
})

test_that("the region-size guard caps a prior as coloc.abf() does", {
    ## 500 snps at p1 = p2 = 0.01 imply five causal variants per trait, which
    ## leaves nothing for "no causal variant". coloc.abf() warns and caps each
    ## at 1/(nsnps + 1), and coloc() now does the same, once the snps are known
    D1 <- coloc_test_data$D1; D2 <- coloc_test_data$D2
    warnings_of <- function(expr) {
        w <- character(0)
        withCallingHandlers(expr, warning = function(c) {
            w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
        w
    }
    new_w <- warnings_of(new <- coloc(D1, D2, method = "single",
                                      prior = coloc_prior(p1 = 0.01, p2 = 0.01, p12 = 1e-4)))
    invisible(utils::capture.output(
        old_w <- warnings_of(old <- coloc.abf(D1, D2, p1 = 0.01, p2 = 0.01, p12 = 1e-4))))
    expect_equal(new_w, old_w)
    expect_match(new_w, "^p[12] \\* nsnps >= 1")
    expect_equal(new$priors, old$priors)
    expect_equal(new$priors[["p1"]], 1 / 501)
    expect_equal(pp(new), old$summary)

    ## and on log BFs, where coloc.bf_bf() never applied it (allowlist AL03)
    big <- matrix(0, 1, 20000, dimnames = list("sig1", paste0("s", 1:20000)))
    w <- warnings_of(res <- coloc(big, big))
    expect_true(any(grepl("p1 * nsnps >= 1", w, fixed = TRUE)))
    expect_true(any(grepl("p2 * nsnps >= 1", w, fixed = TRUE)))
    expect_equal(res$priors[["p1"]], 1 / 20001)
})

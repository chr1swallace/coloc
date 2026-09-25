## ---------------------------------------------------------------------------
## The regression grid: every case compared against frozen 6.0.3 output.
##
## D6 makes 6.0.3 the specification. Layer B replays the engine on the frozen
## bf1/bf2 matrices, so it touches no susieR and no process.dataset() and runs
## in milliseconds. Layer C goes end to end through the 6.0.3 entry points.
## ---------------------------------------------------------------------------

test_that("the goldens are present and pinned by commit, not version", {
    m <- read_manifest()
    expect_equal(m$coloc_sha, "8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45")
    expect_setequal(m$cases, names(coloc_cases))
    for (id in names(coloc_cases))
        expect_true(file.exists(file.path(golden_dir(), paste0(id, ".rds"))),
                    label = paste0("golden for ", id))
})

test_that("layer B: the engine reproduces 6.0.3 on frozen BF matrices", {
    for (case in coloc_cases) {
        g <- read_golden(case$id)
        got <- capture_run(do.call(coloc.bf_bf,
                                   c(list(g$layerA$bf1, g$layerA$bf2), case$args)))
        expect_matches_6_0_3(got, case$id, layer = "B")
    }
})

test_that("layer C: the old entry points reproduce 6.0.3 end to end", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    for (case in coloc_cases) {
        if (is.null(case$entry)) next
        b <- case$build()
        got <- capture_run(case$entry(b$d1, b$d2))
        expect_matches_6_0_3(got, case$id, layer = "C")
    }
})

test_that("layer A: the BF matrices fed to the engine are still reproducible", {
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    for (case in coloc_cases) {
        g <- read_golden(case$id)
        a <- suppressWarnings(suppressMessages(case_layerA(case)))
        expect_equal(a$bf1, g$layerA$bf1, tolerance = golden_tolerance,
                     label = paste0(case$id, " bf1"))
        expect_equal(a$bf2, g$layerA$bf2, tolerance = golden_tolerance,
                     label = paste0(case$id, " bf2"))
    }
})

test_that("the allowlist holds exactly the entries the build has earned", {
    ## Entries are added by the step that causes them, never in advance, so
    ## that this list is a deliberate act rather than a rubber stamp. Every id
    ## here is a behaviour change someone chose and can defend.
    expect_setequal(allowlist_table()$id, c("AL10", "AL11", "AL12"))
    for (r in coloc_allowlist) {
        expect_true(nzchar(r$why),  label = paste0(r$id, " has a reason"))
        expect_true(nzchar(r$ref),  label = paste0(r$id, " names the 6.0.3 site"))
        expect_true(nzchar(r$news), label = paste0(r$id, " has a NEWS bullet"))
        expect_true(is.function(r$applies_to))
    }
})

test_that("tripwire: no fixture reaches susieR's purity subsampling", {
    ## get_purity() is the ONLY RNG in the susie_rss -> susie_get_cs path, and
    ## it only samples when a credible set exceeds n_purity = 100 SNPs. If a
    ## bigger fixture is ever added, the harness should say it has become
    ## non-deterministic rather than going quietly flaky.
    skip_if_not(requireNamespace("susieR", quietly = TRUE), "susieR not installed")
    threshold <- eval(formals(susieR::susie_get_cs)$n_purity)
    for (nm in c("D1", "D2", "D3", "D4")) {
        s <- suppressMessages(runsusie(coloc_test_data[[nm]]))
        sizes <- vapply(s$sets$cs, length, integer(1))
        expect_true(all(sizes < threshold),
                    label = paste0(nm, " largest credible set (", max(sizes), ")"))
    }
})

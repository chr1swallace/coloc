## ---------------------------------------------------------------------------
## Generate the 6.0.3 golden files.
##
## Committed .rds files are the truth; a 6.0.3 install (or the branch point
## itself) is only the machine that produces them. Run this ONCE at the branch
## point, before anything in R/ changes, or from a rebuilt oracle afterwards:
##
##   Rscript data-raw/make_golden.R                    # from a clean branch point
##   COLOC_ORACLE_LIB=/path/to/lib Rscript data-raw/make_golden.R
##
## PIN BY COMMIT, NEVER BY VERSION. Three commits report Version: 6.0.3 and two
## of them differ in behaviour, so packageVersion() is not a sufficient check.
## ---------------------------------------------------------------------------

ORACLE_SHA <- "8f20f0bc5e60ffc99e4c2f787bd55fd30cfe7c45"   ## 6.0.3, "Fix #210"

pkg_root <- normalizePath(".")
stopifnot(file.exists(file.path(pkg_root, "DESCRIPTION")))

oracle_lib <- Sys.getenv("COLOC_ORACLE_LIB", "")
if (nzchar(oracle_lib)) {
    message("oracle: installed library at ", oracle_lib)
    library(coloc, lib.loc = oracle_lib)
    sha <- Sys.getenv("COLOC_ORACLE_SHA", "")
} else {
    message("oracle: the working tree, via pkgload::load_all()")
    sha <- system2("git", c("-C", pkg_root, "rev-parse", "HEAD"), stdout = TRUE)
    ## load_all() is only 6.0.3 while R/ and DESCRIPTION are untouched
    dirty <- system2("git", c("-C", pkg_root, "diff", "--quiet", "HEAD", "--",
                              "R/", "DESCRIPTION"))
    if (dirty != 0)
        stop("R/ or DESCRIPTION differ from HEAD: the working tree is no longer the oracle.\n",
             "Rebuild one with data-raw/build_oracle.sh and set COLOC_ORACLE_LIB.")
    suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))
}

if (!identical(sha, ORACLE_SHA))
    stop("oracle sha mismatch.\n  expected ", ORACLE_SHA, "\n  got      ", sha)

## the case table and capture_run() live with the tests
source(file.path(pkg_root, "tests/testthat/helper-cases.R"))
source(file.path(pkg_root, "tests/testthat/helper-golden.R"))

data(coloc_test_data, package = "coloc", envir = environment())

out_dir <- file.path(pkg_root, "inst/testdata/golden")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ---- generate --------------------------------------------------------------

set.seed(42)          ## belt and braces only; see rewrite-design section 7 on
                      ## why a seed is NOT the determinism mechanism here

written <- character(0)
for (case in coloc_cases) {
    message("case ", case$id)
    built <- case$build()
    layerA <- case_layerA(case)

    ## layer B: what the engine sees and returns
    layerB <- capture_run(do.call(coloc.bf_bf,
                                  c(list(layerA$bf1, layerA$bf2), case$args)))

    ## layer C: end to end through the 6.0.3 entry point, tiers 1 and 2 only
    layerC <- NULL
    if (!is.null(case$entry))
        layerC <- capture_run(case$entry(built$d1, built$d2))

    saveRDS(list(id = case$id, tier = case$tier, args = case$args,
                 layerA = layerA, layerB = layerB, layerC = layerC),
            file.path(out_dir, paste0(case$id, ".rds")), version = 3)
    written <- c(written, case$id)
}

## ---- the manifest ----------------------------------------------------------

## Frozen files capture the environment exactly once. This is what makes
## "why did this move" an answerable question later.
manifest <- list(
    coloc_sha    = ORACLE_SHA,
    coloc_ver    = as.character(utils::packageVersion("coloc")),
    generated    = Sys.time(),
    r_version    = R.version.string,
    platform     = R.version$platform,
    la_library   = La_library(),
    susieR       = tryCatch(as.character(utils::packageVersion("susieR")),
                            error = function(e) NA_character_),
    data.table   = as.character(utils::packageVersion("data.table")),
    session      = utils::capture.output(utils::sessionInfo()),
    generator_md5 = tools::md5sum(file.path(pkg_root, "data-raw/make_golden.R"))[[1]],
    cases        = written)
saveRDS(manifest, file.path(out_dir, "MANIFEST.rds"), version = 3)

message("\nwrote ", length(written), " goldens to ", out_dir)

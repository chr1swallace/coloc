## ---------------------------------------------------------------------------
## DESCRIPTION's Collate: field lists every file in R/, in the order R loads
## them when it installs the package, and is kept by hand (no @include tags).
## If a file in R/ is missing from it, or listed but gone, R CMD INSTALL stops;
## but devtools::load_all(), which devtools::test() uses, loads the files
## anyway, so the mistake would only show once the package is installed or
## checked. This test makes it show in devtools::test().
## ---------------------------------------------------------------------------

test_that("Collate: covers exactly the files in R/", {
    root <- NULL
    for (cand in c(".", "../..", "../../.."))
        if (file.exists(file.path(cand, "DESCRIPTION")) && dir.exists(file.path(cand, "R")))
            { root <- cand; break }
    skip_if(is.null(root), "not running from the source tree")

    d <- read.dcf(file.path(root, "DESCRIPTION"))
    skip_if_not("Collate" %in% colnames(d), "no Collate: field")

    collated <- scan(text = gsub("'", "\"", d[1, "Collate"]),
                     what = character(), quiet = TRUE)
    on_disk <- list.files(file.path(root, "R"), pattern = "\\.R$")

    expect_setequal(collated, on_disk)
})

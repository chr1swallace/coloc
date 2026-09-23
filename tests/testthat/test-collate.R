## ---------------------------------------------------------------------------
## Collate: is hand-ordered with no @include tags, so every file added to R/
## must be entered there by hand or the package will not build. That is the
## kind of small thing that makes a rename commit fail CI and get reverted,
## so it gets a test rather than a code comment.
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

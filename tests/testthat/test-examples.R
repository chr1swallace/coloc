test_that("examples", {
    ## the plotting examples draw on a null device, not into Rplots.pdf here
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    test_examples(path="../..")
})

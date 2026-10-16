library(coloc)
data(coloc_test_data)
D1 <- coloc_test_data$D1
D2 <- coloc_test_data$D2

## a small dataset with one strong signal, so no p value warning
small <- list(snp = paste0("s", 1:5), type = "quant", beta = c(0.1, 0.5, 0.2, 0, 0.05), varbeta = rep(0.01, 5), sdY = 1)

test_that("coloc_dataset() makes a classed list with the elements it was given", {
    D <- do.call(coloc_dataset, small)
    expect_s3_class(D, "coloc_dataset")
    expect_true(is.list(D))
    expect_setequal(names(D), names(small)) # no NULL elements for arguments not given
    expect_identical(D$beta, small$beta)
})

test_that("coloc_dataset() and as_coloc_dataset() agree", {
    ## coloc_dataset() puts the elements in its own order, so compare them by name
    expect_identical(as_coloc_dataset(small)[names(small)],
                     do.call(coloc_dataset, small)[names(small)])
    expect_identical(as_coloc_dataset(D1)$beta, do.call(coloc_dataset, D1)$beta)
})

test_that("mistakes are reported when the dataset is made, with check_dataset()'s messages", {
    expect_error(coloc_dataset(snp = small$snp, beta = small$beta, varbeta = small$varbeta, sdY = 1), "missing required element\\(s\\) type")
    expect_error(coloc_dataset(type = "quant", beta = small$beta, varbeta = small$varbeta, sdY = 1), "missing required element\\(s\\) snp")
    expect_error(coloc_dataset(snp = small$snp, type = "quant", beta = small$beta, sdY = 1), "require p values and MAF")
    bad <- small; bad$beta[2] <- NA
    expect_error(as_coloc_dataset(bad), "beta contains missing values")
    bad <- small; bad$type <- "binary"
    expect_error(as_coloc_dataset(bad), "type must be quant or cc")
})

test_that("as_coloc_dataset() refuses things that are not datasets", {
    expect_error(as_coloc_dataset(1:5), "expected a coloc dataset")
    expect_error(as_coloc_dataset(matrix(0, 2, 2)), "expected a coloc dataset")
    expect_error(as_coloc_dataset(as.data.frame(small[c("snp", "beta", "varbeta")])), "expected a coloc dataset")
    expect_error(as_coloc_dataset(structure(list(), class = "susie")), "expected a coloc dataset")
})

test_that("as_coloc_dataset() keeps elements coloc does not use, and passes a coloc_dataset through", {
    x <- c(small, list(chr = "1"))
    D <- as_coloc_dataset(x)
    expect_identical(D$chr, "1")
    expect_output(print(D), "other elements: chr")
    expect_identical(as_coloc_dataset(D), D)
})

test_that("the p value warning is given once, when the dataset is made", {
    weak <- small; weak$beta <- rep(0.01, 5)
    expect_warning(D <- as_coloc_dataset(weak), "minimum p value")
    expect_warning(check_dataset(D), NA)
    ## a plain list still warns every time, as before
    expect_warning(check_dataset(weak), "minimum p value")
})

test_that("a dataset changed after it was made is checked again when used", {
    D <- as_coloc_dataset(D1)
    D$beta[1] <- NA
    expect_s3_class(D, "coloc_dataset") # $<- keeps the class, so the class alone cannot be trusted
    expect_error(coloc(D, D2), "beta contains missing values")
})

test_that("a coloc_dataset gives the same results as the plain list", {
    a1 <- as_coloc_dataset(D1)
    a2 <- as_coloc_dataset(D2)
    expect_equal(pp(coloc(a1, a2)), pp(coloc(D1, D2)))
    expect_equal(coloc.abf(a1, a2)$summary, coloc.abf(D1, D2)$summary)
})

test_that("a coloc_dataset still behaves as a list for generics with no method for it", {
    D <- do.call(coloc_dataset, small)
    expect_s3_class(data.table::as.data.table(D), "data.table") # plot_extended_data() does this
})

test_that("print() summarises the dataset", {
    D <- do.call(coloc_dataset, small)
    expect_output(print(D), "quantitative trait, 5 snps")
    expect_output(print(D), "sdY = 1")
    cc <- small; cc$type <- "cc"; cc$s <- 0.4; cc$sdY <- NULL
    expect_output(print(as_coloc_dataset(cc)), "case-control trait")
})

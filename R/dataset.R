## A coloc dataset constructor

## The elements a coloc dataset may have. check_dataset() checks these and
## ignores anything else.
.dataset_items <- c("beta", "varbeta", "pvalues", "MAF", "snp", "position", "N", "type", "s", "sdY", "LD")

##' Make a coloc dataset
##'
##' `coloc_dataset()` builds a dataset from its parts and checks it with
##' [check_dataset()], so a mistake is reported where the dataset is made
##' rather than later, inside `coloc()`. The result has class
##' `"coloc_dataset"`, but is still a list with the same elements, so code that
##' reads `D$beta` works unchanged.
##'
##' `as_coloc_dataset()` does the same for a dataset that is already a list, as
##' made for earlier versions of coloc. Elements coloc does not use are kept.
##'
##' Which elements are needed depends on the data, and [check_dataset()] gives
##' the details. In brief: always give `snp` and `type`; give `beta` and
##' `varbeta` if you have them, or `pvalues`, `MAF` and `N` if not; give `s`
##' for a case-control trait, and `sdY` (or `MAF` and `N`, to estimate it) for
##' a quantitative one.
##'
##' A dataset is checked when it is made. If you change its elements
##' afterwards, coloc's functions check it again when they use it, but don't
##' do that.
##' 
##' @param snp character vector of snp ids
##' @param type `"quant"` for a quantitative trait, `"cc"` for
##'     case-control
##' @param beta,varbeta per-snp effect estimates, and their variances
##'     (std error squared)
##' @param pvalues per-snp p values
##' @param MAF per-snp minor allele frequencies
##' @param position per-snp positions, used for plotting
##' @param N number of samples
##' @param s for a case-control trait, the proportion of samples that
##'     are cases
##' @param sdY for a quantitative trait, the standard deviation of the
##'     trait
##' @param LD matrix of LD between snps, with snp ids as row and
##'     column names
##' @param warn.minp warn if no p value is smaller than this
##' @return a list of class `"coloc_dataset"`
##' @seealso [check_dataset()] for what each element must look like
##' @export
##' @examples
##' D <- coloc_dataset(snp = paste0("s", 1:5),
##'                    type = "quant",
##'                    beta = c(0.1, 0.5, 0.2, 0, 0.05),
##'                    varbeta = rep(0.01, 5),
##'                    sdY = 1)
##' D
##'
##' ## a dataset made as a list for an earlier version of coloc
##' data(coloc_test_data)
##' D1 <- as_coloc_dataset(coloc_test_data$D1)
##'
##' ## mistakes are reported when the dataset is made
##' try(coloc_dataset(snp = paste0("s", 1:5),
##'                   type = "quant",
##'                   beta = c(0.1, 0.5, 0.2, 0, 0.05)))
coloc_dataset <- function(snp = NULL, type = NULL, beta = NULL, varbeta = NULL,
                          pvalues = NULL, MAF = NULL,
                          position = NULL, N = NULL, s = NULL, sdY = NULL,
                          LD = NULL, warn.minp = 1e-6) {
    ## every argument defaults to NULL so that a missing snp or type gets
    ## check_dataset()'s message, not R's "argument is missing"
    d <- list(beta = beta, varbeta = varbeta, pvalues = pvalues, MAF = MAF,
              snp = snp, position = position,
              N = N, type = type, s = s, sdY = sdY, LD = LD)
    d <- d[!sapply(d, is.null)]
    .new_coloc_dataset(d, warn.minp = warn.minp)
}

##' @rdname coloc_dataset
##' @param x a list with the elements of a coloc dataset, or a coloc_dataset
##' @export
as_coloc_dataset <- function(x, warn.minp = 1e-6) {
    if (inherits(x, "coloc_dataset")) return(x)
    ## a data.frame and a susie fit are lists too, but neither is a dataset
    if (!is.list(x) || is.data.frame(x) || inherits(x, "susie"))
        stop("expected a coloc dataset: a list, or the result of coloc_dataset(), not a ", class(x)[1])
    .new_coloc_dataset(unclass(x), warn.minp = warn.minp)
}

## The one place the class is set, so nothing gets it without being checked.
.new_coloc_dataset <- function(d, warn.minp = 1e-6) {
    check_dataset(d, warn.minp = warn.minp)
    structure(d, class = c("coloc_dataset", "list"))
}

##' @export
print.coloc_dataset <- function(x, ...) {
    cat("coloc dataset:",
        if (identical(x$type, "cc")) "case-control trait" else "quantitative trait,",
        length(x$snp), "snps\n")
    cat("  per-snp:", paste(intersect(c("beta", "varbeta", "pvalues", "MAF", "position"), names(x)), collapse = ", "), "\n")
    sc <- unlist(x[intersect(c("N", "s", "sdY"), names(x))])
    if (length(sc))
        cat("  ", paste(names(sc), signif(sc, 4), sep = " = ", collapse = ", "),
            "\n", sep = "")
    if (!is.null(x$LD)) cat("  LD:", nrow(x$LD), "x", ncol(x$LD), "\n")
    other <- setdiff(names(x), .dataset_items)
    if (length(other)) cat("  other elements:", paste(other, collapse = ", "), "\n")
    invisible(x)
}

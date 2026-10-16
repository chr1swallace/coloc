q <- data.table(chr = 6,
                pos = c(90279406, 90280331, 90281412, 90276601, 90288603, 90286261),
                snp = c("6_90279406_G_C", "6_90280331_G_A", "6_90281412_A_T", "6_90276601_A_G", "6_90288603_A_G", "6_90286261_C_G"),
                PP = c(0.675, 0.0747, 0.0702, 0.0505, 0.0422, 0.0412))
ql <- list(data.frame(snp="6_90383144_A_T",PP=0.5))
attr(ql, "annot") <- data.frame(id="blah")

test_that("coloc_query", {
    out <- coloc_query(q, ql)
    expect_equal(nrow(out),length(ql))
    expect_equal(out$nsnps1, nrow(q))
    expect_equal(out$nsnps2, nrow(ql[[1]]))
    expect_true(all(colnames(attr(ql,"annot")) %in% colnames(out)))
    expect_error(coloc_query(ql, q))
    attributes(ql) <- NULL
    expect_error(coloc_query(q, ql))
})

test_that("Open Targets file index is sorted and complete", {
    idx <- coloc:::.ot_index()
    expect_equal(nrow(idx), 200)
    expect_true(all(idx$size > 0))
    expect_true(all(is.na(idx$md5) | grepl("^[0-9a-f]{32}$", idx$md5)))
    chrs <- unique(c(rbind(idx$chr_first, idx$chr_last)))
    k1 <- match(idx$chr_first, chrs) * 1e10 + idx$pos_first
    k2 <- match(idx$chr_last, chrs) * 1e10 + idx$pos_last
    expect_true(all(k1 <= k2))
    expect_true(all(head(k2, -1) <= tail(k1, -1))) # files in order and not overlapping
    expect_setequal(chrs, c(1:22, "X", "Y"))
})

test_that("query regions select the right Open Targets files", {
    sel <- function(c, l, u) sub("^(part-[0-9]+)-.*", "\\1", coloc:::.ot_files_for(c, l, u)$file)
    expect_equal(sel(6, 90179406, 90388603), "part-00166")
    expect_equal(sel(1, 155051361, 155315184), "part-00012")
    expect_equal(sel(1, 155028000, 155029000), c("part-00011", "part-00012")) # straddles a file boundary
    expect_equal(sel(10, 1e6, 2e6), "part-00018") # file shared by end of chr1 and start of chr10
    expect_equal(nrow(coloc:::.ot_files_for(23, 1, 1e9)), 0)
})

test_that("query_ot on an uncovered chromosome returns an empty list without downloading", {
    expect_warning(ot <- query_ot(data.frame(chr = 23, pos = 1e6)), "No Open Targets")
    expect_length(ot, 0)
})

test_that("large downloads require a persistent cache dir", {
    old <- options(coloc.cache_dir = NULL, coloc.max_temp_download = NULL)
    on.exit(options(old), add = TRUE)
    expect_error(suppressWarnings(ot_fetch()), "persistent cache directory")
    expect_error(suppressWarnings(eqtl_cat_fetch(c("ge", "majiq"))), "persistent cache directory")
    expect_error(coloc:::.check_download_size(2e9, NULL, "x"), "2 GB")
    expect_null(coloc:::.check_download_size(15e6, NULL, "x"))     # small: temp dir is fine
    expect_null(coloc:::.check_download_size(2e9, "/some/dir", "x")) # dir supplied
    options(coloc.max_temp_download = 5e9)
    expect_null(coloc:::.check_download_size(2e9, NULL, "x"))       # limit raised
})

test_that("downloads are checked, and bad ones removed", {
    src <- tempfile(fileext = ".bin"); writeBin(as.raw(1:100), src)
    url <- paste0("file://", normalizePath(src)); good <- unname(tools::md5sum(src))
    d <- tempfile("dl")
    expect_error(suppressMessages(coloc:::.download_once(url, d, "a.bin", md5 = "0000")), "Checksum mismatch")
    expect_false(file.exists(file.path(d, "a.bin")))
    expect_length(list.files(d, pattern = "\\.part"), 0)
    expect_error(suppressMessages(coloc:::.download_once(url, d, "a.bin", size = 99)), "Size mismatch")
    p <- suppressMessages(coloc:::.download_once(url, d, "a.bin", md5 = good, size = 100))
    expect_equal(file.size(p), 100)
    writeBin(as.raw(1:5), p) # truncated file is replaced when size is known
    p <- suppressMessages(coloc:::.download_once(url, d, "a.bin", size = 100))
    expect_equal(file.size(p), 100)
})

test_that("eQTL Catalogue file indexes are complete and sorted", {
    for (v in coloc:::.eqtl_cat_releases$version) {
        idx <- coloc:::.eqtl_cat_index(v)
        expect_true(all(c("quant_method", "file", "chromosome", "pos_min", "pos_max", "size", "md5") %in% names(idx)), info = v)
        expect_true(all(grepl(paste0("^eQTL_Catalogue_", v, "_"), idx$file)), info = v)
        expect_true(all(grepl("^[0-9a-f]{32}$", idx$md5)), info = v)
        expect_lte(length(unique(idx$file)) + 2, 100) # Zenodo files per record
        ## within each method and chromosome, chunks are in position order
        expect_true(all(tapply(idx$pos_min, paste(idx$quant_method, idx$chromosome), function(x) !is.unsorted(x))), info = v)
    }
    expect_false("majiq" %in% coloc:::.eqtl_cat_index("r7")$quant_method)
    expect_true("majiq" %in% coloc:::.eqtl_cat_index("r8-beta")$quant_method)
})

test_that("eQTL Catalogue releases are cached separately and validated", {
    d <- tempfile("ecat")
    expect_false(eqtl_cat_dir(d, "r7") == eqtl_cat_dir(d, "r8-beta"))
    expect_equal(dirname(eqtl_cat_dir(d, "r7")), dirname(eqtl_cat_dir(d, "r8-beta")))
    expect_error(query_eqtl_cat(q, quant_method = "majiq", version = "r7", dir = d), "not available in eQTL Catalogue r7")
    expect_error(query_eqtl_cat(q, version = "r9", dir = d), "version must be one of")
    expect_warning(ot <- query_eqtl_cat(data.frame(chr = 6, pos = 1e9), dir = d), "No eQTL Catalogue")
    expect_length(ot, 0)
})

test_that("query_ot", {
    skip_on_cran()
    if(is.null(getOption("coloc.cache_dir"))) {
        expect_warning(ot <- query_ot(q))
    } else { 
        (ot <- query_ot(q))
    }
    out <- coloc_query(q, ot)
    expect_equal(length(ot), 484)
    expect_equal(nrow(out), 484)
    expect_equal(sum(out$PP.H4 > 0.5), 44)
})

test_that("query_ot remote = TRUE matches local", {
    skip_on_cran()
    ot_local <- suppressWarnings(query_ot(q))
    ot_remote <- suppressWarnings(query_ot(q, remote = TRUE))
    expect_equal(names(ot_remote), names(ot_local))
    expect_equal(attr(ot_remote, "annot"), attr(ot_local, "annot"))
})

test_that("query_eqtl_cat r7 (slow)", {
    skip_on_cran()
    ecat <- suppressWarnings(query_eqtl_cat(q, quant_method = "exon", version = "r7"))
    expect_true(is.list(ecat))
    expect_true(all(c("dataset_id", "cs_id") %in% names(attr(ecat, "annot"))))
    expect_false("beta" %in% names(attr(ecat, "annot")))            # variant-level columns stay in the credible sets
    expect_equal(attr(ecat, "annot")$cs_id, names(ecat))           # one annotation row per credible set, in order
})

test_that("query_eqtl_cat (slow)", {
    skip_on_cran()
    if(is.null(getOption("coloc.cache_dir"))) {
        expect_warning(ecat <- query_eqtl_cat(q, quant_method="exon"))
    } else {
        (ecat <- query_eqtl_cat(q, quant_method="exon"))
    }
    out <- coloc_query(q, ecat)
    expect_equal(length(ecat), 1)
    expect_equal(length(ecat), nrow(out))
    expect_equal(out$PP.H4, 0)
})

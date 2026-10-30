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

test_that("chr may be numeric or character, so query_ot() credible sets can be used as queries", {
    for (chr in list(6, 6L, "6", factor("6"))) expect_equal(coloc:::.check_q(data.frame(chr = chr, pos = 1)), "6")
    expect_equal(coloc:::.check_q(data.frame(chr = "chr6", pos = 1)), "6") # tolerated, though not documented
    expect_equal(coloc:::.check_q(data.frame(chr = "X", pos = 1)), "X")
    ot_style <- data.frame(chr = c("1", "1"), pos = c(155051361L, 155215184L), snp = c("1_155051361_A_G", "1_155215184_C_T"), PP = c(0.6, 0.3))
    expect_equal(coloc:::.check_q(ot_style), "1")
    expect_error(coloc:::.check_q(data.frame(chr = c("6", "7"), pos = 1:2)), "only have one value")
    expect_error(coloc:::.check_q(data.frame(chr = TRUE, pos = 1)), "chr must be numeric or character")
    ## character chr reaches file selection: a region with no data warns rather than erroring, without downloading
    expect_warning(query_ecat(data.frame(chr = "6", pos = 1e9)), "No eQTL Catalogue")
    expect_warning(query_ot(data.frame(chr = "23", pos = 1e6)), "No Open Targets")
})

test_that("Open Targets variants on alt/patch contigs are kept, with NA position", {
    locus <- data.frame(variantId = c("17_45000000_A_G", "17_KI270857v1_alt_811563_C_G", "17_45001000_T_C"),
                        posteriorProbability = c(0.5, 0.3, 0.15), beta = 0.1)
    expect_no_warning(cs <- coloc:::.otlocus_to_cs(locus))
    expect_equal(cs$pos, c(45000000L, NA, 45001000L))
    expect_equal(cs$chr, rep("17", 3))
    expect_equal(sum(cs$PP), 0.95) # alt contig variant still counts towards the credible set
    expect_equal(coloc:::.check_q(cs), "17") # usable as a query
    expect_error(coloc:::.check_q(data.frame(chr = 17, pos = NA_integer_)), "at least one non-missing pos")
    expect_silent(coloc:::.check_snp_ids(cs$snp)) # accepted as index credible set ids in coloc_query()
    expect_error(coloc:::.check_snp_ids(c("rs123", "chr6_90279406_G_C")), "not in the format")
})

test_that("*_many functions skip failures, collect warnings, and stop on systematic failure", {
    fq <- function(cs, nm) { warning("same warning"); if (identical(cs$chr[1], "bad")) stop("bad credible set"); nm }
    cs3 <- list(a = data.frame(chr = "1"), b = data.frame(chr = "bad"), c = data.frame(chr = "2"))
    ws <- character()
    res <- withCallingHandlers(coloc:::.run_many(cs3, fq, "test", 5), warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
    expect_null(res$b)
    expect_equal(res$c, "c")
    expect_equal(attr(res, "failed"), c(b = "bad credible set"))
    expect_length(ws, 2) # the repeated warning once, plus the failure summary
    expect_match(ws[1], "credible set\\(s\\): a, b, c")
    bad <- setNames(rep(list(data.frame(chr = "bad")), 6), letters[1:6])
    expect_error(suppressWarnings(coloc:::.run_many(bad, fq, "test", 3)), "3 credible sets in a row failed")
    expect_error(suppressWarnings(coloc:::.run_many(bad[1:2], fq, "test", 5)), "all 2 credible sets failed")
    expect_error(query_ot_many(q), "list of credible sets")
})

test_that("coloc_query_many pairs credible sets with their query results by name", {
    cs_list <- list(first = q, second = q, third = q)
    query_list <- list(first = ql, second = NULL, third = ql) # second's query failed
    attr(query_list, "failed") <- c(second = "query failed")
    out <- coloc_query_many(cs_list, query_list)
    expect_equal(out$query_cs, c("first", "third"))
    expect_equal(out[, -1], rbind(coloc_query(q, ql), coloc_query(q, ql)), check.attributes = FALSE)
    expect_equal(names(attr(out, "failed")), "second")
    expect_error(coloc_query_many(cs_list, list(other = ql)), "names not found in cs_list: other")
})

test_that("large downloads require a persistent cache dir", {
    old <- options(coloc.cache_dir = NULL, coloc.max_temp_download = NULL)
    on.exit(options(old), add = TRUE)
    expect_error(suppressWarnings(ot_fetch()), "persistent cache directory")
    expect_error(suppressWarnings(ecat_fetch(c("ge", "majiq"), version = "r8-beta")), "persistent cache directory")
    expect_error(coloc:::.check_download_size(2e9, "x"), "2 GB")
    expect_null(coloc:::.check_download_size(15e6, "x"))          # small: temp dir is fine
    options(coloc.max_temp_download = 5e9)
    expect_null(coloc:::.check_download_size(2e9, "x"))           # limit raised
    options(coloc.cache_dir = tempfile(), coloc.max_temp_download = NULL)
    expect_null(coloc:::.check_download_size(2e9, "x"))           # cache dir set
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
    for (v in coloc:::.ecat_releases$version) {
        idx <- coloc:::.ecat_index(v)
        expect_true(all(c("quant_method", "file", "chromosome", "pos_min", "pos_max", "size", "md5") %in% names(idx)), info = v)
        expect_true(all(grepl(paste0("^eQTL_Catalogue_", v, "_"), idx$file)), info = v)
        expect_true(all(grepl("^[0-9a-f]{32}$", idx$md5)), info = v)
        expect_lte(length(unique(idx$file)) + 2, 100) # Zenodo files per record
        ## within each method and chromosome, chunks are in position order
        expect_true(all(tapply(idx$pos_min, paste(idx$quant_method, idx$chromosome), function(x) !is.unsorted(x))), info = v)
    }
    expect_false("majiq" %in% coloc:::.ecat_index("r7")$quant_method)
    expect_true("majiq" %in% coloc:::.ecat_index("r8-beta")$quant_method)
})

test_that("eQTL Catalogue releases are cached separately and validated", {
    old <- options(coloc.cache_dir = tempfile("ecat"))
    on.exit(options(old), add = TRUE)
    expect_false(ecat_dir("r7") == ecat_dir("r8-beta"))
    expect_equal(dirname(ecat_dir("r7")), dirname(ecat_dir("r8-beta")))
    expect_error(query_ecat(q, quant_method = "majiq", version = "r7"), "not available in eQTL Catalogue r7")
    expect_error(query_ecat(q, version = "r9"), "version must be one of")
    expect_warning(ot <- query_ecat(data.frame(chr = 6, pos = 1e9)), "No eQTL Catalogue")
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
    ## result <- annotate_ot(out)
    ## plot_query(result, label_n=1)
})

test_that("query_ecat r7", {
    skip_on_cran()
    ecat <- suppressWarnings(query_ecat(q, quant_method = "exon", version = "r7"))
    expect_true(is.list(ecat))
    expect_true(all(c("dataset_id", "cs_id") %in% names(attr(ecat, "annot"))))
    expect_false("beta" %in% names(attr(ecat, "annot")))            # variant-level columns stay in the credible sets
    expect_equal(attr(ecat, "annot")$cs_id, names(ecat))           # one annotation row per credible set, in order
    ## result <- coloc_query(q, ecat)
    ## result <- annotate_ecat(result)
    ## plot_query(result)
})

test_that("query_ecat r8-beta", {
    skip_on_cran()
    if(is.null(getOption("coloc.cache_dir"))) {
        expect_warning(ecat <- query_ecat(q, quant_method="exon", version="r8-beta"))
    } else {
        (ecat <- query_ecat(q, quant_method="exon", version="r8-beta"))
    }
    out <- coloc_query(q, ecat)
    expect_equal(length(ecat), 1)
    expect_equal(length(ecat), nrow(out))
    expect_equal(out$PP.H4, 0)
    ## result <- annotate_ecat(out)
    ## plot_query(result)
})

## *_many versions of the three tests above: the same credible set twice, so both sets should give the single-set results
qs <- list(set1 = q, set2 = q)

test_that("query_ot_many", {
    skip_on_cran()
    if(is.null(getOption("coloc.cache_dir"))) {
        expect_warning(ot <- query_ot_many(qs))
    } else {
        (ot <- query_ot_many(qs))
    }
    expect_named(ot, names(qs))
    expect_equal(lengths(ot), c(set1 = 484, set2 = 484))
    out <- coloc_query_many(qs, ot)
    expect_equal(nrow(out), 2 * 484)
    expect_equal(as.vector(table(out$query_cs)), c(484, 484))
    expect_equal(sum(out$PP.H4 > 0.5), 2 * 44)
    expect_equal(out[out$query_cs == "set2", -1], out[out$query_cs == "set1", -1], check.attributes = FALSE)
    result <- annotate_ot(out)
    plot_query(result, label_n=1)
})

test_that("query_ecat_many r7", {
    skip_on_cran()
    ecat <- suppressWarnings(query_ecat_many(qs, quant_method = "exon", version = "r7"))
    expect_named(ecat, names(qs))
    for (e in ecat) {
        expect_true(is.list(e))
        expect_true(all(c("dataset_id", "cs_id") %in% names(attr(e, "annot"))))
        expect_false("beta" %in% names(attr(e, "annot")))            # variant-level columns stay in the credible sets
        expect_equal(attr(e, "annot")$cs_id, names(e))               # one annotation row per credible set, in order
    }
    expect_equal(ecat$set2, ecat$set1)
    out <- query_ecat_many(qs)
    result <- annotate_ecat(coloc_query_many(qs, ecat))
    plot_query(result)
})

test_that("query_ecat_many r8-beta", {
    skip_on_cran()
    if(is.null(getOption("coloc.cache_dir"))) {
        expect_warning(ecat <- query_ecat_many(qs, quant_method="exon", version="r8-beta"))
    } else {
        (ecat <- query_ecat_many(qs, quant_method="exon", version="r8-beta"))
    }
    expect_equal(lengths(ecat), c(set1 = 1, set2 = 1))
    out <- coloc_query_many(qs, ecat)
    expect_equal(nrow(out), 2)
    expect_equal(out$query_cs, names(qs))
    expect_equal(out$PP.H4, c(0, 0))
})

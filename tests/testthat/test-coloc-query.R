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

test_that("query_ot (slow)", {
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

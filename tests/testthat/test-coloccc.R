## ## delete this
## devtools::load_all("~/RP/coloc")

## start here
(load(system.file("testdata","test-coloccc.RData",package="coloc")))
ret.abf <- coloc.abf(dataset1,dataset2,MAF=MAF)
ret <- coloc.cc(dataset1,dataset2,MAF=MAF,n1=n1,n2=n2,n00=n00,n01=0,n02=0,LD=LD)

test_that("return values have correct items", {
    expect_equal(names(ret),c("summary","results","df3","sdf"))
    expect_equal(length(ret$summary),length(ret.abf$summary))
    expect_identical(names(ret$summary),names(ret.abf$summary))
})

test_that("return values are correct (ish)", {
    expect_true(which.max(ret$summary[-1])==which.max(ret.abf$summary[-1]))
    expect_equal(ret$summary[1],ret.abf$summary[1])
})

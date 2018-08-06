## delete this
devtools::load_all("~/RP/coloc")

## start here
(load(system.file("testdata","test-coloccc.RData",package="coloc")))
ret.abf <- coloc.abf(dataset1,dataset2,MAF=MAF)$summary
ret.multi <- coloc.cc(dataset1,dataset2,MAF=MAF,n1=n1,n2=n2,n00=n00,LD=LD,method="multinom")
ret.corr <- coloc.cc(dataset1,dataset2,MAF=MAF,n1=n1,n2=n2,n00=n00,n01=0,n02=0,LD=LD,method="corr")


test_that("return values have correct items", {
    expect_equal(length(ret.multi),length(ret.abf))
    expect_equal(length(ret.corr),length(ret.abf))
    expect_identical(names(ret.multi),names(ret.abf))
    expect_identical(names(ret.corr),names(ret.abf))
})

test_that("return values are correct (ish)", {
    expect_true(which.max(ret.multi[-1])==which.max(ret.abf[-1]))
    expect_true(which.max(ret.corr[-1])==which.max(ret.abf[-1]))
    expect_equal(ret.multi[1],ret.abf[1])
    expect_equal(ret.corr[1],ret.abf[1])
})

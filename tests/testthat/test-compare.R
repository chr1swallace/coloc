library(coloc)
data(coloc_test_data)
attach(coloc_test_data)

result1=compare_datasets(D1,D2)
D0=subset_dataset(D1, c(1:104,106:500))
result2=suppressWarnings(compare_datasets(D0,D2))

test_that("compare_datasets", {
    expect_equal(result1, c(oneintwo=1,twoinone=1))
    expect_equal(result2, c(oneintwo=1,twoinone=0.4360232))
    expect_warning(compare_datasets(D0,D2))
})

detach(coloc_test_data)


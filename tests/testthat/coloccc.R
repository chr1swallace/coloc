## library(snpStats)
## data(testdata)
## X <- Autosomes[,101:201]
## cs <- col.summary(X)
## X <- X[, cs[,"MAF"]>0.05 & cs[,"Call.rate"]>0.9]
## maf <- col.summary(X)[,"MAF"]
## set.seed(42)

## ## case-control trait
## cc1 <- rbinom(nrow(X),1, p=(1+as(X[,8],"numeric"))/4)
## eff<-snp.rhs.estimates(cc1 ~ 1, snp.data=X, family="binomial")
## beta.cc1 <- sapply(eff@.Data, "[[", "beta")
## vbeta.cc1 <- sapply(eff@.Data, "[[", "Var.beta")
## p.cc1 <- pchisq(beta.cc1^2/vbeta.cc1,df=1,lower.tail=FALSE)

## cc2 <- rbinom(nrow(X),1, p=(1+as(X[,8],"numeric"))/4)
## eff<-snp.rhs.estimates(cc2 ~ 1, snp.data=X, family="binomial")
## beta.cc2 <- sapply(eff@.Data, "[[", "beta")
## vbeta.cc2 <- sapply(eff@.Data, "[[", "Var.beta")
## p.cc2 <- pchisq(beta.cc2^2/vbeta.cc2,df=1,lower.tail=FALSE)

## ## general things
## D1 <- list(beta=beta.cc1,
##             varbeta=vbeta.cc1,
##             snp=colnames(X),
##             MAF=maf)
## D2 <- list(beta=beta.cc2,
##             varbeta=vbeta.cc2,
##             snp=colnames(X),
##             MAF=maf)

## RESULTS <- list(naive = coloc.abf(dataset1=c(D1,list(N=nrow(X),s=sum(cc1)/nrow(X),type="cc")),
##                                   dataset2=c(D2,list(N=nrow(X),s=sum(cc2)/nrow(X),type="cc"))),
##                 adj=coloc.shared(dataset1=D1, dataset2=D2,
##                                  n00=nrow(X)-sum(cc1),n01=0,n02=0,n1=sum(cc1),n2=sum(cc2)))
## lapply(RESULTS,"[[","summary")

## ## coloc.abf with coefficients
## test_that("coloc.abf", {
##     expect_error(coloc.abf(dataset1=DQ,dataset2=DCC), NA)
##     result <- coloc.abf(dataset1=DQ,dataset2=DCC)
##     expect_true(which.max(result$summary[-1]) == 5)
##     expect_true(result$summary[1] == ncol(X))
## })


## ## alternative test data
## ## colocdata<- read.table("inst/tests/test.txt", sep="\t", header=T)
## ## N <- 18124
## ## result <- coloc.abf(dataset1=list(beta=colocdata$beta.dataset1,
## ##             varbeta=colocdata$varbeta.dataset1,
## ##             type="quant",
## ##             snp=colocdata$SNP,
## ##             N=N),
## ##           dataset2=list(beta=colocdata$beta.dataset1,
## ##             varbeta=colocdata$varbeta.dataset1,
## ##             type="quant",
## ##             snp=colocdata$SNP,
## ##             N=N),
## ##           MAF=colocdata$MAF)



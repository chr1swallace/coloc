## functions for running coloc with two quant traits measured on the same samples


##' Bayesian colocalisation analysis for case control studies with (partially) shared controls
##'
##' This function calculates posterior probabilities of different
##' causal variant configurations under the assumption of a single
##' causal variant for each trait.
##'
##' Unlike coloc.abf, this specifically requires coefficients and
##' their standard errors in the case of partially shared controls.
##'
##' @title coloc.abf for two case control studies and (partially) shared controls
##' @param dataset1 a list with the following elements
##' \describe{
##' 
##' \item{beta}{regression coefficient for each SNP from dataset 1}
##' 
##' \item{varbeta}{variance of beta}
##' 
##' \item{snp}{a character vector of snp ids. It will be used to merge dataset1 and dataset2.}
##'
##' }
##'
##' @param dataset2 as above, for dataset 2
##' @param n1 number of cases for dataset 1
##' @param n2 number of cases for dataset 2
##' @param n00 number of shared controls
##' @param n01 number of controls unique to dataset 1
##' @param n02 number of controls unique to dataset 2
##' @param MAF minor allele frequency of the variants, named vector
##' @param LD matrix giving correlation between SNPs.  NB - this is r, not rsquared!  Expect a mix of positive and negative values.  This should be estimated from reference data.
##' @param p1 prior probability a SNP is associated with trait 1,
##'     default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2,
##'     default 1e-4
##' @param p12 prior probability a SNP is associated with both traits,
##'     default 1e-5
##' @param method "multonom" or "corr", depending on which
##'     approximation (multinomial or correlation) you wish to use.
##' @return a vector giving the number of SNPs analysed, and the
##'     posterior probabilities of H0 (no causal variant), H1 (causal
##'     variant for trait 1 only), H2 (causal variant for trait 2
##'     only), H3 (two distinct causal variants) and H4 (one common
##'     causal variant)
##' @export
##' @author Chris Wallace
coloc.qq <- function(dataset1,dataset2,
                     n=dataset1$N, VY1, VY2, covY,
                     MAF, LD, 
                     p1=1e-4, p2=1e-4, p12=1e-5) {
    snps <- intersect(dataset1$snp,dataset2$snp)
    df <- data.table(snp=c(snps))
    df <- merge(df,as.data.table(as.data.frame(dataset1[c("snp","beta","varbeta")])),
                by.x="snp",by.y="snp",all.x=TRUE)
    setnames(df,c("beta","varbeta"),c("b1","v1"))
    df <- merge(df,as.data.table(as.data.frame(dataset2[c("beta","varbeta","snp")])),
                by.x="snp",by.y="snp",all.x=TRUE) 
    setnames(df,c("beta","varbeta"),c("b2","v2"))
    fq <- data.table(snp=snps)
    fq[,f0:=MAF[snp]]
    df <- merge(df,fq[,.(snp,f0)],by.x="snp",by.y="snp",all.x=TRUE)
    
    ## H1, H2, H4 all two beta models
    df[,erho:=(covY - 2*f0*(1-f0)*b1*b2)/(2*n*f0*(1+f0))]
    ## lbf for h1,h2,h4 
    df[, lbf1:=bf2(b1,b2,v1,v2,erho,w2=0)]
    df[, lbf2:=bf2(b1,b2,v1,v2,erho,w1=0)]
    df[, lbf4:=bf2(b1,b2,v1,v2,erho)]
    
    ## for h3, we need pairs of snps
    df3 <- as.data.table(as.data.frame(expand.grid(snps,snps)))
    setnames(df3,c("snpA","snpB"))
    df3 <- df3[snpA!=snpB,]
    df3[,isnpA:=match(as.character(snpA),rownames(LD))]
    df3[,isnpB:=match(as.character(snpB),rownames(LD))]
    df3[,r:=LD[cbind(isnpA,isnpB)]]
    df3 <- merge(df3,df[,.(snp,f0,b1,v1)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 1, snp A
    setnames(df3,c("f0","b1","v1"),c("fA","bA1","vA1"))
    df3 <- merge(df3,df[,.(snp,f0,b1,v1)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 1, snp B
    setnames(df3,c("f0","b1","v1"),c("fB","bB1","vB1"))
    df3 <- merge(df3, df[,.(snp,b2,v2)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 2, snp A
    setnames(df3,c("b2","v2"),c("bA2","vA2"))
    df3 <- merge(df3, df[,.(snp,b2,v2)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 2, snp B
    setnames(df3,c("b2","v2"),c("bB2","vB2"))

    ## special case: r==1. This won't allow a model to be fitted at all.  Instead, add in the bf from h4 for one of the two SNPs
    df3.r1 <- df3[abs(r)>=0.99,]
    m <- match(df3.r1$snpA, df$snp)
    df3.r1[, lbf3:=df$lbf4[m]]

    ## todo special case: r=0. Can fit two bivariate normals, which means we can vectorize
    
    df3 <- df3[abs(r)<0.99,]
    ## reconstruct coefficients from bivariate models
    df3[,c("bA1old","bB1old","bA2old","bB2old"):=list(bA1,bB1,bA2,bB2)]
    df3[,c("bA1","bB1"):=bstar_from_b(bA1,bB1,vA1,vB1,r)]
    df3[,c("bA2","bB2"):=bstar_from_b(bA2,bB2,vA2,vB2,r)]
    df3[,VA:=2*fA*(1-fA)]
    df3[,VB:=2*fB*(1-fB)]
    df3[,VAB := r * sqrt(VA) * sqrt(VB)]
    df3[,sumA2:=2*n*fA*(1+fA)] #sum(bvec^2 * rowSums(geno))
    df3[,sumB2:=2*n*fB*(1+fB)] #sum(bvec^2 * colSums(geno))
    df3[,sumAB:=n*( sqrt(VA)*sqrt(VB)*r+4*fA*fB)] # sum( (bvec %*% t(bvec)) * geno )
    df3[,det:=sumA2*sumB2-sumAB^2]
    df3[,EE:=covY - bA1*bA2*VA - bB1*bB2*VB - (bA1*bB2 + bB1*bA2)*VAB]
    df3[,sumA2:=EE*sumA2/det]
    df3[,sumB2:=EE*sumB2/det]
    df3[,sumAB:=EE*sumAB/det]
    df3[,W1:=VY1*0.15^2]
    df3[,W2:=VY2*0.15^2]

    ## wh <- head(which(df3$r < 0.99 & df3$snpA < df3$snpB),50)
    ## microbenchmark({for(i in wh) {
    ##                     tmp[[i]] <- with(df3[i,],
    ##                                 lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
    ##                                            vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
    ##                                            W1=W1,W2=W2,sumA2,sumAB,sumB2,r))})

    for(i in 1:nrow(df3)) { #& df3$snpA < df3$snpB)) {
        ## microbenchmark({tmp <- with(df3[i,],
        ##                             lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
        ##                                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
        ##                                        W1=W1,W2=W2,sumA2,sumAB,sumB2,r))})
        df3[i,lbf3:=lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
                                               vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
                                               W1=W1,W2=W2,sumA2,sumAB,sumB2,r)]
        ## microbenchmark({df3[i,lbf3:=lbf.h3.qq2(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
        ##                                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
        ##                                        W1=W1,W2=W2,sumA2,sumAB,sumB2,r)]})
        ## microbenchmark({df3[i,lbf3:= lbf.h3.qq2(bA1=df3$bA1[i],bB1=df3$bB1[i],bA2=df3$bA2[i],bB2=df3$bB2[i],
        ##                                        vA1=df3$vA1[i],vB1=df3$vB1[i],vA2=df3$vA2[i],vB2=df3$vB2[i],
        ##                                        W1=df3$W1[i],W2=df3$W2[i],df3$sumA2[i],df3$sumAB[i],df3$sumB2[i],df3$r[i])]})
        ## j <- which(df3$snpA==df3$snpB[i] & df3$snpB==df3$snpA[i])
        ## microbenchmark({tmpj <- with(df3[j,],
        ##              lbf.h3.qq(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
        ##                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
        ##                        n,
        ##                        VY1,VY2,covY,fA,fB,r,W=0.04))})
        ## df3[i,lbf3:=tmp$asis]
        ## df3[j,lbf3:=tmp$reversed]
    } 
    ## put it all together
    o<-intersect(names(df3),names(df3.r1))
    df3 <- rbind(df3.r1[,o,with=FALSE],df3[,o,with=FALSE])
    sdf <- df[, lapply(.SD, logsum), .SDcols=grep("lbf",names(df),value=TRUE)]
    sdf3 <- df3[, lapply(.SD, logsum), .SDcols=grep("lbf",names(df3),value=TRUE)]
    sdf <- rbind(melt(sdf,measure.vars=names(sdf)),
                 melt(sdf3,measure.vars=names(sdf3)))
    ## add h0
    sdf0 <- melt(sdf3,measure.vars=names(sdf3))
    sdf0[,value:=0]
    sdf0[,variable:=sub("3","0",variable)]
    sdf <- rbind(sdf,sdf0)
    
    sdf[,h:=as.integer(sub(".*lbf","",variable))]
    sdf[,method:=sub("[0-4]","",variable)]
    sdf[,abf:=0]
    sdf[h==1,abf:=log(p1) + value]
    sdf[h==2,abf:=log(p2) + value]
    sdf[h==3,abf:=log(p1) + log(p2) + value]
    sdf[h==4,abf:=log(p12) + value]
    sdf[,pp:=exp(abf - logsum(abf)),by="method"]
    sdf <- sdf[order(method,h),]
    ## while checking, return this
    ## return(sdf)
    
    ret <- c(nsnps=length(snps), sdf[method=="lbf",]$pp)
    names(ret) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", 
                    "PP.H4.abf")
    ret
}
 




##' Calculate ABF under H3, paired quant traits
##'
##' Assumes SNP A is causal for disease 1, and SNP B is causal for disease 2
##' @title Approximate LBF under H3, using multinomial approximation
##' @param bA1 log OR for SNP A, disease 1 <--- causal for disease 1
##' @param bB1 log OR for SNP B, disease 1
##' @param bA2 log OR for SNP A, disease 2
##' @param bB2 log OR for SNP A, disease 2 <-- causal for disease 2
##' @param n00 number of shared controls 
##' @param n1 number of disease 1 cases
##' @param n2 number of disease 1 cases
##' @param fA MAF of A in controls
##' @param fB MAF of B in controls
##' @param rho correlation (r) between A and B in controls
##' @param W variance on normal prior for effect sizes at causal variants
##' @return lBF for A, B from the prior involving W to a prior assuming true effect size is 0
##' @author Chris Wallace
lbf.h3.qq2 <- function(bA1,bB1,bA2,bB2,vA1,vB1,vA2,vB2,W1,W2,sumA2,sumAB,sumB2,rho) {
    ## geno is matrix giving relative frequency of individuals with
    ## genotype i,j in population, where i is count of SNP 1 alt
    ## allele, and j is count of SNPB alt allele
    ## bA1=bA1star; bB1=bB1star; bA2=bA2star; bB2=bB2star
    ## geno <- estgeno.3(n,fA,fB,rho)
    ## bvec <- c(0,1,2)
    ## sumA2 <- 2*n*fA*(1+fA) #sum(bvec^2 * rowSums(geno))
    ## sumB2 <- 2*n*fB*(1+fB) #sum(bvec^2 * colSums(geno))
    ## ## sumAB <- sum( (bvec %*% t(bvec)) * geno )
    ## sumAB <- n*( sqrt(VA)*sqrt(VB)*rho+4*fA*fB)
    ## VA <- 2*fA*(1-fA)
    ## VB <- 2*fB*(1-fB)
    ## VAB <- rho * sqrt(VA) * sqrt(VB)
    ## EE <- covY - bA1*bA2*VA - bB1*bB2*VB - (bA1*bB2 + bB1*bA2)*VAB

    ## can vectorise above here TODO
    V12 <- matrix(c(sumB2,-sumAB,-sumAB,sumA2),2,2) ## already inverse and x EE
    ## V12 <- EE * solve(WtW) # cov( (bA1,bA2), (bB1,bB2) )
    ## V12 <- EE * WtW # cov( (bA1,bA2), (bB1,bB2) )
    V1 <- vstar_from_v(vA1,vB1,rho)
    V2 <- vstar_from_v(vA2,vB2,rho)
    V <- rbind(cbind(V1,V12), cbind(V12,V2))
    B <- c(bA1,bB1,bA2,bB2)
    Wmat <- diag(c(W1,0,0,W2))
    ## Wmat <- diag(c(0.04,0,0,0.04))
    dmvnorm(B, mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
      dmvnorm(B, mean=rep(0,4), sigma=V, log=TRUE)
    ## Wmat <- diag(c(0,W1,W2,0))
    ## ## Wmat <- diag(c(0,0.04,0.04,0))
    ## ret2 <- dmvnorm(B, mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
    ##   dmvnorm(B, mean=rep(0,4), sigma=V, log=TRUE)
    ## list(asis=ret1,reversed=ret2)
}

lbf.h3.qq <- function(bA1,bB1,bA2,bB2,vA1,vB1,vA2,vB2,n,vY1,vY2,covY,fA,fB,rho,W) {
    ## geno is matrix giving relative frequency of individuals with
    ## genotype i,j in population, where i is count of SNP 1 alt
    ## allele, and j is count of SNPB alt allele
    ## bA1=bA1star; bB1=bB1star; bA2=bA2star; bB2=bB2star
    geno <- estgeno.3(n,fA,fB,rho)
    bvec <- c(0,1,2)
    sumA2 <- 2*n*fA*(1+fA) #sum(bvec^2 * rowSums(geno))
    sumB2 <- 2*n*fB*(1+fB) #sum(bvec^2 * colSums(geno))
    ## sumAB <- sum( (bvec %*% t(bvec)) * geno )
    sumAB <- n*( sqrt(VA)*sqrt(VB)*rho+4*fA*fB)
    VA <- 2*fA*(1-fA)
    VB <- 2*fB*(1-fB)
    VAB <- rho * sqrt(VA) * sqrt(VB)
    EE <- covY - bA1*bA2*VA - bB1*bB2*VB - (bA1*bB2 + bB1*bA2)*VAB

   ## can vectorise above here TODO
    WtW <- matrix(c(sumA2,sumAB,sumAB,sumB2),2,2)
    V12 <- EE * solve(WtW) # cov( (bA1,bA2), (bB1,bB2) )
    V1 <- vstar_from_v(vA1,vB1,rho)
    V2 <- vstar_from_v(vA2,vB2,rho)
    V <- rbind(cbind(V1,V12), cbind(V12,V2))
    B <- c(bA1,bB1,bA2,bB2)
    Wmat <- diag(c(0.15^2*VY1,0,0,0.15^2*VY2))
    ## Wmat <- diag(c(0.04,0,0,0.04))
    ret1 <- dmvnorm(B, mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
      dmvnorm(B, mean=rep(0,4), sigma=V, log=TRUE)
    Wmat <- diag(c(0,0.15^2*VY1,0.15^2*VY2,0))
    ## Wmat <- diag(c(0,0.04,0.04,0))
    ret2 <- dmvnorm(B, mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
      dmvnorm(B, mean=rep(0,4), sigma=V, log=TRUE)
    list(asis=ret1,reversed=ret2)
}

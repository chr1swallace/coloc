bhat.from.bstar <- function(b1star, b2star, v1star, v2star, rho) {
    w1 <- 1/(1-rho^2)
    w2 <- rho/(sqrt(v1star*v2star) * (1-rho^2))
    b1 <- w1 * b1star - w2 * b2star * v1star
    b2 <- w1 * b2star - w2 * b1star * v2star
    v1 <- w1*v1star
    v2 <- w1*v2star
    list(b1,b2,v1,v2)
}
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title recreate bivariate logistic regression from univariate summary stats
##' @param b1 beta_1
##' @param b2 beta_2
##' @param v1 v(beta_1)
##' @param v2 v(beta_2)
##' @param rho cor(snp1,snp2)
##' @param w1 prior variance on b1 (beta_1 ~ N(b1,v); b1 ~ N(0,w1))
##' @param w2 prior variance on b2
##' @return  log ABF of data under alternative vs null (w1=w2=0)
##' @author Chris Wallace
bf2 <- function(b1,b2,v1,v2,rho,w1=0.04,w2=0.04) {
    s1 <- sqrt(v1)
    s2 <- sqrt(v2)
    zdiff <- (-b1 * (s2 ^ 2 + w2) / (rho ^ 2 * s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * w2 - s2 ^ 2 * w1 - w1 * w2) + b2 * s1 * s2 * rho / (rho ^ 2 * s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * w2 - s2 ^ 2 * w1 - w1 * w2)) * b1 + (b1 * s1 * s2 * rho / (rho ^ 2 * s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * w2 - s2 ^ 2 * w1 - w1 * w2) - b2 * (s1 ^ 2 + w1) / (rho ^ 2 * s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * s2 ^ 2 - s1 ^ 2 * w2 - s2 ^ 2 * w1 - w1 * w2)) * b2 - (-b1 / s1 ^ 2 / (rho ^ 2 - 1) + b2 / s1 / s2 * rho / (rho ^ 2 - 1)) * b1 - (b1 / s1 / s2 * rho / (rho ^ 2 - 1) - b2 / s2 ^ 2 / (rho ^ 2 - 1)) * b2
    DR <- (-rho ^ 2 * s1 ^ 2 * s2 ^ 2 + s1 ^ 2 * s2 ^ 2 + s1 ^ 2 * w2 + s2 ^ 2 * w1 + w1 * w2) / (-rho ^ 2 * s1 ^ 2 * s2 ^ 2 + s1 ^ 2 * s2 ^ 2)
    (-zdiff -log(DR))/2
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title recreate bivariate logistic regression from univariate summary stats
##' @param b1 beta_1
##' @param b2 beta_2
##' @param v1 v(beta_1)
##' @param v2 v(beta_2)
##' @param rho cor(snp1,snp2)
##' @return  list of b1, b2 in bivariate regression
##' @author Chris Wallace
bstar_from_b <- function(b1, b2, v1, v2, rho){
bstar1 <- b1/(1 - rho^2) - rho/(1 - rho^2) * sqrt(v1/v2) * b2
bstar2 <- b2/(1 - rho^2) - rho/(1 - rho^2) * sqrt(v2/v1) * b1
## vstar <- 1/(1 - rho^2) * matrix(c(v1, -rho * sqrt(v1 * v2), -rho * sqrt(v1 * v2), v2), nrow = 2, byrow = TRUE)
list(bstar1, bstar2)
}

##' Calculate ABF under H3, multinomial
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
lbf.h3.multinom <- function(bA1,bB1,bA2,bB2,n00,n1,n2,fA,fB,rho,W) {
    ## geno is matrix giving relative frequency of individuals with
    ## genotype i,j in population, where i is count of SNP 1 alt
    ## allele, and j is count of SNP2 alt allele
    geno00 <- estgeno.3(n00,fA,fB,rho)
    geno1 <- estgeno.3(n1,fA,fB,rho,ORA=exp(bA1),ORB=exp(bB1))
    geno2 <- t(estgeno.3(n2,fA=fB,fB=fA,rho,ORA=exp(bB2),ORB=exp(bA2))) # swap snps and t() because estgeno assumes OR for first SNP
    geno <- geno00+geno1+geno2
    ## estimate a1, a2 
    bvec <- c(0,1,2)
    b1m <- matrix(bvec*bA1,3,3) + matrix(bvec*bB1,3,3,byrow=TRUE)
    b2m <- matrix(bvec*bA2,3,3) + matrix(bvec*bB2,3,3,byrow=TRUE)
    G <- geno00+geno1
    fun <- function(a) {
        denom <- 1+exp(a + b1m)
        num1 <- exp(a + b1m)
        abs(sum(G*num1/denom) - n1)
    }
    ans <- optimize(fun,c(-100,100))
    a1 <- ans$minimum
    G <- geno00+geno2
    fun <- function(a) {
        denom <- 1+exp(a + b2m)
        num2 <- exp(a + b2m)
        abs(sum(G*num2/denom) - n2)
    }
    ans <- optimize(fun,c(-100,100))
    a2 <- ans$minimum
    ## use https://czep.net/stat/mlelr.pdf to fix next few lines, have bA1, bB1 etc
    ## eqn 37 i=indiv, j=trait, k=SNP
    num1 <- exp(a1 + b1m)
    num2 <- exp(a2 + b2m)
    denom <- 1 + num1 + num2
    phat1 <- - geno * num1/denom * (1-num1/denom) # have - here so don't need below
    phat2 <- - geno * num2/denom * (1-num2/denom)
    phat12 <- geno * num1/denom * num2/denom
    ka <- matrix(1,3,3)
    kA <- matrix(bvec,3,3)
    kB <- matrix(bvec,3,3,byrow=TRUE)
    H <- calc.hessian3(ka,kA,kB,phat1,phat2,phat12)
    V <- ginv(-H)[3:6,3:6] # var-cov matrix for (a1,a2,b1,b2)
    ## if(rho>0) {
    ##     Wmat <- diag(c(W,rho^2*W,rho^2*W,W))
    ##     Wmat[cbind(c(1,2,3,4),c(2,1,4,3))] <- rho^2*W
    ## } else {
        Wmat <- diag(c(W,0,0,W))
    ## }
    dmvnorm(c(bA1,bA2,bB1,bB2), mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
      dmvnorm(c(bA1,bA2,bB1,bB2), mean=rep(0,4), sigma=V, log=TRUE)
    ## return new v1, v2, rho
    ## c(v1,v2,rho)
    ## list(v.b1=diag(VV)[3],v.b2=diag(VV)[4],v.b1b2=VV[3,4]/sqrt(prod(diag(VV)[3:4])),
    ##      valt.b1=V2[1,1],valt.b2=V2[2,2],valt.b1b2=V2[1,2]/sqrt(prod(diag(V2))))
}
##' Calculate ABF under H3, using correlation
##'
##' Assumes SNP A is causal for disease 1, and SNP B is causal for disease 2
##' @title approximate LBF under H3 using correlation approximation
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
lbf.h3.corr <- function(bA1,bB1,bA2,bB2,vA1,vB1,vA2,vB2,share.factor,rho=1,W) {
    B <- c(bA1,bB1,bA2,bB2)
    v1star <- vstar_from_v(vA1,vB1,rho)
    v2star <- vstar_from_v(vA2,vB2,rho)
    ov1 <- share.factor * sqrt(v1star[1,1] * v2star[1,1]) ## same variant, different traits
    ov2 <- share.factor * sqrt(v1star[2,2] * v2star[2,2])
    ## ov12 <- share.factor *  sqrt(v1star[1,2] * v2star[1,2])
    ## ov21 <- share.factor *  sqrt(v1star[1,2] * v2star[1,2])
    ## ov12 <- share.factor  * sqrt(v1star[1,2] * v2star[1,2]) 
    ## ov21 <- share.factor  * sqrt(v1star[1,2] * v2star[1,2])
    ##colocc4
    ov12 <-  -rho * sqrt(vA1*vB2) * share.factor /(1-rho^2)
    ov21 <-  -rho * sqrt(vA2*vB1) * share.factor /(1-rho^2)
    V <- as.matrix(bdiag(v1star,v2star))
    V[1,3] <- V[3,1] <- ov1
    V[2,4] <- V[4,2] <- ov2
    V[1,4] <- V[4,1] <- ov12
    V[2,3] <- V[3,2] <- ov21
    ## if(rho!=0) {
    ##     Wmat <- diag(c(W,rho^2*W,rho^2*W,W))
    ##     ## Wmat[cbind(c(1,2,3,4),c(2,1,4,3))] <- sign(rho) * rho^2*W
    ##     Wmat[cbind(c(1,2,3,4),c(2,1,4,3))] <- sign(rho) * rho^2*W
    ## } else {
        Wmat <- diag(c(W,0,0,W))
    ## }
                
    ret <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
      dmvnorm(x=B,mean=rep(0,4), sigma=V, log=TRUE)
    ## if(!is.finite(ret))
    ##     ret <-  dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat, log=TRUE)
    ret
}
##' variance-covariance matrix for multinomial regression at a single SNP from coefficients from marginal logistic regressions
##'
##' assume allele freq known in controls
##' @title vcov estimation for H4, multinomial approx
##' @param b1
##' @param b2
##' @param n00
##' @param n1
##' @param n2
##' @param f0
##' @return list with diagnonal and off diagonal elements of vcov matrix
##' @author Chris Wallace
approxv2.multinom <- function(b1,b2,n00,n1,n2,f0) {
    ## geno is matrix giving relative frequency of individuals with
    ## genotype i,j in population, where i is count of SNP 1 alt
    ## allele, and j is count of SNP2 alt allele
    geno00 <- estgeno.4(n00,f0,f0,rho=1)
    geno1 <- estgeno.4(n1,f0,f0,rho=1,OR=exp(b1))
    geno2 <- estgeno.4(n2,f0,f0,rho=1,OR=exp(b2))
    geno <- geno00 + geno1 + geno2
    ## estimate a1, a2 
    bvec <- c(0,1,2)
    b1m <- matrix(bvec*b1,3,3)
    b2m <- matrix(bvec*b2,3,3,byrow=TRUE)
    G <- geno00+geno1
    fun <- function(a) {
        denom <- 1+exp(a + b1m)
        num1 <- exp(a + b1m)
        abs(sum(G*num1/denom) - n1)
    }
    ## microbenchmark(ans <- optimize(fun,c(-100,100)))
    ans <- optimize(fun,c(-100,100))
    a1 <- ans$minimum
    G <- geno00+geno2
    fun <- function(a) {
        denom <- 1+exp(a + b2m)
        num2 <- exp(a + b2m)
        abs(sum(G*num2/denom) - n2)
    }
    ans <- optimize(fun,c(-100,100))
    a2 <- ans$minimum

    ## use these, and b1, b2, to estimate entries in var-covar matrix
    denom <- 1+exp(a1 + b1m) + exp(a2 + b2m)
    num1 <- exp(a1 + b1m)
    num2 <- exp(a2 + b2m)
    phat1 <- - geno * num1/denom * (1-num1/denom) # have - here so don't need below
    phat2 <- - geno * num2/denom * (1-num2/denom)
    phat12 <- geno * num1/denom * num2/denom
    va1 <- sum(phat1) 
    va2 <- sum(phat2)
    va12 <- sum(phat12)
    vb1 <-  t(c(0,1,2)) %*% phat1 %*% c(0,1,2) # j=j'
    vb2 <-  t(c(0,1,2)) %*% phat2 %*% c(0,1,2) # j=j'
    
    ## va1b1 <- t(c(1,1,1)) %*% ( geno * num1/denom * (1-num1/denom) ) %*% c(0,1,2) # j=j'
    va1b1 <- t(c(0,1,2)) %*% phat1 %*% c(1,1,1) # j=j'
    va2b2 <- t(c(1,1,1)) %*% phat2 %*% c(0,1,2) # j=j'
    ## va2b2o <- t(c(0,1,2)) %*% ( geno * num2/denom * (1-num2/denom) ) %*% c(1,1,1) # j=j'
    va1b2 <- t(c(1,1,1)) %*% phat12 %*% c(0,1,2) # j!=j'
    va2b1 <- t(c(0,1,2)) %*% phat12 %*% c(1,1,1) # j!=j'
    vb1b2 <- t(c(0,1,2)) %*% phat12 %*% c(0,1,2) # j!=j'
    H <- matrix(c(va1,va12,va1b1,va1b2,
                  va12,va2,va2b1,va2b2,
                  va1b1,va2b1,vb1,vb1b2,
                  va1b2,va2b2,vb1b2,vb2),4,4)
    V <- ginv(-H)[3:4,3:4] # var-cov matrix for (a1,a2,b1,b2)
    list(v.b1=diag(V)[1],v.b2=diag(V)[2],v.b1b2=V[1,2]/sqrt(prod(diag(V))))
}
calc.hessian3 <- function(ka,kA,kB,phat1,phat2,phat12) {
    k <- list(ka,ka,kA,kA,kB,kB)
    f <- function(k1,k2,phat) sum(k1*k2*phat)
    H <- matrix(0,6,6,dimnames=list(c("a1","a2","A1","A2","B1","B2"),
                                    c("a1","a2","A1","A2","B1","B2")))
    diag(H) <- mapply(f,k,k,list(phat1,phat2,phat1,phat2,phat1,phat2))
    H["a1","a2"] <- f(ka,ka,phat12)
    H["A1","A2"] <- f(kA,kA,phat12)
    H["B1","B2"] <- f(kB,kB,phat12)
    H["a1","B2"] <- H["a2","B1"] <- f(ka,kB,phat12)
    H["a1","A2"] <- H["a2","A1"] <- f(ka,kA,phat12)
    H["A1","B2"] <- H["A2","B1"] <- f(kA,kB,phat12)
    H["a1","A1"] <- f(ka,kA,phat1)
    H["a1","B1"] <- f(ka,kB,phat1)
    H["A1","B1"] <- f(kA,kB,phat1)
    H["a2","A2"] <- f(ka,kA,phat2)
    H["a2","B2"] <- f(ka,kB,phat2)
    H["A2","B2"] <- f(kA,kB,phat2)
    ## make symmetric
    tH <- t(H)
    tH[upper.tri(tH)] <- H[upper.tri(H)]
    tH 
}
vstar_from_v <- function(v1,v2,r) {
    ov <- -r*sqrt(v1*v2)
    matrix(c(v1,ov,ov,v2),2,2) / (1-r^2)
}
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param n total number of individuals
##' @param fA af of SNP A in controls
##' @param fB af of SNP B in controls
##' @param rho corr(A,B) in controls (set rho=1 and fB=fA for a single SNP)
##' @param OR OR for A (set OR=1 for controls)
##' @return 
##' @author Chris Wallace
estgeno.3 <- function(n,fA,fB,rho,ORA=1,ORB=1) {
    ## allele <- match.arg(allele)
    rr <- rho * sqrt(fA * (1 - fA) * fB * (1 - fB))
    condB <-  c((fB * (1-fA) - rr)/(1-fA), (fB * fA + rr)/fA)
    ## geno freq in controls
    H <- matrix(c((1-fA) * (1-condB[1]), (1-fA) * condB[1],
                  fA * (1-condB[2]), fA * condB[2]),
                nrow=2,ncol=2,dimnames=list(c("A0","A1"),c("B0","B1")))
    G <- n*matrix(c(H[1,1]^2,        2*H[1,1]*H[1,2],                   H[1,2]^2,
                    2*H[1,1]*H[2,1], 2*H[1,1]*H[2,2] + 2*H[1,2]*H[2,1], 2*H[1,2]*H[2,2],
                    H[2,1]^2,        2*H[2,1]*H[2,2],                   H[2,2]^2),
                  nrow=3,ncol=3,dimnames=list(c("A0","A1","A2"),c("B0","B1","B2")))
    if(ORA!=1 | ORB!=1) {
        ## adjust overall allele freqs, using both OR
        OR <- outer(c(1,ORA,ORA^2),c(1,ORB,ORB^2))
        terms <- G * OR / G["A0","B0"]
        G <- n * terms/sum(terms)
    }
    G
}
estgeno.4 <- function(n,fA,fB,rho,OR=1) {
    ## allele <- match.arg(allele)
    rr <- rho * sqrt(fA * (1 - fA) * fB * (1 - fB))
    condB <-  c((fB * (1-fA) - rr)/(1-fA), (fB * fA + rr)/fA)
    if(OR==1) {
      H <- matrix(c((1-fA) * (1-condB[1]), (1-fA) * condB[1],
                    fA * (1-condB[2]), fA * condB[2]),
                  nrow=2,ncol=2,dimnames=list(c("A0","A1"),c("B0","B1")))
    }
    ## h11 <- H[1,1]
    ## h12 <- H[1,2]
    ## h21 <- H[2,1]
    ## h22 <- H[2,2]
    if(OR!=1) {
        ## adjust overall allele freqs, using both OR
        ## g0 <-  OR * fA / (OR * fA + 1-fA)
        ## Allele freq in cases at causal SNP
        gA <-  OR * fA / (OR * fA + 1 - fA)
        H <- matrix(c((1-gA) * (1-condB[1]), (1-gA) * condB[1],
                      gA * (1-condB[2]), gA * condB[2]),
                    nrow=2,ncol=2,dimnames=list(c("A0","A1"),c("B0","B1")))
    }
    n*matrix(c(H[1,1]^2,        2*H[1,1]*H[1,2],                   H[1,2]^2,
               2*H[1,1]*H[2,1], 2*H[1,1]*H[2,2] + 2*H[1,2]*H[2,1], 2*H[1,2]*H[2,2],
               H[2,1]^2,        2*H[2,1]*H[2,2],                   H[2,2]^2),
             nrow=3,ncol=3,dimnames=list(c("A0","A1","A2"),c("B0","B1","B2")))
}

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
coloc.cc <- function(dataset1,dataset2,
                     n1=args$ncse,n2=args$ncse,n00=(1-dataset1$s)*dataset1$N, n01=0,n02=0,
                     MAF, LD,
                     p1=1e-4, p2=1e-4, p12=1e-5,
                     method=c("multinom","corr")) {
    method <- match.arg(method)
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
    if(method %in% c("corr","corr5")) { # vectorisable
        share.factor <- n00 * sqrt(n1*n2)/sqrt((n01+n00)*(n02+n00)*(n01+n00+n1)*(n02+n00+n2))
        df[,ev1:=v1]
        df[,ev2:=v2]
        df[,erho:=share.factor]
    } else {
        for(i in 1:nrow(df)) {
            df[i, c("ev1","ev2","erho"):=approxv2.multinom(b1,b2,n00,n1,n2,f0)]
        }
    }
    ## lbf for h1,h2,h4 
    df[, lbf1:=bf2(b1,b2,ev1,ev2,erho,w2=0)]
    df[, lbf2:=bf2(b1,b2,ev1,ev2,erho,w1=0)]
    df[, lbf4:=bf2(b1,b2,ev1,ev2,erho)]

    ## for h3, we need pairs of snps
    df3 <- as.data.table(as.data.frame(expand.grid(snps,snps)))
    setnames(df3,c("snp1","snp2"))
    df3 <- df3[snp1!=snp2,]
    df3[,r:=LD[cbind(as.character(snp1),as.character(snp2))]]
    df3 <- merge(df3,df[,.(snp,f0,b1,v1)],
                by.x="snp1",by.y="snp",all.x=TRUE) # trait 1, snp A
    setnames(df3,c("f0","b1","v1"),c("fA","bA1","vA1"))
    df3 <- merge(df3,df[,.(snp,f0,b2,v2)],
                by.x="snp2",by.y="snp",all.x=TRUE) # trait 2, snp B
    setnames(df3,c("f0","b2","v2"),c("fB","bB2","vB2"))
    df3 <- merge(df3, df[,.(snp,b1,v1)],
                 by.x="snp2",by.y="snp",all.x=TRUE) # trait 1, snp B
    setnames(df3,c("b1","v1"),c("bB1","vB1"))
    df3 <- merge(df3, df[,.(snp,b2,v2)],
                 by.x="snp1",by.y="snp",all.x=TRUE) # trait 2, snp A
    setnames(df3,c("b2","v2"),c("bA2","vA2"))

    ## reconstruct coefficients from bivariate binomial logistics 
    df3[,c("bA1star","bB1star"):=bstar_from_b(bA1,bB1,vA1,vB1,r)]
    df3[,c("bA2star","bB2star"):=bstar_from_b(bA2,bB2,vA2,vB2,r)]
     ## special case: r==1. This won't allow a multinomial to be fitted at all.  Instead, add in the bf from h4 for one of the two SNPs
    if(any(abs(df3$r) >= 0.99)) {
        wh <- which(abs(df3$r)>0.99)
        m <- match(df3$snp1[wh], df$snp)
        df3[wh, lbf3:=df$lbf4[m]]
    }
    if(method=="corr") {
       for(i in which(df3$r < 0.99)) 
            df3[i,lbf3:=lbf.h3.corr(bA1star,bB1star,bA2star,bB2star,vA1,vB1,vA2,vB2,share.factor,r,W=0.04)]
    }  else {
        for(i in which(df3$r < 0.99)) 
            df3[i,lbf3:=lbf.h3.multinom(bA1=bA1star,bB1=bB1star,bA2star,bB2star,n00,n1,n2,fA,fB,r,W=0.04)]
    }
    
    ## put it all together
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

    ret <- c(length(snps), sdf[method=="lbf",]$pp)
    names(ret) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", 
"PP.H4.abf")
    ret
}
 


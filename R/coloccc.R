## bhat.from.bstar <- function(b1star, b2star, v1star, v2star, rho) {
##     w1 <- 1/(1-rho^2)
##     w2 <- rho/(sqrt(v1star*v2star) * (1-rho^2))
##     b1 <- w1 * b1star - w2 * b2star * v1star
##     b2 <- w1 * b2star - w2 * b1star * v2star
##     v1 <- w1*v1star
##     v2 <- w1*v2star
##     list(b1,b2,v1,v2)
## }
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}
##' approximate bayes factor for bivariate normal
##'
##' @title ABF for bivariate normal
##' @inheritParams bstar_from_b
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
##' estimate bivariate logistic regression beta vector from univariate summary stats
##'
##' supply beta and var(beta) for two univariate logistic regressions, and an estimate of the correlation between explanatory variables in controls
##' @title recreate bivariate logistic regression from univariate summary stats
##' @param b1 beta_1
##' @param b2 beta_2
##' @param v1 v(beta_1)
##' @param v2 v(beta_2)
##' @param rho cor(snp1,snp2)
##' @return  list of b1, b2 in bivariate regression
##' @author Chris Wallace
## bstar_from_b <- function(b1, b2, v1, v2, rho){
## bstar1 <- b1/(1 - rho^2) - rho/(1 - rho^2) * sqrt(v1/v2) * b2
## bstar2 <- b2/(1 - rho^2) - rho/(1 - rho^2) * sqrt(v2/v1) * b1
## ## vstar <- 1/(1 - rho^2) * matrix(c(v1, -rho * sqrt(v1 * v2), -rho * sqrt(v1 * v2), v2), nrow = 2, byrow = TRUE)
## list(bstar1, bstar2)
## }

## ##' Calculate ABF under H3, multinomial
## ##'
## ##' Assumes SNP A is causal for disease 1, and SNP B is causal for disease 2
## ##' @title Approximate LBF under H3, using multinomial approximation
## ##' @param bA1 log OR for SNP A, disease 1 <--- causal for disease 1
## ##' @param bB1 log OR for SNP B, disease 1
## ##' @param bA2 log OR for SNP A, disease 2
## ##' @param bB2 log OR for SNP A, disease 2 <-- causal for disease 2
## ##' @param n00 number of shared controls 
## ##' @param n1 number of disease 1 cases
## ##' @param n2 number of disease 1 cases
## ##' @param fA MAF of A in controls
## ##' @param fB MAF of B in controls
## ##' @param rho correlation (r) between A and B in controls
## ##' @param W variance on normal prior for effect sizes at causal variants
## ##' @return lBF for A, B from the prior involving W to a prior assuming true effect size is 0
## ##' @author Chris Wallace
## lbf.h3.multinom <- function(bA1,bB1,bA2,bB2,n00,n1,n2,fA,fB,rho,W) {
##     ## geno is matrix giving relative frequency of individuals with
##     ## genotype i,j in population, where i is count of SNP 1 alt
##     ## allele, and j is count of SNP2 alt allele
##     geno00 <- estgeno.3(n00,fA,fB,rho)
##     geno1 <- estgeno.3(n1,fA,fB,rho,ORA=exp(bA1),ORB=exp(bB1))
##     geno2 <- t(estgeno.3(n2,fA=fB,fB=fA,rho,ORA=exp(bB2),ORB=exp(bA2))) # swap snps and t() because estgeno assumes OR for first SNP
##     geno <- geno00+geno1+geno2
##     ## estimate a1, a2 
##     bvec <- c(0,1,2)
##     b1m <- matrix(bvec*bA1,3,3) + matrix(bvec*bB1,3,3,byrow=TRUE)
##     b2m <- matrix(bvec*bA2,3,3) + matrix(bvec*bB2,3,3,byrow=TRUE)
##     G <- geno00+geno1
##     fun <- function(a) {
##         denom <- 1+exp(a + b1m)
##         num1 <- exp(a + b1m)
##         abs(sum(G*num1/denom) - n1)
##     }
##     ans <- optimize(fun,c(-100,100))
##     a1 <- ans$minimum
##     G <- geno00+geno2
##     fun <- function(a) {
##         denom <- 1+exp(a + b2m)
##         num2 <- exp(a + b2m)
##         abs(sum(G*num2/denom) - n2)
##     }
##     ans <- optimize(fun,c(-100,100))
##     a2 <- ans$minimum
##     ## use https://czep.net/stat/mlelr.pdf to fix next few lines, have bA1, bB1 etc
##     ## eqn 37 i=indiv, j=trait, k=SNP
##     num1 <- exp(a1 + b1m)
##     num2 <- exp(a2 + b2m)
##     denom <- 1 + num1 + num2
##     phat1 <- - geno * num1/denom * (1-num1/denom) # have - here so don't need below
##     phat2 <- - geno * num2/denom * (1-num2/denom)
##     phat12 <- geno * num1/denom * num2/denom
##     ka <- matrix(1,3,3)
##     kA <- matrix(bvec,3,3)
##     kB <- matrix(bvec,3,3,byrow=TRUE)
##     H <- calc.hessian3(ka,kA,kB,phat1,phat2,phat12)
##     V <- ginv(-H)[3:6,3:6] # var-cov matrix for (a1,a2,b1,b2)
##     ## if(rho>0) {
##     ##     Wmat <- diag(c(W,rho^2*W,rho^2*W,W))
##     ##     Wmat[cbind(c(1,2,3,4),c(2,1,4,3))] <- rho^2*W
##     ## } else {
##         Wmat <- diag(c(W,0,0,W))
##     ## }
##     dmvnorm(c(bA1,bA2,bB1,bB2), mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
##       dmvnorm(c(bA1,bA2,bB1,bB2), mean=rep(0,4), sigma=V, log=TRUE)
##     ## return new v1, v2, rho
##     ## c(v1,v2,rho)
##     ## list(v.b1=diag(VV)[3],v.b2=diag(VV)[4],v.b1b2=VV[3,4]/sqrt(prod(diag(VV)[3:4])),
##     ##      valt.b1=V2[1,1],valt.b2=V2[2,2],valt.b1b2=V2[1,2]/sqrt(prod(diag(V2))))
## }
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
## lbf.h3.corr <- function(bA1,bB1,bA2,bB2,vA1,vB1,vA2,vB2,share.factor,rho=1,W) {
##     B <- c(bA1,bB1,bA2,bB2)
##     v1star <- vstar_from_v(vA1,vB1,rho)
##     v2star <- vstar_from_v(vA2,vB2,rho)
##     ov1 <- share.factor * sqrt(v1star[1,1] * v2star[1,1]) ## same variant, different traits
##     ov2 <- share.factor * sqrt(v1star[2,2] * v2star[2,2])
##     ## ov12 <- share.factor *  sqrt(v1star[1,2] * v2star[1,2])
##     ## ov21 <- share.factor *  sqrt(v1star[1,2] * v2star[1,2])
##     ## ov12 <- share.factor  * sqrt(v1star[1,2] * v2star[1,2]) 
##     ## ov21 <- share.factor  * sqrt(v1star[1,2] * v2star[1,2])
##     ##colocc4
##     ov12 <-  -rho * sqrt(vA1*vB2) * share.factor /(1-rho^2)
##     ov21 <-  -rho * sqrt(vA2*vB1) * share.factor /(1-rho^2)
##     V <- as.matrix(bdiag(v1star,v2star))
##     V[1,3] <- V[3,1] <- ov1
##     V[2,4] <- V[4,2] <- ov2
##     V[1,4] <- V[4,1] <- ov12
##     V[2,3] <- V[3,2] <- ov21
##     ## if(rho!=0) {
##     ##     Wmat <- diag(c(W,rho^2*W,rho^2*W,W))
##     ##     ## Wmat[cbind(c(1,2,3,4),c(2,1,4,3))] <- sign(rho) * rho^2*W
##     ##     Wmat[cbind(c(1,2,3,4),c(2,1,4,3))] <- sign(rho) * rho^2*W
##     ## } else {
##         Wmat <- diag(c(W,0,0,W))
##     ## }
                
##     ret <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat, log=TRUE) -
##       dmvnorm(x=B,mean=rep(0,4), sigma=V, log=TRUE)
##     ## if(!is.finite(ret))
##     ##     ret <-  dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat, log=TRUE)
##     ret
## }
##' variance-covariance matrix for multinomial regression at a single SNP from coefficients from marginal logistic regressions
##'
##' assume allele freq known in controls
##' @title vcov estimation for H4, multinomial approx
##' @inheritParams bstar_from_b
##' @inheritParams coloc.cc
##' @param f0 MAF of SNP
##' @return list with diagnonal and off diagonal elements of vcov matrix
##' @author Chris Wallace
## approxv2.multinom <- function(b1,b2,n00,n1,n2,f0) {
##     ## geno is matrix giving relative frequency of individuals with
##     ## genotype i,j in population, where i is count of SNP 1 alt
##     ## allele, and j is count of SNP2 alt allele
##     geno00 <- estgeno.4(n00,f0)
##     geno1 <- estgeno.4(n1,f0,OR=exp(b1))
##     geno2 <- estgeno.4(n2,f0,OR=exp(b2))
##     geno <- geno00 + geno1 + geno2
##     ## estimate a1, a2 
##     bvec <- c(0,1,2)
##     b1m <- matrix(bvec*b1,3,3)
##     b2m <- matrix(bvec*b2,3,3,byrow=TRUE)
##     G <- geno00+geno1
##     fun <- function(a) {
##         denom <- 1+exp(a + b1m)
##         num1 <- exp(a + b1m)
##         abs(sum(G*num1/denom) - n1)
##     }
##     ## microbenchmark(ans <- optimize(fun,c(-100,100)))
##     ans <- optimize(fun,c(-100,100))
##     a1 <- ans$minimum
##     G <- geno00+geno2
##     fun <- function(a) {
##         denom <- 1+exp(a + b2m)
##         num2 <- exp(a + b2m)
##         abs(sum(G*num2/denom) - n2)
##     }
##     ans <- optimize(fun,c(-100,100))
##     a2 <- ans$minimum

##     ## use these, and b1, b2, to estimate entries in var-covar matrix
##     denom <- 1+exp(a1 + b1m) + exp(a2 + b2m)
##     num1 <- exp(a1 + b1m)
##     num2 <- exp(a2 + b2m)
##     phat1 <- - geno * num1/denom * (1-num1/denom) # have - here so don't need below
##     phat2 <- - geno * num2/denom * (1-num2/denom)
##     phat12 <- geno * num1/denom * num2/denom
##     va1 <- sum(phat1) 
##     va2 <- sum(phat2)
##     va12 <- sum(phat12)
##     vb1 <-  t(c(0,1,2)) %*% phat1 %*% c(0,1,2) # j=j'
##     vb2 <-  t(c(0,1,2)) %*% phat2 %*% c(0,1,2) # j=j'
    
##     ## va1b1 <- t(c(1,1,1)) %*% ( geno * num1/denom * (1-num1/denom) ) %*% c(0,1,2) # j=j'
##     va1b1 <- t(c(0,1,2)) %*% phat1 %*% c(1,1,1) # j=j'
##     va2b2 <- t(c(1,1,1)) %*% phat2 %*% c(0,1,2) # j=j'
##     ## va2b2o <- t(c(0,1,2)) %*% ( geno * num2/denom * (1-num2/denom) ) %*% c(1,1,1) # j=j'
##     va1b2 <- t(c(1,1,1)) %*% phat12 %*% c(0,1,2) # j!=j'
##     va2b1 <- t(c(0,1,2)) %*% phat12 %*% c(1,1,1) # j!=j'
##     vb1b2 <- t(c(0,1,2)) %*% phat12 %*% c(0,1,2) # j!=j'
##     H <- matrix(c(va1,va12,va1b1,va1b2,
##                   va12,va2,va2b1,va2b2,
##                   va1b1,va2b1,vb1,vb1b2,
##                   va1b2,va2b2,vb1b2,vb2),4,4)
##     V <- ginv(-H)[3:4,3:4] # var-cov matrix for (a1,a2,b1,b2)
##     list(v.b1=diag(V)[1],v.b2=diag(V)[2],v.b1b2=V[1,2]/sqrt(prod(diag(V))))
## }
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
##' Genotype count matrix estimation for H3
##'
##' @title estimate genotype counts for two SNPs (H3)
##' @param n total number of individuals
##' @param fA af of SNP A in controls
##' @param fB af of SNP B in controls
##' @param rho corr(A,B) in controls (set rho=1 and fB=fA for a single
##'     SNP)
##' @return Estimated genotype count matrix
##' @author Chris Wallace
## estgeno.2 <- function(n,fA,fB,rho) {
##     ## allele <- match.arg(allele)
##     rr <- rho * sqrt(fA * (1 - fA) * fB * (1 - fB))
##     condB <-  c((fB * (1-fA) - rr)/(1-fA), (fB * fA + rr)/fA)
##     ## geno freq in controls
##     H <- matrix(c((1-fA) * (1-condB[1]), (1-fA) * condB[1],
##                   fA * (1-condB[2]), fA * condB[2]),
##                 nrow=2,ncol=2,dimnames=list(c("A0","A1"),c("B0","B1")))
##     n*matrix(c(H[1,1]^2,        2*H[1,1]*H[1,2],                   H[1,2]^2,
##                2*H[1,1]*H[2,1], 2*H[1,1]*H[2,2] + 2*H[1,2]*H[2,1], 2*H[1,2]*H[2,2],
##                H[2,1]^2,        2*H[2,1]*H[2,2],                   H[2,2]^2),
##              nrow=3,ncol=3,dimnames=list(c("A0","A1","A2"),c("B0","B1","B2")))
## }

#estimating bivariate SNP distribution in cases
estgeno.2.cse <- function(coefs, G0) {
    summands <- sapply(1 : 3, FUN = function(i){
        sapply(1 : 3, FUN = function(j){
            G0[i, j] * exp(sum(c(i - 1, j - 1) * coefs))
        })
    })
    g1.00 <- 1/(1 + (sum(summands) - summands[1, 1])/G0[1, 1])
    G1 <- sapply(1 : 3, FUN = function(i){
        sapply(1 : 3, FUN = function(j){
            g1.00/G0[1, 1] * G0[i, j] * exp(sum(c(i - 1, j - 1) * coefs))
        })
    }) %>% t
    G1
}
estgeno.2.ctl <- function(fA, fB, rho){
    r <- rho * sqrt(fA * (1 - fA) * fB * (1 - fB))
    h <- matrix(c((1 - fA) * (1 - fB) + r, fB * (1 - fA) - r,
                  fA * (1 - fB) - r, fA * fB + r),
                nrow = 2, ncol = 2, byrow = TRUE)
    matrix(c(h[1,1]^2, 2*h[1,1]*h[1,2], h[1,2]^2,
                     2*h[1,1]*h[2,1], 2*h[1,1]*h[2,2] + 2*h[1,2] * h[2,1], 2*h[1,2]*h[2,2],
                     h[2,1]^2, 2*h[2,1]*h[2,2], h[2,2]^2),
                   nrow=3,ncol=3, byrow = TRUE)
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param bA log OR for SNP A
##' @param bB log OR for SNP B
##' @param G0 genotype freq in controls
##' @return 
##' @author Stasia Grinberg
bstar_from_bhat <- function(bA,bB,fA,fB,rho) {
    ## G <- estgeno.2.ctl(fA,fB,rho)
    ## sys <- function(logx) {
    ##     x=exp(logx)
    ##     c(f1 = 2 * ((x[1] * G[3,3] + G[2,3] / 2 + G[3,2]) * x[2] ^ 2 + x[2] * G[2,2] / 2 + x[1] * G[3,1] + G[2,1] / 2) * (G[1,1] + G[1,2] + G[1,3] + G[2,1] / 2 + G[2,2] / 2 + G[2,3] / 2) * x[1] / (G[2,1] + G[2,2] + G[2,3] + 2 * G[3,1] + 2 * G[3,2] + 2 * G[3,3]) / ((x[1] * G[2,3] / 2 + G[1,3]) * x[2] ^ 2 + (x[1] * G[2,2] / 2 + G[1,2]) * x[2] + x[1] * G[2,1] / 2 + G[1,1]) - exp(bA),
    ##       f2 = 2 * (x[1] ^ 2 * x[2] * G[3,3] + ((G[2,3] + G[3,2] / 2) * x[2] + G[2,2] / 2) * x[1] + x[2] * G[1,3] + G[1,2] / 2) * x[2] * (G[1,1] + G[2,1] + G[3,1] + G[1,2] / 2 + G[2,2] / 2 + G[3,2] / 2) / (x[1] ^ 2 * G[3,1] + (x[2] ^ 2 * G[3,2] / 2 + x[2] * G[2,2] / 2 + G[2,1]) * x[1] + x[2] * G[1,2] / 2 + G[1,1]) / (G[1,2] + G[2,2] + G[3,2] + 2 * G[1,3] + 2 * G[2,3] + 2 * G[3,3]) - exp(bB))
    ## }
    G0 <- estgeno.2.ctl(fA,fB,rho)
    a1<-2*sum(G0[1, ])+sum(G0[2, ])
    b1<-2*sum(G0[3, ])+sum(G0[2, ])
    a2<-2*sum(G0[, 1])+sum(G0[, 2])
    b2<-2*sum(G0[, 3])+sum(G0[, 2])
    sys <- function(x) {
        c1 <- 2 * sum(G0[1, ] * exp((1 : 3) * x[2])) +
          sum(G0[2, ] * exp(x[1] + (1 : 3) * x[2]))
        d1 <- 2 * sum(G0[3, ] * exp(2 * x[1] + (1 : 3) * x[2])) +
          sum(G0[2, ] * exp(x[1] + (1 : 3) * x[2]))
        c2 <- 2 * sum(G0[, 1] * exp((1 : 3) * x[1])) +
          sum(G0[, 2] * exp((1 : 3) * x[1] + x[2]))
        d2 <- 2 * sum(G0[, 3] * exp((1 : 3) * x[1] + 2 * x[2])) +
          sum(G0[, 2] * exp((1 : 3) * x[1] + x[2]))
        f1 <- log(a1 * d1 /(b1 * c1)) - bA
        f2 <- log(a2 * d2 /(b2 * c2)) - bB
        c(f1 = f1 , f2 = f2)
    }
    ans <- dfsane(par=c(0,0), fn=sys,control=list(trace=FALSE))
    ans$par
    ## optim(c(0,0),sys)
    ## multiroot(f = sys,start = c(0 ,0),verbose=TRUE)$root
}
bstar_from_bhat_pen <- function(bA,bB,fA,fB,rho,lambda=1) {
    ## G <- estgeno.2.ctl(fA,fB,rho)
    ## sys <- function(logx) {
    ##     x=exp(logx)
    ##     c(f1 = 2 * ((x[1] * G[3,3] + G[2,3] / 2 + G[3,2]) * x[2] ^ 2 + x[2] * G[2,2] / 2 + x[1] * G[3,1] + G[2,1] / 2) * (G[1,1] + G[1,2] + G[1,3] + G[2,1] / 2 + G[2,2] / 2 + G[2,3] / 2) * x[1] / (G[2,1] + G[2,2] + G[2,3] + 2 * G[3,1] + 2 * G[3,2] + 2 * G[3,3]) / ((x[1] * G[2,3] / 2 + G[1,3]) * x[2] ^ 2 + (x[1] * G[2,2] / 2 + G[1,2]) * x[2] + x[1] * G[2,1] / 2 + G[1,1]) - exp(bA),
    ##       f2 = 2 * (x[1] ^ 2 * x[2] * G[3,3] + ((G[2,3] + G[3,2] / 2) * x[2] + G[2,2] / 2) * x[1] + x[2] * G[1,3] + G[1,2] / 2) * x[2] * (G[1,1] + G[2,1] + G[3,1] + G[1,2] / 2 + G[2,2] / 2 + G[3,2] / 2) / (x[1] ^ 2 * G[3,1] + (x[2] ^ 2 * G[3,2] / 2 + x[2] * G[2,2] / 2 + G[2,1]) * x[1] + x[2] * G[1,2] / 2 + G[1,1]) / (G[1,2] + G[2,2] + G[3,2] + 2 * G[1,3] + 2 * G[2,3] + 2 * G[3,3]) - exp(bB))
    ## }
    G0 <- estgeno.2.ctl(fA,fB,rho)
    a1<-2*sum(G0[1, ])+sum(G0[2, ])
    b1<-2*sum(G0[3, ])+sum(G0[2, ])
    a2<-2*sum(G0[, 1])+sum(G0[, 2])
    b2<-2*sum(G0[, 3])+sum(G0[, 2])
    sys <- function(x) {
        c1 <- 2 * sum(G0[1, ] * exp((1 : 3) * x[2])) +
          sum(G0[2, ] * exp(x[1] + (1 : 3) * x[2]))
        d1 <- 2 * sum(G0[3, ] * exp(2 * x[1] + (1 : 3) * x[2])) +
          sum(G0[2, ] * exp(x[1] + (1 : 3) * x[2]))
        c2 <- 2 * sum(G0[, 1] * exp((1 : 3) * x[1])) +
          sum(G0[, 2] * exp((1 : 3) * x[1] + x[2]))
        d2 <- 2 * sum(G0[, 3] * exp((1 : 3) * x[1] + 2 * x[2])) +
          sum(G0[, 2] * exp((1 : 3) * x[1] + x[2]))
        f1 <- log(a1 * d1 /(b1 * c1)) - bA 
        f2 <- log(a2 * d2 /(b2 * c2)) - bB 
        ## c(f1 = f1 , f2 = f2)
        f1*f2 + abs(lambda*sum(x^2))
    }
    ans <- spg(par=c(0,0), fn=sys,control=list(trace=FALSE))
    ans$par
    ## optim(c(0,0),sys)
    ## multiroot(f = sys,start = c(0 ,0),verbose=TRUE)$root
}


##function to solve for a
vmat <- function(a,bA,bB,N) {
    bvec <- c(0,1,2)
    bAm <- matrix(bvec*bA,3,3)
    bBm <- matrix(bvec*bB,3,3,byrow=TRUE)
    num <- exp(a + bAm + bBm)
    denom <- 1 + exp(a + bAm + bBm)
    phat <- N * num/denom^2
    xmat <- cbind(1, rep(0:2,3),rep(0:2,each=3))
    D <- diag(as.vector(phat))
    V <- t(xmat) %*% D %*% xmat
}
    
cc.marg2joint <- function(bA,bB, vA,vB, n0,n1, fA,fB,rho) {
    ## est genotype counts in controls
    bvec <- 0:2
    G0 <- estgeno.2.ctl(fA,fB,rho)
    
    ## convert bhat to bstar
    bstar <- bstar_from_bhat(bA,bB,fA,fB,rho)

    solve_for_a <- function(b1, b2, n1, n0, geno) {
        b1m <- matrix(bvec*b1,3,3)
        b2m <- matrix(bvec*b2,3,3,byrow=TRUE)
        fun <- function(a) {
            denom <- 1 + exp(a + b1m + b2m)
            num <- exp(a + b1m + b2m)
            sum(geno*num/denom) - n1/(n1 + n0)
        }
        ## ans <- dfsane(par=0, fn=fun,control=list(trace=FALSE))
        ## ans$par
        uniroot(fun, c(-10, 10),extendInt="yes")$root
    }
    G1 <- estgeno.2.cse(bstar, G0)
    astar <- solve_for_a(bstar[1], bstar[2],
                          n1 = n1, n0 = n0, geno = G1)

    ## within-trait vcov
    ## https://czep.net/stat/mlelr.pdf (16)
    V <- solve(vmat(astar,bstar[1],bstar[2], G0*n0 + G1*n1))
    return(list(coef=c(astar,bstar),V=V))
}

lbf.h3 <- function(bA1,bB1,bA2,bB2,
                   vA1,vB1,vA2,vB2,
                   n00,n01,n02,n1,n2,
                   fA,fB,rho,W=0.04) {

    ## reconstruct joint models
    joint1 <- cc.marg2joint(bA1,bB1,vA1,vB1,n0=n00+n01,n1=n1,fA,fB,rho)
    joint2 <- cc.marg2joint(bA2,bB2,vA2,vB2,n0=n00+n02,n1=n2,fA,fB,rho)
    
    ## ## est genotype counts in controls
    ## bvec <- 0:2
    ## G0 <- estgeno.2.ctl(fA,fB,rho)
    
    ## ## convert bhat to bstar
    ## bstar1 <- bstar_from_bhat(bA1,bB1,fA,fB,rho)
    ## bstar2 <- bstar_from_bhat(bA2,bB2,fA,fB,rho)
    ## solve_for_a <- function(b1, b2, n, n0, geno) {
    ##     b1m <- matrix(bvec*b1,3,3)
    ##     b2m <- matrix(bvec*b2,3,3,byrow=TRUE)
    ##     fun <- function(a) {
    ##         denom <- 1 + exp(a + b1m + b2m)
    ##         num <- exp(a + b1m + b2m)
    ##         abs(sum(geno*num/denom)) - n/(n + n0)
    ##     }
    ##     uniroot(fun, c(-10, 10),extendInt="yes")$root
    ## }
    ## G1 <- estgeno.2.cse(bstar1, G0)
    ## G2 <- estgeno.2.cse(bstar2, G0)
    ## astar1 <- solve_for_a(bstar1[1], bstar1[2],
    ##                       n = n1, n0 = n00 + n01, geno = G1)
    ## astar2 <- solve_for_a(bstar2[1], bstar2[2],
    ##                       n = n2, n0 = n00 + n02, geno = G2)

    ## ## within-trait vcov
    ## ## https://czep.net/stat/mlelr.pdf (16)
    ## V1 <- vmat(astar1,bstar1[1],bstar1[2],G0*(n00 + n01) + G1*(n1))[-1,-1]
    ## V2 <- vmat(astar2,bstar2[1],bstar2[2],G0*(n00 + n02) + G2*(n2))[-1,-1]
    
    sxz <- rho * sqrt(fA * (1 - fA) * fB * (1 - fB)) + 4 * fA * fB
    H2 <- matrix(c(1, 2*fA, 2*fB,
                   2*fA, 2*fA*(1 + fA), sxz,
                   2*fB, sxz, 2*fB*(1 + fB)),
                 nrow = 3, byrow = TRUE)
    H2.inv <- solve(H2)[-1,-1]
    V12 <- n00 * (1 + exp(joint1$coef[1])) * (1 + exp(joint2$coef[1]))/((n00+n01+n1) * (n00+n02+n2)) * H2.inv

    V <- rbind(cbind(joint1$V[-1,-1],V12),cbind(V12,joint2$V[-1,-1]))
    B <- c(joint1$coef[-1],joint2$coef[-2])
    ## abf
    Wmat <- diag(c(W,0,0,W))
    ## dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V + Wmat, log=TRUE) -
    ##   dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V, log=TRUE)
    dalt <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat, log=TRUE) 
    dnull <- dmvnorm(x=B,mean=rep(0,4), sigma=V, log=TRUE)
    list(dalt,dnull,dalt-dnull)
}


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title estimate MAF in cases
##' @param f0 MAF in controls
##' @param b log BF
##' @return MAF in cases
##' @author Stasia Grinberg
get_maf1 <- function(f0, b) {
    f0 * exp(b)/(1 - f0 + f0 * exp(b))
}

estgeno.1 <- function(n,f) {
    n*c((1-f)^2,2*f*(1-f),f^2)
}
    

univariate_a <- function(b,n0,n1,f0,f1) {
     G <-  estgeno.1(n0,f0) +  estgeno.1(n1,f1)
     bvec <- c(0,1,2)*b
     fun <- function(a) {
        num <- exp(a + bvec)
        denom <- 1+exp(a + bvec)
        sum(G*num/denom) - n1
     }
     ## microbenchmark(ans <- optimize(fun,c(-100,100)))
     ## ans <- optimize(fun,c(-100,100))
     ##uniroot(fun, c(-10,10))$root
     ans <- dfsane(par=c(a=0), fn=fun,control=list(trace=FALSE))
     ans$par
}

approx.h124 <- function(b1,b2,n00,n01,n02,n1,n2,f0,f1,f2) {
    ## geno is matrix giving relative frequency of individuals with
    ## genotype i,j in population, where i is count of SNP 1 alt
    ## allele, and j is count of SNP2 alt allele
    a1 <- univariate_a(b1,n00+n01,n1,f0,f1)
    a2 <- univariate_a(b2,n00+n02,n2,f0,f2)
    v12 <- (1+exp(a1)) * (1+exp(a2)) * n00 / ( 2 * f0 * (1-f0) * (n00+n01+n1) * (n00+n02+n2))
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
##' @param method Currently ignored!!
##'     "multonom" or "corr", depending on which
##'     approximation (multinomial or correlation) you wish to use.
##' @return a vector giving the number of SNPs analysed, and the
##'     posterior probabilities of H0 (no causal variant), H1 (causal
##'     variant for trait 1 only), H2 (causal variant for trait 2
##'     only), H3 (two distinct causal variants) and H4 (one common
##'     causal variant)
##' @export
##' @author Chris Wallace
coloc.cc <- function(dataset1,dataset2,
                     n1=dataset1$s*dataset1$N,
                     n2=dataset2$s*dataset2$N,
                     n00=0,
                     n01=round((1-dataset1$s)*dataset1$N - n00),
                     n02=round((1-dataset2$s)*dataset2$N - n00),
                     MAF, LD,
                     p1=1e-4, p2=1e-4, p12=1e-5,
                     method=c("multinom","corr")) {
    ## method <- match.arg(method)
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
    ## get MAF in cases
    df[,f1:=get_maf1(f0, b1)]
    df[,f2:=get_maf1(f0, b2)]

    ## H1, H2, H4 all two beta models
    ## more accurate
    for(i in 1:nrow(df)) {
            df[i, c("v12"):=approx.h124(b1,b2,n00,n01,n02,n1,n2,f0,f1,f2)]
    }
    df[, erho:=v12/sqrt(v1*v2)]
    ## faster
    ## df[,erho:=n00*sqrt(n1*n2)/((n00+n01+n1)*(n00+n02+n2))]
    
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
    
    ## special case: r==1. This won't allow a model to be fitted at all.  Instead, add in the bf from h4 for one of the two SNPs
    df3.r1 <- df3[abs(r)>=0.95,]
    m <- match(df3.r1$snpA, df$snp)
    df3.r1[, lbf3:=df$lbf4[m]]
    
    ## special case: r=0. Can partition as [ H1 for SNP A ] * [ H2 for SNP B ]
    df3.r0 <- df3[abs(r)<0.001,]
    mA <- match(df3.r0$snpA, df$snp)
    mB <- match(df3.r0$snpB, df$snp)
    df3.r0[, lbf3:=df$lbf1[mA] + df$lbf2[mB] ]
    
    ## what's left requires 4 dim Gaussians
    df3 <- df3[abs(r)<0.95 & abs(r)>=0.001,]
    df3 <- merge(df3,df[,.(snp,f0,f1,b1,v1)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 1, snp A
    setnames(df3,c("f0","f1","b1","v1"),c("fA","fA1","bA1","vA1"))
    df3 <- merge(df3,df[,.(snp,f0,f2,b2,v2)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 2, snp B
    setnames(df3,c("f0","f2","b2","v2"),c("fB","fB2","bB2","vB2"))
    df3 <- merge(df3, df[,.(snp,f1,b1,v1)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 1, snp B
    setnames(df3,c("f1","b1","v1"),c("fB1","bB1","vB1"))
    df3 <- merge(df3, df[,.(snp,f2,b2,v2)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 2, snp A
    setnames(df3,c("f2","b2","v2"),c("fA2","bA2","vA2"))
    
    ## reconstruct coefficients from bivariate binomial logistics 
    ## appears not to be vectorisable any more - move into lbf function
    ## df3[,c("bA1star","bB1star"):=bstar_from_b(bA1,bB1,vA1,vB1,r)]
    ## df3[,c("bA2star","bB2star"):=bstar_from_b(bA2,bB2,vA2,vB2,r)]
    for(i in 1:nrow(df3)) 
        df3[i,c("dalt","dnull","lbf3"):=lbf.h3(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
                           vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
                           n00,n01,n02,n1,n2,
                           fA,fB,rho=r,W=0.04)]
    ## attach(df3[isnpA==34 & isnpB==13,])

    ## sometimes approx is very bad - typically extreme r combined with low MAF. Detect these cases as -Inf likelihood under null. replace with no information
    df3[!is.finite(dnull),lbf3:=0]
    
    ## put it all together
    o<-intersect(names(df3),names(df3.r1))
    df3 <- rbind(df3.r1[,o,with=FALSE],df3.r0[,o,with=FALSE],df3[,o,with=FALSE])
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
    print(ret)
    return(list(summary=ret,results=df,df3=df3,sdf=sdf))
}
    ## ## put it all together
    ##     sdf <- df[, lapply(.SD, logsum), .SDcols=grep("lbf",names(df),value=TRUE)]
##     sdf3 <- df3[, lapply(.SD, logsum), .SDcols=grep("lbf",names(df3),value=TRUE)]
##     sdf <- rbind(melt(sdf,measure.vars=names(sdf)),
##                  melt(sdf3,measure.vars=names(sdf3)))
##     ## add h0
##     sdf0 <- melt(sdf3,measure.vars=names(sdf3))
##     sdf0[,value:=0]
##     sdf0[,variable:=sub("3","0",variable)]
##     sdf <- rbind(sdf,sdf0)
    
##     sdf[,h:=as.integer(sub(".*lbf","",variable))]
##     sdf[,method:=sub("[0-4]","",variable)]
##     sdf[,abf:=0]
##     sdf[h==1,abf:=log(p1) + value]
##     sdf[h==2,abf:=log(p2) + value]
##     sdf[h==3,abf:=log(p1) + log(p2) + value]
##     sdf[h==4,abf:=log(p12) + value]
##     sdf[,pp:=exp(abf - logsum(abf)),by="method"]
##     sdf <- sdf[order(method,h),]
##     ## while checking, return this
##     ## return(sdf)

##     ret <- c(length(snps), sdf[method=="lbf",]$pp)
##     names(ret) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", 
## "PP.H4.abf")
##     ret

 


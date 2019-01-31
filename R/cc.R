
## helper functions, not exported
cor2 <- function (x) {
    1/(NROW(x) - 1) * crossprod( scale(x, TRUE, TRUE) )
}
logsum2 <- function(x,y) {
    my.max <- pmax(x,y)
    my.max + log(exp(x-my.max)+exp(y-my.max))
}
lhood.norm2<-function(b1,b2,v11,v22,v12) {
    dt=v11*v22-v12^2
    top=(b1 ^ 2 * v22 - 2 * b1 * b2 * v12 + b2 ^ 2 * v11) / dt
    (-top - log(dt))/2
}

lhood.norm4<-function(bA1,bB1,bA2,bB2,
                      v11,v22,v12,
                      v33,v44,v34,
                      v13,v24,v14) {
   top <- ((-2 * bA1 * bB2 - 2 * bB1 * bA2) * v14 ^ 3 + (bA1 ^ 2 * v44 + 2 * bA1 * bB1 * v34 + 2 * bA1 * bA2 * v24 + bB1 ^ 2 * v33 + 2 * bB1 * bB2 * v13 + bA2 ^ 2 * v22 + 2 * bA2 * bB2 * v12 + bB2 ^ 2 * v11) * v14 ^ 2 + ((-2 * bB2 ^ 2 * v13 - 2 * bA2 ^ 2 * v24 + (2 * bA1 * bB2 + 2 * bB1 * bA2) * v34 - 2 * bA1 * v44 * bA2 - 2 * v33 * bB1 * bB2) * v12 + ((2 * bA1 * bB2 + 2 * bB1 * bA2) * v24 - 2 * bA1 * v44 * bB1 - 2 * bB1 ^ 2 * v34 - 2 * bA2 * bB2 * v22) * v13 + (-2 * bA1 ^ 2 * v34 - 2 * bA1 * bB1 * v33 - 2 * bA2 * bB2 * v11) * v24 + (-2 * bA1 * bA2 * v22 - 2 * bB1 * bB2 * v11) * v34 + 2 * bA1 * v33 * bB2 * v22 + 2 * v44 * bB1 * bA2 * v11) * v14 + (bA2 ^ 2 * v44 - 2 * bA2 * bB2 * v34 + bB2 ^ 2 * v33) * v12 ^ 2 + ((2 * v24 * bA2 * bB2 + 2 * bB1 * (-bA2 * v44 + bB2 * v34)) * v13 - 2 * ((-bA2 * v34 + bB2 * v33) * v24 - bB1 * (v33 * v44 - v34 ^ 2)) * bA1) * v12 + (bB1 ^ 2 * v44 - 2 * bB1 * bB2 * v24 + bB2 ^ 2 * v22) * v13 ^ 2 - 2 * bA1 * (v24 ^ 2 * bA2 - v34 * bB1 * v24 + v22 * (-bA2 * v44 + bB2 * v34)) * v13 + (bA1 ^ 2 * v33 + bA2 ^ 2 * v11) * v24 ^ 2 + 2 * bB1 * v11 * (-bA2 * v34 + bB2 * v33) * v24 + (bA1 ^ 2 * v22 + bB1 ^ 2 * v11) * v34 ^ 2 + 2 * bA2 * bB2 * v11 * v22 * v34 - bA1 ^ 2 * v22 * v33 * v44 - (v33 * v44 * bB1 ^ 2 + v22 * (bA2 ^ 2 * v44 + bB2 ^ 2 * v33)) * v11) / (-v14 ^ 4 + (v11 * v44 + 2 * v12 * v34 + 2 * v13 * v24 + v22 * v33) * v14 ^ 2 + ((-2 * v13 * v44 - 2 * v24 * v33) * v12 - 2 * v34 * (v11 * v24 + v13 * v22)) * v14 + (v33 * v44 - v34 ^ 2) * v12 ^ 2 + 2 * v12 * v13 * v24 * v34 + (v22 * v44 - v24 ^ 2) * v13 ^ 2 + (v33 * v24 ^ 2 - v22 * (v33 * v44 - v34 ^ 2)) * v11)
   dt <- -v11 * v14 ^ 2 * v44 + 2 * v11 * v14 * v24 * v34 + v11 * v22 * v33 * v44 - v11 * v22 * v34 ^ 2 - v11 * v24 ^ 2 * v33 - v12 ^ 2 * v33 * v44 + v12 ^ 2 * v34 ^ 2 + 2 * v12 * v13 * v14 * v44 - 2 * v12 * v13 * v24 * v34 - 2 * v12 * v14 ^ 2 * v34 + 2 * v12 * v14 * v24 * v33 - v13 ^ 2 * v22 * v44 + v13 ^ 2 * v24 ^ 2 - 2 * v13 * v14 ^ 2 * v24 + 2 * v13 * v14 * v22 * v34 + v14 ^ 4 - v14 ^ 2 * v22 * v33
   (-top - log(dt))/2
}


##' Estimate single snp frequency distibutions
##'
##' @title estgeno1
##' @return relative frequency of genotypes 0, 1, 2
##' @export
##' @author Chris Wallace
##' @param f MAF
##' @rdname estgeno1
##' @seealso estgeno2
##' @examples
##' estgeno.1.ctl(f=0.5)
estgeno.1.ctl <- function(f) {
    c((1-f)^2,2*f*(1-f),f^2)
}

##' @param G0 single snp frequency in controls (vector of length 3) - obtained from estgeno.1.ctl
##' @param b log odds ratio
##' @export
##' @rdname estgeno1
##' @examples
##' G0=estgeno.1.ctl(f=0.5)
##' estgeno.1.cse (G0=G0, b=log(2))
estgeno.1.cse <- function(G0,b) {
    g0 <- 1
    g1 <- exp( b - log(G0[1]/G0[2]))
    g2 <- exp( 2*b - log(G0[1]/G0[3]))
    c(g0,g1,g2)/(g0+g1+g2)
}
    
##' Estimate two snp frequency distibutions
##'
##' @title estgeno2
##' @param fA MAF of SNP A
##' @param fB MAF of SNP B
##' @param rho correlation between SNPs A, B (between -1 and 1)
##' @return matrix giving relative frequency of each genotype. Element (i,j) gives frequency of A=i-1, B=j-1.
##' @author Chris Wallace
##' @rdname estgeno2
##' @export
##' @seealso estgeno1
##' @examples
##' ## two unlinked SNPs
##' estgeno.2.ctl(fA=0.5,fB=0.1,rho=0)
##' ## two perfectly correlated SNPs
##' estgeno.2.ctl(fA=0.5,fB=0.5,rho=1)
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

##' @param coefs vector of marginal log OR for SNP A and SNP B
##' @param G0 matrix of relative frequencies in controls, often derived from estgeno.2.ctl
##' @rdname estgeno2
##' @export
##' @examples
##' ## two unlinked SNPs, only one of which (A) is associated
##' G0=estgeno.2.ctl(fA=0.5,fB=0.5,rho=0)
##' estgeno.2.cse(coefs=c(log(2),log(1)), G0=G0)
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
    }) 
    t(G1)
}

################################################################################

##' Estimate the coefficients of the joint model (bstar) from those of the marginal model (bhat)
##'
##' The marginal model fits
##'
##' Y ~ A
##'
##' Y ~ B
##'
##' and returns separate estimates bhat for A and B.  This is the
##' usual practice in GWAS.  We need bstar, the estimates from the joint
##' model
##'
##' Y ~ A + B
##' @title Estimate beta-star
##' @param bA log OR for SNP A
##' @param bB log OR for SNP B
##' @inheritParams estgeno.2.ctl
##' @return vector of length 2, bstar
##' @author Stasia Grinberg
bstar_from_bhat <- function(bA,bB,fA,fB,rho) {
    G0 <- estgeno.2.ctl(fA,fB,rho)
    a1<-2*sum(G0[1, ])+sum(G0[2, ])
    b1<-2*sum(G0[3, ])+sum(G0[2, ])
    a2<-2*sum(G0[, 1])+sum(G0[, 2])
    b2<-2*sum(G0[, 3])+sum(G0[, 2])
    sys <- function(x) {
        c1 <- 2 * sum(G0[1, ] * exp((0:2) * x[2])) +
          sum(G0[2, ] * exp(x[1] + (0:2) * x[2]))
        d1 <- 2 * sum(G0[3, ] * exp(2 * x[1] + (0:2) * x[2])) +
          sum(G0[2, ] * exp(x[1] + (0:2) * x[2]))
        c2 <- 2 * sum(G0[, 1] * exp((0:2) * x[1])) +
          sum(G0[, 2] * exp((0:2) * x[1] + x[2]))
        d2 <- 2 * sum(G0[, 3] * exp((0:2) * x[1] + 2 * x[2])) +
          sum(G0[, 2] * exp((0:2) * x[1] + x[2]))
        f1 <- log(a1 * d1 /(b1 * c1)) - bA
        f2 <- log(a2 * d2 /(b2 * c2)) - bB
        c(f1 = f1 , f2 = f2)
    }
    ans <- dfsane(par=c(0,0), fn=sys,control=list(trace=FALSE))
    ans$par
    ## optim(c(0,0),sys)
    ## multiroot(f = sys,start = c(0 ,0),verbose=TRUE)$root
}

##function to solve for a
.vmat2<- function(a,bA,bB,N) {
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
    
.cfv2 <- function(bA,bB, vA,vB, n0,n1, fA,fB,rho) {
    ## est genotype counts in controls
    bvec <- 0:2
    G0 <- estgeno.2.ctl(fA,fB,rho)
    
    ## convert bhat to bstar
    bstar <- bstar_from_bhat(bA,bB,fA,fB,rho)
    G1 <- estgeno.2.cse(bstar, G0)

    ## solve_for_a <- function(b1, b2, n1, n0, geno) {
    b1m <- matrix(bvec*bstar[1],3,3)
    b2m <- matrix(bvec*bstar[2],3,3,byrow=TRUE)
    fun <- function(a) {
        denom <- 1 + exp(a + b1m + b2m)
        num <- exp(a + b1m + b2m)
        sum(G1*num/denom) - n1/(n1 + n0)
    }
        ## ans <- dfsane(par=0, fn=fun,control=list(trace=FALSE))
        ## ans$par
    astar <- uniroot(fun, c(-10, 10),extendInt="yes")$root
    ## within-trait vcov
    ## https://czep.net/stat/mlelr.pdf (16)
    V <- solve(.vmat2(astar,bstar[1],bstar[2], G0*n0 + G1*n1))
    return(list(coef=c(astar,bstar),V=V))
}
.makecov1 <- function(G, cf1, cf2) {
    x1 <- 0:2
    pred1 <- exp(cf1[1] + cf1[2]*x1)
    pred2 <- exp(cf2[1] + cf2[2]*x1)
    pred <- c(G) * pred1 /(1+pred1) * pred2/(1+pred2)
    matrix(c(sum(pred), sum(x1 * pred),
             sum(x1 * pred), sum(x1^2 * pred)), 2,2)
}
.makecov2 <- function(G, cf1, cf2) {
    x1 <- c(0,1,2,0,1,2,0,1,2)
    x2 <- c(0,0,0,1,1,1,2,2,2)
    pred1 <- exp(cf1[1] + cf1[2]*x1 + cf1[3]*x2)
    pred2 <- exp(cf2[1] + cf2[2]*x1 + cf2[3]*x2)
    pred <- c(G) * pred1 /(1+pred1) * pred2/(1+pred2)
    matrix(c(sum(pred), sum(x1 * pred), sum(x2 * pred),
             sum(x1 * pred), sum(x1^2 * pred), sum(x1*x2*pred),
             sum(x2 * pred), sum(x1*x2*pred), sum(x2^2*pred)),
           3,3)
}


.V.h124 <- function(b1,b2,
                         v1,v2,
                         n00,n01,n02,n1,n2,
                         f0,W=0.04,method="by4") {

    ## reconstruct joint models - list of coef(a, b1, b2) and V=vcov(coef)
        joint1 <- .cfv1(b1,n0=n00+n01,n1=n1,f0=f0)
        joint2 <- .cfv1(b2,n0=n00+n02,n1=n2,f0=f0)
    ## print(joint1)
    ## print(joint2)
    ##   cr <- sqrt(exp(joint1$coef[1] + joint2$coef[2]))*n00/sqrt((n00+n01+n1)*(n00+n02+n2))
    ## VV <- cr * sqrt(joint1$V * joint2$V)

 
## ## est genotype counts in controls
    bvec <- 0:2
    G0 <- estgeno.1.ctl(f0)
    ## G1 <- estgeno.1.cse(G0,b1)
    ## G2 <- estgeno.1.cse(G0,b2)
     ## I matrices depend on all subjects
    ## I1 <- makeI2(G0*(n00+n01) + G1*n1, joint1$coef)
    
    ## ## p <- matrix(G[,"pred"],3,3)
    ## ## I1 <- XX * p/(1+p)^2
    ## I2 <- makeI2(G0*(n00+n02) + G2*n2, joint2$coef)
    
    ## ## ## cov depends only on shared controls
    cv <- .makecov1(G0*n00,joint1$coef,joint2$coef)
        
    ## ## cross trait vcov
    ## VV <- ( solve(I1) %*% cv %*% solve(I2) )
    VV <- ( joint1$V[,2] %*% cv %*% joint2$V[2,] )

    list(v11=joint1$V[2,2], v12=VV, v22=joint2$V[2,2])
}

    ## ## B <- c(joint1$coef,joint2$coef)
    ## ## ## abf
    ## ## ## dalt <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V + Wmat, log=TRUE) 
    ## ## ## dnull <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V, log=TRUE)

    ## ## ensure VV is symmetric
    ## ## VV[1,2] <- VV[2,1]  <- mean(VV[1,2],VV[2,1])

    ## Wmat0 <- diag(c(0,0))
    ## Wmat1 <- diag(c(W,0))
    ## Wmat2 <- diag(c(0,W))
    ## Wmat4 <- diag(c(W,W))
    ## V <- matrix(c(joint1$V[2,2],VV,VV,joint2$V[2,2]),2,2)
    ## ## print(V)
    ## ## print(Wmat1)
    ## B <- c(joint1$coef[-1],joint2$coef[-1])
    ##     dalt1 <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat1), log=TRUE) 
    ##     dalt2 <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat2), log=TRUE) 
    ##     dalt4 <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat4), log=TRUE) 
    ##     dnull <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat0), log=TRUE)
    ##     ## dalt1 <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat1, log=TRUE) 
    ##     ## dalt2 <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat2, log=TRUE) 
    ##     ## dalt4 <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat4, log=TRUE) 
    ##     ## dnull <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat0, log=TRUE)
    ## list(dalt1-dnull, dalt2-dnull, dalt4-dnull)
    
## }

.lbf.h124 <- function(b1,b2,
                         v1,v2,
                         n00,n01,n02,n1,n2,
                         f0,W=0.04,method="by4") {

    ## reconstruct joint models - list of coef(a, b1, b2) and V=vcov(coef)
        joint1 <- .cfv1(b1,n0=n00+n01,n1=n1,f0=f0)
        joint2 <- .cfv1(b2,n0=n00+n02,n1=n2,f0=f0)
    ##   cr <- sqrt(exp(joint1$coef[1] + joint2$coef[2]))*n00/sqrt((n00+n01+n1)*(n00+n02+n2))
    ## VV <- cr * sqrt(joint1$V * joint2$V)

 
## ## est genotype counts in controls
    bvec <- 0:2
    G0 <- estgeno.1.ctl(f0)
    ## G1 <- estgeno.1.cse(G0,b1)
    ## G2 <- estgeno.1.cse(G0,b2)
     ## I matrices depend on all subjects
    ## I1 <- makeI2(G0*(n00+n01) + G1*n1, joint1$coef)
    
    ## ## p <- matrix(G[,"pred"],3,3)
    ## ## I1 <- XX * p/(1+p)^2
    ## I2 <- makeI2(G0*(n00+n02) + G2*n2, joint2$coef)
    
    ## ## ## cov depends only on shared controls
    cv <- .makecov1(G0*n00,joint1$coef,joint2$coef)
        
    ## ## cross trait vcov
    ## VV <- ( solve(I1) %*% cv %*% solve(I2) )
    VV <- ( joint1$V[,2] %*% cv %*% joint2$V[2,] )
        
    ## B <- c(joint1$coef,joint2$coef)
    ## ## abf
    ## ## dalt <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V + Wmat, log=TRUE) 
    ## ## dnull <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V, log=TRUE)

    ## ensure VV is symmetric
    ## VV[1,2] <- VV[2,1]  <- mean(VV[1,2],VV[2,1])

    Wmat0 <- diag(c(0,0))
    Wmat1 <- diag(c(W,0))
    Wmat2 <- diag(c(0,W))
    Wmat4 <- diag(c(W,W))
    V <- matrix(c(joint1$V[2,2],VV,VV,joint2$V[2,2]),2,2)
    ## print(V)
    ## print(Wmat1)
    B <- c(joint1$coef[-1],joint2$coef[-1])
        dalt1 <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat1), log=TRUE) 
        dalt2 <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat2), log=TRUE) 
        dalt4 <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat4), log=TRUE) 
        dnull <- dmvnorm(x=B,mean=rep(0,2), sigma=(V + Wmat0), log=TRUE)
        ## dalt1 <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat1, log=TRUE) 
        ## dalt2 <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat2, log=TRUE) 
        ## dalt4 <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat4, log=TRUE) 
        ## dnull <- dmvnorm(x=B,mean=rep(0,4), sigma=V + Wmat0, log=TRUE)
    list(dalt1-dnull, dalt2-dnull, dalt4-dnull)
    
}

.lbf.h3 <- function(bA1,bB1,bA2,bB2,
                   vA1,vB1,vA2,vB2,
                   n00,n01,n02,n1,n2,
                   fA,fB,rho,W=0.04) {

    ## reconstruct joint models - list of coef(a, b1, b2) and V=vcov(coef)
        joint1 <- .cfv2(bA1,bB1,vA1,vB1,n0=n00+n01,n1=n1,fA,fB,rho)
        joint2 <- .cfv2(bA2,bB2,vA2,vB2,n0=n00+n02,n1=n2,fA,fB,rho)
    
## ## est genotype counts in controls
    ## bvec <- 0:2
    G0 <- estgeno.2.ctl(fA,fB,rho)
    ## G1 <- estgeno.2.cse(joint1$coef[-1], G0)
    ##     G2 <- estgeno.2.cse(joint2$coef[-1], G0)

    ## approximate - we don't use this, but could if we needed - if we didn't know propn of shared controls - better to est no. of shared controls
    ## cr <- sqrt(exp(joint1$coef[1] + joint2$coef[2]))*n00/sqrt((n00+n01+n1)*(n00+n02+n2))
    ## VV <- cr * sqrt(joint1$V * joint2$V)
    
    ## cov depends only on shared controls
    cv <- .makecov2(G0*n00,joint1$coef,joint2$coef)
        
    ## cross trait vcov
    VV <- ( joint1$V %*% cv %*% joint2$V )
        
    B <- c(joint1$coef,joint2$coef)
   ##  ## abf
   ##  ## dalt <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V + Wmat, log=TRUE) 
   ##  ## dnull <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V, log=TRUE)

    ## ensure VV is symmetric
    VV[1,2] <- VV[2,1]  <- mean(VV[1,2],VV[2,1])
    VV[3,2] <- VV[2,3]  <- mean(VV[3,2],VV[2,3])
    VV[1,3] <- VV[3,1]  <- mean(VV[1,3],VV[3,1])

    WA <- 10 # prior variance on intercept - BIG because we want to approx a flat prior
        Wmat0 <- diag(c(WA,0,0,WA,0,0))
        Wmat1 <- diag(c(WA,W,0,WA,0,W))
        Wmat2 <- diag(c(WA,0,W,WA,W,0))
        V <- rbind(matrix(c(joint1$V,VV),3,6),matrix(c(VV,joint2$V),3,6))
        dalt <- dmvnorm(x=B[c(2:3,5:6)],mean=rep(0,4), sigma=(V + Wmat1)[c(2:3,5:6),c(2:3,5:6)], log=TRUE) 
        dalt2 <- dmvnorm(x=B[c(2:3,5:6)],mean=rep(0,4), sigma=(V + Wmat2)[c(2:3,5:6),c(2:3,5:6)], log=TRUE) 
        dnull <- dmvnorm(x=B[c(2:3,5:6)],mean=rep(0,4), sigma=(V + Wmat0)[c(2:3,5:6),c(2:3,5:6)], log=TRUE)
        ## dalt <- dmvnorm(x=B,mean=rep(0,6), sigma=V + Wmat1, log=TRUE) 
        ## dalt2 <- dmvnorm(x=B,mean=rep(0,6), sigma=V + Wmat2, log=TRUE) 
        ## dnull <- dmvnorm(x=B,mean=rep(0,6), sigma=V + Wmat0, log=TRUE)
    list(logsum(c(dalt2-dnull, dalt-dnull)))
    
}

.V.h3 <- function(bA1,bB1,bA2,bB2,
                   vA1,vB1,vA2,vB2,
                   n00,n01,n02,n1,n2,
                   fA,fB,rho,W=0.04) {

    ## reconstruct joint models - list of coef(a, b1, b2) and V=vcov(coef)
        joint1 <- .cfv2(bA1,bB1,vA1,vB1,n0=n00+n01,n1=n1,fA,fB,rho)
        joint2 <- .cfv2(bA2,bB2,vA2,vB2,n0=n00+n02,n1=n2,fA,fB,rho)
    
## ## est genotype counts in controls
    ## bvec <- 0:2
    G0 <- estgeno.2.ctl(fA,fB,rho)
    ## G1 <- estgeno.2.cse(joint1$coef[-1], G0)
    ##     G2 <- estgeno.2.cse(joint2$coef[-1], G0)

    ## approximate - we don't use this, but could if we needed - if we didn't know propn of shared controls - better to est no. of shared controls
    ## cr <- sqrt(exp(joint1$coef[1] + joint2$coef[2]))*n00/sqrt((n00+n01+n1)*(n00+n02+n2))
    ## VV <- cr * sqrt(joint1$V * joint2$V)
    
    ## cov depends only on shared controls
    cv <- .makecov2(G0*n00,joint1$coef,joint2$coef)
        
    ## cross trait vcov
    VV <- ( joint1$V %*% cv %*% joint2$V )

    list(bA1=joint1$coef[2],bB1=joint1$coef[3],bA2=joint1$coef[2],bB2=joint2$coef[3],
         V11=joint1$V[2,2], V22=joint1$V[3,3], V12=joint1$V[2,3],
         V33=joint2$V[2,2], V44=joint2$V[3,3], V34=joint1$V[2,3],
         V13=VV[2,2], V24=VV[3,3], V14=VV[2,3])
         
         
   ##  B <- c(joint1$coef,joint2$coef)
   ## ##  ## abf
   ## ##  ## dalt <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V + Wmat, log=TRUE) 
   ## ##  ## dnull <- dmvt(x=B,df=min(n00+n01,n00+n02)+min(n1,n2), sigma=V, log=TRUE)

   ##  ## ensure VV is symmetric
   ##  VV[1,2] <- VV[2,1]  <- mean(VV[1,2],VV[2,1])
   ##  VV[3,2] <- VV[2,3]  <- mean(VV[3,2],VV[2,3])
   ##  VV[1,3] <- VV[3,1]  <- mean(VV[1,3],VV[3,1])

   ##  WA <- 10 # prior variance on intercept - BIG because we want to approx a flat prior
   ##      Wmat0 <- diag(c(WA,0,0,WA,0,0))
   ##      Wmat1 <- diag(c(WA,W,0,WA,0,W))
   ##      Wmat2 <- diag(c(WA,0,W,WA,W,0))
   ##      V <- rbind(matrix(c(joint1$V,VV),3,6),matrix(c(VV,joint2$V),3,6))
   ##      dalt <- dmvnorm(x=B[c(2:3,5:6)],mean=rep(0,4), sigma=(V + Wmat1)[c(2:3,5:6),c(2:3,5:6)], log=TRUE) 
   ##      dalt2 <- dmvnorm(x=B[c(2:3,5:6)],mean=rep(0,4), sigma=(V + Wmat2)[c(2:3,5:6),c(2:3,5:6)], log=TRUE) 
   ##      dnull <- dmvnorm(x=B[c(2:3,5:6)],mean=rep(0,4), sigma=(V + Wmat0)[c(2:3,5:6),c(2:3,5:6)], log=TRUE)
   ##      ## dalt <- dmvnorm(x=B,mean=rep(0,6), sigma=V + Wmat1, log=TRUE) 
   ##      ## dalt2 <- dmvnorm(x=B,mean=rep(0,6), sigma=V + Wmat2, log=TRUE) 
   ##      ## dnull <- dmvnorm(x=B,mean=rep(0,6), sigma=V + Wmat0, log=TRUE)
   ##  list(logsum(c(dalt2-dnull, dalt-dnull)))
    
}


## ##' .. content for \description{} (no empty lines) ..
## ##'
## ##' .. content for \details{} ..
## ##' @title estimate MAF in cases
## ##' @param f0 MAF in controls
## ##' @param b log BF
## ##' @return MAF in cases
## ##' @author Stasia Grinberg
## get_maf1 <- function(f0, b) {
##     f0 * exp(b)/(1 - f0 + f0 * exp(b))
## }


.vmat1 <- function(a,b,N) {
    bvec <- c(0,1,2)
    num <- exp(a + bvec*b)
    denom <- 1 + exp(a + bvec*b)
    phat <- N * num/denom^2
    xmat <- cbind(1, bvec)
    D <- diag(as.vector(phat))
    V <- t(xmat) %*% D %*% xmat
}
.cfv1 <- function(b,n0,n1,f0) {
    G0 <- estgeno.1.ctl(f0)
    G1 <- estgeno.1.cse(G0,b)
    G <-  n0*G0 + n1*G1
    bvec <- c(0,1,2)*b
    fun <- function(a) {
        num <- exp(a + bvec)
        denom <- 1+exp(a + bvec)
        sum(G1*num/denom) - n1/(n1+n0)
    }
    ## microbenchmark(ans <- optimize(fun,c(-100,100)))
    ## ans <- optimize(fun,c(-100,100))
    ##uniroot(fun, c(-10,10))$root
    astar <- uniroot(fun, c(-10, 10),extendInt="yes")$root
    ## ans <- dfsane(par=c(a=0), fn=fun,control=list(trace=FALSE))
    ## astar <- ans$par
    list(coef=c(astar,b),V=solve(.vmat1(astar,b,G)))
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
##' @param dataset2 as above, for dataset 2
##' @param n1 number of cases for dataset 1
##' @param n2 number of cases for dataset 2
##' @param n00 number of shared controls
##' @param n01 number of controls unique to dataset 1
##' @param n02 number of controls unique to dataset 2
##' @param MAF minor allele frequency of the variants, named vector
##' @param LD matrix giving correlation between SNPs.  NB - this is r,
##'     not rsquared!  Expect a mix of positive and negative values.
##'     This should be estimated from reference data.
##' @param p1 prior probability a SNP is associated with trait 1,
##'     default 1e-4
##' @param p2 prior probability a SNP is associated with trait 2,
##'     default 1e-4
##' @param p12 prior probability a SNP is associated with both traits,
##'     default 1e-5
##' @param rlow to save computational time, we approximate the
##'     correlation between SNPs as r=0 where r<rlow
##' @param rhigh to save computational time, we approximate the
##'     correlation between SNPs as r=1 where r>rhigh
##' @param w prior variance on the log OR.  w=0.04 corresponds to a 1%
##'     chance that the absolute log OR will exceed 0.1, or that the
##'     OR will exceed 1.1.  This seems reasonable, as learnt from
##'     existing GWAS literature, but may not be appropriate for very
##'     large or small sample GWAS.  Change this with care!
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
                     rlow=0.05,rhigh=0.95,
                     w=0.04) {
  ## workaround for "no visible binding for global variable/function" Notes re data.table
    . <- function() {}
        abf <- b1 <- b2 <- bA1 <- bA2 <- bB1 <- bB2 <- f0 <- fA <- fB <- h <- isnpA <- isnpB <- lbf3 <- 
  method <- pp <- r <- snp <- snpA <- snpB <- v1 <- v2 <- vA1 <- vA2 <- value <- variable <- vB1 <- vB2 <- NA

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
    ## df[,f1:=get_maf1(f0, b1)]
    ## df[,f2:=get_maf1(f0, b2)]
    ## ## intercept
    ## uncommented from here
    for(i in 1:nrow(df)) {
    ## commented from here
        ## df[i, c("V11","V22"):=.V.h124(b1,b2,v1,v2,n00,n01,n02,n1,n2,f0)]
        df[i, c("lbf1","lbf2","lbf4"):=.lbf.h124(b1,b2,v1,v2,n00,n01,n02,n1,n2,f0)]
    ## to here
        ## df[i,a1:=univariate_a(b1,n00+n01,n1,f0,f1)]
        ## df[i,a2:=univariate_a(b2,n00+n02,n2,f0,f2)]
        ## df[i,d0:=dmvnorm(c(b1,b2),sigma=matrix(c(V11,V12,V12,V22),2,2),log=TRUE)]
        ## df[i,lbf1:=dmvnorm(c(b1,b2),sigma=matrix(c(V11+w,V12,V12,V22),2,2),log=TRUE) - d0]
        ## df[i,lbf2:=dmvnorm(c(b1,b2),sigma=matrix(c(V11,V12,V12,V22+w),2,2),log=TRUE) - d0]
        ## df[i,lbf4:=dmvnorm(c(b1,b2),sigma=matrix(c(V11+w,V12,V12,V22+w),2,2),log=TRUE) - d0]
        ## df[,lbf1:=lhood.norm2(b1,b2,V11+w,V22,V12) - d0]
        ## df[,lbf2:=lhood.norm2(b1,b2,V11,V22+w,V12) - d0]
        ## df[,lbf4:=lhood.norm2(b1,b2,V11+w,V22+w,V12) - d0]
    }
    ## df[,d0:=lhood.norm2(b1,b2,V11,V22,V12)]
    ## df[,lbf1:=lhood.norm2(b1,b2,V11+w,V22,V12) - d0]
    ## df[,lbf2:=lhood.norm2(b1,b2,V11,V22+w,V12) - d0]
    ## df[,lbf4:=lhood.norm2(b1,b2,V11+w,V22+w,V12) - d0]

    ## ## H1, H2, H4 all two beta models
    ## ## more accurate
    ## ## for(i in 1:nrow(df)) {
    ##         df[i, c("v12"):=approx.h124(b1,b2,n00,n01,n02,n1,n2,f0,f1,f2)]
    ## }
    ## ## df[,v12:=(1+exp(a1)) * (1+exp(a2)) * n00 / ( 2 * f0 * (1-f0) * (n00+n01+n1) * (n00+n02+n2))]
    ## df[, erho:=v12/sqrt(v1*v2)]
    ## df[erho > 1, erho:=0.99]
    ## df[erho < -1, erho:=-0.99]
    ## ## faster
    ## ## df[,erho:=n00*sqrt(n1*n2)/((n00+n01+n1)*(n00+n02+n2))]
    
    ## ## lbf for h1,h2,h4
    ## df[, lbf1:=bf2(b1,b2,v1,v2,erho,w2=0)]
    ## df[, lbf2:=bf2(b1,b2,v1,v2,erho,w1=0)]
    ## df[, lbf4:=bf2(b1,b2,v1,v2,erho)]
    ## to here
    
    ## for h3, we need pairs of snps
    ## df3 <- as.data.table(as.data.frame(expand_half_grid(snps)))
    df3 <- as.data.table(as.data.frame(coloc:::expand_r(LD)))
    setnames(df3,c("isnpA","isnpB","r"))
    df3[,snpA:=rownames(LD)[isnpA]]
    df3[,snpB:=rownames(LD)[isnpB]]
    ## setnames(df3,c("snpA","snpB"))
    ## df3[,isnpA:=match(as.character(snpA),rownames(LD))]
    ## df3[,isnpB:=match(as.character(snpB),rownames(LD))]
    ## ## df3 <- df3[isnpA<isnpB,]
    ## df3[,r:=LD[cbind(isnpA,isnpB)]]
    
    ## special case: r==1. This won't allow a model to be fitted at all.  Instead, add in the bf from h4 for one of the two SNPs
    df3.r1 <- df3[abs(r)>=rhigh,] # r=0.95 -> r2=0.9
    mA <- match(df3.r1$snpA, df$snp)
    mB <- match(df3.r1$snpB, df$snp)
    df3.r1[, lbf3:=logsum2(df$lbf4[mA],df$lbf4[mB])]
    
    ## special case: r=0. Can partition as [ H1 for SNP A ] * [ H2 for SNP B ]
    df3.r0 <- df3[abs(r)<=rlow,] # r=0.1 -> r2=0.01
    mA <- match(df3.r0$snpA, df$snp)
    mB <- match(df3.r0$snpB, df$snp)
    df3.r0[, lbf3:=logsum2(df$lbf1[mA]+df$lbf2[mB],
                           df$lbf2[mA]+df$lbf1[mB])]
    
    ## what's left requires 4 dim Gaussians
    df3 <- df3[abs(r)<rhigh & abs(r)>rlow,]
    ## df3 <- df3[abs(r)<rhigh,]
    df3 <- merge(df3,df[,.(snp,f0,b1,v1)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 1, snp A
    setnames(df3,c("f0","b1","v1"),c("fA","bA1","vA1"))
    df3 <- merge(df3,df[,.(snp,f0,b2,v2)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 2, snp B
    setnames(df3,c("f0","b2","v2"),c("fB","bB2","vB2"))
    df3 <- merge(df3, df[,.(snp,b1,v1)],
                 by.x="snpB",by.y="snp",all.x=TRUE) # trait 1, snp B
    setnames(df3,c("b1","v1"),c("bB1","vB1"))
    df3 <- merge(df3, df[,.(snp,b2,v2)],
                 by.x="snpA",by.y="snp",all.x=TRUE) # trait 2, snp A
    setnames(df3,c("b2","v2"),c("bA2","vA2"))
    
    ## set all fA, fB < 0.5 to avoid dividing by 0
    asw <- df3$fA>0.5
    bsw <- df3$fB>0.5
    df3[asw, c("bA1","bA2","fA","r") := list(-bA1, -bA2, 1-fA, -r)]
    df3[bsw, c("bB1","bB2","fB","r") := list(-bB1, -bB2, 1-fB, -r)]
    
    ## reconstruct coefficients from bivariate binomial logistics 
    for(i in 1:nrow(df3)) 
      df3[i,c("lbf3"):=.lbf.h3(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
                                vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
                                n00,n01,n02,n1,n2,
                              fA,fB,rho=r,W=0.04)]
   ##  df3[i,c("bA1","bB1","bA2","bB2",
   ##              "V11","V22","V12",
   ##              "V33","V44","V34",
   ##              "V13","V24","V14"):=.V.h3(bA1=bA1,bB1=bB1,bA2=bA2,bB2=bB2,
   ##                                        vA1=vA1,vB1=vB1,vA2=vA2,vB2=vB2,
   ##                                        n00,n01,n02,n1,n2,
   ##                                        fA,fB,rho=r,W=0.04)]
   ##  df3[,d0:=lhood.norm4(bA1,bB1,bA2,bB2,
   ##                       V11, V22, V12, V33, V44, V34, V13, V24, V14)]
   ## df3[,d3a:=lhood.norm4(bA1,bB1,bA2,bB2,
   ##                       V11 + w, V22, V12, V33, V44 + w, V34, V13, V24, V14)]
   ## df3[,d3b:=lhood.norm4(bA1,bB1,bA2,bB2,
   ##                       V11, V22 + w, V12, V33 + w, V44, V34, V13, V24, V14)]
   ##  df3[,lbf3:=logsum2( d3a-d0, d3b-d0 ) ]
    ## attach(df3[i,]);
    ## rho=r
    
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



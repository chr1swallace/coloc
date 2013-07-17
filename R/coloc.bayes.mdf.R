coloc.bayes <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(c(1,1,1,1,1)),r2.trim=0.99,pp.thr=0.005,nmodsnp=1,quiet=TRUE) {
    #we consider all models which contain at most nmodsnp snps for each trait
    snps <- unique(snps)
    n.orig <- length(snps)
    if(n.orig<2)
        return(1)
    prep <- prepare.df(df1, snps, r2.trim=r2.trim, dataset=1, quiet=quiet)
    df1 <- prep$df
    snps <- prep$snps
    
    if(!quiet)
        cat("Dropped",n.orig - length(snps),"of",n.orig,"SNPs due to LD: r2 >",r2.trim,"\n",length(snps),"SNPs remain.\n")
    
    ## remove any completely predictive SNPs
    f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    while(any(is.na(coefficients(lm1)))) {
        drop <- which(is.na(coefficients(lm1))[-1])
        if(!quiet)
            cat("Dropping",length(drop),"inestimable SNPs (most likely due to colinearity):\n",drop,"\n")
        snps <- snps[-drop]
        f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
        capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    }
    n.clean <- length(snps)
    #remove SNPs with low posterior probabilities in the individual models.
    f1 <- as.formula(paste("Y ~", paste(snps,collapse="+")))
    df.trait1<-df1[which(df1[,1]!=2),c("Y",snps)]
    df.trait2<-df1[which(df1[,1]!=1),c("Y",snps)]
    df.trait2[,1]<-df.trait2[,1]/2
    modelsep<-diag(n.clean)
    mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
    mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1], error="binomial", link="logit",models=modelsep)
    pp.trait1<-mod.trait1$bf$postprob[,2]
    pp.trait2<-mod.trait2$bf$postprob[,2]
    whichsnps<-union(which(pp.trait1>pp.thr),which(pp.trait2>pp.thr))
    cat("We consider ",length(whichsnps), " SNPs in the final analysis. \n")
    snps<-snps[whichsnps]
    
    n.clean <- length(snps)
    #covert to a binomial model so we can run glib
    f1 <- as.formula(paste("Y ~ 1 | ", paste(snps,collapse="+")))
    binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1 (we do not want it in our model matrix)
    binX<-binX[,-1]
    #extract the new reponse
    binY<-binmod[,"Y.star"]
    models<-makebinmod(n.clean,nmodsnp)
    category<-apply(models,1,whichcat)
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    twologB10=mods1$bf$twologB10[,2]
    logbf<-rep(0,5)
    logbf[1]<-logsum(0.5*twologB10[which(category==0)])
    logbf[2]<-logsum(0.5*twologB10[which(category==1)])
    logbf[3]<-logsum(0.5*twologB10[which(category==2)])
    logbf[4]<-logsum(0.5*twologB10[which(category==3)])
    logbf[5]<-logsum(0.5*twologB10[which(category==4)])
    postprobs<-vector("list",length(priors))
    for (ii in 1:length(priors)){
        prior<-priors[[ii]]
        postlogbf<-logbf+log(prior)
        postprob<-(exp(postlogbf))
        postprob<-postprob/sum(postprob)
        cat("for prior: ", prior/sum(prior), "\n")
        cat("we have posterior: ",postprob, "\n")
        cat("PP for shared variant is: ", 100*postprob[5], "% \n")
        cat("--- \n")
        postprobs[[ii]]=postprob
    }
    return(postprobs)
}


coloc.bayes.tag <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(c(1,1,1,1,1)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE) {
    #we consider all models which contain at most 1 snps for each trait, using tagging
    snps <- unique(snps)
    n.orig <- length(snps)
    
    #generate tags using r2.trim
    X<-new("SnpMatrix",as.matrix(df1[,-1]))
    tagkey <- tag(X,tag.threshold=r2.trim)
    tags<-unique(tagkey)
    n.clean=length(tags)
    tagsize<-rep(0,n.clean)
    for (ii in 1:n.clean){
        tagsize[ii]<-length(which(tagkey==tags[ii]))
    }
    
    ## remove any completely predictive tags
    f1 <- as.formula(paste("Y ~", paste(tags,collapse="+")))
    capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    while(any(is.na(coefficients(lm1)))) {
        drop <- which(is.na(coefficients(lm1))[-1])
        if(!quiet)
            cat("Dropping",length(drop),"inestimable tagss (most likely due to colinearity):\n",drop,"\n")
        tags <- tags[-drop]
        f1 <- as.formula(paste("Y ~", paste(tags,collapse="+")))
        capture.output(lm1 <- multinom(f1, data=df1,maxit=1000))
    }
    n.clean <- length(tags)
    
    
    #remove tagss with low posterior probabilities in the individual models.
    f1 <- as.formula(paste("Y ~", paste(tags,collapse="+")))
    df.trait1<-df1[which(df1[,1]!=2),c("Y",tags)]
    df.trait2<-df1[which(df1[,1]!=1),c("Y",tags)]
    df.trait2[,1]<-df.trait2[,1]/2
    modelsep<-diag(n.clean)
    mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
    mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1], error="binomial", link="logit",models=modelsep)
    pp.trait1<-mod.trait1$bf$postprob[,2]
    pp.trait2<-mod.trait2$bf$postprob[,2]
    whichtags<-union(which(pp.trait1>pp.thr),which(pp.trait2>pp.thr))
    cat("We consider ",length(whichtags), " tags in the final analysis. \n")
    tags<-tags[whichtags]
    tagsize<-tagsize[which(unique(tagkey) %in% tags)]
    n.clean <- length(tags)
    #covert to a binomial model so we can run glib
    f1 <- as.formula(paste("Y ~ 1 | ", paste(tags,collapse="+")))
    binmod<-mlogit2logit(f1,data=df1,choices=0:2,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1 (we do not want it in our model matrix)
    binX<-binX[,-1]
    #extract the new reponse
    binY<-binmod[,"Y.star"]
    models<-makebinmod(n.clean,1)
    category<-apply(models,1,whichcat)
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    logB10=0.5*mods1$bf$twologB10[,2]
    logbf <- numeric(5)
    if (n.clean>1){
        n1 <- as.vector( models[,1:n.clean] %*% tagsize )
        n2 <- as.vector( models[,(n.clean+1):(2*n.clean)] %*% tagsize )
    }else{
        n1<-rep(tagsize,length(category))
        n2<-rep(tagsize,length(category))
    }
    logbf[1]<-logsum(logB10[which(category==0)])
    wh1 <- which(category==1)
    logbf[2]<-wlogsum(logB10[wh1], n1[wh1])
    wh2 <- which(category==2)
    logbf[3]<-wlogsum(logB10[wh2], n2[wh2])
    wh3 <- which(category==3)
    wh4 <- which(category==4)
    logbf[4]<-wlogsum(c(logB10[wh3],logB10[wh4]), c(n1[wh3] * n2[wh3],n1[wh4]*(n1[wh4]-1)))    
    logbf[5]<-wlogsum(logB10[wh4], n1[wh4]) ## assumes n1==n2 for cat 4, can't see why this wouldn't be true, but untested
    postprobs<-vector("list",length(priors))
    for (ii in 1:length(priors)){
        prior<-priors[[ii]]
        postlogbf<-logbf-max(logbf)+log(prior)
        postprob<-(exp(postlogbf))
        postprob<-postprob/sum(postprob)
        cat("for prior: ", prior/sum(prior), "\n")
        cat("we have posterior: ",postprob, "\n")
        cat("PP for shared variant is: ", 100*postprob[5], "% \n")
        cat("--- \n")
        postprobs[[ii]]=postprob
    }
    return(postprobs)
}


makebinmod<-function(p,m){
    #returns a model matrix for the binomial equivalent model
    #p=number of snps present
    #m=max number of snps in each model
    if (m>p) {m=p}
    snplist<-vector("list",2*m)
    for (ii in 1:(2*m)){snplist[[ii]]=0:p}
    whichsnps<-expand.grid(snplist)
    if (m >1){
        whichsnps<-unique(t(rbind(apply(whichsnps[,1:m],1,sort),apply(whichsnps[,(m+1):(2*m)],1,sort))))
    }
    whichsnps.t1<-as.matrix(whichsnps[,1:m])
    whichsnps.t2<-as.matrix(whichsnps[,(1+m):(2*m)])
    nummod<-nrow(whichsnps)
    models.t1<-matrix(0,nummod,p)
    models.t2<-matrix(0,nummod,p)
    for (ii in 1:nummod){
        for (jj in 1:m){
            if (whichsnps.t1[ii,jj]>0){
                models.t1[ii,whichsnps.t1[ii,jj]]=1
            }
            if (whichsnps.t2[ii,jj]>0){
                models.t2[ii,whichsnps.t2[ii,jj]]=1
            }
        }
    }
    models<-cbind(models.t1,models.t2,rep(1,nummod))
    if (m>1){
        models<-unique(models)
    }
    return(models)
}
whichcat<-function(line){
    #puts the model given in line into one of the five categories
    #assumes m=1
    p<-(length(line)-1)/2
    t1<-line[1:p]
    t2<-line[(p+1):(2*p)]
    if (sum(t1)==0 & sum(t2)==0){
        return(0)
    }else if (sum(t2)==0){
        return(1)
    }else if (sum(t1)==0){
        return(2)
    }else if (sum(abs(t1-t2))==0){
        return(4)
    }else{
        return(3)
    }
}

wlogsum <- function(x, w=NULL) {
    if (length(x)==1){
        if(is.null(w)) {
            return(x)
        }else{
            return(x*w)
        }
    }
    my.max <- max(x) ##take out the maximum value in log form
    if (my.max == -Inf){return(-Inf)}
    if(is.null(w)) {
        my.max + log(sum(exp(x - my.max )))
    } else {
        my.max + log(sum(exp(x - my.max )*w))
    }
}

myr2 <- function(X) {
  r2 <- ld(X,
           depth=ncol(X)-1,
           symmetric=TRUE,
           stat="R.squared")
if(any(is.na(r2))) {
    r2.na <- as(is.na(r2),"matrix")
    use <- rowSums(r2.na)>0
## work around for r2=NA bug.
    r2.cor <- as(cor(as(X[,use,drop=FALSE],"numeric"), use="pairwise.complete.obs")^2,"Matrix")
    r2[ which(r2.na) ] <- r2.cor[ which(r2.na[use,use]) ]
}
  diag(r2) <- 1
return(r2)
}
##' Derive tag SNPs for a SnpMatrix object using heirarchical clustering
##'
##' Uses complete linkage and the \code{\link{hclust}} function to define clusters, then cuts the tree at 1-tag.threshold
##' @title tag
##' @param snps colnames of the SnpMatrix object to be used
##' @param tag.threshold threshold to cut tree, default=0.99
##' @param samples optional, subset of samples to use
##' @return character vector, names are \code{snps}, values are the tag for each SNP
##' @author Chris Wallace
##' @export
tag <- function(X,tag.threshold=0.99, snps=NULL, samples=NULL) {
if(!is.null(snps) || !is.null(samples))
    X <- X[samples,snps]
  r2 <- myr2(X)
   D <- as.dist(1-r2)
   hc <- hclust(D, method="complete")
   clusters <- cutree(hc, h=1-tag.threshold)
   snps.use <- names(clusters)[!duplicated(clusters)]
   r2.use <- r2[snps.use, colnames(X), drop=FALSE]
   tags <- rownames(r2.use)[apply(r2.use,2,which.max)]
   names(tags) <- colnames(r2.use)
return(tags)
}

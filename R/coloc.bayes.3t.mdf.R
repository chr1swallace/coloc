coloc.bayes.3t <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(rep(1,15)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE) {
    #we consider all models which contain at most 1 snp for each of the three traits
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
    df.trait1<-df1[which(df1[,1] %in% c(0,1)),c("Y",snps)]
    df.trait2<-df1[which(df1[,1] %in% c(0,2)),c("Y",snps)]
    df.trait3<-df1[which(df1[,1] %in% c(0,3)),c("Y",snps)]
    df.trait2[,1]<-df.trait2[,1]/2
    df.trait3[,1]<-df.trait3[,1]/3
    modelsep<-diag(n.clean)
    mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
    mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1], error="binomial", link="logit",models=modelsep)
    mod.trait3<-glib(x=df.trait3[,-1], y=df.trait3[,1], error="binomial", link="logit",models=modelsep)
    pp.trait1<-mod.trait1$bf$postprob[,2]
    pp.trait2<-mod.trait2$bf$postprob[,2]
    pp.trait3<-mod.trait3$bf$postprob[,2]
    whichsnps<-union(union(which(pp.trait1>pp.thr),which(pp.trait2>pp.thr)),which(pp.trait3>pp.thr))
    cat("We consider ",length(whichsnps), " SNPs in the final analysis. \n")
    snps<-snps[whichsnps]
    
    n.clean <- length(snps)
    #covert to a binomial model so we can run glib
    f1 <- as.formula(paste("Y ~ 1 | ", paste(snps,collapse="+")))
    binmod<-mlogit2logit(f1,data=df1,choices=0:3,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1 (we do not want it in our model matrix)
    binX<-binX[,-1]
    #extract the new reponse
    binY<-binmod[,"Y.star"]
    models<-makebinmod.3t(n.clean)
    category<-apply(models,1,whichcat.3t)
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    twologB10=mods1$bf$twologB10[,2]
    logbf<-rep(0,15)
    for (ii in 1:15){
        logbf[ii]<-my.logsum(0.5*twologB10[which(category==(ii-1))])  
    }
    postprobs<-vector("list",length(priors))
    for (ii in 1:length(priors)){
        prior<-priors[[ii]]
        postlogbf<-logbf+log(prior)-max(logbf)
        postprob<-(exp(postlogbf))
        postprob<-postprob/sum(postprob)
        cat("for prior: ", prior/sum(prior), "\n")
        cat("we have posterior: ",postprob, "\n")
        cat("PP for shared variant is: ", 100*postprob[15], "% \n")
        cat("--- \n")
        postprobs[[ii]]=postprob
    }
    return(postprobs)
}



##########################

coloc.bayes.3t.tag <- function(df1,snps=setdiff(colnames(df1),response),response="Y",priors=list(rep(1,15)),r2.trim=0.99,pp.thr=0.005,quiet=TRUE) {
    #we consider all models which contain at most 1 snp for each of the three traits
    snps <- unique(snps)
    n.orig <- length(snps)
    
    #generate tags using r2.trim
    invisible(capture.output(X<-new("SnpMatrix",as.matrix(df1[,-1]))))
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
 
    #remove SNPs with low posterior probabilities in the individual models.
    f1 <- as.formula(paste("Y ~", paste(tags,collapse="+")))
    df.trait1<-df1[which(df1[,1] %in% c(0,1)),c("Y",tags)]
    df.trait2<-df1[which(df1[,1] %in% c(0,2)),c("Y",tags)]
    df.trait3<-df1[which(df1[,1] %in% c(0,3)),c("Y",tags)]
    df.trait2[,1]<-df.trait2[,1]/2
    df.trait3[,1]<-df.trait3[,1]/3
    modelsep<-diag(n.clean)
    mod.trait1<-glib(x=df.trait1[,-1], y=df.trait1[,1], error="binomial", link="logit",models=modelsep)
    mod.trait2<-glib(x=df.trait2[,-1], y=df.trait2[,1], error="binomial", link="logit",models=modelsep)
    mod.trait3<-glib(x=df.trait3[,-1], y=df.trait3[,1], error="binomial", link="logit",models=modelsep)
    pp.trait1<-mod.trait1$bf$postprob[,2]
    pp.trait2<-mod.trait2$bf$postprob[,2]
    pp.trait3<-mod.trait3$bf$postprob[,2]
    whichtags<-union(union(which(pp.trait1>pp.thr),which(pp.trait2>pp.thr)),which(pp.trait3>pp.thr))
    cat("We consider ",length(whichtags), " tags in the final analysis. \n")
    tags<-tags[whichtags]
    
    n.clean <- length(tags)
    #covert to a binomial model so we can run glib
    f1 <- as.formula(paste("Y ~ 1 | ", paste(tags,collapse="+")))
    binmod<-mlogit2logit(f1,data=df1,choices=0:3,base.choice = 1)$data
    index<-grep("z",colnames(binmod))
    binX<-binmod[,index]
    #remove z_1 (we do not want it in our model matrix)
    binX<-binX[,-1]
    #extract the new reponse
    binY<-binmod[,"Y.star"]
    models<-makebinmod.3t(n.clean)
    category<-apply(models,1,whichcat.3t)
    #run glib
    mods1 <- glib(x=binX, y=binY, error="binomial", link="logit",models=models)
    #-2 log the bf with a flat prior
    logB10=0.5*mods1$bf$twologB10[,2]
    logbf <- numeric(15)
    tagsize<-tagsize[which(unique(tagkey) %in% tags)]
    if (n.clean>1){
        n1 <- as.vector( models[,1:n.clean] %*% tagsize )
        n2 <- as.vector( models[,(n.clean+1):(2*n.clean)] %*% tagsize )
        n3 <- as.vector( models[,(2*n.clean+2):(3*n.clean+1)] %*% tagsize )
    }else{
        n1<-rep(tagsize,length(category))
        n2<-rep(tagsize,length(category))
        n3<-rep(tagsize,length(category))
    }
    wh0<-which(category==0)
    wh1<-which(category==1)
    wh2<-which(category==2)
    wh3<-which(category==3)
    wh4<-which(category==4)
    wh5<-which(category==5)
    wh6<-which(category==6)
    wh7<-which(category==7)
    wh8<-which(category==8)
    wh9<-which(category==9)
    wh10<-which(category==10)
    wh11<-which(category==11)
    wh12<-which(category==12)
    wh13<-which(category==13)
    wh14<-which(category==14)
    if (n.clean >1){
        tagvec1<-rep(tagsize,times=n.clean)[-seq(from=1,by=n.clean+1,length=n.clean)]
        tagvec2<-rep(tagsize,each=(n.clean-1))
    }else{
        tagvec1<-tagsize
        tagvec2<-tagsize
    }
    #no association
    logbf[1]<-my.logsum(logB10[which(category==0)])
    #association with one trait
    logbf[2]<-wlogsum(logB10[wh1], n1[wh1])
    logbf[3]<-wlogsum(logB10[wh2], n2[wh2])
    logbf[4]<-wlogsum(logB10[wh3], n3[wh3])
    #association with two traits, no colocalisation
    logbf[5]<-wlogsum(c(logB10[wh4],logB10[wh7]),c(n1[wh4]*n2[wh4],n1[wh7]*(n1[wh7]-1)))
    logbf[6]<-wlogsum(c(logB10[wh5],logB10[wh8]),c(n1[wh5]*n3[wh5],n1[wh8]*(n1[wh8]-1)))
    logbf[7]<-wlogsum(c(logB10[wh6],logB10[wh9]),c(n2[wh6]*n3[wh6],n2[wh9]*(n2[wh9]-1)))
    #association with two traits, with colocalisation
    logbf[8]<-wlogsum(logB10[wh7], n1[wh7])
    logbf[9]<-wlogsum(logB10[wh8], n1[wh8])
    logbf[10]<-wlogsum(logB10[wh9], n2[wh9])
    #all associated, non colocalised
    logbf[11]<-wlogsum(c(logB10[wh10],logB10[wh11],logB10[wh12],logB10[wh13],logB10[wh14]),c(n1[wh10]*n2[wh10]*n3[wh10],n1[wh11]*(n1[wh11]-1)*n3[wh11],n1[wh12]*(n1[wh12]-1)*n2[wh12],n2[wh13]*(n2[wh13]-1)*n1[wh13],n1[wh14]*(n1[wh14]-1)*(n1[wh14]-2)))
    #all associated, two colocalised
    logbf[12]<-wlogsum(c(logB10[wh11],logB10[wh14]),c(n1[wh11]*(n1[wh11]-1)*n3[wh11],n1[wh14]*(n1[wh14]-1)))
    logbf[13]<-wlogsum(c(logB10[wh12],logB10[wh14]),c(n1[wh12]*(n1[wh12]-1)*n2[wh12],n1[wh14]*(n1[wh14]-1)))
    logbf[14]<-wlogsum(c(logB10[wh13],logB10[wh14]),c(n2[wh13]*(n2[wh13]-1)*n1[wh13],n1[wh14]*(n1[wh14]-1)))
    #all colocalised
    logbf[15]<-wlogsum(logB10[wh14], n1[wh14])
    #
    postprobs<-vector("list",length(priors))
    for (ii in 1:length(priors)){
        prior<-priors[[ii]]
        postlogbf<-logbf-max(logbf)+log(prior)  #normalise the logbf
        postprob<-(exp(postlogbf))
        postprob<-postprob/sum(postprob)
        cat("for prior: ", prior/sum(prior), "\n")
        cat("we have posterior: ",postprob, "\n")
        cat("PP for shared variant is: ", 100*postprob[15], "% \n")
        cat("--- \n")
        postprobs[[ii]]=postprob
    }
    return(postprobs)
}



#########################

makebinmod.3t<-function(p){
    #returns a model matrix for the binomial equivalent model
    #p=number of snps present
    snplist<-vector("list",3)
    for (ii in 1:3){snplist[[ii]]=0:p}
    whichsnps<-expand.grid(snplist)
    whichsnps.t1<-whichsnps[,1]
    whichsnps.t2<-whichsnps[,2]
    whichsnps.t3<-whichsnps[,3]
    nummod<-nrow(whichsnps)
    models.t1<-matrix(0,nummod,p)
    models.t2<-matrix(0,nummod,p)
    models.t3<-matrix(0,nummod,p)
    for (ii in 1:nummod){
        if (whichsnps.t1[ii]>0){
            models.t1[ii,whichsnps.t1[ii]]=1
        }
        if (whichsnps.t2[ii]>0){
            models.t2[ii,whichsnps.t2[ii]]=1
        }
        if (whichsnps.t3[ii]>0){
            models.t3[ii,whichsnps.t3[ii]]=1
        }
    }
    models<-cbind(models.t1,models.t2,rep(1,nummod),models.t3,rep(1,nummod))
    return(models)
}
mywhich <- function(tline) {
    wh <- which(tline==1)
    if(length(wh)==0){
        wh <- 0
    }
    return(wh)
}
 
divide.hyps <- function(t1,t2,t3) {
    if (t1!=t2 & t2!=t3 & t3!=t1){
        return(10)
    }
    if (t1==t2 & t1==t3){
        return(14)
    }
    ## now know there must be exactly one matching pair
    if (t1==t2){
        return(11)
    }
    if (t1==t3){
        return(12)
    }
    if (t2==t3){
        return(13)
    }
}
 
whichcat.3t<-function(line){
#puts the model given in line into one of the twelve categories
#assumes m=1
    p<-(length(line)-2)/3
    tline1<-line[1:p]
    tline2<-line[(p+1):(2*p)]
    tline3<-line[((2*p)+2):((3*p)+1)]
    t1 <- mywhich(tline1)
    t2 <- mywhich(tline2)
    t3 <- mywhich(tline3)
    ## ttt indicates whether t. is 0 or >0 only
    ttt <- paste(min(t1,1), min(t2,1), min(t3,1),sep="")
    ret <- switch(ttt,
    "000"=0,
    "100"=1,
    "010"=2,
    "001"=3,
    "110"= if(t1==t2) { 7 } else { 4 },
    "101"= if(t1==t3) { 8 } else { 5 },
    "011"= if(t2==t3) { 9 } else { 6 },
    "111"= divide.hyps(t1,t2,t3))
    if(is.null(ret)){
        cat("unrecognised category in whichcat.3t")
    }
    return(ret)
}

my.logsum <- function(x) {
    if (length(x)==1){return(x)}
    my.max <- max(x)                              ##take out the maximum value in log form
    if (my.max == -Inf){return(-Inf)}
    my.res <- my.max + log(sum(exp(x - my.max ))) 
    return(my.res)
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

prior.adjust <- function(summ,newp12,p1=1e-4,p2=1e-4,p12=1e-6) {
    if(is.list(summ) && "summary" %in% names(summ))
        summ <- summ$summary
    if(!identical(names(summ), c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")))
        stop("not a coloc summary vector")
    ## back calculate likelihoods
    f <- function(p12)
        prior.snp2hyp(summ["nsnps"],p12=p12,p1=p1,p2=p2)
    pr1 <- f(newp12)
    pr0 <- matrix(f(p12),nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE)
    ## make everything conformable
    ## if(is.matrix(summ) && nrow(summ)==1) summ <- as.vector(summ)
    ## if(nrow(pr1)==1) pr1 <- as.vector(pr1)
    ## if(nrow(pr1)>1) pr1 <- t(pr1)
    newpp <- matrix(summ[-1],nrow=nrow(pr1),ncol=ncol(pr1),byrow=TRUE) * pr1/pr0 # prop to, not equal to
    newpp/rowSums(newpp)
}


prior.snp2hyp <- function(nsnp,p12=1e-6,p1=1e-4,p2=1e-4) {
    if(any(p12<p1*p2) || any(p12 > p1) || any(p12 > p2))
        return(NULL)
    tmp <- cbind(nsnp * p1,
                 nsnp * p2,
                 nsnp * (nsnp-1) * p1 * p2,
                 nsnp * p12)
    tmp <- cbind(1-rowSums(tmp),tmp)
    ## if(nrow(tmp)==1) {
    ##     tmp <- c(tmp)
    ##     names(tmp) <- paste0("H",0:4)
    ## } else 
        colnames(tmp) <- paste0("H",0:4)
    tmp
}

##' Shows how prior and posterior per-hypothesis probabilities change as a function of p12
##'
##' Function is called mainly for plotting side effect.  It draws two
##' plots, showing how prior and posterior probabilities of each coloc
##' hypothesis change with changing p12.  A decision rule sets the
##' values of the posterior probabilities considered acceptable, and
##' is used to shade in green the region of the plot for which the p12
##' prior would give and acceptable result.  The user is encouraged to
##' consider carefully whether some prior values shown within the
##' green shaded region are sensible before accepting the hypothesis.
##' If no shading is shown, then no priors give rise to an accepted
##' result.
##' @title Prior sensitivity for coloc
##' @param obj output of coloc.detail or coloc.process
##' @param rule a decision rule.  This states what values of posterior
##'     probabilities "pass" some threshold.  This is a string which
##'     will be parsed and evaluated, better explained by examples.
##'     "H4 > 0.5" says post prob of H4 > 0.5 is a pass.
##'     "H4 > 0.9 & H4/H3 > 3" says post prob of H4 must be > 0.9 AND
##'     it must be at least 3 times the post prob of H3."
##' @param npoints the number of points over which to evaluate the prior values for p12, equally spaced on a log scale between p1*p2 and min(p1,p2) - these are logical limits on p12, but not scientifically sensible values.
##' @param doplot draw the plot. set to FALSE if you want to just evaluate the prior and posterior matrices and work with them yourself
##' @return list of prior and posterior matrices returned invisibly
##' @author Chris Wallace
sensitivity.coloc <- function(obj,rule="H4 > 0.5",
                       npoints=100,doplot=TRUE,preserve.par=FALSE,...) {
    rule.init <- rule
    rule <- gsub("(H.)","PP.\\1.abf",rule,perl=TRUE)
    pp <- obj$summary
    p12 <- obj$priors["p12"]
    p1 <- obj$priors["p1"]
    p2 <- obj$priors["p2"]
    check <- function(pp) { with(as.list(pp),eval(parse(text=rule))) }
    pass.init <- check(pp)
    message("Results ",if(check(pp)) { "pass" } else { "fail" }, " decision rule ",rule.init)

    testp12 <- 10^seq(log10(p1*p2),log10(min(p1,p1)),length.out=npoints)
    testH <- prior.snp2hyp(pp["nsnps"],p12=testp12,p1=p1,p2=p2)
    testpp <- as.data.frame(prior.adjust(summ=pp,newp12=testp12,p1=p1,p2=p2,p12=p12))
    colnames(testpp) <- gsub("(H.)","PP.\\1.abf",colnames(testpp),perl=TRUE)
    pass <- check(testpp)
    w <- which(pass)

    if(doplot) {
        H <- as.character(0:4)
        ## palette(c("#ffffffff","#000000ff","#666666ff",viridis::viridis(4,alpha=1)[c(2,4)]))
        palette(c("#ffffffff",viridis::viridis(5,alpha=1)[-1]))
        if(!preserve.par)
            par(mfrow=c(2,1))
       par(mar = c(3, 3, 2, 1), # Dist' from plot to side of page
    mgp = c(2, 0.4, 0), # Dist' plot to label
    las = 1, # Rotate y-axis text
    tck = -.01 # Reduce tick length
    )
        m <- list(testH,as.matrix(testpp))
        ti <- list("Prior probabilities",
                   "Posterior probabilities")
        for(i in 1:2) {
            ym <- if(i==1) { max(m[[i]][,-1]) } else { max(m[[i]]) }
        matplot(testp12,m[[i]],log="x",xlab="p12",ylab="Prob",
                type="b",
                bg = 1:5, # Fill colour
                pch = 21, # Shape: circles that can filed
                col="gray20",
                frame.plot = FALSE, # Remove the frame 
                panel.first = abline(h = seq(0, 1, 0.2), col = "grey80"),
                ylim=c(0,ym))
        title(main=ti[[i]],adj=0)
        ## title(main=paste("Acceptance rule (shaded region):",rule.init))
        ## legend("topleft",pch=rep(21,5),pt.bg=1:5,legend=paste0("H",0:4))
        if(i==1)
            legend("left",inset=c(0.1,0),bg="white",pch=rep(21,5),pt.bg=1:5,pt.cex=2,legend=paste0("H",0:4))
        abline(v=p12,lty="dashed",col="gray")
        text(p12,0.5,"results",srt=90,col="gray40")
        if(any(pass))
            rect(xleft=testp12[min(w)],ybottom=0,
                 xright=testp12[max(w)],ytop=1,
                 col=rgb(0,1,0,alpha=0.1), border="green") 
        ## add text showing rule
        mtext(paste("shaded region:",rule.init),side=3,adj=1)
        }
    }

    invisible(cbind(testpp,p12=testp12,pass=pass))
}




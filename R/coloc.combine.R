coloc.combine <- function(L) {
#  L <- lapply(results,"[[",2)  ##  for testing

if(!is.list(L))
    stop("require a list of coloc objects")
  class.ok <- sapply(L,class) == "coloc"
  if(any(!class.ok)) {
    L <- L[class.ok]
    warning("removing ",sum(!class.ok)," elements of L which are not of class coloc: ", sum(class.ok), " remain.")
  }
if(sum(class.ok)==0)
  return(NULL)
  stats <- sapply(L,slot,"result")
  X2 <- stats["chisquare",,drop=TRUE]
  n <- stats["n",,drop=TRUE]
  ns <- sort(unique(n))
  comb.results <- matrix(NA,7,length(ns)+1)
  rownames(comb.results) <- c("n","m","stat","df1","df2","p.F","p.chisq")
comb.results[,1] <- c(1,sum(!class.ok),rep(NA,5))
i <- 2
  for(k in ns) {
    wh <- which(n==k)
    m <- length(wh)
    if(m<3) {
      comb.results[,i] <- c(k,m,rep(NA,5))
      i <- i+1
      next
    }
    Xbar <- mean(X2)
    Xroot <- sqrt(X2)
    r2 <- var(Xroot) * (1+1/m)
    D2 <- (Xbar/(k-1) - r2 * (m+1) / (m-1)) / (1+r2)
    v2=(k-1)^(-3/m) * (m-1) * (1+1/r2)^2
    pF <- pf(D2, k-1, v2,lower.tail=FALSE)
    pchi <- pchisq(D2,k-1,lower.tail=FALSE)    
    comb.results[,i] <- c(k,m,D2,k-1,v2,pF,pchi)
    i <- i+1
  }
return(comb.results)
}

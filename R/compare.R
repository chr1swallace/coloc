##'  Compare two datasets before running coloc
##'
##' Datasets are merged before running coloc on the basis of snp
##' id. If a snp exists in one dataset and not the other, it will not
##' contribute to the coloc analysis. This may not be a problem if it
##' is a snp without much association signal. But if there is a
##' substantial change that one or more missing snps may be causal for
##' one dataset, this can cause misleading coloc
##' results. compare_datasets will allow you to check how much of the
##' fine mapping posterior in dataset1 is captured by snps in dataset2
##' and vice versa. Ideally the result should be two numbers close to
##' 1, but the function will warn if either drops below thr.
##'
##' Note, thr wasn't chosen on any sensible basis, just
##' intuition. This is a function to give the user information to help
##' decide whether to trust their coloc analysis or not, and different
##' values of thr may apply in different situations. Please take some
##' time to think about the right value for your analysis.
##' 
##' @param dataset1 first coloc dataset
##' @param dataset2 second coloc dataset
##' @param thr minimum fine mapping posterior captured in by snps in
##'     the other dataset
##' @return a vector giving the fraction of the fine mapping posterior
##'     in dataset1 captured by snps in dataset2, and the fraction of
##'     the fine mapping posterior in dataset2 captured by snps in
##'     dataset1
##' @export
##' @author Chris Wallace
compare_datasets=function(dataset1,dataset2,thr=0.8) {
   fm1=finemap.abf(dataset1) 
   fm2=finemap.abf(dataset2)
   m=merge(fm1[,c("snp","SNP.PP")], fm2[,c("snp","SNP.PP")], by="snp", all=TRUE, suffixes=c(".1",".2"))
   oneintwo=with(m[ !is.na(m$SNP.PP.2), ], sum(SNP.PP.1,na.rm=TRUE))
   twoinone=with(m[ !is.na(m$SNP.PP.1), ], sum(SNP.PP.2,na.rm=TRUE))
   message("proportion of posterior in dataset1 covered by snps in dataset2=",oneintwo)
   message("proportion of posterior in dataset2 covered by snps in dataset1=",twoinone)
   if(min(oneintwo,twoinone) < thr) warning("SNPs with substantial posterior probability in at least one dataset not found in the other")
   invisible(c(oneintwo=oneintwo,twoinone=twoinone))
}

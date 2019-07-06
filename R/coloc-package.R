#'Colocalisation tests of two genetic traits
#'
#'Performs the colocalisation tests described in Plagnol et al (2009) and
#'Wallace et al (in preparation) and draws some plots.
#'
#'\code{coloc.test()} tests for colocalisation and returns an object of class
#'\code{coloc}.
#'
#'@name coloc-package
#'@docType package
#'@author Chris Wallace <chris.wallace@@cimr.cam.ac.uk>
#'@references Plagnol et al (2009).  Statistical independence of the
#'colocalized association signals for type 1 diabetes and RPS26 gene expression
#'on chromosome 12q13. Biostatistics 10:327-34.
#'
#'\url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
#'
#'Wallace et al (2013).  Statistical Testing of Shared Genetic Control
#' for Potentially Related Traits. Genetic Epidemiology 37:802-813.
#' 
#'\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4158901/}
#'
#' Giambartolomei et al (2014).  Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics.  PLOS Genet e1004383.
#'
#'\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4022491/}
#' 
#'@keywords package
#' @importFrom BMA glib glib.data.frame
#' @importFrom snpStats single.snp.tests col.summary snp.imputation impute.snps
#' @importFrom methods setMethod setClass
#' @import ggplot2
#' @importFrom graphics abline axis box par
#' @importFrom methods as is new slot
#' @importFrom stats as.dist as.formula coef coefficients complete.cases cor cutree glm integrate lm optimize pchisq pf prcomp qnorm sd var vcov hclust
#' @importFrom utils combn
#' @import data.table
#' @importFrom viridis viridis
#' @importFrom graphics layout legend matplot mtext rect text title
#' @importFrom grDevices colorRampPalette palette rgb
#' @importFrom stats pnorm uniroot
utils::globalVariables(c(".","dfsane","dmvnorm","H0","H1","H2","H3","H4","hit1","hit2","lABF.df1","lABF.df2","lABF.h3","lbf1","lbf2","lbf3","lbf4","nsnps","snp","snp1","snp2","varbeta","z"))
NULL

#'Class \code{"colocPCs"}
#'
#'%% ~~ A concise (1-5 lines) description of what the class is. ~~ Class
#'designed to hold objects returned by function \code{\link{pcs.prepare}} which
#'generates a principal component summary of two genotype matrices in a form
#'suitable for use in the function \code{\link{pcs.model}}.
#'
#'
#'@name colocPCs-class
#'@docType class
#'@section Objects from the Class: Objects can be created by calls to the
#'function \code{\link{pcs.prepare}()}. %% ~~ describe objects here ~~
#'@author Chris Wallace.
#'@seealso \code{\link{pcs.prepare}}, \code{\link{pcs.model}}
#'@references Wallace et al (2012).  Statistical colocalisation of monocyte
#'gene expression and genetic risk variants for type 1 diabetes.  Hum Mol Genet
#'21:2815-2824.  \url{http://europepmc.org/abstract/MED/22403184}
#'
#'Plagnol et al (2009).  Statistical independence of the colocalized
#'association signals for type 1 diabetes and RPS26 gene expression on
#'chromosome 12q13. Biostatistics 10:327-34.
#'\url{http://www.ncbi.nlm.nih.gov/pubmed/19039033}
#'@keywords classes
#'@examples
#'
#'showClass("colocPCs")
#'
NULL




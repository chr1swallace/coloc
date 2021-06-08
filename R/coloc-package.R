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
#' @import ggplot2
#' @importFrom graphics abline axis box par
#' @importFrom methods as is new slot
#' @importFrom stats as.dist as.formula coef coefficients complete.cases cor cutree glm integrate lm optimize pchisq pf prcomp qnorm sd var vcov hclust
#' @importFrom utils combn
#' @import data.table
#' @importFrom viridis viridis
#' @importFrom graphics layout legend matplot mtext rect text title hist
#' @importFrom grDevices colorRampPalette palette rgb
#' @importFrom stats pnorm uniroot
#' @importFrom susieR susie_rss susie_get_cs
utils::globalVariables(c(".","dfsane","dmvnorm","H0","H1","H2","H3","H4","hit1","hit2","lABF.df1","lABF.df2","lABF.h3","lbf1","lbf2","lbf3","lbf4","nsnps","snp","snp1","snp2","varbeta","z"))
NULL

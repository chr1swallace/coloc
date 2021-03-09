#'Colocalisation tests of two genetic traits
#'
#'Performs the colocalisation tests described in Giambartolomei (2012) and
#'follow-up papers
#'
#' Previously contained also tests for proportional colocalisation, but these
#' have now been moved to the separate package coloc.prop
#' \url{https://github.com/chr1swallace/coloc.prop}
#'
#'@name coloc-package
#'@docType package
#'@author Chris Wallace <cew54@cam.ac.uk>
#'@references Giambartolomei et al (2014). Bayesian Test for Colocalisation
#'   between Pairs of Genetic Association Studies Using Summary Statistics. PLOS
#'   Genet e1004383.
#'
#'\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4022491/}
#'
#'Wallace C (2020) Eliciting priors and relaxing the single causal variant
#' assumption in colocalisation analyses. PLOS Genetics 16(4): e1008720.
#'
#'\url{https://doi.org/10.1371/journal.pgen.1008720}
#'
#'Wallace C(2021) A more accurate method for colocalisation analysis allowing
#' for multiple causal variants bioRxiv 2021.02.23.432421;
#'
#'\url{https://doi.org/10.1101/2021.02.23.432421}
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
utils::globalVariables(c(".","dfsane","dmvnorm","H0","H1","H2","H3","H4","hit1","hit2","lABF.df1","lABF.df2","lABF.h3","lbf1","lbf2","lbf3","lbf4","nsnps","snp","snp1","snp2","varbeta","z"))
NULL

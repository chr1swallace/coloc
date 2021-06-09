#'Colocalisation tests of two genetic traits
#'
#'Performs the colocalisation tests described in Plagnol et al (2009) and
#'Wallace et al (2020) and draws some plots.
#'
#'@name coloc-package
#'@docType package
#'@author Chris Wallace <cew54@cam.ac.uk>
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

#'Colocalisation tests of two genetic traits
#'
#'Performs the colocalisation tests described in Plagnol et al (2009) and
#'Wallace et al (2020) and draws some plots.
#'
#'@keywords internal
#' @import ggplot2
#' @importFrom graphics abline axis box par
#' @importFrom methods as is new slot
#' @importFrom stats as.dist as.formula coef coefficients complete.cases cor cutree glm integrate lm optimize pchisq pf prcomp qnorm sd var vcov hclust
#' @importFrom utils combn head
#' @import data.table
#' @importFrom viridis viridis
#' @importFrom graphics layout legend matplot mtext rect text title hist
#' @importFrom grDevices colorRampPalette palette rgb
#' @importFrom stats pnorm uniroot
#' @importFrom susieR susie_rss susie_get_cs
#' @importFrom matrixStats logSumExp
"_PACKAGE"

utils::globalVariables(c(".","dfsane","dmvnorm","H0","H1","H2","H3","H4","hit1","hit2","lABF.df1","lABF.df2","lABF.h3","lbf1","lbf2","lbf3","lbf4","nsnps","snp","snp1","snp2","varbeta","z", "logp"))

##' Package options
##'
##' coloc reads the following options, set with [options()], for
##' example in your `.Rprofile`. They mainly control the functions
##' that download public credible sets ([query_ot()], [ot_fetch()],
##' [query_eqtl_cat()], [ensembl_to_symbol()]).
##'
##' \describe{
##'   \item{`coloc.cache_dir`}{Directory where downloaded files are
##'   kept, so they are reused between R sessions. If unset, files go
##'   to a temporary directory that is deleted when R exits, with a
##'   warning. Default: unset.}
##'   \item{`coloc.max_temp_download`}{Largest download, in bytes,
##'   allowed into that temporary directory. Larger downloads (eg
##'   [ot_fetch()], or several eQTL Catalogue quantification methods at
##'   once) stop with a request to set a cache directory. Default:
##'   `1e9` (1 GB).}
##' }
##'
##' @examples
##' \dontrun{
##' ## in ~/.Rprofile, to keep downloads in one place
##' options(coloc.cache_dir = "~/data/coloc-cache")
##' }
##' @name coloc-options
##' @aliases coloc.cache_dir coloc.max_temp_download
NULL

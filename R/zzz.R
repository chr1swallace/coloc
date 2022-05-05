#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("This is coloc version",packageVersion("coloc")))
}

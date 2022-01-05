.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("This is coloc version",utils::packageVersion("coloc")))
}

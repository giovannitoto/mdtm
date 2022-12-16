#' mdtm : Multiple Descriptors Topic Modeling
#'
#' The mdtm package provides function to estimate topic models.
#'
#' @docType package
#' @name mdtm
#' @useDynLib mdtm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the jungle of Topic Modeling!")
}

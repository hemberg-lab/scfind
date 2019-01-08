#' Essential File so module is loaded
#' 
#' @importFrom Rcpp loadModule
#' @useDynLib scfind
#'
#' @include zzz.R
#' 
#' @export
NULL

loadModule('EliasFanoDB', TRUE)

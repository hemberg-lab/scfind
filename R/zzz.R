#' Essential File so module is loaded
#' @importFrom Rcpp loadModule
#' @useDynLib scfind
#'
#' @export
NULL
loadModule('EliasFanoDB', TRUE)


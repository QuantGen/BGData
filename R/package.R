#' BGData: Memory Mapped Matrices and Data-Structures for Genomic Data for R
#'
#' @docType package
#' @name BGData-package
#' @importFrom methods setClass setClassUnion setMethod
#' @importClassesFrom LinkedMatrix LinkedMatrix
#' @importClassesFrom symDMatrix symDMatrix
#' @importClassesFrom bigmemory big.matrix
#' @importMethodsFrom symDMatrix diag
#' @aliases NULL
NULL


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("The BGData package was supported by the National Institutes of Health (Grant: R01GM101219, R01GM099992).")
    packageStartupMessage()
}

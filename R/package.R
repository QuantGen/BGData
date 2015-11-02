#' BGData: Memory Mapped Matrices and Data-Structures for Genomic Data for R
#' 
#' More info in the wiki: \url{https://github.com/QuantGen/BGData/wiki}
#' 
#' @docType package
#' @name BGData
#' @import methods parallel ff
#' @importFrom bit physical physical<-
#' @importClassesFrom LinkedMatrix LinkedMatrix
#' @aliases NULL
NULL


.onAttach <- function(libname, pkgname) {
    packageStartupMessage("The BGData package was supported by the National Institutes of Health (Grant: R01GM101219, R01GM099992).")
}

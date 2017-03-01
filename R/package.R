#' A Suite of Packages for Analysis of Big Genomic Data.
#'
#' Modern genomic datasets are big (large *n*), high-dimensional (large *p*),
#' and multi-layered. The challenges that need to be addressed are memory
#' requirements and computational demands. Our goal is to develop software that
#' will enable researchers to carry out analyses with big genomic data within
#' the R environment.
#'
#' We have identified several approaches to tackle those challanges within R:
#' - Memory mapping: The data is stored in on the hard drive and users can read
#' in smaller chunks when they are needed.
#' - Linked arrays: For very large datasets a single memory-mapped array may
#' not be enough or convenient. A linked array is an array whose content is
#' distributed over multiple memory-mapped nodes.
#' - Multiple dispatch: Methods are presented to users so that they can treat
#' these arrays pretty much as if they were RAM arrays.
#' - Multi-level parallelism: Exploit multi-core and multi-node computing.
#' - Inputs: Users can create these arrays from standard formats (e.g., PLINK
#' .bed).
#'
#' The BGData package is an umbrella package that comprises several packages:
#' [BEDMatrix][BEDMatrix::BEDMatrix-package],
#' [LinkedMatrix][LinkedMatrix::LinkedMatrix-package], and
#' [symDMatrix][symDMatrix::symDMatrix-package]. It features scalable and
#' efficient computational methods for large genomic datasets such as
#' genome-wide association studies (GWAS) or genomic relationship matrices (G
#' matrix). It also contains a data structure called `BGData` that holds
#' genotypes in the `@@geno` slot, phenotypes in the `@@pheno` slot, and
#' additional information in the `@@map` slot.
#'
#'  @seealso [BEDMatrix::BEDMatrix-package],
#' [LinkedMatrix::LinkedMatrix-package], and [symDMatrix::symDMatrix-package]
#' for an introduction to the respective packages.
#' @docType package
#' @name BGData-package
#' @aliases BGData-package
#' @import methods
#' @importClassesFrom BEDMatrix BEDMatrix
#' @importClassesFrom LinkedMatrix LinkedMatrix
#' @importClassesFrom symDMatrix symDMatrix
#' @importClassesFrom bigmemory big.matrix
NULL

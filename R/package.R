#' A Suite of Packages for Analysis of Big Genomic Data.
#'
#' Modern genomic datasets are big (large *n*), high-dimensional (large *p*),
#' and multi-layered. The challenges that need to be addressed are memory
#' requirements and computational demands. Our goal is to develop software that
#' will enable researchers to carry out analyses with big genomic data within
#' the R environment.
#'
#' We have identified several approaches to tackle those challenges within R:
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
#' @section Memory-mapping:
#' Functions with the `bufferSize` parameter work best with memory-mapped
#' matrices such as [BEDMatrix::BEDMatrix-class] objects. To avoid loading the
#' whole, potentially very large matrix into memory, these functions will load
#' chunks of the memory-mapped matrix into memory and perform the operations on
#' one chunk at a time. The size of the chunks is determined by the
#' `bufferSize` parameter. Care must be taken to not set `bufferSize` too high
#' to avoid memory shortage, particularly when combined with parallel
#' computing.
#'
#' @section Multi-level parallelism:
#' Functions with the `nCores`, `i`, and `j` parameters provide
#' capabilities for both parallel and distributed computing.
#'
#' For parallel computing, `nCores` determines the number of cores the code is
#' run on. Memory usage can be an issue for higher values of `nCores` as R is
#' not particularly memory-efficient. As a rule of thumb, at least around
#' `(nCores * object_size(buffer)) + object_size(result)` MB of total memory
#' will be needed for operations on memory-mapped matrices, not including
#' potential copies of your data that might be created (for example
#' [stats::lsfit()] runs `cbind(1, X)`). `i` and `j` can be used to include or
#' exclude certain rows or columns. Internally, the [parallel::mclapply()]
#' function is used and therefore parallel computing will not work on Windows
#' machines.
#'
#' For distributed computing, `i` and `j` determine the subset of the input
#' matrix that the code runs on. In an HPC environment, this can be used not
#' just to include or exclude certain rows or columns, but also to partition
#' the task among many nodes rather than cores. Scheduler-specific code and
#' code to aggregate the results need to be written by the user. It is
#' recommended to set `nCores` to `1` as nodes are often cheaper than cores.
#'
#' @section Example dataset:
#' The `extdata` folder contains example files that were generated from the
#' 250k SNP and phenotype data in [Atwell et al.
#' (2010)](http://www.nature.com/nature/journal/v465/n7298/full/nature08800.html).
#' Only the first 300 SNPs of chromosome 1, 2, and 3 were included to keep the
#' size of the example dataset small.
#' [PLINK](https://www.cog-genomics.org/plink2) was used to convert the data to
#' [.bed](https://www.cog-genomics.org/plink2/input#bed) and
#' [.raw](https://www.cog-genomics.org/plink2/input#raw) files. `FT10` has been
#' chosen as a phenotype and is provided as an [alternate phenotype
#' file](https://www.cog-genomics.org/plink2/input#pheno). The file is
#' intentionally shuffled to demonstrate that the additional phenotypes are put
#' in the same order as the rest of the phenotypes.
#'
#' @seealso [BEDMatrix::BEDMatrix-package],
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

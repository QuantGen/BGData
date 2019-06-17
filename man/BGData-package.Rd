\docType{package}
\name{BGData-package}
\alias{BGData-package}
\title{A Suite of Packages for Analysis of Big Genomic Data}
\description{
    Modern genomic datasets are big (large \emph{n}), high-dimensional (large
    \emph{p}), and multi-layered. The challenges that need to be addressed are
    memory requirements and computational demands. Our goal is to develop
    software that will enable researchers to carry out analyses with big
    genomic data within the R environment.
}
\details{
    We have identified several approaches to tackle those challenges within R:

    \itemize{
        \item File-backed matrices: The data is stored in on the hard drive and
            users can read in smaller chunks when they are needed.
        \item Linked arrays: For very large datasets a single file-backed array
            may not be enough or convenient. A linked array is an array whose
            content is distributed over multiple file-backed nodes.
        \item Multiple dispatch: Methods are presented to users so that they
            can treat these arrays pretty much as if they were RAM arrays.
        \item Multi-level parallelism: Exploit multi-core and multi-node
            computing.
        \item Inputs: Users can create these arrays from standard formats
            (e.g., PLINK .bed).
    }

    The \code{BGData} package is an umbrella package that comprises several
    packages: \code{BEDMatrix}, \code{LinkedMatrix}, and \code{symDMatrix}. It
    features scalable and efficient computational methods for large genomic
    datasets such as genome-wide association studies (GWAS) or genomic
    relationship matrices (G matrix). It also contains a data structure called
    \code{BGData} that holds genotypes in the \code{@geno} slot, phenotypes in
    the \code{@pheno} slot, and additional information in the \code{@map} slot.
}
\section{Example dataset}{
    The \code{extdata} folder contains example files that were generated from
    the 250k SNP and phenotype data in
    \href{http://www.nature.com/nature/journal/v465/n7298/full/nature08800.html}{Atwell
    et al. (2010)}. Only the first 300 SNPs of chromosome 1, 2, and 3 were
    included to keep the size of the example dataset small.
    \href{https://www.cog-genomics.org/plink2}{PLINK} was used to convert the
    data to \href{https://www.cog-genomics.org/plink2/input#bed}{.bed} and
    \href{https://www.cog-genomics.org/plink2/input#raw}{.raw} files.
    \code{FT10} has been chosen as a phenotype and is provided as an
    \href{https://www.cog-genomics.org/plink2/input#pheno}{alternate phenotype
    file}. The file is intentionally shuffled to demonstrate that the
    additional phenotypes are put in the same order as the rest of the
    phenotypes.
}
\seealso{
    \code{\link[BEDMatrix]{BEDMatrix-package}},
    \code{\link[LinkedMatrix]{LinkedMatrix-package}}, and
    \code{\link[symDMatrix]{symDMatrix-package}} for an introduction to the
    respective packages.

    \code{\link{file-backed-matrices}} for more information on file-backed
    matrices. \code{\link{multi-level-parallelism}} for more information on
    multi-level parallelism.
}

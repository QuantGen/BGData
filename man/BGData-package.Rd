\docType{package}
\name{BGData-package}
\alias{BGData-package}
\title{A Suite of Packages for Analysis of Big Genomic Data}
\description{
Modern genomic datasets are big (large \emph{n}), high-dimensional (large \emph{p}),
and multi-layered. The challenges that need to be addressed are memory
requirements and computational demands. Our goal is to develop software that
will enable researchers to carry out analyses with big genomic data within
the R environment.
}
\details{
We have identified several approaches to tackle those challenges within R:
\itemize{
\item File-backed matrices: The data is stored in on the hard drive and users
can read in smaller chunks when they are needed.
\item Linked arrays: For very large datasets a single file-backed array may not
be enough or convenient. A linked array is an array whose content is
distributed over multiple file-backed nodes.
\item Multiple dispatch: Methods are presented to users so that they can treat
these arrays pretty much as if they were RAM arrays.
\item Multi-level parallelism: Exploit multi-core and multi-node computing.
\item Inputs: Users can create these arrays from standard formats (e.g., PLINK
.bed).
}

The BGData package is an umbrella package that comprises several packages:
\link[BEDMatrix:BEDMatrix-package]{BEDMatrix},
\link[LinkedMatrix:LinkedMatrix-package]{LinkedMatrix}, and
\link[symDMatrix:symDMatrix-package]{symDMatrix}. It features scalable and
efficient computational methods for large genomic datasets such as
genome-wide association studies (GWAS) or genomic relationship matrices (G
matrix). It also contains a data structure called \code{BGData} that holds
genotypes in the \code{@geno} slot, phenotypes in the \code{@pheno} slot, and
additional information in the \code{@map} slot.
}
\section{File-backed matrices}{

Functions with the \code{chunkSize} parameter work best with file-backed matrices
such as \link[BEDMatrix:BEDMatrix-class]{BEDMatrix::BEDMatrix} objects. To avoid loading the whole,
potentially very large matrix into memory, these functions will load chunks
of the file-backed matrix into memory and perform the operations on one
chunk at a time. The size of the chunks is determined by the \code{chunkSize}
parameter. Care must be taken to not set \code{chunkSize} too high to avoid
memory shortage, particularly when combined with parallel computing.
}

\section{Multi-level parallelism}{

Functions with the \code{nCores}, \code{i}, and \code{j} parameters provide
capabilities for both parallel and distributed computing.

For parallel computing, \code{nCores} determines the number of cores the code is
run on. Memory usage can be an issue for higher values of \code{nCores} as R is
not particularly memory-efficient. As a rule of thumb, at least around
\code{(nCores * object_size(chunk)) + object_size(result)} MB of total memory
will be needed for operations on file-backed matrices, not including
potential copies of your data that might be created (for example
\code{\link[stats:lsfit]{stats::lsfit()}} runs \code{cbind(1, X)}). \code{i} and \code{j} can be used to include or
exclude certain rows or columns. Internally, the \code{\link[parallel:mclapply]{parallel::mclapply()}}
function is used and therefore parallel computing will not work on Windows
machines.

For distributed computing, \code{i} and \code{j} determine the subset of the input
matrix that the code runs on. In an HPC environment, this can be used not
just to include or exclude certain rows or columns, but also to partition
the task among many nodes rather than cores. Scheduler-specific code and
code to aggregate the results need to be written by the user. It is
recommended to set \code{nCores} to \code{1} as nodes are often cheaper than cores.
}

\section{Example dataset}{

The \code{extdata} folder contains example files that were generated from the
250k SNP and phenotype data in \href{http://www.nature.com/nature/journal/v465/n7298/full/nature08800.html}{Atwell et al. (2010)}.
Only the first 300 SNPs of chromosome 1, 2, and 3 were included to keep the
size of the example dataset small.
\href{https://www.cog-genomics.org/plink2}{PLINK} was used to convert the data to
\href{https://www.cog-genomics.org/plink2/input#bed}{.bed} and
\href{https://www.cog-genomics.org/plink2/input#raw}{.raw} files. \code{FT10} has been
chosen as a phenotype and is provided as an \href{https://www.cog-genomics.org/plink2/input#pheno}{alternate phenotype file}. The file is
intentionally shuffled to demonstrate that the additional phenotypes are put
in the same order as the rest of the phenotypes.
}

\seealso{
\link[BEDMatrix:BEDMatrix-package]{BEDMatrix::BEDMatrix-package},
\link[LinkedMatrix:LinkedMatrix-package]{LinkedMatrix::LinkedMatrix-package}, and \link[symDMatrix:symDMatrix-package]{symDMatrix::symDMatrix-package}
for an introduction to the respective packages.
}

\docType{class}
\name{BGData-class}
\alias{BGData-class}
\title{An S4 Class to Represent Phenotype and Genotype Data}
\description{
    This class is inspired by the phenotype/genotype file format .raw and its
    binary companion (also known as .bed) of
    \href{https://www.cog-genomics.org/plink2}{PLINK}. It is used by several
    functions of this package such as \code{GWAS} for performing a Genome Wide
    Association Study or \code{getG} for calculating a genomic relationship
    matrix.
}
\details{
    There are several ways to create an instance of this class:

    \itemize{
        \item from arbitrary phenotype/genotype data using the \code{BGData}
            constructor function.
        \item from a BED file using \code{as.BGData} and \code{BEDMatrix}.
        \item from a previously saved \code{BGData} object using
            \code{load.BGData}.
        \item from multiple files (even a mixture of different file types)
            using \code{LinkedMatrix}.
        \item from a .raw file (or a .ped-like file) using
            \code{readRAW}, \code{readRAW_matrix}, or
            \code{readRAW_big.matrix}.
    }

    A .ped file can be recoded to a .raw file in
    \href{https://www.cog-genomics.org/plink2}{PLINK} using \code{plink --file
    myfile --recodeA}, or converted to a BED file using \code{plink --file
    myfile --make-bed}. Conversely, a BED file can be transformed back to a
    .ped file using \code{plink --bfile myfile --recode} or to a .raw file
    using \code{plink --bfile myfile --recodeA} without losing information.
}
\section{Slots}{
    \describe{
        \item{\code{geno}}{
            A \code{geno} object that contains genotypes.
        }
        \item{\code{pheno}}{
            A \code{data.frame} that contains sample information (including
            phenotypes).
        }
        \item{\code{map}}{
            A \code{data.frame} that contains variant information.
        }
    }
}
\seealso{
    \code{\link{BGData}}, \code{\link{as.BGData}}, \code{\link{load.BGData}},
    \code{\link{readRAW}} to create \code{BGData} objects.

    \code{\link[LinkedMatrix]{LinkedMatrix-class}} and
    \code{\link[BEDMatrix]{BEDMatrix-class}} for more information on the above
    mentioned classes.
}
\examples{
X <- matrix(data = rnorm(100), nrow = 10, ncol = 10)
Y <- data.frame(y = runif(10))
MAP <- data.frame(means = colMeans(X), freqNA = colMeans(is.na(X)))
DATA <- BGData(geno = X, pheno = Y, map = MAP)

dim(DATA@geno)
head(DATA@pheno)
head(DATA@map)
}

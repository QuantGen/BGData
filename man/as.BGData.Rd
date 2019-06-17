\name{as.BGData}
\alias{as.BGData}
\alias{as.BGData.BEDMatrix}
\alias{as.BGData.ColumnLinkedMatrix}
\alias{as.BGData.RowLinkedMatrix}
\title{Convert Other Objects to BGData Objects}
\description{
    Converts other objects to \code{BGData} objects by loading supplementary
    phenotypes and map files referenced by the object to be used for the
    \code{@pheno} and \code{@map} slot, respectively.

    Currently supported are \code{BEDMatrix} objects, plain or nested in
    \code{ColumnLinkedMatrix} objects.
}
\usage{
as.BGData(x, alternatePhenotypeFile = NULL, ...)

\method{as.BGData}{BEDMatrix}(x, alternatePhenotypeFile = NULL, ...)

\method{as.BGData}{ColumnLinkedMatrix}(x, alternatePhenotypeFile = NULL,
  ...)

\method{as.BGData}{RowLinkedMatrix}(x, alternatePhenotypeFile = NULL,
  ...)
}
\arguments{
    \item{x}{
        An object. Currently supported are \code{BEDMatrix} objects, plain or
        nested in \code{ColumnLinkedMatrix} objects.
    }
    \item{alternatePhenotypeFile}{
        Path to an
        \href{https://www.cog-genomics.org/plink2/input#pheno}{alternate
        phenotype file}.
    }
    \item{...}{
        Additional arguments to the \code{read.table} or \code{fread} call (if
        data.table package is installed) call to parse the alternate pheno
        file.
    }
}
\details{
    The .ped and .raw formats only allows for a single phenotype. If more
    phenotypes are required it is possible to store them in an
    \href{https://www.cog-genomics.org/plink2/input#pheno}{alternate phenotype
    file}. The path to such a file can be provided with
    \code{alternatePhenotypeFile} and will be merged with the data in the
    \code{@pheno} slot. The first and second columns of that file must contain
    family and within-family IDs, respectively.

    For \code{BEDMatrix} objects: If a .fam file (which corresponds to the
    first six columns of a .ped or .raw file) of the same name and in the same
    directory as the BED file exists, the \code{@pheno} slot will be populated
    with the data stored in that file. Otherwise a stub that only contains an
    \code{IID} column populated with the rownames of \code{@geno} will be
    generated. The same will happen for a .bim file for the \code{@map} slot.

    For \code{ColumnLinkedMatrix} objects: See the case for \code{BEDMatrix}
    objects, but only the .fam file of the first node of the
    \code{LinkedMatrix} will be read and used for the \code{@pheno} slot, and
    the .bim files of all nodes will be combined and used for the \code{@map}
    slot.
}
\value{
    A \code{BGData} object.
}
\seealso{
    \code{\link[=readRAW]{readRAW()}} to convert text files to \code{BGData}
    objects. \code{\link{BGData-class}},
    \code{\link[BEDMatrix]{BEDMatrix-class}},
    \code{\link[LinkedMatrix]{ColumnLinkedMatrix-class}} for more information
    on the above mentioned classes. \code{\link[utils]{read.table}} and
    \code{\link[data.table]{fread}} to learn more about extra arguments that
    can be passed via \code{...}.
}
\examples{
# Path to example data
path <- system.file("extdata", package = "BGData")

# Convert a single BEDMatrix object to a BGData object
chr1 <- BEDMatrix::BEDMatrix(paste0(path, "/chr1.bed"))
bg1 <- as.BGData(chr1)

# Convert multiple BEDMatrix objects in a ColumnLinkedMatrix to a BGData object
chr2 <- BEDMatrix::BEDMatrix(paste0(path, "/chr2.bed"))
chr3 <- BEDMatrix::BEDMatrix(paste0(path, "/chr3.bed"))
clm <- ColumnLinkedMatrix(chr1, chr2, chr3)
bg2 <- as.BGData(clm)

# Load additional (alternate) phenotypes
bg3 <- as.BGData(clm, alternatePhenotypeFile = paste0(path, "/pheno.txt"))
}

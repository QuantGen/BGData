\docType{methods}
\name{initialize,BGData-method}
\alias{initialize,BGData-method}
\title{Creates a New BGData Instance}
\description{
    This method is run when a \linkS4class{BGData} object is created using
    \code{BGData(...)} or \code{new("BGData", ...)}.
}
\usage{
\S4method{initialize}{BGData}(.Object, geno, pheno, map)
}
\arguments{
    \item{.Object}{
        The \linkS4class{BGData} instance to be initialized. This argument is passed in by R and can be ignored, but still needs to be documented.
    }
    \item{geno}{
        A \linkS4class{geno} object that contains genotypes. \linkS4class{geno}
        is a class union of several matrix-like types, many of them suitable
        for very large datasets. Currently supported are
        \link[LinkedMatrix:LinkedMatrix-class]{LinkedMatrix::LinkedMatrix},
        \link[BEDMatrix:BEDMatrix-class]{BEDMatrix::BEDMatrix},
        \link[bigmemory:big.matrix-class]{bigmemory::big.matrix},
        \code{ff_matrix}, and \code{matrix}.
    }
    \item{pheno}{
        A \code{data.frame} that contains sample information (including
        phenotypes). A stub that only contains a \code{sample_id} column
        populated with either the rownames of \code{@geno} or a sequence will
        be generated if missing.
    }
    \item{map}{
        A \code{data.frame} that contains variant information. A stub that only
        contains a \code{variant_id} column populated with either the colnames
        of \code{@geno} or a sequence will be generated if missing.
    }
}

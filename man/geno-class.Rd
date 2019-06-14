\docType{class}
\name{geno-class}
\alias{geno-class}
\title{An Abstract S4 Class Union of Matrix-Like Types.}
\description{
\linkS4class{geno} is a class union of several matrix-like types, many of them
suitable for very large datasets. Currently supported are
\link[LinkedMatrix:LinkedMatrix-class]{LinkedMatrix::LinkedMatrix}, \link[BEDMatrix:BEDMatrix-class]{BEDMatrix::BEDMatrix},
\link[bigmemory:big.matrix-class]{bigmemory::big.matrix}, \code{ff_matrix}, and \code{matrix}.
}
\seealso{
The \code{@geno} slot of \linkS4class{BGData} that accepts \linkS4class{geno}
objects.
}

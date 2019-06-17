\docType{class}
\name{geno-class}
\alias{geno-class}
\title{An Abstract S4 Class Union of Matrix-Like Types}
\description{
    \code{geno} is a class union of several matrix-like types, many of
    them suitable for very large datasets.

    Currently supported are \code{LinkedMatrix}, \code{BEDMatrix},
    \code{big.matrix},\code{ff_matrix}, and \code{matrix}.
}
\seealso{
    \code{\link[LinkedMatrix]{LinkedMatrix-class}},
    \code{\link[BEDMatrix]{BEDMatrix-class}},
    \code{\link[bigmemory]{big.matrix-class}}, \code{\link[ff]{ff}}, and
    \code{\link[base]{matrix}} for more information on each matrix-like type.

    \code{\link{BGData-class}} for more information on the \code{BGData} class,
    in particular its \code{@geno} slot that accepts \code{geno} objects.
}

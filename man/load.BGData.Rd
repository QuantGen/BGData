\name{load.BGData}
\alias{load.BGData}
\title{Loads BGData (and Other) Objects from .RData Files}
\description{
    This function is similar to \code{load}, but also initializes the different
    types of objects that the \code{@geno} slot of a \code{BGData} object can
    take.

    Currently supported are \code{ff_matrix}, \code{big.matrix}, and
    \code{BEDMatrix} objects. If the object is of type \code{LinkedMatrix}, all
    nodes will be initialized with their appropriate method.
}
\usage{
load.BGData(file, envir = parent.frame())
}
\arguments{
    \item{file}{
        The name of the .RData file to be loaded.
    }
    \item{envir}{
        The environment where to load the data.
    }
}
\seealso{
    \code{\link{BGData-class}}, \code{\link[ff]{ff}},
    \code{\link[bigmemory]{big.matrix-class}},
    \code{\link[BEDMatrix]{BEDMatrix-class}}, and
    \code{\link[LinkedMatrix]{LinkedMatrix-class}} for more information on the
    above mentioned classes.
}

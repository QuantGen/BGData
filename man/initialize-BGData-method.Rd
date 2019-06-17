\docType{methods}
\name{initialize,BGData-method}
\alias{initialize,BGData-method}
\title{Creates a New BGData Instance}
\description{
    This method is run when a \code{BGData} object is created using
    \code{BGData(...)} or \code{new("BGData", ...)}.
}
\usage{
\S4method{initialize}{BGData}(.Object, geno, pheno, map)
}
\arguments{
    \item{.Object}{
        The \code{BGData} instance to be initialized. This argument is passed
        in by R and can be ignored, but still needs to be documented.
    }
    \item{geno}{
        A \code{geno} object that contains genotypes.
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
\seealso{
     \code{\link{BGData-class}} and \code{\link{geno-class}} for more
     information on the above mentioned classes.
}

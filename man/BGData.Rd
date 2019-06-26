\name{BGData}
\alias{BGData}
\title{Creates a New BGData Instance}
\description{
    This function constructs a new \code{BGData} object.
}
\usage{
BGData(geno, pheno = NULL, map = NULL)
}
\arguments{
    \item{geno}{
        A \code{geno} object that contains genotypes.
    }
    \item{pheno}{
        A \code{data.frame} that contains sample information (including
        phenotypes). A stub that only contains a \code{sample_id} column
        populated with either the rownames of \code{geno} or a sequence will be
        generated if \code{NULL}
    }
    \item{map}{
        A \code{data.frame} that contains variant information. A stub that only
        contains a \code{variant_id} column populated with either the colnames
        of \code{geno} or a sequence will be generated if \code{NULL}
    }
}
\seealso{
     \code{\link{BGData-class}} and \code{\link{geno-class}} for more
     information on the above mentioned classes.
}

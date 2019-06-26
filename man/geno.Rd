\name{geno}
\alias{geno}
\alias{geno<-}
\alias{pheno}
\alias{pheno<-}
\alias{map}
\alias{map<-}
\title{Getting/Setting Genotypes, Sample Information, and Variant Information}
\description{
    A set of generic functions for getting/setting the genotypes, sample
    information, and variant information.
}
\usage{
geno(x)
geno(x) <- value

pheno(x)
pheno(x) <- value

map(x)
map(x) <- value
}
\arguments{
    \item{x}{
        The object from/on which to get/set genotypes, sample information, and
        variant information. Typically a \code{BGData} object.
    }
    \item{value}{
        Typically a \code{geno} object for the \code{geno} setter.

        Typically a \code{data.frame} object for the \code{pheno} setter.

        Typically a \code{data.frame} object for the \code{map} setter.
    }
}
\seealso{
    \itemize{
        \item \code{\link{BGData-class}}
        \item \code{\link{geno-class}}
    }
}
\examples{
# Load example data
bg <- BGData:::loadExample()

# Access genotypes
geno(bg)

# Access sample information
pheno(bg)

# Access variant information
map(bg)
}
\keyword{methods}

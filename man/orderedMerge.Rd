\name{orderedMerge}
\alias{orderedMerge}
\title{Merge Two Data Frames Keeping the Order of the First}
\description{
    This is a simplified version of \code{merge} useful for merging additional
    data into a \code{BGData} object while keeping the order of the data in the
    \code{BGData} object.
}
\usage{
orderedMerge(x, y, by = c(1L, 2L))
}
\arguments{
    \item{x}{
        Data frame
    }
    \item{y}{
        Data frame
    }
    \item{by}{
        Specifications of the columns used for merging. Defaults to the first
        two columns of the data frame, which traditionally has the family ID
        and the individual ID.
    }
}
\value{
    Merged data frame
}
\seealso{
    \code{\link{BGData-class}} for more information on the \code{BGData} class.
}

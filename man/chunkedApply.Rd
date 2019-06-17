\name{chunkedApply}
\alias{chunkedApply}
\title{Applies a Function on Each Row or Column of a File-Backed Matrix}
\description{
    Similar to \code{apply}, but designed for file-backed matrices. The
    function brings chunks of an object into physical memory by taking subsets,
    and applies a function on either the rows or the columns of the chunks
    using an optimized version of \code{apply}. If \code{nCores} is greater
    than 1, the function will be applied in parallel using \code{mclapply}. In
    that case the subsets of the object are taken on the slaves.
}
\usage{
chunkedApply(X, MARGIN, FUN, i = seq_len(nrow(X)),
  j = seq_len(ncol(X)), chunkSize = 5000L,
  nCores = getOption("mc.cores", 2L), verbose = FALSE, ...)
}
\arguments{
    \item{X}{
        A file-backed matrix, typically \code{@geno} of a \code{BGData} object.
    }
    \item{MARGIN}{
        The subscripts which the function will be applied over. 1 indicates
        rows, 2 indicates columns.
    }
    \item{FUN}{
        The function to be applied.
    }
    \item{i}{
        Indicates which rows of \code{X} should be used. Can be integer,
        boolean, or character. By default, all rows are used.
    }
    \item{j}{
        Indicates which columns of \code{X} should be used. Can be integer,
        boolean, or character. By default, all columns are used.
    }
    \item{chunkSize}{
        The number of rows or columns of \code{X} that are brought into
        physical memory for processing per core. If \code{NULL}, all elements
        in \code{i} or \code{j} are used. Defaults to 5000.
    }
    \item{nCores}{
        The number of cores (passed to \code{mclapply}). Defaults to the number
        of cores as detected by \code{detectCores}.
    }
    \item{verbose}{
        Whether progress updates will be posted. Defaults to \code{FALSE}.
    }
    \item{...}{
        Additional arguments to be passed to the \code{apply} like function.
    }
}
\seealso{
    \code{\link{file-backed-matrices}} for more information on file-backed
    matrices. \code{\link{multi-level-parallelism}} for more information on
    multi-level parallelism. \code{\link{BGData-class}} for more information on
    the \code{BGData} class.
}
\examples{
# Restrict number of cores to 1 on Windows
if (.Platform$OS.type == "windows") {
    options(mc.cores = 1)
}

# Load example data
bg <- BGData:::loadExample()

# Compute standard deviation of columns
chunkedApply(X = bg@geno, MARGIN = 2, FUN = sd)
}

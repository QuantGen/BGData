\name{getG_symDMatrix}
\alias{getG_symDMatrix}
\title{Computes a Very Large Genomic Relationship Matrix}
\description{
    Computes a positive semi-definite symmetric genomic relation matrix G=XX'
    offering options for centering and scaling the columns of \code{X}
    beforehand.
}
\usage{
getG_symDMatrix(X, center = TRUE, scale = TRUE, scaleG = TRUE,
  minVar = 1e-05, blockSize = 5000L,
  folderOut = paste0("symDMatrix_", randomString()), vmode = "double",
  i = seq_len(nrow(X)), j = seq_len(ncol(X)), chunkSize = 5000L,
  nCores = getOption("mc.cores", 2L), verbose = FALSE)
}
\arguments{
    \item{X}{
        A matrix-like object, typically \code{@geno} of a \code{BGData} object.
    }
    \item{center}{
        Either a logical value or a numeric vector of length equal to the
        number of columns of \code{X}. If \code{FALSE}, no centering is done.
        Defaults to \code{TRUE}.
    }
    \item{scale}{
        Either a logical value or a numeric vector of length equal to the
        number of columns of \code{X}. If \code{FALSE}, no scaling is done.
        Defaults to \code{TRUE}.
    }
    \item{scaleG}{
        TRUE/FALSE whether xx' must be scaled.
    }
    \item{minVar}{
        Columns with variance lower than this value will not be used in the
        computation (only if \code{scale} is not \code{FALSE}).
    }
    \item{blockSize}{
        The number of rows and columns of each block. If \code{NULL}, a single
        block of the same length as \code{i} will be created. Defaults to 5000.
    }
    \item{folderOut}{
        The path to the folder where to save the \code{symDMatrix} object.
        Defaults to a random string prefixed with "symDMatrix_".
    }
    \item{vmode}{
        vmode of \code{ff} objects.
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
        The number of columns of \code{X} that are brought into physical memory
        for processing per core. If \code{NULL}, all columns of \code{X} are
        used. Defaults to 5000.
    }
    \item{nCores}{
        The number of cores (passed to \code{mclapply}). Defaults to the number
        of cores as detected by \code{detectCores}.
    }
    \item{verbose}{
        Whether progress updates will be posted. Defaults to \code{FALSE}.
    }
}
\details{
    Even very large genomic relationship matrices are supported by partitioning
    \code{X} into blocks and calling \code{getG} on these blocks. This function
    performs the block computations sequentially, which may be slow. In an HPC
    environment, performance can be improved by manually distributing these
    operations to different nodes.
}
\value{
    A \code{symDMatrix} object.
}
\seealso{
    \code{\link{multi-level-parallelism}} for more information on multi-level
    parallelism. \code{\link[symDMatrix]{symDMatrix-class}} and
    \code{\link{BGData-class}} for more information on the \code{BGData} class.
    \code{\link{getG}} to learn more about the underlying method.
}

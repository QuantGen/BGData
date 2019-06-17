\name{getG}
\alias{getG}
\title{Computes a Genomic Relationship Matrix}
\description{
    Computes a positive semi-definite symmetric genomic relation matrix G=XX'
    offering options for centering and scaling the columns of \code{X}
    beforehand.
}
\usage{
getG(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05,
  i = seq_len(nrow(X)), j = seq_len(ncol(X)), i2 = NULL,
  chunkSize = 5000L, nCores = getOption("mc.cores", 2L),
  verbose = FALSE)
}
\arguments{
    \item{X}{
        A matrix-like object, typically \code{@geno} of a \code{BGData} object.
    }
    \item{center}{
        Either a logical value or a numeric vector of length equal to the
        number of columns of \code{X}. Numeric vector required if \code{i2} is
        used. If \code{FALSE}, no centering is done. Defaults to \code{TRUE}.
    }
    \item{scale}{
        Either a logical value or a numeric vector of length equal to the
        number of columns of \code{X}. Numeric vector required if \code{i2} is
        used. If \code{FALSE}, no scaling is done. Defaults to \code{TRUE}.
    }
    \item{scaleG}{
        Whether XX' should be scaled. Defaults to \code{TRUE}.
    }
    \item{minVar}{
        Columns with variance lower than this value will not be used in the
        computation (only if \code{scale} is not \code{FALSE}).
    }
    \item{i}{
        Indicates which rows of \code{X} should be used. Can be integer,
        boolean, or character. By default, all rows are used.
    }
    \item{j}{
        Indicates which columns of \code{X} should be used. Can be integer,
        boolean, or character. By default, all columns are used.
    }
    \item{i2}{
        Indicates which rows should be used to compute a block of the genomic
        relationship matrix. Will compute XY' where X is determined by \code{i}
        and \code{j} and Y by \code{i2} and \code{j}. Can be integer, boolean,
        or character. If \code{NULL}, the whole genomic relationship matrix XX'
        is computed. Defaults to \code{NULL}.
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
    If \code{center = FALSE}, \code{scale = FALSE} and \code{scaleG = FALSE},
    \code{getG} produces the same outcome than \code{tcrossprod}.
}
\value{
    A positive semi-definite symmetric numeric matrix.
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

# Compute a scaled genomic relationship matrix from centered and scaled
# genotypes
g1 <- getG(X = bg@geno)

# Disable scaling of G
g2 <- getG(X = bg@geno, scaleG = FALSE)

# Disable centering of genotypes
g3 <- getG(X = bg@geno, center = FALSE)

# Disable scaling of genotypes
g4 <- getG(X = bg@geno, scale = FALSE)

# Provide own scales
scales <- chunkedApply(X = bg@geno, MARGIN = 2, FUN = sd)
g4 <- getG(X = bg@geno, scale = scales)

# Provide own centers
centers <- chunkedApply(X = bg@geno, MARGIN = 2, FUN = mean)
g5 <- getG(X = bg@geno, center = centers)

# Only use the first 50 individuals (useful to account for population structure)
g6 <- getG(X = bg@geno, i = 1:50)

# Only use the first 100 markers (useful to ignore some markers)
g7 <- getG(X = bg@geno, j = 1:100)

# Compute unscaled G matrix by combining blocks of $XX_{i2}'$ where $X_{i2}$ is
# a horizontal partition of X. This is useful for distributed computing as each
# block can be computed in parallel. Centers and scales need to be precomputed.
block1 <- getG(X = bg@geno, i2 = 1:100, center = centers, scale = scales)
block2 <- getG(X = bg@geno, i2 = 101:199, center = centers, scale = scales)
g8 <- cbind(block1, block2)

# Compute unscaled G matrix by combining blocks of $X_{i}X_{i2}'$ where both
# $X_{i}$ and $X_{i2}$ are horizontal partitions of X. Similarly to the example
# above, this is useful for distributed computing, in particular to compute
# very large G matrices. Centers and scales need to be precomputed. This
# approach is similar to the one taken by the symDMatrix package, but the
# symDMatrix package adds memory-mapped blocks, only stores the upper side of
# the triangular matrix, and provides a type that allows for indexing as if the
# full G matrix is in memory.
block11 <- getG(X = bg@geno, i = 1:100, i2 = 1:100, center = centers, scale = scales)
block12 <- getG(X = bg@geno, i = 1:100, i2 = 101:199, center = centers, scale = scales)
block21 <- getG(X = bg@geno, i = 101:199, i2 = 1:100, center = centers, scale = scales)
block22 <- getG(X = bg@geno, i = 101:199, i2 = 101:199, center = centers, scale = scales)
g9 <- rbind(
    cbind(block11, block12),
    cbind(block21, block22)
)
}

\name{chunkedMap}
\alias{chunkedMap}
\title{Applies a Function on Each Chunk of a File-Backed Matrix}
\usage{
chunkedMap(X, FUN, i = seq_len(nrow(X)), j = seq_len(ncol(X)),
  chunkBy = 2L, chunkSize = 5000L, nCores = getOption("mc.cores",
  2L), verbose = FALSE, ...)
}
\arguments{
\item{X}{A file-backed matrix, typically \code{@geno} of a \linkS4class{BGData}
object.}

\item{FUN}{The function to be applied on each chunk.}

\item{i}{Indicates which rows of \code{X} should be used. Can be integer,
boolean, or character. By default, all rows are used.}

\item{j}{Indicates which columns of \code{X} should be used. Can be integer,
boolean, or character. By default, all columns are used.}

\item{chunkBy}{Whether to extract chunks by rows (1) or by columns (2).
Defaults to columns (2).}

\item{chunkSize}{The number of rows or columns of \code{X} that are brought into
physical memory for processing per core. If \code{NULL}, all elements in \code{i} or
\code{j} are used. Defaults to 5000.}

\item{nCores}{The number of cores (passed to \code{\link[parallel:mclapply]{parallel::mclapply()}}).
Defaults to the number of cores as detected by \code{\link[parallel:detectCores]{parallel::detectCores()}}.}

\item{verbose}{Whether progress updates will be posted. Defaults to \code{FALSE}.}

\item{...}{Additional arguments to be passed to the \code{\link[base:apply]{base::apply()}} like
function.}
}
\description{
Similar to \code{\link[base:lapply]{base::lapply()}}, but designed for file-backed matrices. The
function brings chunks of an object into physical memory by taking subsets,
and applies a function on them. If \code{nCores} is greater than 1, the function
will be applied in parallel using \code{\link[parallel:mclapply]{parallel::mclapply()}}. In that case the
subsets of the object are taken on the slaves.
}
\section{File-backed matrices}{

Functions with the \code{chunkSize} parameter work best with file-backed matrices
such as \link[BEDMatrix:BEDMatrix-class]{BEDMatrix::BEDMatrix} objects. To avoid loading the whole,
potentially very large matrix into memory, these functions will load chunks
of the file-backed matrix into memory and perform the operations on one
chunk at a time. The size of the chunks is determined by the \code{chunkSize}
parameter. Care must be taken to not set \code{chunkSize} too high to avoid
memory shortage, particularly when combined with parallel computing.
}

\section{Multi-level parallelism}{

Functions with the \code{nCores}, \code{i}, and \code{j} parameters provide
capabilities for both parallel and distributed computing.

For parallel computing, \code{nCores} determines the number of cores the code is
run on. Memory usage can be an issue for higher values of \code{nCores} as R is
not particularly memory-efficient. As a rule of thumb, at least around
\code{(nCores * object_size(chunk)) + object_size(result)} MB of total memory
will be needed for operations on file-backed matrices, not including
potential copies of your data that might be created (for example
\code{\link[stats:lsfit]{stats::lsfit()}} runs \code{cbind(1, X)}). \code{i} and \code{j} can be used to include or
exclude certain rows or columns. Internally, the \code{\link[parallel:mclapply]{parallel::mclapply()}}
function is used and therefore parallel computing will not work on Windows
machines.

For distributed computing, \code{i} and \code{j} determine the subset of the input
matrix that the code runs on. In an HPC environment, this can be used not
just to include or exclude certain rows or columns, but also to partition
the task among many nodes rather than cores. Scheduler-specific code and
code to aggregate the results need to be written by the user. It is
recommended to set \code{nCores} to \code{1} as nodes are often cheaper than cores.
}

\examples{
# Restrict number of cores to 1 on Windows
if (.Platform$OS.type == "windows") {
    options(mc.cores = 1)
}

# Load example data
bg <- BGData:::loadExample()

# Compute column sums of each chunk
chunkedMap(X = bg@geno, FUN = colSums)
}

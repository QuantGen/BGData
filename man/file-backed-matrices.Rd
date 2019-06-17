\name{file-backed-matrices}
\alias{file-backed-matrices}
\title{File-Backed Matrices}
\description{
    Functions with the \code{chunkSize} parameter work best with file-backed
    matrices such as \code{BEDMatrix} objects. To avoid loading the whole,
    potentially very large matrix into memory, these functions will load chunks
    of the file-backed matrix into memory and perform the operations on one
    chunk at a time. The size of the chunks is determined by the
    \code{chunkSize} parameter. Care must be taken to not set \code{chunkSize}
    too high to avoid memory shortage, particularly when combined with parallel
    computing.
}
\seealso{
    \code{\link[BEDMatrix]{BEDMatrix-class}} as an example of a file-backed
    matrix.
}

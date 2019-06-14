\name{findRelated}
\alias{findRelated}
\alias{findRelated.matrix}
\alias{findRelated.symDMatrix}
\title{Find related individuals in a relationship matrix}
\usage{
findRelated(x, ...)

\method{findRelated}{matrix}(x, cutoff = 0.03, ...)

\method{findRelated}{symDMatrix}(x, cutoff = 0.03, verbose = FALSE,
  ...)
}
\arguments{
\item{x}{A matrix-like object with dimnames.}

\item{...}{Additional arguments for methods.}

\item{cutoff}{The cutoff between 0 and 1 for related individuals to be included
in the output. Defaults to 0.03.}

\item{verbose}{Whether progress updates will be posted. Defaults to \code{FALSE}.}
}
\value{
A vector of names or indices of related individuals.
}
\description{
Find related individuals in a relationship matrix.
}
\section{Methods (by class)}{
\itemize{
\item \code{matrix}: Find related individuals in matrices

\item \code{symDMatrix}: Find related individuals in symDMatrix objects
}}

\examples{
# Load example data
bg <- BGData:::loadExample()

G <- getG(bg@geno)
findRelated(G)
}

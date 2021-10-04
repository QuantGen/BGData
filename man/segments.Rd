\name{segments}
\alias{segments}
\title{Find non-overlapping segments based on a summary statistic}
\description{
    Given a summary statistic and a threshold, this function computes the
    number of non-overlapping segments, each defined as a discovery (i.e.,
    \code{statistic[i] <= threshold)} +/- a gap, in the same units as \code{bp}
    (often base-pair position).
}
\usage{
segments(statistic, chr, bp, threshold, gap, trim = FALSE, verbose = FALSE)
}
\arguments{
    \item{statistic}{
        A statistic (e.g., BFDR or p-values).
    }
    \item{chr}{
        A vector containing the chromosome for each value of \code{statistic}.
    }
    \item{bp}{
        A vector containing the base-pair positions for each value of
        \code{statistic}.
    }
    \item{threshold}{
        The threshold to determine 'significance' (e.g., \code{1e-5} for
        p-values).
    }
    \item{gap}{
        1/2 of the length of the desired segments.
    }
    \item{trim}{
        Whether to collapse segments that were artifically inflated by
        \code{gap}. Defaults to \code{FALSE}.
    }
    \item{verbose}{
        Whether progress updates will be posted. Defaults to \code{FALSE}.
    }
}
\value{
    A data frame containing the following information:
    \item{chr}{
        Chromosome
    }
    \item{start}{
        Index where segment starts within \code{statistic}.
    }
    \item{end}{
        Index where segment ends within \code{statistic}.
    }
    \item{length}{
        Length of segment.
    }
    \item{bpStart}{
        Base-pair position where segment starts.
    }
    \item{bpEnd}{
        Base-pair position where segment ends.
    }
    \item{bpLength}{
        Length of segment in base-pair positions.
    }
    \item{minValue}{
        Smallest value of \code{statistic} within segment.
    }
    \item{minValuePos}{
        Position of variant with the smallest value of \code{statistic} within
        segment.
    }
}
\examples{
library(BGData)

# Load example data
bg <- BGData:::loadExample()

# Perform GWAS
pValues <- GWAS(
    formula = FT10 ~ 1,
    data = bg,
    method = "rayOLS"
)

# Determine segments within +/- 1MB from a significant variant
segments <- segments(
    statistic = pValues[, 4],
    chr = map(bg)$chromosome,
    bp = map(bg)$base_pair_position,
    threshold = 1e-5,
    gap = 1e6,
    trim = FALSE,
    verbose = FALSE
)
}

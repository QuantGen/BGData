\name{preprocess}
\alias{preprocess}
\title{Center, scale, and impute data}
\description{
    A faster version of \code{\link[base]{scale}} with a similar interface that
    also allows for imputation. The main difference is that this version scales
    by the standard deviation regardless of whether centering is enabled or
    not. If centering is enabled, missing values are imputed by 0, otherwise by
    the mean of the column that contains the value.
}
\usage{
preprocess(X, center = FALSE, scale = FALSE, impute = FALSE)
}
\arguments{
    \item{X}{
        A numeric matrix.
    }
    \item{center}{
        Either a logical value or numeric vector of length equal to the number
        of columns of \code{X}.
    }
    \item{scale}{
        Either a logical value or numeric vector of length equal to the number
        of columns of \code{X}.
    }
    \item{impute}{
        Indicates whether missing values should be imputed.
    }
}
\value{
    The centered, scaled, and imputed matrix.
}
\seealso{
    \code{\link[base]{scale}}, which this function tries to improve upon.
}
\examples{
# Load example data
bg <- BGData:::loadExample()

# Center and scale genotypes
W <- preprocess(as.matrix(geno(bg)), center = TRUE, scale = TRUE)
}

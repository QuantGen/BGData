\name{FWD}
\alias{FWD}
\title{Performs Forward Regressions}
\description{
    Performs forward regression of \code{y} on the columns of \code{X}.
    Predictors are added, one at a time, each time adding the one that produces
    the largest reduction in the residual sum of squares (RSS). The function
    returns estimates and summaries for the entire forward search. This
    function performs a similar search than that of \code{step(,
        direction='forward')}, however, \code{FWD()} is optimized for
    computational speed for linear models with very large sample size. To
    achieve fast computations, the software first computes the sufficient
    statistics X'X and X'y. At each step, the function first finds the
    predictor that produces the largest reduction in the sum of squares (this
    can be derived from X'X, X'y and the current solution of effects), and then
    updates the estimates of effects for the resulting model using Gauss Seidel
    iterations performed on the linear system (X'X)b=X'y, iterating only over
    the elements of b that are active in the model.
}
\usage{
FWD(y, X, df = 20, tol = 1e-7, maxIter = 1000, centerImpute = TRUE,
    verbose = TRUE)
}
\arguments{
    \item{y}{
        The response vector (numeric nx1).
    }
    \item{X}{
        An (nxp) numeric matrix. Columns are the features (aka predictors)
        considered in the forward search. The rows of \code{X} must be matched
        to the entries of \code{y}.
    }
    \item{df}{
        Defines the maximum number of predictors to be included in the model.
        For complete forward search, set \code{df = ncol(X)}.
    }
    \item{tol}{
        A tolerance parameter to control when to stop the Gauss Seidel
        algorithm.
    }
    \item{maxIter}{
        The maximum number of iterations for the Gauss Seidel algorithm (only
        used when the algorithm is not stopped by the tolerance parameter).
    }
    \item{centerImpute}{
        Whether to center the columns of \code{X} and impute the missing values
        with the column means.
    }
    \item{verbose}{
        Use \code{verbose = TRUE} to print summaries of the forward search.
    }
}
\value{
    A list with two entries:
    \itemize{
        \item \code{B}: (pxdf+1) includes the estimated effects for each
        predictor (rows) at each step of the forward search (df, in columns).
        \item \code{path}: A data frame providing the order in which variables
        were added to the model (\code{variable}) and statistics for each step
        of the forward search (\code{RSS}, \code{LogLik}, \code{VARE} (the
        residual variance), \code{DF}, \code{AIC}, and \code{BIC}).
    }
}

preprocess <- function(X, center = FALSE, scale = FALSE, impute = FALSE, nCores = getOption("mc.cores", 2L)) {
    if (!(is.numeric(X) && length(dim(X)) == 2)) {
        stop("'X' needs to be a numeric matrix")
    }
    if (!(is.logical(center) && length(center) == 1L) && !(is.numeric(center) && length(center) == ncol(X))) {
        stop("'center' needs to be either a logical vector of size 1 or a numeric vector of size 'ncol(X)'")
    }
    if (!(is.logical(scale) && length(scale) == 1L) && !(is.numeric(scale) && length(scale) == ncol(X))) {
        stop("'scale' needs to be either a logical vector of size 1 or a numeric vector of size 'ncol(X)'")
    }
    if (!(is.logical(impute) && length(impute) == 1L)) {
        stop("'impute' needs to be a logical vector of size 1")
    }
    if (!(is.numeric(nCores) && nCores > 0L)) {
        stop("'nCores' needs to be a positive number")
    }
    .Call(C_preprocess, X, center, scale, impute, as.integer(nCores))
}

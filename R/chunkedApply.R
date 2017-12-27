#' Applies a Function on Each Row or Column of a Memory-Mapped Matrix-Like
#' Object.
#'
#' Similar to [base::apply()], but designed for memory-mapped matrix-like
#' objects. The function brings chunks of an object into physical memory by
#' taking subsets, and applies a function on either the rows or the columns of
#' the chunks using an optimized version of [base::apply()]. If `nCores` is
#' greater than 1, the function will be applied in parallel using
#' [parallel::mclapply()]. In that case the subsets of the object are taken on
#' the slaves.
#'
#' @inheritSection BGData-package Memory-mapping
#' @inheritSection BGData-package Multi-level parallelism
#' @param X A memory-mapped matrix-like object, typically `@@geno` of a
#' [BGData-class] object.
#' @param MARGIN The subscripts which the function will be applied over. 1
#' indicates rows, 2 indicates columns.
#' @param FUN The function to be applied.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param bufferSize The number of rows or columns of `X` that are brought into
#' physical memory for processing per core. If `NULL`, all elements in `i` or `j` are
#' used. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @param ... Additional arguments to be passed to the [base::apply()] like
#' function.
#' @example man/examples/chunkedApply.R
#' @export
chunkedApply <- function(X, MARGIN, FUN, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {
    if (!length(dim(X))) {
        stop("dim(X) must have a positive length")
    }
    i <- convertIndexTypes(i, rownames(X))
    j <- convertIndexTypes(j, colnames(X))
    dimX <- c(length(i), length(j))
    if (is.null(bufferSize)) {
        bufferSize <- dimX[MARGIN]
        nBuffers <- 1L
    } else {
        nBuffers <- ceiling(dimX[MARGIN] / bufferSize)
    }
    bufferRanges <- LinkedMatrix:::chunkRanges(dimX[MARGIN], nBuffers)
    bufferApply <- function(curBuffer, ...) {
        if (verbose) {
            message("Buffer ", curBuffer, " of ", nBuffers, " ...")
        }
        if (MARGIN == 2L) {
            buffer <- X[i, j[seq(bufferRanges[1L, curBuffer], bufferRanges[2L, curBuffer])], drop = FALSE]
        } else {
            buffer <- X[i[seq(bufferRanges[1L, curBuffer], bufferRanges[2L, curBuffer])], j, drop = FALSE]
        }
        apply2(X = buffer, MARGIN = MARGIN, FUN = FUN, ...)
    }
    if (nCores == 1L) {
        res <- lapply(X = seq_len(nBuffers), FUN = bufferApply, ...)
    } else {
        res <- parallel::mclapply(X = seq_len(nBuffers), FUN = bufferApply, ..., mc.preschedule = TRUE, mc.cores = nCores)
    }
    simplifyList(res)
}


# A more memory-efficient version of apply.
#
# apply always makes a copy of the data.
apply2 <- function(X, MARGIN, FUN, ...) {
    d <- dim(X)
    if (MARGIN == 1L) {
        subset <- X[1L, ]
    } else {
        subset <- X[, 1L]
    }
    sample <- FUN(subset, ...)
    if (is.table(sample)) {
        stop("tables are not supported.")
    } else if (is.list(sample)) {
        # List
        OUT <- vector(mode = "list", length = d[MARGIN])
        names(OUT) <- dimnames(X)[[MARGIN]]
        OUT[[1L]] <- sample
        if (d[MARGIN] > 1L) {
            for (i in seq(2L, d[MARGIN])) {
                if (MARGIN == 1L) {
                    subset <- X[i, ]
                } else {
                    subset <- X[, i]
                }
                OUT[[i]] <- FUN(subset, ...)
            }
        }
    } else {
        if (length(sample) > 1L) {
            # Matrix or atomic vector of length > 1
            OUT <- matrix(data = normalizeType(typeof(sample)), nrow = length(sample), ncol = d[MARGIN])
            if (!is.matrix(sample) && !is.null(names(sample))) {
                if (MARGIN == 1L) {
                    dimnames(OUT) <- list(NULL, names(sample))
                } else {
                    dimnames(OUT) <- list(names(sample), NULL)
                }
            }
            OUT[, 1L] <- sample
            if (d[MARGIN] > 1L) {
                for (i in seq(2L, d[MARGIN])) {
                    if (MARGIN == 1L) {
                        subset <- X[i, ]
                    } else {
                        subset <- X[, i]
                    }
                    OUT[, i] <- FUN(subset, ...)
                }
            }
        } else {
            # Atomic vector of length 1
            OUT <- vector(mode = typeof(sample), length = d[MARGIN])
            names(OUT) <- dimnames(X)[[MARGIN]]
            OUT[1L] <- sample
            if (d[MARGIN] > 1L) {
                for (i in seq(2L, d[MARGIN])) {
                    if (MARGIN == 1L) {
                        subset <- X[i, ]
                    } else {
                        subset <- X[, i]
                    }
                    OUT[i] <- FUN(subset, ...)
                }
            }
        }
    }
    return(OUT)
}


simplifyList <- function(x) {
    sample <- x[[1L]]
    if (is.matrix(sample)) {
        x <- matrix(data = unlist(x), nrow = nrow(sample), byrow = FALSE)
        rownames(x) <- rownames(sample)
    } else {
        x <- unlist(x)
    }
    return(x)
}

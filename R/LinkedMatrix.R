apply.chunk <- function(W, fn, MARGIN, chunks, chunk, ...) {
    W <- W[, chunks == chunk, drop = FALSE]
    ans <- apply(FUN = fn, X = W, MARGIN = MARGIN, ...)
    return(ans)
}


apply.parallel <- function(X, FUN, MARGIN, nChunks = detectCores(), mc.cores = detectCores(), ...) {
    n <- ifelse(MARGIN == 1, nrow(X), ncol(X))
    chunks <- rep(1:nChunks, each = ceiling(n/nChunks))[1:n]
    TMP <- mclapply(X = 1:nChunks, FUN = apply.chunk, W = X, MARGIN = MARGIN, fn = FUN, chunks = chunks, mc.cores = mc.cores, ...)
    if (is.atomic(TMP[[1]])) {
        if (is.vector(TMP[[1]])) {
            ANS <- unlist(TMP)
            return(ANS)
        }
        if (is.matrix(TMP[[1]])) {
            ANS <- TMP[[1]]
            for (i in 2:length(TMP)) {
                ANS <- cbind(ANS, TMP[[i]])
            }
            return(ANS)
        }
    } else {
        return(TMP)
    }
}


apply.LinkedMatrix <- function(X, MARGIN, FUN, chunkSize = 1000, verbose = FALSE, ...) {
    FUN <- match.fun(FUN)
    n <- ifelse(MARGIN == 1, nrow(X), ncol(X))
    # Apply function on first element to check output type.
    x <- ifelse(MARGIN == 1, X[1, ], X[, 1])
    outputType <- FUN(x, ...)
    if (is.atomic(outputType)) {
        ANS <- matrix(nrow = length(outputType), ncol = n, NA)
        rownames(ANS) <- names(outputType)
        if (MARGIN == 1) {
            colnames(ANS) <- rownames(X)
        } else {
            colnames(ANS) <- colnames(X)
        }
        nChunks <- ceiling(n/chunkSize)
        end <- 0
        for (i in 1:nChunks) {
            if (verbose) {
                cat(i, " out of ", nChunks, " \n")
            }
            ini <- end + 1
            end <- min(ini + chunkSize - 1, n)
            if (MARGIN == 1) {
                Z <- X[ini:end, ]
            } else {
                Z <- X[, ini:end]
            }
            ANS[, ini:end] <- apply(Z, MARGIN, FUN, ...)
        }
    } else {
        ANS <- vector("list", n)
        names(ANS) <- ifelse(MARGIN == 1, rownames(X), colnames(X))
        end <- 0
        for (i in 1:n) {
            if (verbose) {
                cat(i, " out of ", n, " \n")
            }
            if (MARGIN == 1) {
                ANS[[i]] <- FUN(X[i, ], ...)
            } else {
                ANS[[i]] <- FUN(X[, i], ...)
            }
        }
    }
    return(ANS[, , drop = TRUE])
}


colMeans.LinkedMatrix <- function(x, na.rm = TRUE, chunkSize = 1000, ...) {
    if (na.rm) {
        warning("Ignoring missing values")
    }
    ANS <- apply.LinkedMatrix(X = x, MARGIN = 2, FUN = mean, chunkSize = chunkSize, na.rm = na.rm, ...)
    return(ANS)
}


colSums.LinkedMatrix <- function(x, na.rm = TRUE, chunkSize = 1000, ...) {
    if (na.rm) {
        warning("Ignoring missing values")
    }
    ANS <- apply.LinkedMatrix(X = x, MARGIN = 2, FUN = sum, chunkSize = chunkSize, na.rm = na.rm, ...)
    return(ANS)
}


rowMeans.LinkedMatrix <- function(x, na.rm = TRUE, chunkSize = 1000, ...) {
    if (na.rm) {
        warning("Ignoring missing values")
    }
    ANS <- apply.LinkedMatrix(X = x, MARGIN = 1, FUN = mean, chunkSize = chunkSize, na.rm = na.rm, ...)
    return(ANS)
}


rowSums.LinkedMatrix <- function(x, na.rm = TRUE, chunkSize = 1000, ...) {
    if (na.rm) {
        warning("Ignoring missing values")
    }
    ANS <- apply.LinkedMatrix(X = x, MARGIN = 1, FUN = sum, chunkSize = chunkSize, na.rm = na.rm, ...)
    return(ANS)
}


summary.num <- function(x) {
    out <- c(range(x, na.rm = T), mean(x, na.rm = T), stats::sd(x, na.rm = T), mean(is.na(x)))
    names(out) <- c("min", "max", "mean", "sd", "prop NAs")
    return(out)
}


summary.char <- function(x) {
    out <- table(x, useNA = "always")
    out <- out/length(x)
    return(out)
}


summary.LinkedMatrix <- function(object, MARGIN = 2, chunkSize = 1000, ...) {
    sample <- object[1, 1]
    if (is.numeric(sample)) {
        fun <- summary.num
    } else if (is.character(sample) | is.logical(sample)) {
        fun <- summary.char
    } else {
        fun <- summary
    }
    apply.LinkedMatrix(X = object, MARGIN = MARGIN, FUN = fun, chunkSize = chunkSize, ...)
}


#' Apply function for \code{\linkS4class{ColumnLinkedMatrix}} or 
#' \code{\linkS4class{RowLinkedMatrix}} objects.
#' 
#' This function brings chunks (of size \code{chunkSize}) of rows (if 
#' \code{MARGIN} is 1) or columns (if \code{MARGIN} is 2) of the 
#' \code{\linkS4class{LinkedMatrix}} instance into RAM as \code{matrix} objects
#' and calls the \code{apply} function of the base package for each chunk.
#' Results from all the chunks are collected and returned.
#' 
#' @param X Either a \code{\linkS4class{ColumnLinkedMatrix}} or a 
#'   \code{\linkS4class{RowLinkedMatrix}} object.
#' @param MARGIN Use 1 to apply function over rows or 2 to apply function over
#'   columns.
#' @param FUN The function to be applied.
#' @param chunkSize The number of columns or rows that are processed at a time 
#'   (see Details).
#' @param verbose Whether to print additional information.
#' @param ... Optional arguments to FUN.
#' @return Returns a \code{matrix} or a \code{list} with results from FUN.
#' @export
setMethod("apply", signature("LinkedMatrix"), apply.LinkedMatrix)


#' Form column means.
#' 
#' @inheritParams base::colMeans
#' @param chunkSize The number of columns that are processed at a time.
#' @param ... Optional arguments to \code{mean}.
#' @export
setMethod("colMeans", signature("LinkedMatrix"), colMeans.LinkedMatrix)


#' Form column sums.
#' 
#' @inheritParams base::colSums
#' @param chunkSize The number of columns that are processed at a time.
#' @param ... Optional arguments to \code{mean}.
#' @export
setMethod("colSums", signature("LinkedMatrix"), colSums.LinkedMatrix)


#' Form row means.
#' 
#' @inheritParams base::rowMeans
#' @param chunkSize The number of rows that are processed at a time.
#' @param ... Optional arguments to \code{mean}.
#' @export
setMethod("rowMeans", signature("LinkedMatrix"), rowMeans.LinkedMatrix)


#' Form row sums.
#' 
#' @inheritParams base::rowSums
#' @param chunkSize The number of rows that are processed at a time.
#' @param ... Optional arguments to \code{sum}.
#' @export
setMethod("rowSums", signature("LinkedMatrix"), rowSums.LinkedMatrix)


#' Summary function for \code{\linkS4class{ColumnLinkedMatrix}} or 
#' \code{\linkS4class{RowLinkedMatrix}} objects.
#' 
#' This function brings chunks (of size \code{chunkSize}) of rows (if 
#' \code{MARGIN} is 1) or columns (if \code{MARGIN} is 2) of the 
#' \code{\linkS4class{LinkedMatrix}} instance into RAM as \code{matrix} objects
#' and calls an appropriate summary function based on the type of the matrix for
#' each chunk. Results from all the chunks are collected and returned.
#' 
#' @param object Either a \code{\linkS4class{ColumnLinkedMatrix}} or a 
#'   \code{\linkS4class{RowLinkedMatrix}} object.
#' @param MARGIN Use 1 to obtain row summaries or 2 to obtain column summaries.
#' @param chunkSize The number of rows or columns that are processed at a time 
#'   (see Details).
#' @param ... Optional arguments to summary functions.
#' @return Returns a \code{matrix} of summaries.
#' @export
setMethod("summary", signature("LinkedMatrix"), summary.LinkedMatrix)

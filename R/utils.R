# A more memory-efficient version of apply.
#
# apply always makes a copy of the data.
apply2 <- function(X, MARGIN, FUN, ...) {
    d <- dim(X)
    if (MARGIN == 1) {
        subset <- X[1, ]
    } else {
        subset <- X[, 1]
    }
    sample <- FUN(subset, ...)
    if (is.table(sample)) {
        stop("tables are not supported.")
    } else if (is.list(sample)) {
        # List
        OUT <- vector(mode = "list", length = d[MARGIN])
        names(OUT) <- dimnames(X)[[MARGIN]]
        OUT[[1]] <- sample
        if (d[MARGIN] > 1) {
            for (i in seq(2, d[MARGIN])) {
                if (MARGIN == 1) {
                    subset <- X[i, ]
                } else {
                    subset <- X[, i]
                }
                OUT[[i]] <- FUN(subset, ...)
            }
        }
    } else {
        if (length(sample) > 1) {
            # Matrix or atomic vector of length > 1
            OUT <- matrix(data = normalizeType(typeof(sample)), nrow = length(sample), ncol = d[MARGIN])
            if (!is.matrix(sample) && !is.null(names(sample))) {
                if (MARGIN == 1) {
                    dimnames(OUT) <- list(NULL, names(sample))
                } else {
                    dimnames(OUT) <- list(names(sample), NULL)
                }
            }
            OUT[, 1] <- sample
            if (d[MARGIN] > 1) {
                for (i in seq(2, d[MARGIN])) {
                    if (MARGIN == 1) {
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
            OUT[1] <- sample
            if (d[MARGIN] > 1) {
                for (i in seq(2, d[MARGIN])) {
                    if (MARGIN == 1) {
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


#' Applies a Function on Each Row or Column of a Matrix in Parallel.
#'
#' Similar to [base::apply()], but designed to carry out operations in
#' parallel.  The input matrix `X` is broken into `nTasks` chunks and passed to
#' [parallel::mclapply()] which applies `FUN` on either the rows or the columns
#' of each chunk.
#'
#' `nTasks` has to be chosen carefully to avoid running out of memory. As a
#' rule of thumb, at least around \code{object_size(X) + (nCores *
#' (object_size(X) / nTasks)) + object_size(result)} MB of total memory will be
#' needed, not including potential copies of your data that might be created
#' (for example [stats::lsfit()] runs `cbind(1, X)`). Therefore, for 20 nodes
#' and 20 tasks you will need at least `2 * object_size(X)` MB, for 20 nodes
#' and 40 tasks `1.5 * object_size(X)` MB, etc.
#'
#' If `nTasks` is `1`, [base::apply()] will be called directly without
#' parallelism.
#'
#' @param X A matrix or matrix-like object.
#' @param MARGIN The subscripts which the function will be applied over. `1`
#' indicates rows, `2` indicates columns.
#' @param FUN The function to be applied.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param ... Additional arguments to be passed to [base::apply()].
#' @seealso [chunkedApply()] if `X` is a memory-mapped matrix and too large to
#' hold in memory.
#' @export
parallelApply <- function(X, MARGIN, FUN, nTasks = nCores, nCores = parallel::detectCores(), ...) {
    d <- dim(X)
    if (!length(d)) {
        stop("dim(X) must have a positive length")
    }
    nTasks <- as.integer(nTasks)
    if (is.na(nTasks) || nTasks < 1) {
        stop("nTasks has to be greater than 0")
    }
    if (nTasks == 1) {
        apply2(X = X, MARGIN = MARGIN, FUN = FUN, ...)
    } else {
        res <- parallel::mclapply(X = seq_len(nTasks), FUN = function(i, ...) {
            range <- LinkedMatrix:::chunkRanges(d[MARGIN], nTasks, i)
            if (MARGIN == 2) {
                subset <- X[, seq(range[1], range[2]), drop = FALSE]
            } else {
                subset <- X[seq(range[1], range[2]), , drop = FALSE]
            }
            apply2(X = subset, MARGIN = MARGIN, FUN = FUN, ...)
        }, ..., mc.preschedule = FALSE, mc.cores = nCores)
        simplifyList(res)
    }
}


#' Reads Chunks of Data from a Memory-Mapped File into Memory and Applies a
#' Function on Each Row or Column of a Matrix in Parallel.
#'
#' Similar to [base::apply()], but designed to bring chunks of data into memory
#' and carry out operations on them in parallel. `nBufferSize` rows or columns
#' of the input matrix `X` are read into memory and handed over to
#' [parallelApply()]. This function is only useful for memory-mapped files. For
#' data that is already in memory, use [parallelApply()] directly.
#'
#' `bufferSize` and `nTasks` have to be chosen carefully to avoid running out
#' of memory. As a rule of thumb, at least around \code{object_size(buffer) +
#' (nCores * (object_size(buffer) / nTasks)) + object_size(result)} MB of total
#' memory will be needed, not including potential copies of your data that
#' might be created (for example [stats::lsfit()] runs `cbind(1, X)`).
#' Therefore, for 20 nodes and 20 tasks you will need at least `2 *
#' object_size(buffer)` MB, for 20 nodes and 40 tasks `1.5 *
#' object_size(buffer)` MB, etc.
#'
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param MARGIN The subscripts which the function will be applied over. 1
#' indicates rows, 2 indicates columns.
#' @param FUN The function to be applied.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param bufferSize The number of rows or columns of `X` that are brought into
#' RAM for processing. Overwrites `nBuffers`. If both parameters are `NULL`,
#' all elements in `i` or `j` are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the rows or columns of `X` that
#' are brought into RAM for processing. Is overwritten by `bufferSize`. If both
#' parameters are `NULL`, all elements in `i` or `j` are used. Defaults to
#' `NULL`.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @param ... Additional arguments to be passed to [parallelApply()].
#' @seealso [parallelApply()] if `X` is not a memory-mapped matrix or can be
#' held in memory.
#' @export
chunkedApply <- function(X, MARGIN, FUN, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = FALSE, ...) {
    if (!length(dim(X))) {
        stop("dim(X) must have a positive length")
    }
    # Convert index types
    if (is.logical(i)) {
        i <- which(i)
    } else if (is.character(i)) {
        i <- match(i, rownames(X))
    }
    if (is.logical(j)) {
        j <- which(j)
    } else if (is.character(j)) {
        j <- match(j, colnames(X))
    }
    d <- c(length(i), length(j))
    if (is.null(bufferSize) && is.null(nBuffers)) {
        bufferSize <- d[MARGIN]
        nBuffers <- 1
    } else if (is.null(bufferSize) && !is.null(nBuffers)) {
        bufferSize <- ceiling(d[MARGIN] / nBuffers)
    } else {
        nBuffers <- ceiling(d[MARGIN] / bufferSize)
    }
    nBuffers <- ceiling(d[MARGIN] / bufferSize)
    ranges <- LinkedMatrix:::chunkRanges(d[MARGIN], nBuffers)
    res <- lapply(seq_len(nBuffers), function(k) {
        if (verbose) {
            message("Processing chunk ", k, " of ", nBuffers, " (", round(k / nBuffers * 100, 3), "%) ...")
        }
        if (MARGIN == 2) {
            subset <- X[i, j[seq(ranges[1, k], ranges[2, k])], drop = FALSE]
        } else {
            subset <- X[i[seq(ranges[1, k], ranges[2, k])], j, drop = FALSE]
        }
        parallelApply(X = subset, MARGIN = MARGIN, FUN = FUN, nTasks = nTasks, nCores = nCores, ...)
    })
    simplifyList(res)
}


crossprods <- function(x, y = NULL, use_tcrossprod = FALSE, nTasks = nCores, nCores = parallel::detectCores()) {
    dx <- dim(x)
    if (!is.null(y)) {
        y <- as.matrix(y)
        dy <- dim(y)
        if (use_tcrossprod) {
            if (dx[2] != ncol(y)) {
                stop("Error in tcrossprod_parallel: non-conformable arguments.")
            }
        } else {
            if (dx[1] != dy[1]) {
                stop("Error in crossprod_parallel: non-conformable arguments.")
            }
        }
    }
    if (nTasks == 1) {
        if (use_tcrossprod) {
            Xy <- tcrossprod(x, y)
        } else {
            Xy <- crossprod(x, y)
        }
    } else {
        nX <- ifelse(use_tcrossprod, dx[2], dx[1])
        if (!is.null(y)) {
            nY <- ifelse(use_tcrossprod, dy[2], dy[1])
        }
        chunks <- parallel::mclapply(X = seq_len(nTasks), FUN = function(i, ...) {
            if (!is.null(y)) {
                ranges <- LinkedMatrix:::chunkRanges(nY, nTasks, i)
                if (use_tcrossprod) {
                    Y <- y[, seq(ranges[1], ranges[2]), drop = FALSE]
                } else {
                    Y <- y[seq(ranges[1], ranges[2]), , drop = FALSE]
                }
            } else {
                Y <- NULL
            }
            ranges <- LinkedMatrix:::chunkRanges(nX, nTasks, i)
            if (use_tcrossprod) {
                X <- x[, seq(ranges[1], ranges[2]), drop = FALSE]
                Xy <- tcrossprod(X, Y)
            } else {
                X <- x[seq(ranges[1], ranges[2]), , drop = FALSE]
                Xy <- crossprod(X, Y)
            }
            return(Xy)
        }, mc.preschedule = FALSE, mc.cores = nCores)
        # We now need to add up chunks sequentially
        Xy <- chunks[[1]]
        if (length(chunks) > 1) {
            for (i in 2:length(chunks)) {
                Xy <- Xy + chunks[[i]]
            }
        }
    }
    return(Xy)
}


#' Computes crossprod (x'x or x'y) or tcrossprod (xx' or xy') in Parallel.
#'
#' Similar to [base::crossprod()] and [base::tcrossprod()], but designed to
#' carry out operations in parallel. The input matrix `x` (and `y` if not
#' `NULL`) is broken into `nTasks` chunks and passed to [parallel::mclapply()]
#' which performs [base::crossprod()] or [base::tcrossprod()] on each chunk.
#' The results are added up and returned.
#'
#' If `nTasks` is `1`, [base::crossprod()] or [base::tcrossprod()] will be
#' called directly without parallelism.
#'
#' @param x A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param y vector or matrix-like object. `NULL` by default.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @return x'x or x'y (`crossprod_parallel`), or xx' or xy'
#' (`tcrossprod_parallel`), depending on whether `y` is provided.
#' @seealso [getG()] to compute a genomic relationship matrix.
#' @export
crossprod_parallel <- function(x, y = NULL, nTasks = nCores, nCores = parallel::detectCores()) {
    crossprods(x = x, y = y, use_tcrossprod = FALSE, nTasks = nTasks, nCores = nCores)
}


#' @rdname crossprod_parallel
#' @export
tcrossprod_parallel <- function(x, y = NULL, nTasks = nCores, nCores = parallel::detectCores()) {
    crossprods(x = x, y = y, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
}



#' Computes a Genomic Relationship Matrix.
#'
#' Computes a positive semi-definite symmetric genomic relation matrix G=XX'
#' offering options for centering and scaling the columns of `X` beforehand.
#'
#' If `center = FALSE`, `scale = FALSE` and `scaleG = FALSE`, [getG()] produces
#' the same outcome than [base::tcrossprod()].
#'
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param center Either a logical value or a numeric vector of length equal to
#' the number of columns of `X`. Numeric vector required if `i2` is used. If
#' `FALSE`, no centering is done. Defaults to `TRUE`.
#' @param scale Either a logical value or a numeric vector of length equal to
#' the number of columns of `X`. Numeric vector required if `i2` is used. If
#' `FALSE`, no scaling is done.  Defaults to `TRUE`.
#' @param scaleG Whether XX' should be scaled. Defaults to `TRUE`.
#' @param minVar Columns with variance lower than this value will not be used
#' in the computation (only if `scale` is not `FALSE`).
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param i2 Indicates which rows should be used to compute a block of the
#' genomic relationship matrix. Will compute XY' where X is determined by `i`
#' and `j` and Y by `i2` and `j`. Can be integer, boolean, or character. If
#' `NULL`, the whole genomic relationship matrix XX' is computed.  Defaults to
#' `NULL`.
#' @param bufferSize The number of columns of `X` that are brought into RAM for
#' processing. Overwrites `nBuffers`. If both parameters are `NULL`, all
#' columns of `X` are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the columns of `X` that are
#' brought into RAM for processing. Is overwritten by `bufferSize`. If both
#' parameters are `NULL`, all columns of `X` are used. Defaults to `NULL`.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `TRUE`.
#' @return A positive semi-definite symmetric numeric matrix.
#' @export
getG <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05, i = seq_len(nrow(X)), j = seq_len(ncol(X)), i2 = NULL, bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = TRUE) {

    # compute XY' rather than XX'
    hasY <- !is.null(i2)

    if (hasY) {
        if (!is.numeric(center)) {
            stop("center need to be precomputed.")
        }
        if (!is.numeric(scale)) {
            stop("scale need to be precomputed.")
        }
    }

    if (is.null(bufferSize) && is.null(nBuffers)) {
        bufferSize <- length(j)
        nBuffers <- 1
    } else if (is.null(bufferSize) && !is.null(nBuffers)) {
        bufferSize <- ceiling(length(j) / nBuffers)
    } else {
        nBuffers <- ceiling(length(j) / bufferSize)
    }

    # Convert index types
    if (is.logical(i)) {
        i <- which(i)
    } else if (is.character(i)) {
        i <- match(i, rownames(X))
    }
    if (is.logical(j)) {
        j <- which(j)
    } else if (is.character(j)) {
        j <- match(j, colnames(X))
    }
    if (hasY) {
        if (is.logical(i2)) {
            i2 <- which(i2)
        } else if (is.character(i2)) {
            i2 <- match(i2, rownames(X))
        }
    }

    nX <- nrow(X)
    pX <- ncol(X)

    if ((min(i) < 1) || (max(i) > nX)) {
        stop("Index out of bounds")
    }
    if ((min(j) < 1) || (max(j) > pX)) {
        stop("Index out of bounds")
    }
    if (hasY) {
        if ((min(i2) < 1) || (max(i2) > nX)) {
            stop("Index out of bounds")
        }
    }

    n <- length(i)
    p <- length(j)
    if (hasY) {
        n2 <- length(i2)
    }

    if (!hasY) {
        G <- matrix(data = 0, nrow = n, ncol = n, dimnames = list(rownames(X)[i], rownames(X)[i]))
    } else {
        G <- matrix(data = 0, nrow = n, ncol = n2, dimnames = list(rownames(X)[i], rownames(X)[i2]))
        K <- 0
    }

    end <- 0
    for (k in seq_len(nBuffers)) {
        ini <- end + 1
        end <- min(p, ini + bufferSize - 1)

        if (verbose) {
            message("Chunk: ", k, " (markers ", ini, ":", end, " ~", round(100 * end / p, 1), "% done)")
            message("  => Acquiring genotypes...")
        }

        # subset
        localColIndex <- j[ini:end]
        X1 <- X[i, localColIndex, drop = FALSE]
        if (hasY) {
            X2 <- X[i2, localColIndex, drop = FALSE]
        }

        # compute centers
        if (is.logical(center) && center == TRUE) {
            center.chunk <- colMeans(X1, na.rm = TRUE)
        } else if (is.numeric(center)) {
            center.chunk <- center[localColIndex]
        } else {
            center.chunk = FALSE
        }

        # compute scales
        if (is.logical(scale) && scale == TRUE) {
            scale.chunk <- apply(X = X1, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
        } else if (is.numeric(scale)) {
            scale.chunk <- scale[localColIndex]
        } else {
            scale.chunk <- FALSE
        }

        # remove constant columns
        if (!(is.logical(scale) && scale == FALSE)) {
            removeCols <- which(scale.chunk < minVar)
            if (length(removeCols) > 0) {
                X1 <- X1[, -removeCols]
                if (hasY) {
                    X2 <- X2[, -removeCols]
                }
                scale.chunk <- scale.chunk[-removeCols]
                center.chunk <- center.chunk[-removeCols]
            }
        }

        # compute XX'
        if (ncol(X1) > 0) {

            if (verbose) {
              message("  => Computing...")
            }

            # scale and impute X
            X1 <- scale(X1, center = center.chunk, scale = scale.chunk)
            X1[is.na(X1)] <- 0
            if (hasY) {
                X2 <- scale(X2, center = center.chunk, scale = scale.chunk)
                X2[is.na(X2)] <- 0
            }

            if (nTasks > 1) {
                if (!hasY) {
                    G_chunk <- crossprods(x = X1, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
                } else {
                    G_chunk <- crossprods(x = X1, y = X2, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
                }
            } else {
                if (!hasY) {
                    G_chunk <- tcrossprod(X1)
                } else {
                    G_chunk <- tcrossprod(x = X1, y = X2)
                }
            }

            G[] <- G + G_chunk

            if (hasY && scaleG) {
                if (is.logical(scale) && scale == FALSE) {
                    K <- K + sum(scale.chunk^2)
                } else {
                    K <- K + ncol(X1)
                }
            }

        }

    }

    if (!hasY && scaleG) {
        # Use seq instead of diag to avoid copy as it does not increase ref count
        K <- mean(G[seq(from = 1, to = n * n, by = n + 1)])
    }

    if (scaleG) {
        G[] <- G / K
    }

    return(G)

}


#' Computes a Genomic Relationship Matrix G=xx' without Ever Loading G in RAM
#' by Creating a symDMatrix.
#'
#' Offers options for centering and scaling the columns of x before computing
#' xx'.
#'
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param nBlocks The number of blocks.
#' @param blockSize The number of columns of a block (if NULL inferred from block).
#' @param centers Precomputed centers.
#' @param scales Precomputed scales.
#' @param centerCol TRUE/FALSE whether columns must be centered before
#' computing xx'.
#' @param scaleCol TRUE/FALSE whether columns must be scaled before computing
#' xx'.
#' @param scaleG TRUE/FALSE whether xx' must be scaled.
#' @param folderOut The path to the folder where to save the
#' [symDMatrix::symDMatrix-class] object. Defaults to a random string prefixed
#' with "symDMatrix_".
#' @param vmode vmode of `ff` objects.
#' @param saveRData Whether to save an RData file to easily reload the
#' [symDMatrix::symDMatrix-class] object.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `TRUE`.
#' @return A positive semi-definite symmetric numeric matrix.
#' @export
getG.symDMatrix <- function(X, nBlocks = 5, blockSize = NULL, centers = NULL, scales = NULL, centerCol = TRUE, scaleCol = TRUE, scaleG = TRUE, folderOut = paste0("symDMatrix_", randomString()), vmode = "double", saveRData = TRUE, i = seq_len(nrow(X)), j = seq_len(ncol(X)), nTasks = nCores, nCores = parallel::detectCores(), verbose = TRUE) {

    nX <- nrow(X)
    pX <- ncol(X)

    # Convert index types
    if (is.logical(i)) {
        i <- which(i)
    } else if (is.character(i)) {
        i <- match(i, rownames(X))
    }
    if (is.logical(j)) {
        j <- which(j)
    } else if (is.character(j)) {
        j <- match(j, colnames(X))
    }

    n <- length(i)
    p <- length(j)

    if ((min(i) < 1) || (max(i) > nX)) {
        stop("Index out of bounds")
    }
    if ((min(j) < 1) || (max(j) > pX)) {
        stop("Index out of bounds")
    }

    if ((centerCol || scaleCol) && (is.null(centers) || is.null(scales))) {
        if (is.null(centers) && is.null(scales)) {
            centers <- vector(mode = "double", length = p)
            scales <- vector(mode = "double", length = p)
            for (k in seq_len(p)) {
                xi <- X[i, j[k]]
                scales[k] <- stats::sd(xi, na.rm = TRUE) * sqrt((nX - 1) / nX)
                centers[k] <- mean(xi, na.rm = TRUE)
            }
        } else if ((!is.null(centers)) && (is.null(scales))) {
            scales <- vector(mode = "double", length = p)
            for (k in seq_len(p)) {
                xi <- X[i, j[k]]
                scales[k] <- stats::sd(xi, na.rm = TRUE) * sqrt((nX - 1) / nX)
            }
        } else if ((is.null(centers)) && (!is.null(scales))) {
            centers <- vector(mode = "double", length = p)
            for (k in seq_len(p)) {
                xi <- X[i, j[k]]
                centers[k] <- mean(xi, na.rm = TRUE)
            }
        }
    }
    if (!centerCol) {
        centers <- rep(0, p)
    }
    if (!scaleCol) {
        scales <- rep(1, p)
    }

    if (is.null(blockSize)) {
        blockSize <- ceiling(n / nBlocks)
    }
    blockIndex <- cbind(i, ceiling(seq_len(n) / blockSize))

    nFiles <- nBlocks * (nBlocks + 1) / 2
    DATA <- vector(mode = "list", length = nBlocks)

    if (file.exists(folderOut)) {
        stop(folderOut, " already exists")
    }
    curDir <- getwd()
    dir.create(folderOut)
    setwd(folderOut)

    counter <- 1
    for (r in seq_len(nBlocks)) {

        DATA[[r]] <- vector(mode = "list", length = nBlocks - r)

        rowIndex_r <- blockIndex[which(blockIndex[, 2] == r), 1]
        Xi <- X[rowIndex_r, j, drop = FALSE]

        # centering/scaling
        for (k in seq_len(p)) {
            xik <- Xi[, k]
            xik <- (xik - centers[j[k]]) / scales[j[k]]
            xik[is.na(xik)] <- 0
            Xi[, k] <- xik
        }

        for (s in r:nBlocks) {

            if (verbose) {
                message("Working on block ", r, "-", s, " (", round(100 * counter / (nBlocks * (nBlocks + 1) / 2)), "%)")
            }

            rowIndex_s <- blockIndex[which(blockIndex[, 2] == s), 1]
            Xj <- X[rowIndex_s, j, drop = FALSE]

            # centering/scaling
            for (k in seq_len(p)) {
                xjk <- Xj[, k]
                xjk <- (xjk - centers[j[k]]) / scales[j[k]]
                xjk[is.na(xjk)] <- 0
                Xj[, k] <- xjk
            }

            Gij <- tcrossprod_parallel(x = Xi, y = Xj, nTasks = nTasks, nCores = nCores)

            blockName <- paste0("data_", padDigits(r, nBlocks), "_", padDigits(s, nBlocks), ".bin")
            block <- ff::ff(dim = dim(Gij), vmode = vmode, initdata = as.vector(Gij), filename = blockName, dimnames = list(rownames(X)[rowIndex_r], rownames(X)[rowIndex_s]))
            # Change ff path to a relative one
            bit::physical(block)$pattern <- "ff"
            bit::physical(block)$filename <- blockName
            DATA[[r]][[s - r + 1]] <- block

            counter <- counter + 1

            if (verbose) {
                message("  => Done")
            }
        }
    }

    names(centers) <- colnames(X)[j]
    names(scales) <- colnames(X)[j]

    G <- new("symDMatrix", data = DATA, centers = centers, scales = scales)

    if (scaleG) {
        K <- mean(diag(G))
        for (r in seq_len(length(G@data))) {
            for (s in seq_len(length(G@data[[r]]))) {
                G@data[[r]][[s]][] <- G@data[[r]][[s]][] / K
            }
        }
    }

    if (saveRData) {
        save(G, file = "G.RData")
    }

    setwd(curDir)

    return(G)
}


#' Performs Single Marker Regressions Using BGData Objects.
#'
#' Implements single marker regressions. The regression model includes all the
#' covariates specified in the right-hand-side of the `formula` plus one column
#' of `@@geno` at a time. The data from the association tests is obtained from
#' a [BGData-class] object.
#'
#' @param formula The formula for the GWAS model without including the marker,
#' e.g. `y ~ 1` or `y ~ factor(sex) + age`. The variables included in the
#' formula must be in the `@@pheno` object of the [BGData-class].
#' @param data A [BGData-class] object.
#' @param method The regression method to be used. Currently, the following
#' methods are implemented: [stats::lm()], [stats::lm.fit()], [stats::lsfit()],
#' [stats::glm()], [lme4::lmer()], and [SKAT::SKAT()]. Defaults to `lsfit`.
#' @param i Indicates which rows of `@@geno` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `@@geno` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param bufferSize The number of columns of `@@geno` that are brought into
#' RAM for processing. Overwrites `nBuffers`. If both parameters are `NULL`,
#' all elements in `j` are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the columns of `@@geno` that are
#' brought into RAM for processing. Is overwritten by `bufferSize`. If both
#' parameters are `NULL`, all elements in `j` are used. Defaults to `NULL`.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @param ... Additional arguments for chunkedApply and regression method.
#' @return The same matrix that would be returned by `coef(summary(model))`.
#' @export
GWAS <- function(formula, data, method = "lsfit", i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = FALSE, ...) {

    if (class(data) != "BGData") {
        stop("data must BGData")
    }

    if (!method %in% c("lm", "lm.fit", "lsfit", "glm", "lmer", "SKAT")) {
        stop("Only lm, lm.fit, lsfit, glm, lmer and SKAT have been implemented so far.")
    }

    if (method == "lsfit") {
        OUT <- GWAS.lsfit(formula = formula, data = data, i = i, j = j, bufferSize = bufferSize, nTasks = nTasks, nCores = nCores, verbose = verbose, ...)
    } else if (method == "SKAT") {
        OUT <- GWAS.SKAT(formula = formula, data = data, i = i, j = j, verbose = verbose, ...)
    } else {
        if (method == "lmer") {
            if (!requireNamespace("lme4", quietly = TRUE)) {
                stop("lme4 needed for this function to work. Please install it.", call. = FALSE)
            }
            FUN <- lme4::lmer
        } else {
            FUN <- match.fun(method)
        }
        pheno <- data@pheno
        GWAS.model <- stats::update(stats::as.formula(formula), ".~z+.")
        OUT <- chunkedApply(X = data@geno, MARGIN = 2, FUN = function(col, ...) {
            pheno$z <- col
            fm <- FUN(GWAS.model, data = pheno, ...)
            getCoefficients(fm)
        }, i = i, j = j, bufferSize = bufferSize, nBuffers = nBuffers, nTasks = nTasks, nCores = nCores, verbose = verbose, ...)
        colnames(OUT) <- colnames(data@geno)[j]
        OUT <- t(OUT)
    }

    return(OUT)
}


# OLS for the regression y=xb+e (data is assumed to be pre-adjusted by non-genetic effects
rayOLS <- function(y, x, n = length(y)){
    tmp <- !(is.na(x) | is.na(y))
    x <- x[tmp]
    y <- y[tmp]
    rhs <- sum(x * y)
    XtX <- sum(x^2)
    sol <- rhs / XtX
    error <- y - x * sol
    vE <- sum(error^2) / (n - 1)
    SE <- sqrt(vE / XtX)
    z_stat <- sol / SE
    return(c(sol, SE, z_stat, stats::pt(q = abs(z_stat), df = n - 1, lower.tail = FALSE) * 2))
}


GWAS.lsfit <- function(formula, data, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = FALSE, ...) {

    # subset of model.frame has bizarre scoping issues
    frame <- stats::model.frame(formula = formula, data = data@pheno)[i, , drop = FALSE]
    model <- stats::model.matrix(formula, frame)
    model <- cbind(1, model) # Reserve space for marker column

    y <- data@pheno[i, as.character(stats::terms(formula)[[2]]), drop = TRUE]

    res <- chunkedApply(X = data@geno, MARGIN = 2, FUN = function(col, ...) {
        model[, 1] <- col
        fm <- stats::lsfit(x = model, y = y, intercept = FALSE)
        stats::ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
    }, i = i, j = j, bufferSize = bufferSize, nBuffers = nBuffers, nTasks = nTasks, nCores = nCores, verbose = verbose, ...)
    colnames(res) <- colnames(data@geno)[j]
    res <- t(res)

    return(res)
}


# formula: the formula for the GWAS model without including the markers, e.g.
# y~1 or y~factor(sex)+age
# all the variables in the formula must be in data@pheno data (BGData)
# containing slots @pheno and @geno
# groups: a vector mapping markers into groups (can be integer, character or
# factor)
GWAS.SKAT <- function(formula, data, groups, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), verbose = FALSE, ...) {

    if (!requireNamespace("SKAT", quietly = TRUE)) {
        stop("SKAT needed for this function to work. Please install it.", call. = FALSE)
    }

    uniqueGroups <- unique(groups)

    OUT <- matrix(data = double(), nrow = length(uniqueGroups), ncol = 2)
    colnames(OUT) <- c("nMrk", "p-value")
    rownames(OUT) <- uniqueGroups

    H0 <- SKAT::SKAT_Null_Model(formula, data = data@pheno[i, , drop = FALSE], ...)

    for (group in seq_len(length(uniqueGroups))) {
        Z <- data@geno[i, groups == uniqueGroups[group], drop = FALSE]
        fm <- SKAT::SKAT(Z = Z, obj = H0, ...)
        OUT[group, ] <- c(ncol(Z), fm$p.value)
        if (verbose) {
            message("Group ", group, " of ", length(uniqueGroups), " (", round(group / length(uniqueGroups) * 100, 3), "% done)")
        }
    }

    return(OUT)
}


getCoefficients <- function(x) {
    UseMethod("getCoefficients")
}


getCoefficients.lm <- function(x) {
    summary(x)$coef[2, ]
}


getCoefficients.glm <- function(x) {
    summary(x)$coef[2, ]
}


getCoefficients.lmerMod <- function(x) {
    ans <- summary(x)$coef[2, ]
    ans <- c(ans, c(1 - stats::pnorm(ans[3])))
    return(ans)
}


#' Generates Various Summary Statistics.
#'
#' Computes the frequency of missing values, the (minor) allele frequency, and
#' standard deviation of each column of `X`.
#'
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param bufferSize The number of columns of `X` that are brought into RAM for
#' processing. Overwrites `nBuffers`. If both parameters are `NULL`, all
#' elements in `j` are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the columns of `X` that are
#' brought into RAM for processing. Is overwritten by `bufferSize`. If both
#' parameters are `NULL`, all elements in `j` are used. Defaults to `NULL`.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A `data.frame` with three columns: `freq_na` for frequencies of
#' missing values, `allele_freq` for (minor) allele frequencies, and `sd` for
#' standard deviations.
#' @export
summarize <- function(X, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000, nTasks = nCores, nBuffers = NULL, nCores = parallel::detectCores(), verbose = FALSE) {
    m <- chunkedApply(X = X, MARGIN = 2, FUN = function(col) {
        freqNA <- mean(is.na(col))
        alleleFreq <- mean(col, na.rm = TRUE) / 2
        sd <- stats::sd(col, na.rm = TRUE)
        cbind(freqNA, alleleFreq, sd)
    }, i = i, j = j, bufferSize = bufferSize, nBuffers = nBuffers, nTasks = nTasks, nCores = nCores, verbose = verbose)
    df <- data.frame(
        freq_na = m[1, ],
        allele_freq = m[2, ],
        sd = m[3, ],
        stringsAsFactors = FALSE
    )
    rownames(df) <- colnames(X)[j]
    return(df)
}


getLineCount <- function(path, header) {
    file <- file(path, open = "r")
    n <- 0
    while (length(readLines(file, n = 1)) > 0) {
        n <- n + 1
    }
    if (header) {
        n <- n - 1
    }
    close(file)
    return(n)
}


getFileHeader <- function(path, sep = "") {
    file <- file(path, open = "r")
    header <- scan(file, nlines = 1, what = character(), sep = sep, quiet = TRUE)
    close(file)
    return(header)
}


getColumnCount <- function(path, sep = "") {
    header <- getFileHeader(path, sep)
    p <- length(header)
    return(p)
}


randomString <- function() {
    paste(sample(c(0:9, letters, LETTERS), size = 5, replace = TRUE), collapse = "")
}


normalizeType <- function(val) {
    type <- typeof(val)
    # detect strings
    if (type == "character" && length(val) > 0) {
        # convert to type if type and value match
        convert <- try(vector(mode = val), silent = TRUE)
        if (class(convert) == "try-error") {
            # return a character type if conversion failed
            warning("could no convert type, using character instead")
            character()
        } else {
            # return conversion result otherwise
            convert
        }
        # value doesn't contain type information and can be handled by typeof
    } else {
        val
    }
}


simplifyList <- function(x) {
    sample <- x[[1]]
    if (is.matrix(sample)) {
        x <- matrix(data = unlist(x), nrow = nrow(sample), byrow = FALSE)
        rownames(x) <- rownames(sample)
    } else {
        x <- unlist(x)
    }
    return(x)
}


padDigits <- function(x, total) {
    formatC(x, width = as.integer(log10(total) + 1), format = "d", flag = "0")
}


loadExample <- function() {
    path <- system.file("extdata", package = "BGData")
    message("Loading chromosomes as BED files...")
    m <- do.call("ColumnLinkedMatrix", lapply(c("chr1", "chr2", "chr3"), function(chr) {
        suppressMessages(BEDMatrix::BEDMatrix(paste0(path, "/", chr)))
    }))
    as.BGData(m)
}

padDigits <- function(x, total) {
    formatC(x, width = as.integer(log10(total) + 1L), format = "d", flag = "0")
}


#' Computes a Genomic Relationship Matrix.
#'
#' Computes a positive semi-definite symmetric genomic relation matrix G=XX'
#' offering options for centering and scaling the columns of `X` beforehand.
#'
#' If `center = FALSE`, `scale = FALSE` and `scaleG = FALSE`, [getG()] produces
#' the same outcome than [base::tcrossprod()].
#'
#' @inheritSection BGData-package Memory-mapping
#' @inheritSection BGData-package Multi-level parallelism
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param center Either a logical value or a numeric vector of length equal to
#' the number of columns of `X`. Numeric vector required if `i2` is used. If
#' `FALSE`, no centering is done. Defaults to `TRUE`.
#' @param scale Either a logical value or a numeric vector of length equal to
#' the number of columns of `X`. Numeric vector required if `i2` is used. If
#' `FALSE`, no scaling is done. Defaults to `TRUE`.
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
#' `NULL`, the whole genomic relationship matrix XX' is computed. Defaults to
#' `NULL`.
#' @param bufferSize The number of columns of `X` that are brought into
#' physical memory for processing per core. If `NULL`, all columns of `X` are
#' used. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A positive semi-definite symmetric numeric matrix.
#' @example man/examples/getG.R
#' @export
getG <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05, i = seq_len(nrow(X)), j = seq_len(ncol(X)), i2 = NULL, bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {

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

    i <- convertIndexTypes(i, rownames(X))
    j <- convertIndexTypes(j, colnames(X))
    if (hasY) {
        i2 <- convertIndexTypes(i2, rownames(X))
    }

    nX <- nrow(X)
    pX <- ncol(X)

    if (min(i) < 1L || max(i) > nX) {
        stop("Index out of bounds")
    }
    if (min(j) < 1L || max(j) > pX) {
        stop("Index out of bounds")
    }
    if (hasY) {
        if (min(i2) < 1L || max(i2) > nX) {
            stop("Index out of bounds")
        }
    }

    n <- length(i)
    p <- length(j)
    if (hasY) {
        n2 <- length(i2)
    }

    if (is.null(bufferSize)) {
        bufferSize <- p
        nBuffers <- 1L
    } else {
        nBuffers <- ceiling(p / bufferSize)
    }

    if (hasY) {
        G <- bigmemory::big.matrix(nrow = n, ncol = n2, type = "double", init = 0.0, dimnames = list(rownames(X)[i], rownames(X)[i2]))
    } else {
        G <- bigmemory::big.matrix(nrow = n, ncol = n, type = "double", init = 0.0, dimnames = list(rownames(X)[i], rownames(X)[i]))
    }

    mutex <- synchronicity::boost.mutex()

    bufferRanges <- LinkedMatrix:::chunkRanges(p, nBuffers)
    bufferApply <- function(curBuffer) {

        if (verbose) {
            message("Buffer ", curBuffer, " of ", nBuffers, " ...")
        }

        # subset
        localColIndex <- j[seq(bufferRanges[1L, curBuffer], bufferRanges[2L, curBuffer])]
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
            scale.chunk <- apply(X = X1, MARGIN = 2L, FUN = stats::sd, na.rm = TRUE)
        } else if (is.numeric(scale)) {
            scale.chunk <- scale[localColIndex]
        } else {
            scale.chunk <- FALSE
        }

        # remove constant columns
        if (!(is.logical(scale) && scale == FALSE)) {
            removeCols <- which(scale.chunk < minVar)
            if (length(removeCols) > 0L) {
                X1 <- X1[, -removeCols]
                if (hasY) {
                    X2 <- X2[, -removeCols]
                }
                scale.chunk <- scale.chunk[-removeCols]
                center.chunk <- center.chunk[-removeCols]
            }
        }

        p <- ncol(X1)

        # compute XX'
        if (p > 0L) {

            # scale and impute X
            X1 <- scale(X1, center = center.chunk, scale = scale.chunk)
            X1[is.na(X1)] <- 0L
            if (hasY) {
                X2 <- scale(X2, center = center.chunk, scale = scale.chunk)
                X2[is.na(X2)] <- 0L
            }

            if (hasY) {
                G_chunk <- tcrossprod(x = X1, y = X2)
            } else {
                G_chunk <- tcrossprod(X1)
            }

            synchronicity::lock(mutex)
            G[] <- G[] + G_chunk
            synchronicity::unlock(mutex)

        }

        return(p)

    }

    if (nCores == 1L) {
        res <- lapply(X = seq_len(nBuffers), FUN = bufferApply)
    } else {
        res <- parallel::mclapply(seq_len(nBuffers), bufferApply)
    }

    # Convert big.matrix to matrix
    G <- G[]

    if (scaleG) {
        if (hasY) {
            K <- do.call(base::sum, res)
        } else {
            # Use seq instead of diag to avoid copy as it does not increase ref count
            K <- mean(G[seq(from = 1L, to = n * n, by = n + 1L)])
        }
        G[] <- G / K
    }

    return(G)

}


#' Computes a Very Large Genomic Relationship Matrix.
#'
#' Computes a positive semi-definite symmetric genomic relation matrix G=XX'
#' offering options for centering and scaling the columns of `X` beforehand.
#'
#' Even very large genomic relationship matrices are supported by partitioning
#' `X` into blocks and calling [getG()] on these blocks. This function performs
#' the block computations sequentially, which may be slow. In an HPC
#' environment, performance can be improved by manually distributing these
#' operations to different nodes.
#'
#' @inheritSection BGData-package Multi-level parallelism
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param center Either a logical value or a numeric vector of length equal to
#' the number of columns of `X`. If `FALSE`, no centering is done. Defaults to
#' `TRUE`.
#' @param scale Either a logical value or a numeric vector of length equal to
#' the number of columns of `X`. If `FALSE`, no scaling is done. Defaults to
#' `TRUE`.
#' @param scaleG TRUE/FALSE whether xx' must be scaled.
#' @param minVar Columns with variance lower than this value will not be used
#' in the computation (only if `scale` is not `FALSE`).
#' @param folderOut The path to the folder where to save the
#' [symDMatrix::symDMatrix-class] object. Defaults to a random string prefixed
#' with "symDMatrix_".
#' @param vmode vmode of `ff` objects.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param blockSize The number of rows and columns of each block. If `NULL`, a
#' single block of the same length as `i` will be created. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A [symDMatrix::symDMatrix-class] object.
#' @export
getG_symDMatrix <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05, folderOut = paste0("symDMatrix_", randomString()), vmode = "double", i = seq_len(nrow(X)), j = seq_len(ncol(X)), blockSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {

    i <- convertIndexTypes(i, rownames(X))
    j <- convertIndexTypes(j, colnames(X))

    nX <- nrow(X)
    pX <- ncol(X)

    if (min(i) < 1L || max(i) > nX) {
        stop("Index out of bounds")
    }
    if (min(j) < 1L || max(j) > pX) {
        stop("Index out of bounds")
    }

    n <- length(i)
    p <- length(j)

    if (is.null(blockSize)) {
        blockSize <- n
        nBlocks <- 1L
    } else {
        nBlocks <- ceiling(n / blockSize)
    }

    if (is.logical(center) && center == TRUE) {
        center <- chunkedApply(X = X, MARGIN = 2L, FUN = mean, i = i, j = j, bufferSize = blockSize, nCores = nCores, verbose = FALSE, na.rm = TRUE)
    } else if (is.logical(center) && center == FALSE) {
        center <- rep(0L, p)
    }
    names(center) <- colnames(X)[j]

    if (is.logical(scale) && scale == TRUE) {
        scale <- chunkedApply(X = X, MARGIN = 2L, FUN  = stats::sd, i = i, j = j, bufferSize = blockSize, nCores = nCores, verbose = FALSE, na.rm = TRUE)
        scale <- scale * sqrt((nX - 1L) / nX) # to avoid NaN
    } else if (is.logical(scale) && scale == FALSE) {
        scale <- rep(1L, p)
    }
    names(scale) <- colnames(X)[j]

    if (file.exists(folderOut)) {
        stop(folderOut, " already exists")
    }
    dir.create(folderOut)

    blockIndices <- split(i, ceiling(i / blockSize))
    args <- vector(mode = "list", length = nBlocks)
    counter <- 1L
    for (rowIndex in 1L:nBlocks) {
        rowArgs <- vector(mode = "list", length = nBlocks)
        for (colIndex in 1L:nBlocks) {
            if (verbose) {
                message("Block ", rowIndex, "-", colIndex, " ...")
            }
            if (colIndex >= rowIndex) {
                blockName <- paste0("data_", padDigits(rowIndex, nBlocks), "_", padDigits(colIndex, nBlocks), ".bin")
                block <- ff::as.ff(getG(X, center = center, scale = scale, scaleG = FALSE, minVar = minVar, i = blockIndices[[rowIndex]], j = j, i2 = blockIndices[[colIndex]], bufferSize = blockSize, nCores = nCores, verbose = FALSE), filename = paste0(folderOut, "/", blockName), vmode = vmode)
                # Change ff path to a relative one
                bit::physical(block)$filename <- blockName
                rowArgs[[colIndex]] <- block
                counter <- counter + 1L
            } else {
                rowArgs[[colIndex]] <- ff::vt(args[[colIndex]][[rowIndex]])
            }
        }
        args[[rowIndex]] <- do.call(LinkedMatrix::ColumnLinkedMatrix, rowArgs)
    }

    G <- do.call(symDMatrix::symDMatrix, args)

    # Keep centers and scales
    attr(G, "centers") <- center
    attr(G, "scales") <- scale

    if (scaleG) {
        K <- mean(diag(G))
        for (rowIndex in seq_len(nBlocks)) {
            for (colIndex in seq(rowIndex, nBlocks)) {
                G[[rowIndex]][[colIndex]][] <- G[[rowIndex]][[colIndex]][] / K
            }
        }
    }

    save(G, file = paste0(folderOut, "/G.RData"))

    return(G)

}

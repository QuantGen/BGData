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


#' Applies a function on each row or column of a matrix in parallel.
#'
#' Similar to \code{apply}, designed to carry out operations in parallel.
#'
#' The input matrix \code{X} is broken into \code{nTasks} chunks and passed to
#' \code{\link[parallel]{mclapply}}. The number of cores can be configured using
#' \code{nCores}.
#'
#' \code{nTasks} has to be chosen carefully to avoid running out of memory. As a
#' rule of thumb, at least around \code{object_size(X) + (nCores *
#' (object_size(X) / nTasks)) + object_size(result)} MB of total memory will be
#' needed, not including potential copies of your data that might be created
#' (for example \code{lsfit} runs \code{cbind(1, X)}). Therefore, for 20 nodes
#' and 20 tasks you will need at least \code{2 * object_size(X)} MB, for 20
#' nodes and 40 tasks \code{1.5 * object_size(X)} MB, etc.
#'
#' If \code{nTasks} equals 1, the regular \code{apply} function will be called
#' to preserve memory.
#'
#' @param X A matrix.
#' @param MARGIN The subscripts which the function will be applied over. 1
#'   indicates rows, 2 indicates columns.
#' @param FUN The function to be applied.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @param ... Additional arguments to be passed to \code{apply}.
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


#' Reads chunks of data from a memory-mapped file into memory and applies a
#' function on each row or column of a matrix in parallel.
#'
#' Similar to \code{apply}, designed to bring chunks of data into memory and
#' carry out operations on them in parallel.
#'
#' \code{bufferSize} and \code{nTasks} have to be chosen carefully to avoid
#' running out of memory. As a rule of thumb, at least around
#' \code{object_size(buffer) + (nCores * (object_size(buffer) / nTasks)) +
#' object_size(result)} MB of total memory will be needed, not including
#' potential copies of your data that might be created (for example \code{lsfit}
#' runs \code{cbind(1, X)}). Therefore, for 20 nodes and 20 tasks you will need
#' at least \code{2 * object_size(buffer)} MB, for 20 nodes and 40 tasks
#' \code{1.5 * object_size(buffer)} MB, etc.
#'
#' This function is only useful for memory-mapped files. For data that is
#' already in memory, use \code{\link{parallelApply}} directly.
#'
#' @param X A matrix-like object, typically \code{@@geno} of a
#'   \code{\link[=BGData-class]{BGData}} object.
#' @param MARGIN The subscripts which the function will be applied over. 1
#'   indicates rows, 2 indicates columns.
#' @param FUN The function to be applied.
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be
#'   used. By default, all columns are used.
#' @param bufferSize The number of rows or columns of \code{X} that are brought
#'   into RAM for processing. Overwrites \code{nBuffers}. If both parameters are
#'   \code{NULL}, all elements in \code{i} or \code{j} are used. Defaults to
#'   5000.
#' @param nBuffers The number of partitions of the rows or columns of \code{X}
#'   that are brought into RAM for processing. Is overwritten by
#'   \code{bufferSize}. If both parameters are \code{NULL}, all elements in
#'   \code{i} or \code{j} are used.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{FALSE}.
#' @param ... Additional arguments to be passed to \code{parallelApply}.
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


# Computes crossprod(x,y) or tcrossprod(x,y)
crossprods <- function(x, y = NULL, use_tcrossprod = FALSE, nTasks = nCores, nCores = parallel::detectCores()) {
    dx <- dim(x)
    if (!is.null(y)) {
        y <- as.matrix(y)
        dy <- dim(y)
        if (use_tcrossprod) {
            if (dx[2] != ncol(y)) {
                stop("Error in tcrossprod.parallel: non-conformable arguments.")
            }
        } else {
            if (dx[1] != dy[1]) {
                stop("Error in crossprod.parallel: non-conformable arguments.")
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


#' Computes crossprod (x'y or x'x) in parallel.
#'
#' @param x A matrix-like object, typically \code{@@geno} of a
#'   \code{\link[=BGData-class]{BGData}} object.
#' @param y vector or matrix-like object. NULL by default.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @return x'y' or x'x depending on whether y is provided.
#' @export
crossprod.parallel <- function(x, y = NULL, nTasks = nCores, nCores = parallel::detectCores()) {
    crossprods(x = x, y = y, use_tcrossprod = FALSE, nTasks = nTasks, nCores = nCores)
}


#' Computes tcrossprod (xy' or xx') in parallel.
#'
#' @param x A matrix-like object, typically \code{@@geno} of a
#'   \code{\link[=BGData-class]{BGData}} object.
#' @param y vector or matrix-like object. NULL by default.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @return xy' or xx' depending on whether y is provided.
#' @export
tcrossprod.parallel <- function(x, y = NULL, nTasks = nCores, nCores = parallel::detectCores()) {
    crossprods(x = x, y = y, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
}



#' Computes a genomic relationship matrix G=xx'.
#'
#' Offers options for centering and scaling the columns of x before computing
#' xx'. If \code{centerCol=FALSE}, \code{scaleCol=FALSE} and
#' \code{scaleG=FALSE}, \code{getG} produces the same outcome than
#' \code{tcrossprod}.
#'
#' @param x A matrix-like object, typically \code{@@geno} of a
#'   \code{\link[=BGData-class]{BGData}} object.
#' @param scaleCol TRUE/FALSE whether columns must be scaled before computing
#'   xx'.
#' @param scales Precomputed scales if i2 is used.
#' @param centerCol TRUE/FALSE whether columns must be centered before computing
#'   xx'.
#' @param centers Precomputed centers if i2 is used.
#' @param scaleG TRUE/FALSE whether xx' must be scaled.
#' @param minVar Columns with variance lower than this value will not be used in
#'   the computation (only if \code{scaleCol} is set).
#' @param saveG Whether to save genomic relationship matrix into file.
#' @param saveType File format to save genomic relationship matrix in. Either
#'   \code{RData} or \code{ff}.
#' @param saveName Name without extension to save genomic relationship matrix
#'   with.
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be
#'   used. By default, all columns are used.
#' @param i2 (integer, boolean or character) Indicates which rows should be used
#'   to divide matrix into blocks.
#' @param bufferSize The number of columns of \code{x} that are brought into
#'   RAM for processing. Overwrites \code{nBuffers}. If both parameters are
#'   \code{NULL}, all columns of \code{x} are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the columns of \code{x} that are
#'   brought into RAM for processing. Is overwritten by \code{bufferSize}. If
#'   both parameters are \code{NULL}, all columns of \code{x} are used.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{TRUE}.
#' @return A positive semi-definite symmetric numeric matrix.
#' @export
getG <- function(x, scaleCol = TRUE, scales = NULL, centerCol = TRUE, centers = NULL, scaleG = TRUE, minVar = 1e-05, saveG = FALSE, saveType = "RData", saveName = "Gij", i = seq_len(nrow(x)), j = seq_len(ncol(x)), i2 = NULL, bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = TRUE) {
    if (is.null(bufferSize) && is.null(nBuffers)) {
        bufferSize <- length(j)
        nBuffers <- 1
    } else if (is.null(bufferSize) && !is.null(nBuffers)) {
        bufferSize <- ceiling(length(j) / nBuffers)
    } else {
        nBuffers <- ceiling(length(j) / bufferSize)
    }
    if (is.null(i2)) {
        G <- getGi(x = x, scales = scales, centers = centers, scaleCol = scaleCol, centerCol = centerCol, scaleG = scaleG, minVar = minVar, i = i, j = j, bufferSize = bufferSize, nBuffers = nBuffers, nTasks = nTasks, nCores = nCores, verbose = verbose)
    } else {
        G <- getGij(x = x, scales = scales, centers = centers, scaleCol = scaleCol, centerCol = centerCol, scaleG = scaleG, minVar = minVar, i = i, i2 = i2, j = j, bufferSize = bufferSize, nBuffers = nBuffers, nTasks = nTasks, nCores = nCores, verbose = verbose)
    }
    if (saveG) {
        if (saveType == "RData") {
            save(G, file = paste0(saveName, ".RData"))
        }
        if (saveType == "ff") {
            Gij <- ff::as.ff(G, file = paste0(saveName, ".bin"))
            save(Gij, file = paste0(saveName, ".ff"))
        }
    }
    return(G)
}


getGi <- function(x, scaleCol = TRUE, scales = NULL, centerCol = FALSE, centers = NULL, scaleG = TRUE, minVar = 1e-05, i = seq_len(nrow(x)), j = seq_len(ncol(x)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = TRUE) {
    nX <- nrow(x)
    pX <- ncol(x)

    # Convert index types
    if (is.logical(i)) {
        i <- which(i)
    } else if (is.character(i)) {
        i <- match(i, rownames(x))
    }
    if (is.logical(j)) {
        j <- which(j)
    } else if (is.character(j)) {
        j <- match(j, colnames(x))
    }

    n <- length(i)
    p <- length(j)

    if (n > nX | p > pX) {
        stop("Index out of bounds")
    }
    if (is.numeric(i)) {
        if ((min(i) < 1) | (max(i) > nX)) {
            stop("Index out of bounds")
        }
    }
    if (is.numeric(j)) {
        if ((min(j) < 1) | (max(j) > pX)) {
            stop("Index out of bounds")
        }
    }

    G <- matrix(data = 0, nrow = n, ncol = n, dimnames = list(rownames(x)[i], rownames(x)[i]))

    bufferSize <- ceiling(p / nBuffers)

    end <- 0
    for (k in seq_len(nBuffers)) {
        ini <- end + 1
        if (ini <= p) {
            end <- min(p, ini + bufferSize - 1)
            if (verbose) {
                message("Chunk: ", k, " (markers ", ini, ":", end, " ~", round(100 * end / p, 1), "% done)")
                message("  => Acquiring genotypes...")
            }

            # subset
            localColIndex <- j[ini:end]
            X <- x[i, localColIndex, drop = FALSE]

            # compute centers
            if (centerCol) {
                if (is.null(centers)) {
                    centers.chunk <- colMeans(X, na.rm = TRUE)
                } else {
                    centers.chunk <- centers[localColIndex]
                }
            } else {
                centers.chunk = FALSE
            }

            # compute scales
            if (scaleCol) {
                if (is.null(scales)) {
                    scales.chunk <- apply(X = X, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
                } else {
                    scales.chunk <- scales[localColIndex]
                }
                removeCols <- which(scales.chunk < minVar)
                if (length(removeCols) > 0) {
                  X <- X[, -removeCols]
                  scales.chunk <- scales.chunk[-removeCols]
                  centers.chunk <- centers.chunk[-removeCols]
                }
            } else {
                scales.chunk <- FALSE
            }

            if (ncol(X) > 0) {
                if (verbose) {
                  message("  =>Computing...")
                }
                X <- scale(X, center = centers.chunk, scale = scales.chunk)
                X[is.na(X)] <- 0

                if (nTasks > 1) {
                  G_chunk <- crossprods(x = X, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
                } else {
                  G_chunk <- tcrossprod(X)
                }
                G[] <- G + G_chunk
            }
        }
    }
    if (scaleG) {
        # Use seq instead of diag to avoid copy as it does not increase ref count
        G[] <- G / mean(G[seq(from = 1, to = n * n, by = n + 1)])
    }

    return(G)
}


getGij <- function(x, scaleCol = TRUE, scales, centerCol = FALSE, centers, scaleG = TRUE, minVar = 1e-05, i, i2, j = seq_len(ncol(x)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = TRUE) {

    if (scaleCol && is.null(scales)) {
        stop("scales need to be precomputed.")
    }
    if (centerCol && is.null(centers)) {
        stop("centers need to be precomputed.")
    }

    nX <- nrow(x)
    pX <- ncol(x)

    # Convert index types
    if (is.logical(i)) {
        i <- which(i)
    } else if (is.character(i)) {
        i <- match(i, rownames(x))
    }
    if (is.logical(i2)) {
        i2 <- which(i2)
    } else if (is.character(i2)) {
        i2 <- match(i2, rownames(x))
    }
    if (is.logical(j)) {
        j <- which(j)
    } else if (is.character(j)) {
        j <- match(j, colnames(x))
    }

    n1 <- length(i)
    p <- length(j)
    n2 <- length(i2)

    if ((min(i) < 1) | (max(i) > nX)) {
        stop("Index out of bounds")
    }
    if ((min(i2) < 1) | (max(i2) > nX)) {
        stop("Index out of bounds")
    }
    if ((min(j) < 1) | (max(j) > pX)) {
        stop("Index out of bounds")
    }

    K <- 0

    G <- matrix(data = 0, nrow = n1, ncol = n2, dimnames = list(rownames(x)[i], rownames(x)[i2]))

    bufferSize <- ceiling(p / nBuffers)

    end <- 0
    for (k in seq_len(nBuffers)) {
        ini <- end + 1
        if (ini <= p) {
            end <- min(p, ini + bufferSize - 1)
            if (verbose) {
                message("Working on chunk: ", k, " (markers ", ini, ":", end, " ~", round(100 * ini / p, 1), "% done)")
                message("  => Acquiring genotypes...")
            }

            # subset
            localColIndex <- j[ini:end]
            X1 <- x[i, localColIndex, drop = FALSE]
            X2 <- x[i2, localColIndex, drop = FALSE]
            if (centerCol) {
                centers.chunk <- centers[localColIndex]
            } else {
                centers.chunk <- FALSE
            }
            if (scaleCol) {
                scales.chunk <- scales[localColIndex]
            } else {
                scales.chunk <- FALSE
            }

            if (scaleCol) {
                removeCols <- which(scales.chunk < sqrt(minVar))
                if (length(removeCols) > 0) {
                  X1 <- X1[, -removeCols]
                  X2 <- X2[, -removeCols]
                  scales.chunk <- scales.chunk[-removeCols]
                  centers.chunk <- centers.chunk[-removeCols]
                }
            }

            if (ncol(X1) > 0) {
                if (verbose) {
                    message("  => Computing...")
                }
                X1 <- scale(X1, center = centers.chunk, scale = scales.chunk)
                X1[is.na(X1)] <- 0
                X2 <- scale(X2, center = centers.chunk, scale = scales.chunk)
                X2[is.na(X2)] <- 0
                G_chunk <- tcrossprod.parallel(x = X1, y = X2, nTasks = nTasks, nCores = nCores)
                G[] <- G + G_chunk
            }
            if (scaleG) {
                if (scaleCol) {
                    K <- K + ncol(X1)
                } else {
                    K <- K + sum(scales.chunk^2)
                }
            }
        }
    }
    if (scaleG) {
        G[] <- G / K
    }

    return(G)
}


#' Computes a genomic relationship matrix G=xx' without ever loading G in RAM by
#' creating a \code{\link[=symDMatrix-class]{symDMatrix}}.
#'
#' Offers options for centering and scaling the columns of x before computing
#' xx'.
#'
#' @param X A matrix-like object, typically \code{@@geno} of a
#'   \code{\link[=BGData-class]{BGData}} object.
#' @param nBlocks The number of blocks.
#' @param blockSize The number of columns of a block (if NULL inferred from block).
#' @param centers Precomputed centers.
#' @param scales Precomputed scales.
#' @param centerCol TRUE/FALSE whether columns must be centered before computing
#'   xx'.
#' @param scaleCol TRUE/FALSE whether columns must be scaled before computing
#'   xx'.
#' @param scaleG TRUE/FALSE whether xx' must be scaled.
#' @param folderOut The path to the folder where to save the
#'   \code{\link[=symDMatrix-class]{symDMatrix}} object. Defaults to a random
#'   string prefixed with "symDMatrix_".
#' @param vmode vmode of \code{ff} objects.
#' @param saveRData Whether to save an RData file to easily reload
#'   \code{\link[=symDMatrix-class]{symDMatrix}}
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be
#'   used. By default, all columns are used.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{TRUE}.
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

    if (n > nX | p > pX) {
        stop("Index out of bounds")
    }
    if (is.numeric(i)) {
        if ((min(i) < 1) | (max(i) > nX)) {
            stop("Index out of bounds")
        }
    }
    if (is.numeric(j)) {
        if ((min(j) < 1) | (max(j) > pX)) {
            stop("Index out of bounds")
        }
    }

    if ((centerCol | scaleCol) & (is.null(centers) | is.null(scales))) {
        if (is.null(centers) & is.null(scales)) {
            centers <- vector(mode = "double", length = p)
            scales <- vector(mode = "double", length = p)
            for (k in seq_len(p)) {
                xi <- X[i, j[k]]
                scales[k] <- stats::sd(xi, na.rm = TRUE) * sqrt((nX - 1) / nX)
                centers[k] <- mean(xi, na.rm = TRUE)
            }
        }
        if ((!is.null(centers)) & (is.null(scales))) {
            scales <- vector(mode = "double", length = p)
            for (k in seq_len(p)) {
                xi <- X[i, j[k]]
                scales[k] <- stats::sd(xi, na.rm = TRUE) * sqrt((nX - 1) / nX)
            }
        }
        if ((is.null(centers)) & (!is.null(scales))) {
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

            Gij <- tcrossprod.parallel(x = Xi, y = Xj, nTasks = nTasks, nCores = nCores)

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


#' Performs single marker regressions using a
#' \code{\link[=BGData-class]{BGData}} object.
#'
#' Implements single marker regressions. The regression model includes all the
#' covariates specified in the right-hand-side of the \code{formula} plus one
#' column of \code{@@geno}, one column at a time. The data from the association
#' tests is obtained from a \code{\link[=BGData-class]{BGData}} object.
#'
#' @param formula A formula (e.g. weight~sex+age) with the response on the
#'   left-hand side and predictors (all the covariates except the markers) on
#'   the right-hand side. The variables included in the formula must be in the
#'   \code{@@pheno} object of the \code{\link[=BGData-class]{BGData}}.
#' @param data A \code{\link[=BGData-class]{BGData}} object.
#' @param method The regression method to be used. Currently, the following
#'   methods are implemented: \code{\link{lm}}, \code{\link{lm.fit}},
#'   \code{\link{lsfit}}, \code{\link{glm}}, \code{\link[lme4]{lmer}},
#'   and \code{\link[SKAT]{SKAT}}.
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be
#'   used. By default, all columns are used.
#' @param bufferSize The number of columns of \code{@@geno} that are brought into
#'   RAM for processing. Overwrites \code{nBuffers}. If both parameters are
#'   \code{NULL}, all elements in \code{j} are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the columns of \code{@@geno}
#'   that are brought into RAM for processing. Is overwritten by
#'   \code{bufferSize}. If both parameters are \code{NULL}, all elements in
#'   \code{j} are used.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{FALSE}.
#' @param ... Additional arguments for chunkedApply and regression method.
#' @return Returns a matrix with estimates, SE, p-value, etc.
#' @export
GWAS <- function(formula, data, method, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = FALSE, ...) {

    if (class(data) != "BGData") {
        stop("data must BGData")
    }

    if (!method %in% c("lm", "lm.fit", "lsfit", "glm", "lmer", "SKAT")) {
        stop("Only lm, glm, lmer and SKAT have been implemented so far.")
    }

    # We can have specialized methods, for instance for OLS it is better to use
    # lsfit (that is what GWAS.ols does)
    if (method %in% c("lm", "lm.fit", "lsfit")) {
        OUT <- GWAS.ols(formula = formula, data = data, i = i, j = j, bufferSize = bufferSize, nTasks = nTasks, nCores = nCores, verbose = verbose, ...)
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

## OLS for the regression y=xb+e (data is assumed to be pre-adjusted by non-genetic effects
    
rayOLS=function(y,x,n=length(y)){
	tmp=!(is.na(x)|is.na(y))
	x=x[tmp]
	y=y[tmp]
	
	rhs=sum(x*y)
	XtX=sum(x^2)
	sol=rhs/XtX
	error=y-x*sol
	vE=sum(error^2)/(n-1)
	SE=sqrt(vE/XtX)
	z_stat=sol/SE
	return(c(sol,SE,z_stat,pt(q=abs(z_stat),df=n-1,lower.tail=FALSE)*2))
}

    
# GWAS 'Ordinary Least Squares' (e.g., lsfit, lm.fit, lm)
# formula: the formula for the GWAS model without including the marker, e.g.
# y~1 or y~factor(sex)+age
# all the variables in the formula must be in data@pheno data (BGData)
# containing slots @pheno and @geno    
GWAS.ols <- function(formula, data, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000, nBuffers = NULL, nTasks = nCores, nCores = parallel::detectCores(), verbose = FALSE, ...) {

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


#' Calculates frequencies of missing values and alleles.
#'
#' @param X A matrix-like object, typically \code{@@geno} of a
#'   \code{\link[=BGData-class]{BGData}} object.
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be
#'   used. By default, all columns are used.
#' @param bufferSize The number of columns of \code{X} that are brought into
#'   RAM for processing. Overwrites \code{nBuffers}. If both parameters are
#'   \code{NULL}, all elements in \code{j} are used. Defaults to 5000.
#' @param nBuffers The number of partitions of the columns of \code{X} that
#'   are brought into RAM for processing. Is overwritten by \code{bufferSize}. If
#'   both parameters are \code{NULL}, all elements in \code{j} are used.
#' @param nTasks The number of tasks the problem should be broken into to be
#'   distributed among \code{nCores} cores. Defaults to \code{nCores}.
#' @param nCores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}). Defaults to the number of cores as
#'   detected by \code{\link[parallel]{detectCores}}).
#' @param verbose Whether progress updates will be posted. Defaults to
#'   \code{FALSE}.
#' @export
summarize <- function(X, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000, nTasks = nCores, nBuffers = NULL, nCores = parallel::detectCores(), verbose = FALSE) {
    res <- chunkedApply(X = X, MARGIN = 2, FUN = function(col) {
        freqNA <- mean(is.na(col))
        alleleFreq <- mean(col, na.rm = TRUE) / 2
        sd <- stats::sd(col, na.rm = TRUE)
        cbind(freqNA, alleleFreq, sd)
    }, i = i, j = j, bufferSize = bufferSize, nBuffers = nBuffers, nTasks = nTasks, nCores = nCores, verbose = verbose)
    rownames(res) <- c("freq_na", "allele_freq", "sd")
    colnames(res) <- colnames(X)[j]
    t(res)
}


getLineCount <- function(path, header) {
    # gzfile and readLines throw some warnings, but since it works, let's disable
    # warnings for this block
    warnLevel <- unlist(options("warn"))
    options(warn = -1)
    file <- gzfile(path, open = "r")
    n <- 0
    while (length(readLines(file, n = 1)) > 0) {
        n <- n + 1
    }
    if (header) {
        n <- n - 1
    }
    close(file)
    # restore previous warning level
    options(warn = warnLevel)
    return(n)
}


getFileHeader <- function(path, sep = "") {
    file <- gzfile(path, open = "r")
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

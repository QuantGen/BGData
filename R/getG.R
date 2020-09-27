padDigits <- function(x, total) {
    formatC(x, width = as.integer(log10(total) + 1L), format = "d", flag = "0")
}

getG <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05, i = seq_len(nrow(X)), j = seq_len(ncol(X)), i2 = NULL, chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {

    # compute XY' rather than XX'
    hasY <- !is.null(i2)

    if (hasY) {
        if (is.logical(center) && center == TRUE) {
            stop("centers need to be precomputed.")
        }
        if (is.logical(scale) && scale == TRUE) {
            stop("scales need to be precomputed.")
        }
    }

    i <- convertIndex(X, i, "i")
    j <- convertIndex(X, j, "j")
    if (hasY) {
        i2 <- convertIndex(X, i2, "i")
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

    if (is.null(chunkSize)) {
        chunkSize <- p
        nChunks <- 1L
    } else {
        nChunks <- ceiling(p / chunkSize)
    }

    if (hasY) {
        G <- big.matrix(nrow = n, ncol = n2, type = "double", init = 0.0, dimnames = list(rownames(X)[i], rownames(X)[i2]))
    } else {
        G <- big.matrix(nrow = n, ncol = n, type = "double", init = 0.0, dimnames = list(rownames(X)[i], rownames(X)[i]))
    }

    mutex <- boost.mutex()

    chunkApply <- function(curChunk) {

        if (verbose) {
            if (nCores > 1) {
                message("Process ", Sys.getpid(), ": Chunk ", curChunk, " of ", nChunks, " ...")
            } else {
                message("Chunk ", curChunk, " of ", nChunks, " ...")
            }
        }

        # subset
        range <- seq(
            ((curChunk - 1L) * chunkSize) + 1L,
            min(curChunk * chunkSize, p)
        )
        X1 <- X[i, j[range], drop = FALSE]
        if (hasY) {
            X2 <- X[i2, j[range], drop = FALSE]
        }

        # compute centers
        if (is.logical(center) && center == TRUE) {
            center.chunk <- colMeans(X1, na.rm = TRUE)
        } else if (is.numeric(center)) {
            center.chunk <- center[j[range]]
        } else {
            center.chunk = FALSE
        }

        # compute scales
        if (is.logical(scale) && scale == TRUE) {
            scale.chunk <- apply(X = X1, MARGIN = 2L, FUN = sd, na.rm = TRUE)
        } else if (is.numeric(scale)) {
            scale.chunk <- scale[j[range]]
        } else {
            scale.chunk <- FALSE
        }

        # remove constant columns
        if (is.numeric(scale.chunk)) {
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

            # center, scale and impute without duplications
            # set nCores to 1 here because section is already parallelized
            X1 <- preprocess(X1, center = center.chunk, scale = scale.chunk, impute = TRUE, nCores = 1)
            if (hasY) {
                X2 <- preprocess(X2, center = center.chunk, scale = scale.chunk, impute = TRUE, nCores = 1)
            }

            if (hasY) {
                G_chunk <- tcrossprod(x = X1, y = X2)
            } else {
                G_chunk <- tcrossprod(X1)
            }

            lock(mutex)
            G[] <- G[] + G_chunk
            unlock(mutex)

        }

        return(p)

    }

    if (nCores == 1L) {
        res <- lapply(X = seq_len(nChunks), FUN = chunkApply)
    } else {
        res <- mclapply(X = seq_len(nChunks), FUN = chunkApply, mc.cores = nCores)
    }

    # Convert big.matrix to matrix
    G <- G[]

    if (scaleG) {
        if (hasY) {
            K <- do.call(sum, res)
        } else {
            # Use seq instead of diag to avoid copy as it does not increase ref count
            K <- mean(G[seq(from = 1L, to = n * n, by = n + 1L)])
        }
        G[] <- G / K
    }

    return(G)

}

getG_symDMatrix <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05, blockSize = 5000L, folderOut = paste0("symDMatrix_", randomString()), vmode = "double", i = seq_len(nrow(X)), j = seq_len(ncol(X)), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {

    i <- convertIndex(X, i, "i")
    j <- convertIndex(X, j, "j")

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

    if (is.null(chunkSize)) {
        chunkSize <- p
        nChunks <- 1L
    } else {
        nChunks <- ceiling(p / chunkSize)
    }

    if (is.logical(center) && center == TRUE) {
        if (verbose) {
            message("Computing centers ...")
        }
        center <- rep(0, pX)
        names(center) <- colnames(X)
        center[j] <- chunkedApply(X = X, MARGIN = 2L, FUN = mean, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = FALSE, na.rm = TRUE)
    }

    if (is.logical(scale) && scale == TRUE) {
        if (verbose) {
            message("Computing scales ...")
        }
        scale <- rep(1, pX)
        names(scale) <- colnames(X)
        scale[j] <- chunkedApply(X = X, MARGIN = 2L, FUN = sd, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = FALSE, na.rm = TRUE)
    }

    if (file.exists(folderOut)) {
        stop(folderOut, " already exists")
    }
    dir.create(folderOut)

    if (is.null(blockSize)) {
        blockSize <- n
        nBlocks <- 1L
    } else {
        nBlocks <- ceiling(n / blockSize)
    }

    blockIndices <- split(i, ceiling(seq_along(i) / blockSize))
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
                block <- as.ff(getG(X, center = center, scale = scale, scaleG = FALSE, minVar = minVar, i = blockIndices[[rowIndex]], j = j, i2 = blockIndices[[colIndex]], chunkSize = chunkSize, nCores = nCores, verbose = FALSE), filename = paste0(folderOut, "/", blockName), vmode = vmode)
                # Change ff path to a relative one
                physical(block)[["filename"]] <- blockName
                rowArgs[[colIndex]] <- block
                counter <- counter + 1L
            } else {
                rowArgs[[colIndex]] <- vt(args[[colIndex]][[rowIndex]])
            }
        }
        args[[rowIndex]] <- do.call(ColumnLinkedMatrix, rowArgs)
    }

    G <- do.call(symDMatrix, args)

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

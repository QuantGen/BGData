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
        res <- parallel::mclapply(X = seq_len(nBuffers), FUN = bufferApply, ..., mc.preschedule = FALSE, mc.cores = nCores)
    }
    simplifyList(res)
}


crossprods <- function(x, y = NULL, use_tcrossprod = FALSE, nTasks = nCores, nCores = getOption("mc.cores", 2L)) {
    dx <- dim(x)
    if (!is.null(y)) {
        y <- as.matrix(y)
        dy <- dim(y)
        if (use_tcrossprod) {
            if (dx[2L] != ncol(y)) {
                stop("Error in tcrossprod_parallel: non-conformable arguments.")
            }
        } else {
            if (dx[1L] != dy[1L]) {
                stop("Error in crossprod_parallel: non-conformable arguments.")
            }
        }
    }
    if (nTasks == 1L) {
        x <- as.matrix(x)
        if (!is.null(y)) {
            y <- as.matrix(y)
        }
        if (use_tcrossprod) {
            Xy <- tcrossprod(x, y)
        } else {
            Xy <- crossprod(x, y)
        }
    } else {
        nX <- ifelse(use_tcrossprod, dx[2L], dx[1L])
        if (!is.null(y)) {
            nY <- ifelse(use_tcrossprod, dy[2L], dy[1L])
        }
        chunks <- parallel::mclapply(X = seq_len(nTasks), FUN = function(i, ...) {
            if (!is.null(y)) {
                ranges <- LinkedMatrix:::chunkRanges(nY, nTasks, i)
                if (use_tcrossprod) {
                    Y <- y[, seq(ranges[1L], ranges[2L]), drop = FALSE]
                } else {
                    Y <- y[seq(ranges[1L], ranges[2L]), , drop = FALSE]
                }
            } else {
                Y <- NULL
            }
            ranges <- LinkedMatrix:::chunkRanges(nX, nTasks, i)
            if (use_tcrossprod) {
                X <- x[, seq(ranges[1L], ranges[2L]), drop = FALSE]
                Xy <- tcrossprod(X, Y)
            } else {
                X <- x[seq(ranges[1L], ranges[2L]), , drop = FALSE]
                Xy <- crossprod(X, Y)
            }
            return(Xy)
        }, mc.preschedule = FALSE, mc.cores = nCores)
        # We now need to add up chunks sequentially
        Xy <- chunks[[1L]]
        if (length(chunks) > 1L) {
            for (i in 2L:length(chunks)) {
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
#' @inheritSection BGData-package Multi-level parallelism
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
#' @example man/examples/crossprod_parallel.R
#' @export
crossprod_parallel <- function(x, y = NULL, nTasks = nCores, nCores = getOption("mc.cores", 2L)) {
    crossprods(x = x, y = y, use_tcrossprod = FALSE, nTasks = nTasks, nCores = nCores)
}


#' @rdname crossprod_parallel
#' @export
tcrossprod_parallel <- function(x, y = NULL, nTasks = nCores, nCores = getOption("mc.cores", 2L)) {
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
#' @param bufferSize The number of columns of `X` that are brought into RAM for
#' processing. If `NULL`, all columns of `X` are used. Defaults to 5000.
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param gpu Whether to use the GPU to perform the computations. By setting
#' this parameter, `nCores` and `nTasks` will be set to 1.
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A positive semi-definite symmetric numeric matrix.
#' @example man/examples/getG.R
#' @export
getG <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, minVar = 1e-05, i = seq_len(nrow(X)), j = seq_len(ncol(X)), i2 = NULL, bufferSize = 5000L, nTasks = nCores, nCores = getOption("mc.cores", 2L), gpu = FALSE, verbose = FALSE) {

    if (gpu) {
        if (!requireNamespace("gpuR", quietly = TRUE)) {
            stop("gpuR needed for this function to work. Please install it.", call. = FALSE)
        }
        nCores <- 1L
        nTasks <- 1L
    }

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

    if (!hasY) {
        G <- matrix(data = 0.0, nrow = n, ncol = n, dimnames = list(rownames(X)[i], rownames(X)[i]))
    } else {
        G <- matrix(data = 0.0, nrow = n, ncol = n2, dimnames = list(rownames(X)[i], rownames(X)[i2]))
        K <- 0.0
    }

    end <- 0L
    for (whichBuffer in seq_len(nBuffers)) {
        ini <- end + 1L
        end <- min(p, ini + bufferSize - 1L)

        if (verbose) {
            message("Buffer ", whichBuffer, " of ", nBuffers, " ...")
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

        # compute XX'
        if (ncol(X1) > 0L) {

            # scale and impute X
            X1 <- scale(X1, center = center.chunk, scale = scale.chunk)
            X1[is.na(X1)] <- 0L
            if (hasY) {
                X2 <- scale(X2, center = center.chunk, scale = scale.chunk)
                X2[is.na(X2)] <- 0L
            }

            if (nTasks > 1L) {
                if (!hasY) {
                    G_chunk <- crossprods(x = X1, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
                } else {
                    G_chunk <- crossprods(x = X1, y = X2, use_tcrossprod = TRUE, nTasks = nTasks, nCores = nCores)
                }
            } else {
                if (gpu) {
                    # Currently, subsets have to be converted to double as gpuR
                    # does not support integer matrices
                    storage.mode(X1) <- "double"
                    if (!hasY) {
                        G_chunk <- gpuR::tcrossprod(gpuR::gpuMatrix(X1))
                    } else {
                        storage.mode(X2) <- "double"
                        G_chunk <- gpuR::tcrossprod(x = gpuR::gpuMatrix(X1), y = gpuR::gpuMatrix(X2))
                    }
                    G_chunk <- as.matrix(G_chunk)
                } else {
                    if (!hasY) {
                        G_chunk <- tcrossprod(X1)
                    } else {
                        G_chunk <- tcrossprod(x = X1, y = X2)
                    }
                }
            }

            G[] <- G + G_chunk

            if (hasY && scaleG) {
                K <- K + ncol(X1)
            }

        }

    }

    if (!hasY && scaleG) {
        # Use seq instead of diag to avoid copy as it does not increase ref count
        K <- mean(G[seq(from = 1L, to = n * n, by = n + 1L)])
    }

    if (scaleG) {
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
#' @param nTasks The number of tasks the problem should be broken into to be
#' distributed among `nCores` cores. Defaults to `nCores`.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param gpu Whether to use the GPU to perform the computations. By setting
#' this parameter, `nCores` and `nTasks` will be set to 1.
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A [symDMatrix::symDMatrix-class] object.
#' @export
getG_symDMatrix <- function(X, center = TRUE, scale = TRUE, scaleG = TRUE, folderOut = paste0("symDMatrix_", randomString()), vmode = "double", i = seq_len(nrow(X)), j = seq_len(ncol(X)), blockSize = 5000L, nTasks = nCores, nCores = getOption("mc.cores", 2L), gpu = FALSE, verbose = FALSE) {

    if (gpu) {
        if (!requireNamespace("gpuR", quietly = TRUE)) {
            stop("gpuR needed for this function to work. Please install it.", call. = FALSE)
        }
        nCores <- 1L
        nTasks <- 1L
    }

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
    blocks <- vector(mode = "list", length = nBlocks)
    counter <- 1L
    for (r in 1L:nBlocks) {
        blocks[[r]] <- vector(mode = "list", length = nBlocks - r + 1L)
        for (s in r:nBlocks) {
            if (verbose) {
                message("Block ", r, "-", s, " ...")
            }
            blockName <- paste0("data_", padDigits(r, nBlocks), "_", padDigits(s, nBlocks), ".bin")
            block <- ff::as.ff(getG(X, center = center, scale = scale, scaleG = FALSE, i = blockIndices[[r]], j = j, i2 = blockIndices[[s]], bufferSize = blockSize, nTasks = nTasks, nCores = nCores, gpu = gpu, verbose = FALSE), filename = paste0(folderOut, "/", blockName), vmode = vmode)
            # Change ff path to a relative one
            bit::physical(block)$filename <- blockName
            blocks[[r]][[s - r + 1L]] <- block
            counter <- counter + 1L
        }
    }

    G <- symDMatrix::symDMatrix(blocks, center, scale)

    if (scaleG) {
        K <- mean(diag(G))
        for (r in seq_len(length(G@data))) {
            for (s in seq_len(length(G@data[[r]]))) {
                G@data[[r]][[s]][] <- G@data[[r]][[s]][] / K
            }
        }
    }

    save(G, file = paste0(folderOut, "/G.RData"))

    return(G)

}


#' Performs Single Marker Regressions Using BGData Objects.
#'
#' Implements single marker regressions. The regression model includes all the
#' covariates specified in the right-hand-side of the `formula` plus one column
#' of `@@geno` at a time. The data from the association tests is obtained from
#' a [BGData-class] object.
#'
#' @inheritSection BGData-package Memory-mapping
#' @inheritSection BGData-package Multi-level parallelism
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
#' physical memory for processing per core. If `NULL`, all elements in `j` are
#' used. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @param ... Additional arguments for chunkedApply and regression method.
#' @return The same matrix that would be returned by `coef(summary(model))`.
#' @example man/examples/GWAS.R
#' @export
GWAS <- function(formula, data, method = "lsfit", i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {

    if (class(data) != "BGData") {
        stop("data must BGData")
    }

    if (!method %in% c("lm", "lm.fit", "lsfit", "glm", "lmer", "SKAT", "rayOLS")) {
        stop("Only lm, lm.fit, lsfit, glm, lmer, SKAT, and rayOLS have been implemented so far.")
    }

    i <- convertIndexTypes(i, rownames(data@geno))
    j <- convertIndexTypes(j, colnames(data@geno))

    if (method == "lsfit") {
        OUT <- GWAS.lsfit(formula = formula, data = data, i = i, j = j, bufferSize = bufferSize, nCores = nCores, verbose = verbose, ...)
    } else if (method == "rayOLS") {
        if (length(attr(stats::terms(formula), "term.labels")) > 0L) {
            stop("method rayOLS can only be used with y~1 formula, if you want to add covariates pre-adjust your phenotype.")
        }
        OUT <- GWAS.rayOLS(formula = formula, data = data, i = i, j = j, bufferSize = bufferSize, nCores = nCores, verbose = verbose, ...)
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
        OUT <- chunkedApply(X = data@geno, MARGIN = 2L, FUN = function(col, ...) {
            pheno$z <- col
            fm <- FUN(GWAS.model, data = pheno, ...)
            getCoefficients(fm)
        }, i = i, j = j, bufferSize = bufferSize, nCores = nCores, verbose = verbose, ...)
        colnames(OUT) <- colnames(data@geno)[j]
        OUT <- t(OUT)
    }

    return(OUT)
}


# OLS for the regression y=xb+e (data is assumed to be pre-adjusted by non-genetic effects
# rayOLS <- function(y, x){
#    tmp <- !(is.na(x) | is.na(y))
#    x <- x[tmp]
#    y <- y[tmp]
#    x <- x - mean(x)
#    n <- length(y)
#    rhs <- sum(x * y)
#    XtX <- sum(x^2L)
#    sol <- rhs / XtX
#    error <- y - x * sol
#    vE <- sum(error^2L) / (n - 2L)
#    SE <- sqrt(vE / XtX)
#    z_stat <- sol / SE
#    p_val<-stats::pt(q = abs(z_stat), df = n - 2L, lower.tail = FALSE) * 2L
#    return(c(sol, SE, z_stat, p_val))
#}

## An alternative that is faster than rayOLS
rayOLS <- function(y, x,SSy,Int,n,isNAY){
	
    isNAX<-which(is.na(x))
    isNAXY<-unique(c(isNAX,isNAY))
    SSy<-SSy-sum(y[isNAX]^2,na.rm=T)
    
    if(length(isNAXY)>0){
      y[isNAXY]=0
      x[isNAXY]=0  
    }
    n<- n-length(isNAXY)
     
    # crossprodUCTS
    sX<-crossprod(Int,x)
    sY<-crossprod(Int,y)
    XX<-crossprod(x)-sX*sX/n
    Xy<-crossprod(x,y)-sX*sY/n 
    
    # solution and SE
    sol <- Xy/XX
    RSS<-(SSy)-XX*(sol^2)
    SE <- sqrt(RSS/(n-2L)/ XX)
    z_stat <- sol / SE
    p_val<-stats::pt(q = abs(z_stat), df = n - 2L, lower.tail = FALSE) * 2L
    return(c(sol, SE, z_stat, p_val))
}


# the GWAS method for rayOLS
GWAS.rayOLS <- function(formula, data, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {
    y <- data@pheno[i, as.character(stats::terms(formula)[[2L]]), drop = TRUE]
    y <- y - mean(y,na.rm=TRUE)
    n<-length(y)
    Int<-rep(1,n)
    SSy<-sum(y^2,na.rm=TRUE)
    isNAY<-which(is.na(y))
    res <- chunkedApply(X = data@geno, MARGIN = 2L, FUN = rayOLS, y = y , i = i, j = j,Int=Int,SSy=SSy,n=n,isNAY=isNAY,
		       bufferSize = bufferSize, nTasks = nTasks, nCores = nCores, verbose = verbose, ...)
    return(t(res))
}


GWAS.lsfit <- function(formula, data, i = seq_len(nrow(data@geno)), j = seq_len(ncol(data@geno)), bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE, ...) {

    # subset of model.frame has bizarre scoping issues
    frame <- stats::model.frame(formula = formula, data = data@pheno)[i, , drop = FALSE]
    model <- stats::model.matrix(formula, frame)
    model <- cbind(1L, model) # Reserve space for marker column

    y <- data@pheno[i, as.character(stats::terms(formula)[[2L]]), drop = TRUE]

    res <- chunkedApply(X = data@geno, MARGIN = 2L, FUN = function(col, ...) {
        model[, 1L] <- col
        fm <- stats::lsfit(x = model, y = y, intercept = FALSE)
        stats::ls.print(fm, print.it = FALSE)$coef.table[[1L]][1L, ]
    }, i = i, j = j, bufferSize = bufferSize, nCores = nCores, verbose = verbose, ...)
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

    OUT <- matrix(data = double(), nrow = length(uniqueGroups), ncol = 2L)
    colnames(OUT) <- c("nMrk", "p-value")
    rownames(OUT) <- uniqueGroups

    H0 <- SKAT::SKAT_Null_Model(formula, data = data@pheno[i, , drop = FALSE], ...)

    for (group in seq_len(length(uniqueGroups))) {
        Z <- data@geno[i, groups == uniqueGroups[group], drop = FALSE]
        fm <- SKAT::SKAT(Z = Z, obj = H0, ...)
        OUT[group, ] <- c(ncol(Z), fm$p.value)
        if (verbose) {
            message("Group ", group, " of ", length(uniqueGroups), " ...")
        }
    }

    return(OUT)
}


getCoefficients <- function(x) {
    UseMethod("getCoefficients")
}


getCoefficients.lm <- function(x) {
    summary(x)$coef[2L, ]
}


getCoefficients.glm <- function(x) {
    summary(x)$coef[2L, ]
}


getCoefficients.lmerMod <- function(x) {
    ans <- summary(x)$coef[2L, ]
    ans <- c(ans, c(1L - stats::pnorm(ans[3L])))
    return(ans)
}


#' Generates Various Summary Statistics.
#'
#' Computes the frequency of missing values, the (minor) allele frequency, and
#' standard deviation of each column of `X`.
#'
#' @inheritSection BGData-package Memory-mapping
#' @inheritSection BGData-package Multi-level parallelism
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param i Indicates which rows of `X` should be used. Can be integer,
#' boolean, or character. By default, all rows are used.
#' @param j Indicates which columns of `X` should be used. Can be integer,
#' boolean, or character. By default, all columns are used.
#' @param bufferSize The number of columns of `X` that are brought into
#' physical memory for processing per core. If `NULL`, all elements in `j` are
#' used. Defaults to 5000.
#' @param nCores The number of cores (passed to [parallel::mclapply()]).
#' Defaults to the number of cores as detected by [parallel::detectCores()].
#' @param verbose Whether progress updates will be posted. Defaults to `FALSE`.
#' @return A `data.frame` with three columns: `freq_na` for frequencies of
#' missing values, `allele_freq` for (minor) allele frequencies, and `sd` for
#' standard deviations.
#' @example man/examples/summarize.R
#' @export
summarize <- function(X, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {
    i <- convertIndexTypes(i, rownames(X))
    j <- convertIndexTypes(j, colnames(X))
    m <- chunkedApply(X = X, MARGIN = 2L, FUN = function(col) {
        freqNA <- mean(is.na(col))
        alleleFreq <- mean(col, na.rm = TRUE) / 2L
        sd <- stats::sd(col, na.rm = TRUE)
        cbind(freqNA, alleleFreq, sd)
    }, i = i, j = j, bufferSize = bufferSize, nCores = nCores, verbose = verbose)
    df <- data.frame(
        freq_na = m[1L, ],
        allele_freq = m[2L, ],
        sd = m[3L, ],
        stringsAsFactors = FALSE
    )
    rownames(df) <- colnames(X)[j]
    return(df)
}


getLineCount <- function(path, header) {
    file <- file(path, open = "r")
    n <- 0L
    while (length(readLines(file, n = 1L)) > 0L) {
        n <- n + 1L
    }
    if (header) {
        n <- n - 1L
    }
    close(file)
    return(n)
}


getFileHeader <- function(path, sep = "") {
    file <- file(path, open = "r")
    header <- scan(file, nlines = 1L, what = character(), sep = sep, quiet = TRUE)
    close(file)
    return(header)
}


getColumnCount <- function(path, sep = "") {
    header <- getFileHeader(path, sep)
    p <- length(header)
    return(p)
}


randomString <- function() {
    paste(sample(c(0L:9L, letters, LETTERS), size = 5L, replace = TRUE), collapse = "")
}


normalizeType <- function(val) {
    type <- typeof(val)
    # detect strings
    if (type == "character" && length(val) > 0L) {
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
    sample <- x[[1L]]
    if (is.matrix(sample)) {
        x <- matrix(data = unlist(x), nrow = nrow(sample), byrow = FALSE)
        rownames(x) <- rownames(sample)
    } else {
        x <- unlist(x)
    }
    return(x)
}


padDigits <- function(x, total) {
    formatC(x, width = as.integer(log10(total) + 1L), format = "d", flag = "0")
}


loadExample <- function() {
    path <- system.file("extdata", package = "BGData")
    message("Loading chromosomes as BED files...")
    m <- do.call("ColumnLinkedMatrix", lapply(c("chr1", "chr2", "chr3"), function(chr) {
        suppressMessages(BEDMatrix::BEDMatrix(paste0(path, "/", chr)))
    }))
    as.BGData(m, alternatePhenotypeFile = paste0(path, "/pheno.txt"))
}


convertIndexTypes <- function(index, names) {
    if (is.logical(index)) {
        index <- which(index)
    } else if (is.character(index)) {
        index <- match(index, names)
    }
    return(index)
}

chunkRanges <- function(a, n, i = NULL) {
    k <- as.integer(a / n)
    r <- as.integer(a %% n)
    range <- function(i, k, r) {
        c((i - 1) * k + min(i - 1, r) + 1, i * k + min(i, r))
    }
    if (!is.null(i)) {
        range(i, k, r)
    } else {
        sapply(seq_len(n), range, k, r)
    }
}


simplifyList <- function(x) {
    sample <- x[[1]]
    if (is.vector(sample)) {
        x <- unlist(x)
    } else if (is.matrix(sample)) {
        x <- matrix(unlist(x), nrow = nrow(sample), byrow = FALSE)
        rownames(x) <- rownames(sample)
    }
    return(x)
}


#' Applies a function on each row or column of a matrix in parallel.
#' 
#' The input matrix \code{X} is broken into \code{nTasks} chunks and passed to
#' \code{\link[parallel]{mclapply}}. The number of cores can be configured
#' using \code{mc.cores}. Uses \code{apply} from base internally.
#' 
#' @param X A matrix.
#' @param MARGIN The subscripts which the function will be applied over. 1
#'   indicates rows, 2 indicates columns.
#' @param FUN The function to be applied.
#' @param nTasks The number of submatrices of \code{X} to be processed in
#'   parallel.
#' @param mc.cores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}).
#' @param ... Additional arguments to be passed to \code{apply}.
#' @export
parallelApply <- function(X, MARGIN, FUN, nTasks = detectCores(), mc.cores = detectCores(), ...) {
    d <- dim(X)
    if (!length(d)) {
        stop("dim(X) must have a positive length")
    }
    res <- parallel::mclapply(X = seq_len(nTasks), FUN = function(i, ...) {
        range <- chunkRanges(ifelse(MARGIN == 2, ncol(X), nrow(X)), nTasks, i)
        if (MARGIN == 2) {
            subset <- X[, seq(range[1], range[2]), drop = FALSE]
        } else {
            subset <- X[seq(range[1], range[2]), , drop = FALSE]
        }
        base::apply(subset, MARGIN, FUN, ...)
    }, ..., mc.cores = mc.cores)
    simplifyList(res)
}


#' Reads chunks of data into memory and applies a function on each row or
#' column of a matrix.
#' 
#' \code{chunkedApply} uses \code{parallelApply} internally, so \code{nTasks}
#' and \code{mc.cores} can be passed as \code{...}.
#' 
#' @param X A matrix.
#' @param MARGIN The subscripts which the function will be applied over. 1
#'   indicates rows, 2 indicates columns.
#' @param FUN The function to be applied.
#' @param bufferSize The number of rows or columns of \code{X} that are brought
#'   into memory for processing.
#' @param verbose Whether to print additional information.
#' @param ... Additional arguments to be passed to \code{parallelApply}.
#' @export
chunkedApply <- function(X, MARGIN, FUN, bufferSize, verbose = FALSE, ...) {
    d <- dim(X)
    if (!length(d)) {
        stop("dim(X) must have a positive length")
    }
    nChunks <- ceiling(d[MARGIN] / bufferSize)
    ranges <- chunkRanges(d[MARGIN], nChunks)
    res <- lapply(seq_len(nChunks), function(i) {
        if (verbose) {
            cat("Processing chunk ", i, " of ", nChunks, " (", round(i / nChunks * 100, 3), "%) ...", "\n", sep = "")
        }
        if (MARGIN == 2) {
            subset <- X[, seq(ranges[1, i], ranges[2, i]), drop = FALSE]
        } else {
            subset <- X[seq(ranges[1, i], ranges[2, i]), , drop = FALSE]
        }
        parallelApply(X = subset, MARGIN = MARGIN, FUN = FUN, ...)
    })
    simplifyList(res)
}


# Performs crossprod() or tcrossprod() for a chunk (set of columns or sets of
# rows) of x.
crossprods.chunk <- function(chunk, x, y = NULL, nChunks, use_tcrossprod = FALSE) {
    if (!is.null(y)) {
        y <- as.matrix(y)
        nY <- ifelse(use_tcrossprod, ncol(y), nrow(y))
        ranges <- chunkRanges(nY, nChunks, chunk)
        if (use_tcrossprod) {
            y <- y[, seq(ranges[1], ranges[2]), drop = FALSE]
        } else {
            y <- y[seq(ranges[1], ranges[2]), , drop = FALSE]
        }
    }
    nX <- ifelse(use_tcrossprod, ncol(x), nrow(x))
    ranges <- chunkRanges(nX, nChunks, chunk)
    if (use_tcrossprod) {
        X <- x[, seq(ranges[1], ranges[2]), drop = FALSE]
        Xy <- tcrossprod(X, y)
    } else {
        X <- x[seq(ranges[1], ranges[2]), , drop = FALSE]
        Xy <- crossprod(X, y)
    }
    return(Xy)
}


crossprods <- function(x, y = NULL, nChunks = detectCores(), use_tcrossprod = FALSE, mc.cores = detectCores()) {
    if (!is.null(y)) {
        if (use_tcrossprod) {
            if (ncol(x) != ncol(y)) {
                stop("Error in tcrossprod.parallel: non-conformable arguments.")
            }
        } else {
            if (nrow(x) != nrow(y)) {
                stop("Error in crossprod.parallel: non-conformable arguments.")
            }
        }
    }

    # Computes crossprod(x,y) or tcrossprod(x,y)
    if (nChunks == 1) {
        if (use_tcrossprod) {
            Xy <- tcrossprod(x, y)
        } else {
            Xy <- crossprod(x, y)
        }
    } else {
        TMP <- mclapply(X = 1:nChunks, FUN = crossprods.chunk, x = x, y = y, nChunks = nChunks, use_tcrossprod = use_tcrossprod, mc.cores = mc.cores)
        ## We now need to add up chunks sequentially
        Xy <- TMP[[1]]
        if (length(TMP) > 1) {
            for (i in 2:length(TMP)) {
                Xy <- Xy + TMP[[i]]
            }
        }
    }
    return(Xy)
}


#' Computes crossprod (x'y or x'x) in parallel.
#' 
#' @param x matrix, ff_matrix, RowLinkedMatrix or ColumnLinkedMatrix.
#' @param y vector, matrix, ff_matrix, RowLinkedMatrix or ColumnLinkedMatrix. NULL by default.
#' @param nChunks The number of chunks used when X and y are partitioned.
#' @param mc.cores The number of cores (passed to \code{\link[parallel]{mclapply}}).
#' @return x'y' or x'x depending on whether y is provided.
#' @export
crossprod.parallel <- function(x, y = NULL, nChunks = detectCores(), mc.cores = detectCores()) {
    ans <- crossprods(x = x, y = y, nChunks = nChunks, mc.cores = mc.cores, use_tcrossprod = FALSE)
    return(ans)
}


#' Computes tcrossprod (xy' or xx') in parallel.
#' 
#' @param x matrix, ff_matrix, RowLinkedMatrix or ColumnLinkedMatrix.
#' @param y vector, matrix, ff_matrix, RowLinkedMatrix or ColumnLinkedMatrix. NULL by default.
#' @param nChunks The number of chunks used when X and y are partitioned.
#' @param mc.cores The number of cores (passed to \code{\link[parallel]{mclapply}}).
#' @return xy' or xx' depending on whether y is provided.
#' @export
tcrossprod.parallel <- function(x, y = NULL, nChunks = detectCores(), mc.cores = detectCores()) {
    ans <- crossprods(x = x, y = y, nChunks = nChunks, mc.cores = mc.cores, use_tcrossprod = TRUE)
    return(ans)
}



#' Computes a genomic relationship matrix G=xx'.
#' 
#' Offers options for centering and scaling the columns of x before computing
#' xx'. If \code{centerCol=FALSE}, \code{scaleCol=FALSE} and
#' \code{scaleG=FALSE}, \code{getG} produces the same outcome than
#' \code{tcrossprod}.
#' 
#' @param x matrix, ff_matrix, RowLinkedMatrix or ColumnLinkedMatrix
#' @param nChunks The number of columns that are processed at a time.
#' @param scaleCol TRUE/FALSE whether columns must be scaled before computing
#'   xx'.
#' @param scaleG TRUE/FALSE whether columns must be scaled before computing xx'.
#' @param verbose If TRUE more messages are printed.
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be 
#'   used. By default, all columns are used.
#' @param i2 (integer, boolean or character) Indicates which rows should be used 
#'   to divide matrix into blocks.
#' @param minVar Columns with variance lower than this value will not be used 
#'   in the computation (only if \code{scaleCol} is set).
#' @param nChunks2 The number of chunks that each chunk is split into for 
#'   processing in parallel.
#' @param scales Precomputed scales if i2 is used.
#' @param centers Precomputed centers if i2 is used.
#' @param impute Whether to impute data if i2 is used.
#' @param saveG Whether to save genomic relationship matrix into file.
#' @param saveType File format to save genomic relationship matrix in. Either 
#'   \code{RData} or \code{ff}.
#' @param saveName Name without extension to save genomic relationship matrix with.
#' @param mc.cores The number of cores (passed to \code{\link[parallel]{mclapply}}).
#' @return A positive semi-definite symmetric numeric matrix.
#' @export
getG <- function(x, nChunks = ceiling(ncol(x) / 10000), scaleCol = TRUE, scaleG = TRUE, verbose = TRUE, i = 1:nrow(x), j = 1:ncol(x), i2 = NULL, minVar = 1e-05, nChunks2 = detectCores(), scales = NULL, centers = NULL, impute = TRUE, saveG = FALSE, saveType = "RData", saveName = "Gij", mc.cores = detectCores()) {
    if (is.null(i2)) {
        G <- getGi(x = x, nChunks = nChunks, scaleCol = scaleCol, scaleG = scaleG, verbose = verbose, i = i, j = j, minVar = minVar, nChunks2 = nChunks2, mc.cores = mc.cores)
    } else {
        if (is.null(scales) || is.null(centers)) stop("scales and centers need to be precomputed.")
        G <- getGij(x = x, i1 = i, i2 = i2, scales = scales, centers = centers, scaleCol = scaleCol, scaleG = scaleG, verbose = verbose, nChunks = nChunks, j = j, minVar = minVar, nChunks2 = nChunks2, impute = impute, mc.cores = mc.cores)
    }
    if (saveG) {
        if (saveType == "RData") {
            save(G, file = paste0(saveName, ".RData"))
        }
        if (saveType == "ff") {
            Gij <- as.ff(G, file = paste0(saveName, ".bin"))
            save(Gij, file = paste0(saveName, ".ff"))
        }
    }
    return(G)
}


getGi <- function(x, nChunks = ceiling(ncol(x) / 10000), scaleCol = TRUE, scaleG = TRUE, verbose = TRUE, i = 1:nrow(x), j = 1:ncol(x), minVar = 1e-05, nChunks2 = detectCores(), mc.cores = detectCores()) {
    nX <- nrow(x)
    pX <- ncol(x)

    # converting boolean to integer index (it leads to a more efficient subsetting
    # than booleans)
    if (is.logical(i)) {
        i <- which(i)
    }
    if (is.logical(j)) {
        j <- which(j)
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

    tmp <- x[i, 1:2]
    n <- nrow(tmp)

    G <- matrix(0, nrow = n, ncol = n)
    rownames(G) <- rownames(tmp)
    colnames(G) <- rownames(G)

    end <- 0
    delta <- ceiling(p / nChunks)

    for (k in 1:nChunks) {
        ini <- end + 1
        if (ini <= p) {
            end <- min(p, ini + delta - 1)
            if (verbose) {
                cat("Chunk: ", k, " (markers ", ini, ":", end, " ~", round(100 * end / p, 1), "% done)\n", sep = "")
                cat("  =>Acquiring genotypes...\n")
            }

            # subset
            tmp <- j[ini:end]
            X <- x[i, tmp, drop = FALSE]

            if (scaleCol) {
                VAR <- apply(X = X, FUN = var, MARGIN = 2, na.rm = TRUE)
                tmp <- which(VAR < minVar)
                if (length(tmp) > 0) {
                  X <- X[, -tmp]
                  VAR <- VAR[-tmp]
                }
            }

            if (ncol(X) > 0) {
                if (verbose) {
                  cat("  =>Computing...\n")
                }
                X <- scale(X, center = TRUE, scale = scaleCol)
                TMP <- is.na(X)
                if (any(TMP)) {
                  X <- ifelse(TMP, 0, X)
                }

                if (nChunks2 > 1) {
                  TMP <- crossprods(x = X, use_tcrossprod = TRUE, nChunks = nChunks2, mc.cores = mc.cores)
                } else {
                  TMP <- tcrossprod(X)
                }
                G <- G + TMP
            }
        }
    }
    if (scaleG) {
        tmp <- mean(diag(G))
        G <- G / tmp
    }

    return(G)
}


getGij <- function(x, i1, i2, scales, centers, scaleCol = TRUE, scaleG = TRUE, verbose = TRUE, nChunks = ceiling(ncol(x) / 10000), j = 1:ncol(x), minVar = 1e-05, nChunks2 = detectCores(), impute = TRUE, mc.cores = detectCores()) {

    nX <- nrow(x)
    pX <- ncol(x)
    K <- 0

    # need to make this more general, convert character to boolean, booleand to
    # integer
    if (is.logical(i1)) {
        i1 <- which(i1)
    }
    if (is.logical(i2)) {
        i1 <- which(i2)
    }
    if (is.logical(j)) {
        j <- which(j)
    }

    n1 <- length(i1)
    p <- length(j)
    n2 <- length(i2)

    if ((min(i1) < 1) | (max(i1) > nX)) {
        stop("Index out of bounds")
    }
    if ((min(i2) < 1) | (max(i2) > nX)) {
        stop("Index out of bounds")
    }
    if ((min(j) < 1) | (max(j) > pX)) {
        stop("Index out of bounds")
    }

    G <- matrix(nrow = n1, ncol = n2, 0)
    tmp <- rownames(x)
    rownames(G) <- tmp[i1]
    colnames(G) <- tmp[i2]

    end <- 0
    delta <- ceiling(p / nChunks)

    for (k in 1:nChunks) {
        ini <- end + 1
        if (ini <= p) {
            end <- min(p, ini + delta - 1)
            if (verbose) {
                cat("Working with chunk: ", k, " (markers ", ini, ":", end, " ~", round(100 * ini / p, 1), "% done)\n", sep = "")
                cat("  =>Acquiring genotypes...\n")
            }

            # subset
            tmpCol <- j[ini:end]
            # K<-K+length(tmpCol)
            X1 <- x[i1, tmpCol, drop = FALSE]
            X2 <- x[i2, tmpCol, drop = FALSE]
            centers.chunk <- centers[tmpCol]
            scales.chunk <- scales[tmpCol]

            if (scaleCol) {
                tmp <- which(scales.chunk < sqrt(minVar))

                if (length(tmp) > 0) {
                  # K<-K-length(tmp)
                  X1 <- X1[, -tmp]
                  X2 <- X2[, -tmp]
                  scales.chunk <- scales.chunk[-tmp]
                  centers.chunk <- centers.chunk[-tmp]
                }
            }

            if (ncol(X1) > 0) {
                if (!scaleCol) {
                  scales.chunk <- FALSE
                }
                if (verbose) {
                  cat("  =>Computing...\n")
                }
                X1 <- scale(X1, center = centers.chunk, scale = scales.chunk)
                TMP <- is.na(X1)
                if (any(TMP)) {
                  X1 <- ifelse(TMP, 0, X1)
                }
                X2 <- scale(X2, center = centers.chunk, scale = scales.chunk)
                TMP <- is.na(X2)
                if (any(TMP)) {
                  X2 <- ifelse(TMP, 0, X2)
                }
                TMP <- tcrossprod.parallel(x = X1, y = X2, mc.cores = mc.cores, nChunks = nChunks2)
                G <- G + TMP
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
        G <- G / K
    }

    return(G)
}


#' @export
getG.symDMatrix <- function(X, nChunks = 5, chunkSize = NULL, centers = NULL, scales = NULL,
    centerCol = T, scaleCol = T, nChunks2 = 1, folder = randomString(5), vmode = "double",
    verbose = TRUE, saveRData = TRUE, mc.cores = 1, scaleG = T) {

    timeIn <- proc.time()[3]
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(chunkSize)) {
        chunkSize <- ceiling(n/nChunks)
    }

    if ((centerCol | scaleCol) & (is.null(centers) | is.null(scales))) {
        if (is.null(centers) & is.null(scales)) {
            centers <- rep(NA, p)
            scales <- rep(NA, p)
            for (i in 1:p) {
                xi <- X[, i]
                scales[i] <- sd(xi, na.rm = TRUE) * sqrt((n - 1)/n)
                centers[i] <- mean(xi, na.rm = TRUE)
            }
        }
        if ((!is.null(centers)) & (is.null(scales))) {
            scales <- rep(NA, p)
            for (i in 1:p) {
                xi <- X[, i]
                scales[i] <- sd(xi, na.rm = TRUE) * sqrt((n - 1)/n)
            }
        }
        if ((is.null(centers)) & (!is.null(scales))) {
            centers <- rep(NA, p)
            for (i in 1:p) {
                xi <- X[, i]
                centers[i] <- mean(xi, na.rm = TRUE)
            }
        }
    }

    if (!centerCol){ centers <- rep(0, p) }
    if (!scaleCol){ scales <- rep(1, p) }

    chunkID <- ceiling(1:n/chunkSize)
    nChunks <- max(chunkID)
    nFiles <- nChunks * (nChunks + 1)/2
    DATA <- list()
    counter <- 1

    tmpDir <- getwd()
    dir.create(folder)
    setwd(folder)

    for (i in 1:nChunks) {
    
    if (verbose) {
                cat(" Working pair ", i, "-", j, " (", round(100 * counter/(nChunks *
                  (nChunks + 1)/2)), "% ", round(proc.time()[3] - timeIn, 3), " seconds).\n",
                  sep = "")
    }
        DATA[[i]] <- list()
        rowIndex_i <- which(chunkID == i)
        Xi <- X[rowIndex_i, ]

        # centering/scaling
        for (k in 1:p) {
            xik <- Xi[, k]
            xik <- (xik - centers[k])/scales[k]
            xik[is.na(xik)] <- 0
            Xi[, k] <- xik
        }

        for (j in i:nChunks) {
            rowIndex_j <- which(chunkID == j)
            Xj <- X[rowIndex_j, ]

            # centering/scaling
            for (k in 1:p) {
                xjk <- Xj[, k]
                xjk <- (xjk - centers[k])/scales[k]
                xjk[is.na(xjk)] <- 0
                Xj[, k] <- xjk
            }

            Gij <- tcrossprod.parallel(x = Xi, y = Xj, mc.cores = mc.cores, nChunks = nChunks2)

            DATA[[i]][[j - i + 1]] <- ff(dim = dim(Gij), vmode = vmode, initdata = as.vector(Gij),
                                         filename = paste0("data_", i, "_", j, ".bin"))
            colnames(DATA[[i]][[j - i + 1]]) <- colnames(X)[rowIndex_j]
            rownames(DATA[[i]][[j - i + 1]]) <- rownames(X)[rowIndex_i]
            counter <- counter + 1
            physical(DATA[[i]][[j - i + 1]])$pattern <- "ff"
            physical(DATA[[i]][[j - i + 1]])$filename <- paste0("data_", i, "_",
                j, ".bin")

            if (verbose) {
                cat("  =>Done ")
            }
        }
    }
    if (is.null(rownames(X)))
        rownames(X) <- 1:n
    names(centers) <- colnames(X)
    names(scales) <- colnames(X)
    G <- new("symDMatrix", names = rownames(X), data = DATA, centers = centers, scales = scales)
    if (scaleG) {
        K <- mean(diag(G))
        for (i in 1:length(G@data)) {
            for (j in 1:length(G@data[[i]])) {
                G@data[[i]][[j]][] <- G@data[[i]][[j]][]/K
            }
        }
    }
    if (saveRData) {
        save(G, file = "G.RData")
    }
    setwd(tmpDir)
    return(G)
}


#' Performs single marker regressions using a \code{\linkS4class{BGData}} 
#' object.
#' 
#' Implements single marker regressions. The regression model includes all the
#' covariates specified in the right-hand-side of the \code{formula} plus one
#' column of \code{@@geno}, one column at a time. The data from the association
#' tests is obtained from a \code{\linkS4class{BGData}} object.
#' 
#' @param formula A formula (e.g. weight~sex+age) with the response on the 
#'   left-hand side and predictors (all the covariates except the markers) on 
#'   the right-hand side. The variables included in the formula must be in the 
#'   \code{@@pheno} object of the \code{\linkS4class{BGData}}.
#' @param data A \code{\linkS4class{BGData}} object.
#' @param method The regression method to be used. Currently, the following 
#'   methods are implemented: \code{\link{lm}}, \code{\link{lm.fit}}, 
#'   \code{\link{lsfit}}, \code{\link{glm}} and \code{\link[lme4]{lmer}}.
#' @param plot If TRUE a Manhattan plot is produced and filled with points as 
#'   the association tests are run.
#' @param verbose If TRUE more messages are printed.
#' @param min.pValue Numeric, the minimum p-value expected, used to determine 
#'   the limits of the vertical axis of the Manhattan plot.
#' @param chunkSize Represents the number of columns of \code{@@geno} that are 
#'   brought into RAM for processing (5000 by default).
#' @param ... Optional arguments for chunkedApply and regression method.
#' @return Returns a matrix with estimates, SE, p-value, etc.
#' @export
GWAS <- function(formula, data, method, plot = FALSE, verbose = FALSE, min.pValue = 1e-10, chunkSize = 5000, ...) {
    if (class(data) != "BGData") {
        stop("data must BGData")
    }

    if (!method %in% c("lm", "lm.fit", "lsfit", "glm", "lmer", "SKAT")) {
        stop("Only lm, glm, lmer and SKAT have been implemented so far.")
    }
    ## We can have specialized methods, for instance for OLS it is better to use
    ## lsfit (that is what GWAS.ols does)
    if (method %in% c("lm", "lm.fit", "lsfit", "SKAT")) {
        if (method %in% c("lm", "lm.fit", "lsfit")) {
            OUT <- GWAS.ols(formula = formula, data = data, plot = plot, verbose = verbose, min.pValue = min.pValue, chunkSize = chunkSize, ...)
        }
        if (method == "SKAT") {
            OUT <- GWAS.SKAT(formula = formula, data = data, plot = plot, verbose = verbose, min.pValue = min.pValue, ...)
        }
    } else {
        if (method == "lmer" && !requireNamespace("lme4", quietly = TRUE)) {
            stop("lme4 needed for this function to work. Please install it.", call. = FALSE)
        }
        if (method == "lmer") {
            FUN <- lme4::lmer
        } else {
            FUN <- match.fun(method)
        }
        # could subset based on NAs so that subsetting does not take place in each
        # iteration of the GWAS loop
        pheno <- data@pheno
        fm <- FUN(formula, data = pheno, ...)
        tmp <- getCoefficients(fm)
        p <- ncol(data@geno)
        OUT <- matrix(nrow = p, ncol = length(tmp), NA)
        rownames(OUT) <- colnames(data@geno)
        colnames(OUT) <- colnames(tmp)
        GWAS.model <- update(as.formula(formula), ".~z+.")
        if (plot) {
            tmp <- paste(as.character(GWAS.model[2]), as.character(GWAS.model[3]), sep = "~")
            plot(numeric() ~ numeric(), xlim = c(0, p), ylim = c(0, -log(min.pValue, base = 10)), ylab = "-log(p-value)", xlab = "Marker", main = tmp)
        }
        nChunks <- ceiling(p / chunkSize)
        end <- 0
        tmpRow <- 0

        for (i in 1:nChunks) {
            time.in <- proc.time()[3]
            ini <- end + 1
            end <- min(ini + chunkSize - 1, p)
            Z <- data@geno[, ini:end, drop = FALSE]

            for (j in 1:(end - ini + 1)) {
                pheno$z <- Z[, j]
                fm <- FUN(GWAS.model, data = pheno, ...)
                tmp <- getCoefficients(fm)
                tmpRow <- tmpRow + 1
                OUT[tmpRow, ] <- tmp
                if (plot) {
                  x <- c(tmpRow - 1, tmpRow)
                  y <- -log(OUT[c(tmpRow - 1, tmpRow), 4], base = 10)
                  if (tmpRow > 1) {
                    lines(x = x, y = y, col = 8, lwd = 0.5)
                  }
                  points(y = -log(tmp[4], base = 10), col = 2, cex = 0.5, x = tmpRow)
                }
            }
            if (verbose) {
                cat(sep = "", "Chunk ", i, " of ", nChunks, " (", round(proc.time()[3] - time.in, 2), " seconds / chunk, ", round(i / nChunks * 100, 3), "% done )\n")
            }
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
    ans <- c(ans, c(1 - pnorm(ans[3])))
    return(ans)
}


# GWAS 'Ordinary Least Squares' (e.g., lsfit, lm.fit, lm)
# formula: the formula for the GWAS model without including the marker, e.g.
# y~1 or y~factor(sex)+age
# all the variables in the formula must be in data@pheno data (BGData)
# containing slots @pheno and @geno
GWAS.ols <- function(formula, data, plot = FALSE, verbose = FALSE, min.pValue = 1e-10, chunkSize = 10, ...) {

    X <- model.matrix(formula, data@pheno)
    X <- X[match(rownames(data@pheno), rownames(X)), ]
    X <- cbind(0, X) # Reserve space for marker column

    y <- data@pheno[, as.character(terms(formula)[[2]])]

    res <- chunkedApply(data@geno, 2, function(col) {
        X[, 1] <- col
        fm <- lsfit(x = X, y = y, intercept = FALSE)
        ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
    }, bufferSize = chunkSize, verbose = verbose, ...)
    colnames(res) <- colnames(data@geno)
    res <- t(res)

    if (plot) {
        title <- paste(as.character(formula[2]), as.character(formula[3]), sep = "~")
        plot(numeric() ~ numeric(), xlim = c(0, ncol(data@geno)), ylim = c(0, -log(min.pValue, base = 10)), ylab = "-log(p-value)", xlab = "Marker", main = title)
        for (i in seq_len(nrow(res))) {
            row <- res[i, ]
            x <- c(i - 1, i)
            y <- -log(res[c(i - 1, i), 4], base = 10)
            if (i > 1) {
                lines(x = x, y = y, col = 8, lwd = 0.5)
            }
            points(y = -log(row[4], base = 10), col = 2, cex = 0.5, x = i)
        }
    }

    return(res)
}


# formula: the formula for the GWAS model without including the markers, e.g.
# y~1 or y~factor(sex)+age
# all the variables in the formula must be in data@pheno data (BGData)
# containing slots @pheno and @geno
# groups: a vector mapping markers into groups (can be integer, character or
# factor)
GWAS.SKAT <- function(formula, data, groups, plot = FALSE, verbose = FALSE, min.pValue = 1e-10, ...) {

    if (!requireNamespace("SKAT", quietly = TRUE)) {
        stop("SKAT needed for this function to work. Please install it.", call. = FALSE)
    }

    p <- length(unique(groups))

    OUT <- matrix(nrow = p, ncol = 2, NA)
    colnames(OUT) <- c("nMrk", "p-value")
    levels <- unique(groups)
    rownames(OUT) <- levels

    H0 <- SKAT::SKAT_Null_Model(formula, data = data@pheno, ...)

    if (plot) {
        tmp <- paste(as.character(formula[2]), as.character(formula[3]), sep = "~")
        plot(numeric() ~ numeric(), xlim = c(0, p), ylim = c(0, -log(min.pValue, base = 10)), ylab = "-log(p-value)", xlab = "Marker", main = tmp)
    }

    for (i in 1:p) {
        time.in <- proc.time()[3]
        Z <- data@geno[, groups == levels[i], drop = FALSE]
        fm <- SKAT::SKAT(Z = Z, obj = H0, ...)
        OUT[i, ] <- c(ncol(Z), fm$p.value)

        if (plot) {
            tmp.x <- c(i - 1, i)
            tmp.y <- -log(OUT[tmp.x, 2], base = 10)
            if (i > 1) {
                lines(x = tmp.x, y = tmp.y, col = 8, lwd = 0.5)
            }
            points(y = tmp.y[2], col = 2, cex = 0.5, x = i)
        }
        if (verbose) {
            cat(sep = "", "Group ", i, " of ", p, " (", round(proc.time()[3] - time.in, 2), " seconds / chunk, ", round(i / p * 100, 3), "% done )\n")
        }
    }
    return(OUT)
}


#' Calculate frequencies of missing values and alleles.
#' 
#' @param X matrix, ff_matrix, RowLinkedMatrix or ColumnLinkedMatrix
#' @param verbose If TRUE more messages are printed.
#' @param bufferSize Represents the number of columns of \code{@@geno} that are 
#'   brought into RAM for processing (5000 by default).
#' @param nTasks Represents the number of parallel tasks each buffer is split
#'   into.
#' @param mc.cores The number of cores (passed to
#'   \code{\link[parallel]{mclapply}}).
#' @export
summarize <- function(X, verbose = FALSE, bufferSize = 5000, nTasks = detectCores(), mc.cores = detectCores()) {
    res <- chunkedApply(X, 2, function(col) {
        freqNA <- mean(is.na(col))
        allFreq <- mean(col, na.rm = TRUE) / 2
        cbind(freqNA, allFreq)
    }, bufferSize = bufferSize, verbose = verbose, nTasks = nTasks, mc.cores = mc.cores)
    rownames(res) <- c("freq_na", "freq_all")
    colnames(res) <- colnames(X)
    t(res)
}


#' Generate and store a simulated plaintext raw PED file (see \code{--recodeA}
#' in PLINK) or PED-like file for testing purposes.
#' 
#' @param filename The path where to save the generated file.
#' @param n The number of observations to generate.
#' @param p The number of markers to generate.
#' @param genoChars The alphabet used to generate the genotypes.
#' @param na.string The symbol used to denote missing values.
#' @param propNA The probability of generating NAs.
#' @param returnGenos Whether to return the genotypes from the function.
#' @export
simPED <- function(filename, n, p, genoChars = 0:2, na.string = NA, propNA = 0.02, returnGenos = FALSE) {
    if (file.exists(filename)) {
        stop(paste("File", filename, "already exists. Please move it or pick a different name."))
    }
    markerNames <- paste0("mrk_", 1:p)
    subjectNames <- paste0("id_", 1:n)
    if (returnGenos) {
        OUT <- matrix(nrow = n, ncol = p, NA)
        colnames(OUT) <- markerNames
        rownames(OUT) <- subjectNames
    }
    fileOut <- file(filename, open = "w")
    pedP <- 6 + p
    header <- c(c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"), markerNames)
    write(header, ncolumns = pedP, append = TRUE, file = fileOut)
    for (i in 1:n) {
        geno <- sample(genoChars, size = p, replace = TRUE)
        geno[runif(p) < propNA] <- na.string
        pheno <- c(0, subjectNames[i], rep(NA, 4))
        x <- c(pheno, geno)
        write(x, ncolumns = pedP, append = TRUE, file = fileOut)
        if (returnGenos) {
            OUT[i, ] <- geno
        }
    }
    close(fileOut)
    if (returnGenos) {
        return(OUT)
    }
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


getFileHeader <- function(path) {
    file <- gzfile(path, open = "r")
    header <- scan(file, nlines = 1, what = character(), quiet = TRUE)
    close(file)
    return(header)
}


getColumnCount <- function(path) {
    header <- getFileHeader(path)
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

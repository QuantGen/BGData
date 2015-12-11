library(parallel)

context("utils")

hasCores <- function(numCores) {
    if (Sys.getenv("_R_CHECK_LIMIT_CORES_") == TRUE || numCores > detectCores()) {
        skip("Not enough cores or number of cores capped for CRAN submission checks.")
    }
}

for (nCores in seq_len(2)) {

    test_that(paste("crossprod.parallel", "on", nCores, "cores"), {

        hasCores(nCores)

        W <- matrix(nrow = 10, ncol = 20, rnorm(200))
        Z <- matrix(nrow = 10, ncol = 2, rnorm(20))

        # Testing X'X
        TMP <- crossprod(W)

        TMP2 <- crossprod.parallel(W, nChunks = 3, mc.cores = nCores)
        expect_equal(TMP, TMP2)

        TMP2 <- crossprod.parallel(W, nChunks = 1, mc.cores = nCores)
        expect_equal(TMP, TMP2)

        # Testing X'y
        TMP <- crossprod(W, y = Z)

        TMP2 <- crossprod.parallel(W, y = Z, nChunks = 3, mc.cores = nCores)
        expect_equal(TMP, TMP2)

        TMP2 <- crossprod.parallel(W, y = Z, nChunks = 1, mc.cores = nCores)
        expect_equal(TMP, TMP2)

    })


    test_that(paste("tcrossprod.parallel", "on", nCores, "cores"), {

        hasCores(nCores)

        W <- matrix(nrow = 10, ncol = 20, rnorm(200))
        Z <- matrix(nrow = 5, ncol = 20, rnorm(100))

        # Testing XX'
        TMP <- tcrossprod(W)

        TMP2 <- tcrossprod.parallel(W, nChunks = 3, mc.cores = nCores)
        expect_equal(TMP, TMP2)

        TMP2 <- tcrossprod.parallel(W, nChunks = 1, mc.cores = nCores)
        expect_equal(TMP, TMP2)

        # Testing XY'
        TMP <- tcrossprod(W, y = Z)

        TMP2 <- tcrossprod.parallel(W, y = Z, nChunks = 3, mc.cores = nCores)
        expect_equal(TMP, TMP2)

        TMP2 <- tcrossprod.parallel(W, y = Z, nChunks = 1, mc.cores = nCores)
        expect_equal(TMP, TMP2)

    })

    test_that(paste("getGi", "on", nCores, "cores"), {

        hasCores(nCores)

        n <- 10
        p <- 100
        X <- matrix(nrow = n, ncol = p, data = rnorm(n * p))

        for (nChunks in c(1, 3)) {
            for (nChunks2 in c(1, 3)) {

                # all scalings
                G <- tcrossprod(scale(X))
                G <- G/mean(diag(G))
                G2 <- getG(x = X, scaleG = T, scaleCol = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G, G2)

                # without scaling to average diagonal=1
                G <- tcrossprod(scale(X))
                G2 <- getG(x = X, scaleG = F, scaleCol = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G, G2)

                # without scaling columns, but scaling average diagonal =1
                G <- tcrossprod(scale(X, center = T, scale = F))
                G <- G/mean(diag(G))
                G2 <- getG(x = X, scaleG = T, scaleCol = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G, G2)

                # no scaling at all
                G <- tcrossprod(scale(X, center = T, scale = F))
                G2 <- getG(x = X, scaleG = F, scaleCol = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G, G2)

            }
        }

        X[sample(1:length(X), size = 20)] <- NA
        G <- getG(X, verbose = F)
        expect_true(!any(is.na(G)))

    })


    test_that(paste("getGij", "on", nCores, "cores"), {

        hasCores(nCores)

        n <- 10
        p <- 100
        X <- matrix(nrow = n, ncol = p, data = rnorm(n * p))

        for (nChunks in c(1, 3)) {
            for (nChunks2 in c(1, 3)) {

                i <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)

                centers <- colMeans(X)
                scales <- apply(X, 2, sd) * sqrt((n - 1)/n)

                # all scalings
                G <- tcrossprod(scale(X))
                G <- G/mean(diag(G))
                G_12 <- getG(x = X, i = i, i2 = i2, centers = centers, scales = scales, scaleG = T, verbose = F, scaleCol = TRUE, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i2], G_12)

                G_12 <- getG(x = X, i = i, i2 = i, centers = centers, scales = scales, scaleG = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i], G_12)

                # without scaling to average diagonal = 1
                G <- tcrossprod(scale(X) * sqrt(n/(n - 1)))
                G_12 <- getG(x = X, i = i, i2 = i2, centers = centers, scales = scales, scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i2], G_12)

                G_12 <- getG(x = X, i = i, i2 = i, centers = centers, scales = scales, scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i], G_12)

                # without scaling columns, but scaling average diagonal = 1
                scales <- rep(1, ncol(X))

                G <- tcrossprod(scale(X, center = T, scale = F))
                G <- G/ncol(X)
                G_12 <- getG(x = X, i = i, i2 = i2, centers = centers, scales = scales, scaleG = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i2], G_12)

                G_12 <- getG(x = X, i = i, i2 = i, centers = centers, scales = scales, scaleG = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i], G_12)

                # no scaling at all
                G <- tcrossprod(scale(X, center = T, scale = F))
                G_12 <- getG(x = X, i = i, i2 = i2, centers = centers, scales = scales, scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i2], G_12)

                G_12 <- getG(x = X, i = i, i2 = i, centers = centers, scales = scales, scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i, i], G_12)

            }
        }
    })


    test_that(paste("summarize", "on", nCores, "cores"), {

        hasCores(nCores)

        genotypes <- matrix(nrow = 3, ncol = 6, c(0, 0, 1, 0, 2, 2, 1, 2, 0, 1, 2, 0, 0, 1, 2, 0, NA, 0))

        dummy <- matrix(nrow = ncol(genotypes), ncol = 2, NA)
        colnames(dummy) <- c("freq_na", "freq_all")
        for (col in seq_len(ncol(genotypes))) {
            Z <- genotypes[, col]
            NAs <- sum(is.na(Z))
            dummy[col, 1] <- NAs / length(Z)
            dummy[col, 2] <- sum(Z, na.rm = TRUE) / ((length(Z) - NAs) * 2)
        }

        for (bufferSize in c(3, 6)) {
            for (nTasks in c(1, 3)) {
                expect_equal(summarize(genotypes, bufferSize = bufferSize, nTasks = nTasks, mc.cores = nCores), dummy)
                expect_equal(summarize(genotypes, bufferSize = bufferSize, nTasks = nTasks, mc.cores = nCores), dummy)
                expect_equal(summarize(genotypes, bufferSize = bufferSize, nTasks = nTasks, mc.cores = nCores), dummy)
                expect_equal(summarize(genotypes, mc.cores = nCores), dummy)
            }
        }

    })

}


test_that("normalizeType", {

    expect_equal(typeof(normalizeType("double")), "double")
    expect_equal(typeof(normalizeType(double())), "double")
    expect_equal(typeof(normalizeType("integer")), "integer")
    expect_equal(typeof(normalizeType(integer())), "integer")
    expect_equal(typeof(normalizeType("character")), "character")
    expect_equal(typeof(normalizeType(character())), "character")
    expect_equal(typeof(normalizeType("complex")), "complex")
    expect_equal(typeof(normalizeType(complex())), "complex")
    expect_warning(normalizeType("test"))
    expect_equal(suppressWarnings(typeof(normalizeType("test"))), "character")
    expect_equal(typeof(normalizeType(1)), "double")
    expect_equal(typeof(normalizeType(1L)), "integer")

})

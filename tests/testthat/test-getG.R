context("getG")

for (nCores in seq_len(2)) {

    test_that(paste("getGi", "on", nCores, "cores"), {

        hasCores(nCores)

        n <- 10
        p <- 100
        X <- matrix(data = rnorm(n * p), nrow = n, ncol = p)

        for (chunkSize in c(NULL, p, ceiling(p / 3))) {

            # both scalings
            G <- tcrossprod(scale(X))
            G <- G / mean(diag(G))
            G2 <- getG(X = X, scale = TRUE, scaleG = TRUE, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G, G2)

            # without scaling to average diagonal = 1 (scaleG)
            G <- tcrossprod(scale(X))
            G2 <- getG(X = X, scale = TRUE, scaleG = FALSE, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G, G2)

            # without scaling columns, but scaling average diagonal = 1 (scaleG)
            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G <- G / mean(diag(G))
            G2 <- getG(X = X, scale = FALSE, scaleG = TRUE, chunkSize = chunkSize, nCores = nCores)

            expect_equivalent(G, G2)

            # no scaling at all
            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G2 <- getG(X = X, scale = FALSE, scaleG = FALSE, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G, G2)

            # neither scaling nor centering
            G <- tcrossprod(X)
            G2 <- getG(X = X, center = FALSE, scale = FALSE, scaleG = FALSE, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G, G2)

        }

        X[sample(1:length(X), size = 20)] <- NA
        G <- getG(X, nCores = nCores)
        expect_true(!any(is.na(G)))

    })

    test_that(paste("getGij", "on", nCores, "cores"), {

        hasCores(nCores)

        n <- 10
        p <- 100
        X <- matrix(data = rnorm(n * p), nrow = n, ncol = p)

        for (chunkSize in c(NULL, p, ceiling(p / 3))) {

            i <- sample(1:nrow(X), size = 3)
            i2 <- sample(1:nrow(X), size = 4)

            centers <- colMeans(X)
            scales <- apply(X, 2, sd) * sqrt((n - 1)/n)

            # all scalings
            G <- tcrossprod(scale(X))
            G <- G / mean(diag(G))
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i2, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i2], G_12)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i], G_12)

            # without scaling to average diagonal = 1
            G <- tcrossprod(scale(X) * sqrt(n/(n - 1)))
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i2, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i2], G_12)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i], G_12)

            # without scaling columns, but scaling average diagonal = 1
            scales <- rep(1, ncol(X))

            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G <- G / ncol(X)
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i2, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i2], G_12)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = TRUE, i = i, i2 = i, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i], G_12)

            # no scaling at all
            G <- tcrossprod(scale(X, center = TRUE, scale = FALSE))
            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i2, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i2], G_12)

            G_12 <- getG(X = X, center = centers, scale = scales, scaleG = FALSE, i = i, i2 = i, chunkSize = chunkSize, nCores = nCores)
            expect_equivalent(G[i, i], G_12)

        }
    })

    test_that(paste("getG_symDMatrix", "on", nCores, "cores"), {

        hasCores(nCores)

        W <- matrix(data = rnorm(200), nrow = 10, ncol = 20)
        G1 <- tcrossprod(scale(W))
        G1 <- G1 / mean(diag(G1))

        G2 <- getG_symDMatrix(X = W, blockSize = ceiling(nrow(W) / 3), folderOut = testDir(), nCores = nCores)
        expect_equivalent(G2[], G1) # use equivalent to correct slight difference in NULL dimnames handling

    })

}

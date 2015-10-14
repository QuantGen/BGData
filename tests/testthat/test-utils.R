library(BGData)

context("utils")


test_that("crossprod.parallel", {
    
    W <- matrix(nrow = 10, ncol = 20, rnorm(200))
    Z <- matrix(nrow = 10, ncol = 2, rnorm(20))
    
    # TESTING X'X;
    TMP <- crossprod(W)
    
    TMP2 <- crossprod.parallel(W)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, nChunks = 3, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, nChunks = 1, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, nChunks = 3, mc.cores = 2)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, nChunks = 1, mc.cores = 1)
    expect_equal(TMP, TMP2)
    
    # Testing X'y
    TMP <- crossprod(W, y = Z)
    
    TMP2 <- crossprod.parallel(W, y = Z)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, y = Z, nChunks = 3, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, y = Z, nChunks = 1, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, y = Z, nChunks = 3, mc.cores = 2)
    expect_equal(TMP, TMP2)
    
    TMP2 <- crossprod.parallel(W, y = Z, nChunks = 1, mc.cores = 1)
    expect_equal(TMP, TMP2)
    
})


test_that("tcrossprod.parallel", {
    
    W <- matrix(nrow = 10, ncol = 20, rnorm(200))
    Z <- matrix(nrow = 10, ncol = 2, rnorm(20))
    
    # TESTING XX';
    TMP <- tcrossprod(W)
    
    TMP2 <- tcrossprod.parallel(W)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, nChunks = 3, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, nChunks = 1, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, nChunks = 3, mc.cores = 2)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, nChunks = 1, mc.cores = 1)
    expect_equal(TMP, TMP2)
    
    # Testing XY'
    
    W <- matrix(nrow = 10, ncol = 20, rnorm(200))
    Z <- matrix(nrow = 5, ncol = 20, rnorm(100))
    
    TMP <- tcrossprod(W, y = Z)
    
    TMP2 <- tcrossprod.parallel(W, y = Z)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, y = Z, nChunks = 3, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, y = Z, nChunks = 1, mc.cores = 3)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, y = Z, nChunks = 3, mc.cores = 2)
    expect_equal(TMP, TMP2)
    
    TMP2 <- tcrossprod.parallel(W, y = Z, nChunks = 1, mc.cores = 1)
    expect_equal(TMP, TMP2)
    
})


test_that("getG", {
    
    n <- 155
    p <- 1237
    X <- matrix(nrow = n, ncol = p, data = rnorm(n * p))
    
    nChunks <- c(1, 2, 3)
    mc.cores <- nChunks
    nChunks2 <- nChunks
    
    for (i in 1:3) {
        for (j in 1:3) {
            for (k in 1:3) {
                
                # all scalings
                G <- tcrossprod(scale(X))
                G <- G/mean(diag(G))
                G2 <- getG(x = X, scaleG = T, scaleCol = T, verbose = F, nChunks = i, mc.cores = j, nChunks2 = k)
                expect_equal(G, G2)
                
                # without scaling to average diagonal=1
                G <- tcrossprod(scale(X))
                G2 <- getG(x = X, scaleG = F, scaleCol = T, verbose = F, nChunks = i, mc.cores = j, nChunks2 = k)
                expect_equal(G, G2)
                
                # without scaling columns, but scaling average diagonal =1
                G <- tcrossprod(scale(X, center = T, scale = F))
                G <- G/mean(diag(G))
                G2 <- getG(x = X, scaleG = T, scaleCol = F, verbose = F, nChunks = i, mc.cores = j, nChunks2 = k)
                expect_equal(G, G2)
                
                # no scaling at all
                G <- tcrossprod(scale(X, center = T, scale = F))
                G2 <- getG(x = X, scaleG = F, scaleCol = F, verbose = F, nChunks = i, mc.cores = j, nChunks2 = k)
                expect_equal(G, G2)
                
            }
        }
    }
    
    X[sample(1:100, size = 20)] <- NA
    G <- getG(X, verbose = F)
    expect_true(!any(is.na(G)))
    
})


test_that("getGij", {
    
    n <- 155
    p <- 1237
    X <- matrix(nrow = n, ncol = p, data = rnorm(n * p))
    
    centers <- colMeans(X)
    scales <- apply(X = X, MARGIN = 2, FUN = sd)
    
    nChunks <- c(1, 2, 3)
    mc.cores <- nChunks
    nChunks2 <- nChunks
    
    for (i in 1:3) {
        for (j in 1:3) {
            for (k in 1:3) {
                
                # all scalings
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                G <- tcrossprod(scale(X))
                G <- G/mean(diag(G))
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = T, verbose = F, scaleCol = TRUE)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = T, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
                # without scaling to average diagonal=1
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                
                G <- tcrossprod(scale(X) * sqrt(n/(n - 1)))
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = F, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = F, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
                # without scaling columns, but scaling average diagonal =1
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                
                G <- tcrossprod(scale(X, center = T, scale = F))
                G <- G/ncol(X)
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = T, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = T, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
                # no scaling at all
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                
                G <- tcrossprod(scale(X, center = T, scale = F))
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = F, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = F, verbose = F)
                expect_equal(G[i1, i2], G_12)
                
            }
        }
    }
})


test_that("normalizeType", {
    
    expect_equal(typeof(normalizeType("double")), "double")
    expect_equal(typeof(normalizeType(double())), "double")
    expect_equal(typeof(normalizeType("integer")), "integer")
    expect_equal(typeof(normalizeType(integer())), "integer")
    expect_equal(typeof(normalizeType("character")), "character")
    expect_equal(typeof(normalizeType(character())), "character")
    expect_equal(typeof(normalizeType("complex")), "complex")
    expect_equal(typeof(normalizeType(complex())), "complex")
    expect_equal(typeof(normalizeType("test")), "character")
    expect_equal(typeof(normalizeType(1)), "double")
    expect_equal(typeof(normalizeType(1L)), "integer")
    
})

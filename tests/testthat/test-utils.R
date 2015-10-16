library(parallel)

context("utils")

hasCores <- function(numCores) {
    if (Sys.getenv("_R_CHECK_LIMIT_CORES_") == TRUE || numCores > detectCores()) {
        skip("Not enough cores or number of cores capped for CRAN submission checks.")
    }
}

for (nCores in seq_len(4)) {
    
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
    
}


for (nCores in seq_len(4)) {

    test_that(paste("getG", "on", nCores, "cores"), {
        
        hasCores(nCores)
        
        n <- 155
        p <- 1237
        X <- matrix(nrow = n, ncol = p, data = rnorm(n * p))
        
        for (nChunks in 1:3) {
            for (nChunks2 in 1:3) {
                
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
        
        X[sample(1:100, size = 20)] <- NA
        G <- getG(X, verbose = F)
        expect_true(!any(is.na(G)))
        
    })
    
    
    test_that(paste("getGij", "on", nCores, "cores"), {
        
        hasCores(nCores)
        
        n <- 155
        p <- 1237
        X <- matrix(nrow = n, ncol = p, data = rnorm(n * p))
        
        centers <- colMeans(X)
        scales <- apply(X = X, MARGIN = 2, FUN = sd)
        
        for (nChunks in 1:3) {
            for (nChunks2 in 1:3) {
                
                # all scalings
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                G <- tcrossprod(scale(X))
                G <- G/mean(diag(G))
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = T, verbose = F, scaleCol = TRUE, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                # without scaling to average diagonal=1
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                
                G <- tcrossprod(scale(X) * sqrt(n/(n - 1)))
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = apply(FUN = sd, X = X, MARGIN = 2) * sqrt((n - 1)/n), scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                # without scaling columns, but scaling average diagonal =1
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                
                G <- tcrossprod(scale(X, center = T, scale = F))
                G <- G/ncol(X)
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = T, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                # no scaling at all
                i1 <- sample(1:nrow(X), size = 3)
                i2 <- sample(1:nrow(X), size = 4)
                
                G <- tcrossprod(scale(X, center = T, scale = F))
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
                i2 <- i1
                G_12 <- getGij(x = X, i1 = i1, i2 = i2, centers = colMeans(X), scales = rep(1, ncol(X)), scaleG = F, verbose = F, nChunks = nChunks, nChunks2 = nChunks2, mc.cores = nCores)
                expect_equal(G[i1, i2], G_12)
                
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
    expect_equal(typeof(normalizeType("test")), "character")
    expect_equal(typeof(normalizeType(1)), "double")
    expect_equal(typeof(normalizeType(1L)), "integer")
    
})

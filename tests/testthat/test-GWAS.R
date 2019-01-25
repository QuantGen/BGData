context("GWAS")

set.seed(1)

nRows <- 15
nCols <- 50
percentNA <- 0.1

lm_test <- function(X, y, covariates = NULL) {
    res <- apply(X, 2, function(x) {
        data <- data.frame(
            y = y,
            x = x
        )
        if (!is.null(covariates)) {
            data <- cbind(data, covariates)
        }
        fm <- lm(y ~ ., data = data)
        coefficients(summary(fm))[2, ]
    })
    res <- t(res)
    rownames(res) <- colnames(X)
    return(res)
}

lsfit_test <- function(X, y, covariates = NULL) {
    res <- apply(X, 2, function(x) {
        fm <- lsfit(x = cbind(x, covariates), y = y)
        ls.print(fm, print.it = FALSE)$coef.table[[1]][2, ]
    })
    res <- t(res)
    rownames(res) <- colnames(X)
    return(res)
}

test_that("GWAS without covariates", {

    for (mode in c("integer", "double")) {

        X <- matrix(data = rnorm(nRows * nCols, sd = 100), nrow = nRows, ncol = nCols)
        X[sample(seq_along(X), size = ceiling(length(X) * percentNA))] <- NA
        storage.mode(X) <- mode

        y <- rnorm(nRows, sd = 100)
        y[sample(seq_along(y), size = ceiling(length(y) * percentNA))] <- NA

        lsfit_res <- suppressWarnings(lsfit_test(X, y))
        lm_res <- lm_test(X, y)

        DATA <- BGData(geno = X, pheno = data.frame(
            y = y
        ))

        for (method in c("rayOLS", "lsfit", "lm")) {

            for (nCores in seq_len(2)) {

                hasCores(nCores)

                GWAS_res <- suppressWarnings(GWAS(formula = y ~ 1, data = DATA, method = method, nCores = nCores))

                expect_equivalent(GWAS_res, lsfit_res)
                expect_equivalent(GWAS_res, lm_res)

            }

        }

    }

})

test_that("GWAS with covariates", {

    for (mode in c("integer", "double")) {

        X <- matrix(data = rnorm(nRows * nCols, sd = 100), nrow = nRows, ncol = nCols)

        PCs <- svd(X, nu = 2, nv = 0)$u
        colnames(PCs) <- c("pc1", "pc2")
        PCs[sample(seq_along(PCs), size = ceiling(length(PCs) * percentNA))] <- NA

        X[sample(seq_along(X), size = ceiling(length(X) * percentNA))] <- NA
        storage.mode(X) <- mode

        y <- rnorm(nRows, sd = 100)
        y[sample(seq_along(y), size = ceiling(length(y) * percentNA))] <- NA

        lsfit_res <- suppressWarnings(lsfit_test(X, y, PCs))
        lm_res <- lm_test(X, y, PCs)

        DATA <- BGData(geno = X, pheno = data.frame(
            y = y,
            pc1 = PCs[, 1],
            pc2 = PCs[, 2]
        ))

        for (method in c("lsfit", "lm")) {

            for (nCores in seq_len(2)) {

                hasCores(nCores)

                GWAS_res <- suppressWarnings(GWAS(formula = y ~ pc1 + pc2, data = DATA, method = method, nCores = nCores))

                expect_equivalent(GWAS_res, lsfit_res)
                expect_equivalent(GWAS_res, lm_res)

            }

        }

    }

})

context("GWAS")

set.seed(1)

nRows <- 5
nCols <- 10
nNAs <- 2

lsfit_R <- function(X, y) {
    res <- apply(X, 2, function(x) {
        fm <- lsfit(x = cbind(x, 1), y = y, intercept = FALSE)
        ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
    })
    res <- t(res)
    rownames(res) <- colnames(X)
    return(res)
}

test_that("lsfit", {

    for (mode in c("integer", "double")) {

        X <- matrix(data = rnorm(nRows * nCols, sd = 100), nrow = nRows, ncol = nCols)
        X[sample(seq_along(X), size = nNAs)] <- NA
        storage.mode(X) <- mode

        y <- data.frame(y = rnorm(nRows, sd = 100))
        y$y[sample(seq_along(y$y), size = nNAs)] <- NA

        DATA <- BGData(geno = X, pheno = y)

        expect_equal(
            GWAS(formula = y ~ 1, data = DATA, method = "lsfit"),
            suppressWarnings(lsfit_R(DATA@geno, DATA@pheno))
        )

    }

})

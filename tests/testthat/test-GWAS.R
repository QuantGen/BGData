context("GWAS")

set.seed(1)

nRows <- 5
nCols <- 10
percentNA <- 0.15

lsfit_R <- function(X, y) {
    res <- apply(X, 2, function(x) {
        fm <- lsfit(x = cbind(x, 1), y = y, intercept = FALSE)
        ls.print(fm, print.it = FALSE)$coef.table[[1]][1, ]
    })
    res <- t(res)
    rownames(res) <- colnames(X)
    return(res)
}

test_that("GWAS", {

    for (mode in c("integer", "double")) {

        X <- matrix(data = rnorm(nRows * nCols, sd = 100), nrow = nRows, ncol = nCols)
        X[sample(seq_along(X), size = as.integer(length(X) * percentNA))] <- NA
        storage.mode(X) <- mode

        y <- rnorm(nRows, sd = 100)
        y[sample(seq_along(y), size = as.integer(length(y) * percentNA))] <- NA

        DATA <- BGData(geno = X, pheno = data.frame(
            y = y
        ))

        for (method in c("rayOLS", "lsfit")) {

            expect_equal(
                GWAS(formula = y ~ 1, data = DATA, method = method),
                suppressWarnings(lsfit_R(DATA@geno, DATA@pheno))
            )

        }

    }

})

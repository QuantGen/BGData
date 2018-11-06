context("summarize")

for (nCores in seq_len(2)) {

    test_that(paste("summarize", "on", nCores, "cores"), {

        hasCores(nCores)

        genotypes <- matrix(data = c(0, 0, 1, 0, 2, 2, 1, 2, 0, 1, 2, 0, 0, 1, 2, 0, NA, 0), nrow = 3, ncol = 6)

        computeDummy <- function(i = seq_len(nrow(genotypes)), j = seq_len(ncol(genotypes))) {
            dummy <- data.frame(
                freq_na = vector(mode = "double", length = length(j)),
                allele_freq = vector(mode = "double", length = length(j)),
                sd = vector(mode = "double", length = length(j))
            )
            for (col in seq_along(j)) {
                Z <- genotypes[i, j[col]]
                NAs <- sum(is.na(Z))
                dummy$freq_na[col] <- NAs / length(Z)
                dummy$allele_freq[col] <- sum(Z, na.rm = TRUE) / ((length(Z) - NAs) * 2)
                dummy$sd[col] <- sd(Z, na.rm = TRUE)
            }
            return(dummy)
        }

        for (chunkSize in c(3, 6)) {
            expect_equal(summarize(X = genotypes, chunkSize = chunkSize, nCores = nCores), computeDummy())
            expect_equal(summarize(X = genotypes, i = c(1, 3), chunkSize = chunkSize, nCores = nCores), computeDummy(i = c(1, 3)))
            expect_equal(summarize(X = genotypes, j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), computeDummy(j =  c(1, 3, 5)))
            expect_equal(summarize(X = genotypes, i = c(1, 3), j = c(1, 3, 5), chunkSize = chunkSize, nCores = nCores), computeDummy(i = c(1, 3), j = c(1, 3, 5)))
            expect_equal(summarize(X = genotypes, nCores = nCores), computeDummy())
        }

    })

}

library(testthat)

context("dMatrix")

randomString <- function () {
    paste(sample(c(0:9, letters, LETTERS), size = 5, replace = TRUE), collapse = "")
}

# Create temporary directory
tmpPath <- paste0("/tmp/dMatrix-", randomString(), "/")
dir.create(tmpPath)

# Create example PED file
pedPath <- paste0(tmpPath, 'ped-', randomString(), '.txt')
nRows <- 3
nCols <- 3
genotypes <- simPED(filename = pedPath, n = nRows, p = nCols, returnGenos = TRUE)

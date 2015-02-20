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

test_that("setGenData complains if folderOut already exists", {
    dirExistsPath <- paste0(tmpPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(
        setGenData(fileIn = pedPath, dataType = 'integer', header = TRUE,
                   n = nRows, folderOut = dirExistsPath)
    )
})

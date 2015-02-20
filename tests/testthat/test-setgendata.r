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


context("setGenData")

test_that("setGenData complains if folderOut already exists", {
    dirExistsPath <- paste0(tmpPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(
        setGenData(fileIn = pedPath, dataType = 'integer', header = TRUE,
                   n = nRows, folderOut = dirExistsPath)
    )
})


#####
context("rDMatrix")

rGenData <- setGenData(fileIn = pedPath, dataType = 'integer',
                      header = TRUE, n = nRows,
                      folderOut = paste0(tmpPath, randomString()))

test_that("subsetting", {

    expect_equal(all.equal(rGenData@geno[], genotypes), TRUE)

    # expect_equal(all.equal(rGenData@geno[1], genotypes[1]), TRUE) Not implemented yet
    expect_equal(all.equal(rGenData@geno[1, ], genotypes[1, ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 1], genotypes[, 1]), TRUE)
    expect_equal(all.equal(rGenData@geno[1, 1], genotypes[1, 1]), TRUE)
    expect_equal(all.equal(rGenData@geno[1, , drop = FALSE], genotypes[1, , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 1, drop = FALSE], genotypes[, 1, drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[1, 1, drop = FALSE], genotypes[1, 1, drop = FALSE]), TRUE)

    # expect_equal(all.equal(rGenData@geno[1:2], genotypes[1:2]), TRUE)  Not implemented yet
    expect_equal(all.equal(rGenData@geno[1:2, ], genotypes[1:2, ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 1:2], genotypes[, 1:2]), TRUE)
    expect_equal(all.equal(rGenData@geno[1:2, 1:2], genotypes[1:2, 1:2]), TRUE)
    expect_equal(all.equal(rGenData@geno[1:2, , drop = FALSE], genotypes[1:2, , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 1:2, drop = FALSE], genotypes[, 1:2, drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[1:2, 1:2, drop = FALSE], genotypes[1:2, 1:2, drop = FALSE]), TRUE)

    # expect_equal(all.equal(rGenData@geno[2:1], genotypes[2:1]), TRUE) Not implemented yet
    expect_equal(all.equal(rGenData@geno[2:1, ], genotypes[2:1, ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 2:1], genotypes[, 2:1]), TRUE)
    expect_equal(all.equal(rGenData@geno[2:1, 2:1], genotypes[2:1, 2:1]), TRUE)
    expect_equal(all.equal(rGenData@geno[2:1, , drop = FALSE], genotypes[2:1, , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 2:1, drop = FALSE], genotypes[, 2:1, drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[2:1, 2:1, drop = FALSE], genotypes[2:1, 2:1, drop = FALSE]), TRUE)

    # expect_equal(all.equal(rGenData@geno[c(3, 1)], genotypes[c(3, 1)]), TRUE) Not implemented yet
    expect_equal(all.equal(rGenData@geno[c(3, 1), ], genotypes[c(3, 1), ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c(3, 1)], genotypes[, c(3, 1)]), TRUE)
    expect_equal(all.equal(rGenData@geno[c(3, 1), c(3, 1)], genotypes[c(3, 1), c(3, 1)]), TRUE)
    expect_equal(all.equal(rGenData@geno[c(3, 1), , drop = FALSE], genotypes[c(3, 1), , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c(3, 1), drop = FALSE], genotypes[, c(3, 1), drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[c(3, 1), c(3, 1), drop = FALSE], genotypes[c(3, 1), c(3, 1), drop = FALSE]), TRUE)

})


#####
context("cDMatrix")

cGenData <- setGenData(fileIn = pedPath, dataType = 'integer',
                       header = TRUE, n = nRows, distributed.by = "columns",
                       folderOut = paste0(tmpPath, randomString()))

test_that("subsetting", {

    expect_equal(all.equal(cGenData@geno[], genotypes), TRUE)

    # expect_equal(all.equal(cGenData@geno[1], genotypes[1]), TRUE) Not implemented yet
    expect_equal(all.equal(cGenData@geno[1, ], genotypes[1, ]), TRUE)
    expect_equal(all.equal(cGenData@geno[, 1], genotypes[, 1]), TRUE)
    expect_equal(all.equal(cGenData@geno[1, 1], genotypes[1, 1]), TRUE)
    expect_equal(all.equal(cGenData@geno[1, , drop = FALSE], genotypes[1, , drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[, 1, drop = FALSE], genotypes[, 1, drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[1, 1, drop = FALSE], genotypes[1, 1, drop = FALSE]), TRUE)

    # expect_equal(all.equal(cGenData@geno[1:2], genotypes[1:2]), TRUE)  Not implemented yet
    expect_equal(all.equal(cGenData@geno[1:2, ], genotypes[1:2, ]), TRUE)
    expect_equal(all.equal(cGenData@geno[, 1:2], genotypes[, 1:2]), TRUE)
    expect_equal(all.equal(cGenData@geno[1:2, 1:2], genotypes[1:2, 1:2]), TRUE)
    expect_equal(all.equal(cGenData@geno[1:2, , drop = FALSE], genotypes[1:2, , drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[, 1:2, drop = FALSE], genotypes[, 1:2, drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[1:2, 1:2, drop = FALSE], genotypes[1:2, 1:2, drop = FALSE]), TRUE)

    # expect_equal(all.equal(cGenData@geno[2:1], genotypes[2:1]), TRUE) Not implemented yet
    expect_equal(all.equal(cGenData@geno[2:1, ], genotypes[2:1, ]), TRUE)
    expect_equal(all.equal(cGenData@geno[, 2:1], genotypes[, 2:1]), TRUE)
    expect_equal(all.equal(cGenData@geno[2:1, 2:1], genotypes[2:1, 2:1]), TRUE)
    expect_equal(all.equal(cGenData@geno[2:1, , drop = FALSE], genotypes[2:1, , drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[, 2:1, drop = FALSE], genotypes[, 2:1, drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[2:1, 2:1, drop = FALSE], genotypes[2:1, 2:1, drop = FALSE]), TRUE)

    # expect_equal(all.equal(cGenData@geno[c(3, 1)], genotypes[c(3, 1)]), TRUE) Not implemented yet
    expect_equal(all.equal(cGenData@geno[c(3, 1), ], genotypes[c(3, 1), ]), TRUE)
    expect_equal(all.equal(cGenData@geno[, c(3, 1)], genotypes[, c(3, 1)]), TRUE)
    expect_equal(all.equal(cGenData@geno[c(3, 1), c(3, 1)], genotypes[c(3, 1), c(3, 1)]), TRUE)
    expect_equal(all.equal(cGenData@geno[c(3, 1), , drop = FALSE], genotypes[c(3, 1), , drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[, c(3, 1), drop = FALSE], genotypes[, c(3, 1), drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[c(3, 1), c(3, 1), drop = FALSE], genotypes[c(3, 1), c(3, 1), drop = FALSE]), TRUE)
})
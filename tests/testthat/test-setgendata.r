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

    # expect_equal(all.equal(rGenData@geno[genotypes > 1], genotypes[genotypes > 1]), TRUE) Not implemented yet
    expect_equal(all.equal(rGenData@geno[c(TRUE, FALSE, TRUE), ], genotypes[c(TRUE, FALSE, TRUE), ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c(TRUE, FALSE, TRUE)], genotypes[, c(TRUE, FALSE, TRUE)]), TRUE)
    expect_equal(all.equal(rGenData@geno[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], genotypes[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)]), TRUE)
    expect_equal(all.equal(rGenData@geno[c(TRUE, FALSE, TRUE), , drop = FALSE], genotypes[c(TRUE, FALSE, TRUE), , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c(TRUE, FALSE, TRUE), drop = FALSE], genotypes[, c(TRUE, FALSE, TRUE), drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], genotypes[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE]), TRUE)

})

test_that("replacement", {

    # Generate new genotypes for replacement
    replacmentPath <- paste0(tmpPath, 'ped-', randomString(), '.txt')
    replacement <- simPED(filename = replacmentPath, n = nRows, p = nCols, returnGenos = TRUE)
    comparison <- genotypes

    testAndRestore <- function () {
        expect_equal(all.equal(rGenData@geno[], comparison), TRUE)
        rGenData@geno[] <- genotypes # no environment change necessary
        assign('comparison', genotypes, parent.frame())
    }

    rGenData@geno[] <- replacement
    comparison[] <- replacement
    testAndRestore()

    rGenData@geno[1, ] <- replacement[1, ]
    comparison[1, ] <- replacement[1, ]
    testAndRestore()
    rGenData@geno[, 1] <- replacement[, 1]
    comparison[, 1] <- replacement[, 1]
    testAndRestore()
    rGenData@geno[1, 1] <- replacement[1, 1]
    comparison[1, 1] <- replacement[1, 1]
    testAndRestore()

    rGenData@geno[1:2, ] <- replacement[1:2, ]
    comparison[1:2, ] <- replacement[1:2, ]
    testAndRestore()
    rGenData@geno[, 1:2] <- replacement[, 1:2]
    comparison[, 1:2] <- replacement[, 1:2]
    testAndRestore()
    rGenData@geno[1:2, 1:2] <- replacement[1:2, 1:2]
    comparison[1:2, 1:2] <- replacement[1:2, 1:2]
    testAndRestore()

    rGenData@geno[2:1, ] <- replacement[2:1, ]
    comparison[2:1, ] <- replacement[2:1, ]
    testAndRestore()
    rGenData@geno[, 2:1] <- replacement[, 2:1]
    comparison[, 2:1] <- replacement[, 2:1]
    testAndRestore()
    rGenData@geno[2:1, 2:1] <- replacement[2:1, 2:1]
    comparison[2:1, 2:1] <- replacement[2:1, 2:1]
    testAndRestore()

    rGenData@geno[c(3, 1), ] <- replacement[c(3, 1), ]
    comparison[c(3, 1), ] <- replacement[c(3, 1), ]
    testAndRestore()
    rGenData@geno[, c(3, 1)] <- replacement[, c(3, 1)]
    comparison[, c(3, 1)] <- replacement[, c(3, 1)]
    testAndRestore()
    rGenData@geno[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
    comparison[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
    testAndRestore()

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

    # expect_equal(all.equal(cGenData@geno[genotypes > 1], genotypes[genotypes > 1]), TRUE) Not implemented yet
    expect_equal(all.equal(cGenData@geno[c(TRUE, FALSE, TRUE), ], genotypes[c(TRUE, FALSE, TRUE), ]), TRUE)
    expect_equal(all.equal(cGenData@geno[, c(TRUE, FALSE, TRUE)], genotypes[, c(TRUE, FALSE, TRUE)]), TRUE)
    expect_equal(all.equal(cGenData@geno[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)], genotypes[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE)]), TRUE)
    expect_equal(all.equal(cGenData@geno[c(TRUE, FALSE, TRUE), , drop = FALSE], genotypes[c(TRUE, FALSE, TRUE), , drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[, c(TRUE, FALSE, TRUE), drop = FALSE], genotypes[, c(TRUE, FALSE, TRUE), drop = FALSE]), TRUE)
    expect_equal(all.equal(cGenData@geno[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE], genotypes[c(TRUE, FALSE, TRUE), c(TRUE, FALSE, TRUE), drop = FALSE]), TRUE)

})

test_that("replacement", {

    # Generate new genotypes for replacement
    replacmentPath <- paste0(tmpPath, 'ped-', randomString(), '.txt')
    replacement <- simPED(filename = replacmentPath, n = nRows, p = nCols, returnGenos = TRUE)
    comparison <- genotypes

    testAndRestore <- function () {
        expect_equal(all.equal(cGenData@geno[], comparison), TRUE)
        cGenData@geno[] <- genotypes # no environment change necessary
        assign('comparison', genotypes, parent.frame())
    }

    cGenData@geno[] <- replacement
    comparison[] <- replacement
    testAndRestore()

    cGenData@geno[1, ] <- replacement[1, ]
    comparison[1, ] <- replacement[1, ]
    testAndRestore()
    cGenData@geno[, 1] <- replacement[, 1]
    comparison[, 1] <- replacement[, 1]
    testAndRestore()
    cGenData@geno[1, 1] <- replacement[1, 1]
    comparison[1, 1] <- replacement[1, 1]
    testAndRestore()

    cGenData@geno[1:2, ] <- replacement[1:2, ]
    comparison[1:2, ] <- replacement[1:2, ]
    testAndRestore()
    cGenData@geno[, 1:2] <- replacement[, 1:2]
    comparison[, 1:2] <- replacement[, 1:2]
    testAndRestore()
    cGenData@geno[1:2, 1:2] <- replacement[1:2, 1:2]
    comparison[1:2, 1:2] <- replacement[1:2, 1:2]
    testAndRestore()

    cGenData@geno[2:1, ] <- replacement[2:1, ]
    comparison[2:1, ] <- replacement[2:1, ]
    testAndRestore()
    cGenData@geno[, 2:1] <- replacement[, 2:1]
    comparison[, 2:1] <- replacement[, 2:1]
    testAndRestore()
    cGenData@geno[2:1, 2:1] <- replacement[2:1, 2:1]
    comparison[2:1, 2:1] <- replacement[2:1, 2:1]
    testAndRestore()

    cGenData@geno[c(3, 1), ] <- replacement[c(3, 1), ]
    comparison[c(3, 1), ] <- replacement[c(3, 1), ]
    testAndRestore()
    cGenData@geno[, c(3, 1)] <- replacement[, c(3, 1)]
    comparison[, c(3, 1)] <- replacement[, c(3, 1)]
    testAndRestore()
    cGenData@geno[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
    comparison[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
    testAndRestore()

})
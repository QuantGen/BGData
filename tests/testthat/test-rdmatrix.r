library(dMatrix)

context("rDMatrix")

# Create temporary directory
tmpPath <- paste0("/tmp/dMatrix-", randomString(), "/")
dir.create(tmpPath)

# Create example PED file
pedPath <- paste0(tmpPath, 'ped-', randomString(), '.txt')
nRows <- 3
nCols <- 3
genotypes <- simPED(filename = pedPath, n = nRows, p = nCols, returnGenos = TRUE)

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
    
    expect_equal(all.equal(rGenData@geno['id_1', ], genotypes['id_1', ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 'mrk_1'], genotypes[, 'mrk_1']), TRUE)
    expect_equal(all.equal(rGenData@geno['id_1', 'mrk_1'], genotypes['id_1', 'mrk_1']), TRUE)
    expect_equal(all.equal(rGenData@geno['id_1', , drop = FALSE], genotypes['id_1', , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, 'mrk_1', drop = FALSE], genotypes[, 'mrk_1', drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno['id_1', 'mrk_1', drop = FALSE], genotypes['id_1', 'mrk_1', drop = FALSE]), TRUE)
    
    expect_equal(all.equal(rGenData@geno[c('id_1', 'id_2'), ], genotypes[c('id_1', 'id_2'), ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c('mrk_1', 'mrk_2')], genotypes[, c('mrk_1', 'mrk_2')]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_1', 'id_2'), c('mrk_1', 'mrk_2')], genotypes[c('id_1', 'id_2'), c('mrk_1', 'mrk_2')]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_1', 'id_2'), , drop = FALSE], genotypes[c('id_1', 'id_2'), , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c('mrk_1', 'mrk_2'), drop = FALSE], genotypes[, c('mrk_1', 'mrk_2'), drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_1', 'id_2'), c('mrk_1', 'mrk_2'), drop = FALSE], genotypes[c('id_1', 'id_2'), c('mrk_1', 'mrk_2'), drop = FALSE]), TRUE)
    
    expect_equal(all.equal(rGenData@geno[c('id_2', 'id_1'), ], genotypes[c('id_2', 'id_1'), ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c('mrk_2', 'mrk_1')], genotypes[, c('mrk_2', 'mrk_1')]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_2', 'id_1'), c('mrk_2', 'mrk_1')], genotypes[c('id_2', 'id_1'), c('mrk_2', 'mrk_1')]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_2', 'id_1'), , drop = FALSE], genotypes[c('id_2', 'id_1'), , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c('mrk_2', 'mrk_1'), drop = FALSE], genotypes[, c('mrk_2', 'mrk_1'), drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_2', 'id_1'), c('mrk_2', 'mrk_1'), drop = FALSE], genotypes[c('id_2', 'id_1'), c('mrk_2', 'mrk_1'), drop = FALSE]), TRUE)
    
    expect_equal(all.equal(rGenData@geno[c('id_3', 'id_1'), ], genotypes[c('id_3', 'id_1'), ]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c('mrk_3', 'mrk_1')], genotypes[, c('mrk_3', 'mrk_1')]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_3', 'id_1'), c('mrk_3', 'mrk_1')], genotypes[c('id_3', 'id_1'), c('mrk_3', 'mrk_1')]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_3', 'id_1'), , drop = FALSE], genotypes[c('id_3', 'id_1'), , drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[, c('mrk_3', 'mrk_1'), drop = FALSE], genotypes[, c('mrk_3', 'mrk_1'), drop = FALSE]), TRUE)
    expect_equal(all.equal(rGenData@geno[c('id_3', 'id_1'), c('mrk_3', 'mrk_1'), drop = FALSE], genotypes[c('id_3', 'id_1'), c('mrk_3', 'mrk_1'), drop = FALSE]), TRUE)
    
})

test_that("replacement", {
    
    # Generate new genotypes for replacement
    replacmentPath <- paste0(tmpPath, 'ped-', randomString(), '.txt')
    replacement <- simPED(filename = replacmentPath, n = nRows, p = nCols, returnGenos = TRUE)
    comparison <- genotypes
    
    testAndRestore <- function (label) {
        expect_equal(all.equal(rGenData@geno[], comparison), TRUE, label=label)
        rGenData@geno[] <- genotypes # no environment change necessary
        assign('comparison', genotypes, parent.frame())
    }
    
    rGenData@geno[] <- replacement
    comparison[] <- replacement
    testAndRestore('[]')
    
    rGenData@geno[1, ] <- replacement[1, ]
    comparison[1, ] <- replacement[1, ]
    testAndRestore('[1, ]')
    rGenData@geno[, 1] <- replacement[, 1]
    comparison[, 1] <- replacement[, 1]
    testAndRestore('[, 1]')
    rGenData@geno[1, 1] <- replacement[1, 1]
    comparison[1, 1] <- replacement[1, 1]
    testAndRestore('[1, 1]')
    
    rGenData@geno[1:2, ] <- replacement[1:2, ]
    comparison[1:2, ] <- replacement[1:2, ]
    testAndRestore('[1:2, ]')
    rGenData@geno[, 1:2] <- replacement[, 1:2]
    comparison[, 1:2] <- replacement[, 1:2]
    testAndRestore('[, 1:2]')
    rGenData@geno[1:2, 1:2] <- replacement[1:2, 1:2]
    comparison[1:2, 1:2] <- replacement[1:2, 1:2]
    testAndRestore('[1:2, 1:2]')
    
    rGenData@geno[2:1, ] <- replacement[2:1, ]
    comparison[2:1, ] <- replacement[2:1, ]
    testAndRestore('[2:1, ]')
    rGenData@geno[, 2:1] <- replacement[, 2:1]
    comparison[, 2:1] <- replacement[, 2:1]
    testAndRestore('[, 2:1]')
    rGenData@geno[2:1, 2:1] <- replacement[2:1, 2:1]
    comparison[2:1, 2:1] <- replacement[2:1, 2:1]
    testAndRestore('[2:1, 2:1]')
    
    rGenData@geno[c(3, 1), ] <- replacement[c(3, 1), ]
    comparison[c(3, 1), ] <- replacement[c(3, 1), ]
    testAndRestore('[c(3, 1), ]')
    rGenData@geno[, c(3, 1)] <- replacement[, c(3, 1)]
    comparison[, c(3, 1)] <- replacement[, c(3, 1)]
    testAndRestore('[, c(3, 1)]')
    rGenData@geno[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
    comparison[c(3, 1), c(3, 1)] <- replacement[c(3, 1), c(3, 1)]
    testAndRestore('[c(3, 1), c(3, 1)]')
    
})

test_that("apply", {
    
    expect_equal(all.equal(colMeans(rGenData@geno), colMeans(genotypes)), TRUE)
    expect_equal(all.equal(colSums(rGenData@geno), colSums(genotypes)), TRUE)
    expect_equal(all.equal(rowMeans(rGenData@geno), rowMeans(genotypes)), TRUE)
    expect_equal(all.equal(rowSums(rGenData@geno), rowSums(genotypes)), TRUE)
    
    # Introduce NA
    genotypes_na <- genotypes
    genotypes_na[1,1] <- NA
    rGenData@geno[1,1] <- NA
    
    expect_warning(colMeans(rGenData@geno))
    expect_warning(colSums(rGenData@geno))
    expect_warning(rowMeans(rGenData@geno))
    expect_warning(rowSums(rGenData@geno))
    
    expect_equal(all.equal(colMeans(rGenData@geno, na.rm=FALSE), colMeans(genotypes_na, na.rm=FALSE)), TRUE)
    expect_equal(all.equal(colSums(rGenData@geno, na.rm=FALSE), colSums(genotypes_na, na.rm=FALSE)), TRUE)
    expect_equal(all.equal(rowMeans(rGenData@geno, na.rm=FALSE), rowMeans(genotypes_na, na.rm=FALSE)), TRUE)
    expect_equal(all.equal(rowSums(rGenData@geno, na.rm=FALSE), rowSums(genotypes_na, na.rm=FALSE)), TRUE)
    
    expect_equal(all.equal(colMeans(rGenData@geno, na.rm=TRUE), colMeans(genotypes_na, na.rm=TRUE)), TRUE)
    expect_equal(all.equal(colSums(rGenData@geno, na.rm=TRUE), colSums(genotypes_na, na.rm=TRUE)), TRUE)
    expect_equal(all.equal(rowMeans(rGenData@geno, na.rm=TRUE), rowMeans(genotypes_na, na.rm=TRUE)), TRUE)
    expect_equal(all.equal(rowSums(rGenData@geno, na.rm=TRUE), rowSums(genotypes_na, na.rm=TRUE)), TRUE)
    
    # Revert NA
    rGenData@geno[] <- genotypes
    
})

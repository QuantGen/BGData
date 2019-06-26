context("BGData")

# Create dummy path
testPath <- paste0(tempdir(), "/BGData-", BGData:::randomString(), "/")
dir.create(testPath)

restoreGenotypes <- function() {
    genotypes <- matrix(c(4, 4, 4, 3, 2, 3, 1, 2, 1), nrow = nRows, ncol = nCols)
    colnames(genotypes) <- paste0("mrk_", 1:3)
    rownames(genotypes) <- paste0("1_", 1:3)
    return(genotypes)
}

# Create example .raw files
pedPath <- paste0(testPath, "ped-", BGData:::randomString(), ".txt")
nRows <- 3
nCols <- 3
phenotypes <- data.frame(FID = c("1", "1", "1"), IID = c("1", "2", "3"), 
    PAT = c("NA", "NA", "NA"), MAT = c("NA", "NA", "NA"), SEX = c("NA", "NA", "NA"), 
    PHENOTYPE = c("NA", "NA", "NA"), stringsAsFactors = FALSE)
phenotypes[] <- lapply(phenotypes, type.convert, as.is = TRUE)
rownames(phenotypes) <- paste0("1_", 1:3)
genotypes <- restoreGenotypes()
ped <- cbind(phenotypes, genotypes)
outFile <- file(pedPath, "w")
write.table(ped, file = outFile, quote = FALSE, row.names = FALSE)
close(outFile)


context("initialize")

test_that("it requires at least geno", {
    expect_error(BGData())
})

test_that("it checks if pheno is a data.frame", {
    expect_error(BGData(geno = genotypes, pheno = rownames(genotypes)))
})

test_that("it checks if map is a data.frame", {
    expect_error(BGData(geno = genotypes, map = colnames(genotypes)))
})

test_that("it checks if the number of rows of geno match with the number of rows of pheno", {
    expect_error(BGData(geno = genotypes, pheno = phenotypes[-1, ]))
})

test_that("it checks if the number of rows of geno match with the number of rows of pheno", {
    map <- data.frame(mrk = colnames(genotypes))
    expect_error(BGData(geno = genotypes, map = map[-1, ]))
})

test_that("it warns if the row names of pheno do not match the row names of geno", {
    expect_warning(BGData(geno = genotypes, pheno = phenotypes[nrow(phenotypes):1, ]))
})

test_that("it warns if the row names of map do not match the columns names of geno", {
    map <- data.frame(mrk = rev(colnames(genotypes)))
    expect_warning(BGData(geno = genotypes, map = map))
})


context("readRAW")

test_that("it complains if folderOut already exists", {
    dirExistsPath <- paste0(testPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(readRAW(fileIn = pedPath, n = nRows, folderOut = dirExistsPath))
})


test_that("it reads .raw files into BGData objects", {

    # With minimum number of parameters (with exception of folderOut)
    BGData <- readRAW(fileIn = pedPath, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equivalent(geno(BGData)[], genotypes)

    # With n
    BGData <- readRAW(fileIn = pedPath, n = nRows, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equivalent(geno(BGData)[], genotypes)

    # With p
    BGData <- readRAW(fileIn = pedPath, p = nCols, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equivalent(geno(BGData)[], genotypes)

    # With both n and p
    BGData <- readRAW(fileIn = pedPath, n = nRows, p = nCols, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equivalent(geno(BGData)[], genotypes)

    # As integer
    class(genotypes) <- "integer"
    BGData <- readRAW(fileIn = pedPath, dataType = integer(), folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equivalent(geno(BGData)[], genotypes)
    BGData <- readRAW(fileIn = pedPath, dataType = "integer", folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equivalent(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

    # As double
    class(genotypes) <- "double"
    BGData <- readRAW(fileIn = pedPath, dataType = double(), folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equivalent(geno(BGData)[], genotypes)
    BGData <- readRAW(fileIn = pedPath, dataType = "double", folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equivalent(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

    # As character
    expect_error(readRAW(fileIn = pedPath, dataType = character(), folderOut = paste0(testPath, "test-", BGData:::randomString())))
    expect_error(readRAW(fileIn = pedPath, dataType = "character", folderOut = paste0(testPath, "test-", BGData:::randomString())))

})


context("readRAW_matrix")

test_that("it reads a .raw file into a matrix object", {

    # With minimum number of parameters (with exception of folderOut)
    BGData <- readRAW_matrix(fileIn = pedPath)
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # With n
    BGData <- readRAW_matrix(fileIn = pedPath, n = nRows)
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # With p
    BGData <- readRAW_matrix(fileIn = pedPath, p = nCols)
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # With both n and p
    BGData <- readRAW_matrix(fileIn = pedPath, n = nRows, p = nCols)
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # As integer
    class(genotypes) <- "integer"
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = integer())
    expect_equal(geno(BGData)[], genotypes)
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = "integer")
    expect_equal(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

    # As double
    class(genotypes) <- "double"
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = double())
    expect_equal(geno(BGData)[], genotypes)
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = "double")
    expect_equal(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

    # As character
    class(genotypes) <- "character"
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = character())
    expect_equal(geno(BGData)[], genotypes)
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = "character")
    expect_equal(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

})

context("readRAW_big.matrix")

test_that("it reads a .raw file into a big.matrix object", {

    # With minimum number of parameters (with exception of folderOut)
    BGData <- readRAW_big.matrix(fileIn = pedPath, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # With n
    BGData <- readRAW_big.matrix(fileIn = pedPath, n = nRows, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # With p
    BGData <- readRAW_big.matrix(fileIn = pedPath, p = nCols, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # With both n and p
    BGData <- readRAW_big.matrix(fileIn = pedPath, n = nRows, p = nCols, folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(pheno(BGData), phenotypes)
    expect_equal(geno(BGData)[], genotypes)

    # As integer
    class(genotypes) <- "integer"
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = integer(), folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(geno(BGData)[], genotypes)
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = "integer", folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

    # As double
    class(genotypes) <- "double"
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = double(), folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(geno(BGData)[], genotypes)
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = "double", folderOut = paste0(testPath, "test-", BGData:::randomString()))
    expect_equal(geno(BGData)[], genotypes)
    genotypes <- restoreGenotypes()

    # As character
    expect_error(readRAW(fileIn = pedPath, dataType = character(), folderOut = paste0(testPath, "test-", BGData:::randomString())))
    expect_error(readRAW(fileIn = pedPath, dataType = "character", folderOut = paste0(testPath, "test-", BGData:::randomString())))

})

context("load.BGData")

test_that("it loads BGData objects created by readRAW", {

    # Create dummy BGData object without returning data
    path <- paste0(testPath, "test-", BGData:::randomString())
    readRAW(fileIn = pedPath, folderOut = path)
    expect_true(!("BGData" %in% ls()))

    # Append BGData.RData to path
    path <- paste0(path, "/", "BGData.RData")

    # Load BGData object and test if all nodes have been opened
    load.BGData(path)
    expect_true("BGData" %in% ls())
    for (node in seq_len(LinkedMatrix::nNodes(geno(BGData)))) {
        expect_true(ff::is.open(geno(BGData)[[node]]))
    }
    expect_equal(dim(geno(BGData)), c(nRows, nCols))

})

test_that("it loads BGData objects created by readRAW_matrix", {

    # Create dummy BGData object
    path <- paste0(testPath, "test-", BGData:::randomString(), "/", "BGData.RData")
    dir.create(dirname(path))
    BGData <- readRAW_matrix(fileIn = pedPath)
    save(BGData, file = path)
    rm(BGData)
    expect_true(!("BGData" %in% ls()))

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(geno(BGData)), c(nRows, nCols))

})

test_that("it loads BGData objects created by readRAW_big.matrix", {

    # Create dummy BGData object
    path <- paste0(testPath, "test-", BGData:::randomString())
    readRAW_big.matrix(fileIn = pedPath, dataType = integer(), folderOut = path)
    expect_true(!("BGData" %in% ls()))

    # Append BGData.RData to path
    path <- paste0(path, "/", "BGData.RData")

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(geno(BGData)), c(nRows, nCols))

})

test_that("it loads BGData objects containing a BEDMatrix object", {

    # Create dummy objects
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"))
    bedDims <- dim(bedMatrix)
    bedDNames <- dimnames(bedMatrix)
    bedRow <- bedMatrix[1, ]
    BGData <- BGData(geno = bedMatrix)

    # Save BGData object
    path <- paste0(testPath, "test-", BGData:::randomString(), "/", "BGData.RData")
    dir.create(dirname(path))
    save(BGData, file = path)
    rm(BGData)
    expect_true(!("BGData" %in% ls()))

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(geno(BGData)), bedDims)
    expect_equal(dimnames(geno(BGData)), bedDNames)
    expect_equal(geno(BGData)[1, ], bedRow)

})

context("as.BGData")

test_that("it converts a regular BEDMatrix object to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"))
    bgData <- as.BGData(bedMatrix)
    expect_is(bgData, "BGData")
    expect_equal(dim(geno(bgData)), dim(bedMatrix))
    expect_equal(nrow(pheno(bgData)), nrow(bedMatrix))
    expect_equal(rownames(pheno(bgData)), rownames(bedMatrix))
    expect_equal(nrow(map(bgData)), ncol(bedMatrix))
    expect_equal(rownames(map(bgData)), colnames(bedMatrix))
})

test_that("it converts a BEDMatrix object created with the n parameter to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"), n = 199)
    bgData <- as.BGData(bedMatrix)
    expect_is(bgData, "BGData")
    expect_equal(dim(geno(bgData)), dim(bedMatrix))
    expect_equal(nrow(pheno(bgData)), nrow(bedMatrix))
    expect_equal(nrow(map(bgData)), ncol(bedMatrix))
    expect_equal(rownames(map(bgData)), colnames(bedMatrix))
})

test_that("it converts a BEDMatrix object created with the p parameter to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"), p = 300)
    bgData <- as.BGData(bedMatrix)
    expect_is(bgData, "BGData")
    expect_equal(dim(geno(bgData)), dim(bedMatrix))
    expect_equal(nrow(pheno(bgData)), nrow(bedMatrix))
    expect_equal(rownames(pheno(bgData)), rownames(bedMatrix))
    expect_equal(nrow(map(bgData)), ncol(bedMatrix))
})

test_that("it converts a BEDMatrix object created with the n and p parameters to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"), n = 199, p = 300)
    bgData <- as.BGData(bedMatrix)
    expect_is(bgData, "BGData")
    expect_equal(dim(geno(bgData)), dim(bedMatrix))
    expect_equal(nrow(pheno(bgData)), nrow(bedMatrix))
    expect_equal(nrow(map(bgData)), ncol(bedMatrix))
})

test_that("it throws an error if an alternate phenotype file does not exist when converting a BEDMatrix object to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"))
    expect_error(as.BGData(bedMatrix, alternatePhenotypeFile = "NOT_FOUND"))
})

test_that("it reads an alternate phenotype file when converting a BEDMatrix object to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"))
    bgData <- as.BGData(bedMatrix, alternatePhenotypeFile = system.file("extdata", "pheno.txt", package = "BGData"))
    expect_is(bgData, "BGData")
    # Test if pheno has an extra column for the phenotype
    expect_equal(ncol(pheno(bgData)), 7)
    # Test merging and NA handling
    expect_equal(pheno(bgData)[1, 7], 57.0)
    expect_equal(nrow(pheno(bgData)), nrow(geno(bgData)))
    expect_true(all(is.na(pheno(bgData)[c(178, 180, 189, 190, 196), 7])))
    # Test if rownames are retained
    expect_equal(rownames(pheno(bgData)), rownames(bedMatrix))
})

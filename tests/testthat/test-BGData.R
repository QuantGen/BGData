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


context("readRAW")

test_that("it complains if folderOut already exists", {
    dirExistsPath <- paste0(testPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(readRAW(fileIn = pedPath, n = nRows, folderOut = dirExistsPath))
})


test_that("it reads .raw files into BGData objects", {

    # With minimum number of parameters (with exception of folderOut)
    BGData <- readRAW(fileIn = pedPath, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)

    # With n
    BGData <- readRAW(fileIn = pedPath, n = nRows, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)

    # With p
    BGData <- readRAW(fileIn = pedPath, p = nCols, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)

    # With both n and p
    BGData <- readRAW(fileIn = pedPath, n = nRows, p = nCols, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)

    # As integer
    class(genotypes) <- "integer"
    BGData <- readRAW(fileIn = pedPath, dataType = integer(), folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)
    BGData <- readRAW(fileIn = pedPath, dataType = "integer", folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)
    genotypes <- restoreGenotypes()

    # As double
    class(genotypes) <- "double"
    BGData <- readRAW(fileIn = pedPath, dataType = double(), folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)
    BGData <- readRAW(fileIn = pedPath, dataType = "double", folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes, check.attributes = FALSE)
    genotypes <- restoreGenotypes()

    # As character
    expect_error(readRAW(fileIn = pedPath, dataType = character(), folderOut = paste0(testPath, "test-", randomString())))
    expect_error(readRAW(fileIn = pedPath, dataType = "character", folderOut = paste0(testPath, "test-", randomString())))

})


context("readRAW_matrix")

test_that("it reads a .raw file into a matrix object", {

    # With minimum number of parameters (with exception of folderOut)
    BGData <- readRAW_matrix(fileIn = pedPath)
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # With n
    BGData <- readRAW_matrix(fileIn = pedPath, n = nRows)
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # With p
    BGData <- readRAW_matrix(fileIn = pedPath, p = nCols)
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # With both n and p
    BGData <- readRAW_matrix(fileIn = pedPath, n = nRows, p = nCols)
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # As integer
    class(genotypes) <- "integer"
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = integer())
    expect_equal(BGData@geno[], genotypes)
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = "integer")
    expect_equal(BGData@geno[], genotypes)
    genotypes <- restoreGenotypes()

    # As double
    class(genotypes) <- "double"
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = double())
    expect_equal(BGData@geno[], genotypes)
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = "double")
    expect_equal(BGData@geno[], genotypes)
    genotypes <- restoreGenotypes()

    # As character
    class(genotypes) <- "character"
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = character())
    expect_equal(BGData@geno[], genotypes)
    BGData <- readRAW_matrix(fileIn = pedPath, dataType = "character")
    expect_equal(BGData@geno[], genotypes)
    genotypes <- restoreGenotypes()

})

context("readRAW_big.matrix")

test_that("it reads a .raw file into a big.matrix object", {

    # With minimum number of parameters (with exception of folderOut)
    BGData <- readRAW_big.matrix(fileIn = pedPath, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # With n
    BGData <- readRAW_big.matrix(fileIn = pedPath, n = nRows, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # With p
    BGData <- readRAW_big.matrix(fileIn = pedPath, p = nCols, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # With both n and p
    BGData <- readRAW_big.matrix(fileIn = pedPath, n = nRows, p = nCols, folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@pheno, phenotypes)
    expect_equal(BGData@geno[], genotypes)

    # As integer
    class(genotypes) <- "integer"
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = integer(), folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes)
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = "integer", folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes)
    genotypes <- restoreGenotypes()

    # As double
    class(genotypes) <- "double"
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = double(), folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes)
    BGData <- readRAW_big.matrix(fileIn = pedPath, dataType = "double", folderOut = paste0(testPath, "test-", randomString()))
    expect_equal(BGData@geno[], genotypes)
    genotypes <- restoreGenotypes()

    # As character
    expect_error(readRAW(fileIn = pedPath, dataType = character(), folderOut = paste0(testPath, "test-", randomString())))
    expect_error(readRAW(fileIn = pedPath, dataType = "character", folderOut = paste0(testPath, "test-", randomString())))

})

context("load.BGData")

test_that("it loads BGData objects created by readRAW", {

    # Create dummy BGData object without returning data
    path <- paste0(testPath, "test-", randomString())
    readRAW(fileIn = pedPath, folderOut = path)
    expect_true(!("BGData" %in% ls()))

    # Append BGData.RData to path
    path <- paste0(path, "/", "BGData.RData")

    # Load BGData object and test if all nodes have been opened
    load.BGData(path)
    expect_true("BGData" %in% ls())
    for (node in seq_len(LinkedMatrix::nNodes(BGData@geno))) {
        expect_true(ff::is.open(BGData@geno[[node]]))
    }
    expect_equal(dim(BGData@geno), c(nRows, nCols))

})

test_that("it loads BGData objects created by readRAW_matrix", {

    # Create dummy BGData object
    path <- paste0(testPath, "test-", randomString(), "/", "BGData.RData")
    dir.create(dirname(path))
    BGData <- readRAW_matrix(fileIn = pedPath)
    save(BGData, file = path)
    rm(BGData)
    expect_true(!("BGData" %in% ls()))

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(BGData@geno), c(nRows, nCols))

})

test_that("it loads BGData objects created by readRAW_big.matrix", {

    # Create dummy BGData object
    path <- paste0(testPath, "test-", randomString())
    readRAW_big.matrix(fileIn = pedPath, dataType = integer(), folderOut = path)
    expect_true(!("BGData" %in% ls()))

    # Append BGData.RData to path
    path <- paste0(path, "/", "BGData.RData")

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(BGData@geno), c(nRows, nCols))

})

test_that("it loads BGData objects containing a BEDMatrix object", {

    # Create dummy objects
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"))
    bedDims <- dim(bedMatrix)
    bedDNames <- dimnames(bedMatrix)
    bedRow <- bedMatrix[1, ]
    BGData <- BGData(geno = bedMatrix)

    # Save BGData object
    path <- paste0(testPath, "test-", randomString(), "/", "BGData.RData")
    dir.create(dirname(path))
    save(BGData, file = path)
    rm(BGData)
    expect_true(!("BGData" %in% ls()))

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(BGData@geno), bedDims)
    expect_equal(dimnames(BGData@geno), bedDNames)
    expect_equal(BGData@geno[1, ], bedRow)

})

context("as.BGData")

test_that("it converts a BEDMatrix object to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "chr1.bed", package = "BGData"))
    bgData <- as.BGData(bedMatrix)
    expect_is(bgData, "BGData")
    expect_equal(dim(bgData@geno), dim(bedMatrix))
    expect_equal(nrow(bgData@pheno), nrow(bedMatrix))
    expect_equal(rownames(bgData@pheno), rownames(bedMatrix))
    expect_equal(paste0(bgData@pheno$FID, "_", bgData@pheno$IID), rownames(bedMatrix))
    expect_equal(nrow(bgData@map), ncol(bedMatrix))
    expect_equal(rownames(bgData@map), colnames(bedMatrix))
    expect_equal(paste0(bgData@map$snp_id, "_", bgData@map$allele_1), colnames(bedMatrix))
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
    expect_equal(ncol(bgData@pheno), 7)
    # Test merging and NA handling
    expect_equal(bgData@pheno[1, 7], 57.0)
    expect_equal(nrow(bgData@pheno), nrow(bgData@geno))
    expect_true(all(is.na(bgData@pheno[c(178, 180, 189, 190, 196), 7])))
})

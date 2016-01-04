context("BGData")

# Create temporary directory
tmpPath <- paste0("/tmp/BGData-", randomString(), "/")
dir.create(tmpPath)

restoreGenotypes <- function() {
    genotypes <- matrix(c(4, 4, 4, 3, 2, 3, 1, 2, 1), nrow = nRows, ncol = nCols)
    colnames(genotypes) <- paste0("mrk_", 1:3)
    rownames(genotypes) <- paste0("1_", 1:3)
    return(genotypes)
}

# Create example PED files
pedPath <- paste0(tmpPath, "ped-", randomString(), ".txt")
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
outGzFile <- gzfile(paste0(pedPath, ".gz"), "w")
write.table(ped, file = outFile, quote = FALSE, row.names = FALSE)
write.table(ped, file = outGzFile, quote = FALSE, row.names = FALSE)
close(outFile)
close(outGzFile)


context("readPED")

test_that("it complains if folderOut already exists", {
    dirExistsPath <- paste0(tmpPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(readPED(fileIn = pedPath, header = TRUE, dataType = "integer", n = nRows, folderOut = dirExistsPath))
})


test_that("it reads PED files into BGData objects", {

    for (ext in c("", ".gz")) {

        adjPedPath <- paste0(pedPath, ext)

        # With minimum number of parameters (with exception of folderOut)
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = integer(), folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With n
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = integer(), n = nRows, folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With p
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = integer(), p = nCols, folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With both n and p
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = integer(), n = nRows, p = nCols, folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # As integer
        class(genotypes) <- "integer"
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = integer(), folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = "integer", folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

        # As double
        class(genotypes) <- "double"
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = double(), folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED(fileIn = adjPedPath, header = TRUE, dataType = "double", folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

        # As character
        expect_error(readPED(fileIn = adjPedPath, header = TRUE, dataType = character(), folderOut = paste0(tmpPath, "test-", randomString())))
        expect_error(readPED(fileIn = adjPedPath, header = TRUE, dataType = "character", folderOut = paste0(tmpPath, "test-", randomString())))
    }

})


context("readPED.matrix")

test_that("it reads a PED file into a matrix object", {

    for (ext in c("", ".gz")) {

        adjPedPath <- paste0(pedPath, ext)

        # With minimum number of parameters (with exception of folderOut)
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer())
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With n
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), n = nRows)
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With p
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), p = nCols)
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With both n and p
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), n = nRows, p = nCols)
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # As integer
        class(genotypes) <- "integer"
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer())
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = "integer")
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

        # As double
        class(genotypes) <- "double"
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = double())
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = "double")
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

        # As character
        class(genotypes) <- "character"
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = character())
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED.matrix(fileIn = adjPedPath, header = TRUE, dataType = "character")
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

    }

})

context("readPED.big.matrix")

test_that("it reads a PED file into a big.matrix object", {

    for (ext in c("", ".gz")) {

        adjPedPath <- paste0(pedPath, ext)

        # With minimum number of parameters (with exception of folderOut)
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With n
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), n = nRows, folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With p
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), p = nCols, folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # With both n and p
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), n = nRows, p = nCols, folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@pheno, phenotypes)
        expect_equal(BGData@geno[], genotypes)

        # As integer
        class(genotypes) <- "integer"
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = integer(), folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = "integer", folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

        # As double
        class(genotypes) <- "double"
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = double(), folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        BGData <- readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = "double", folderOut = paste0(tmpPath, "test-", randomString()))
        expect_equal(BGData@geno[], genotypes)
        genotypes <- restoreGenotypes()

        # As character
        class(genotypes) <- "character"
        expect_error(readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = character(), folderOut = paste0(tmpPath, "test-", randomString())))
        expect_error(readPED.big.matrix(fileIn = adjPedPath, header = TRUE, dataType = "character", folderOut = paste0(tmpPath, "test-", randomString())))
        genotypes <- restoreGenotypes()

    }

})

context("load.BGData")

test_that("it loads BGData objects created by readPED", {

    # Create dummy BGData object without returning data
    path <- paste0(tmpPath, "test-", randomString())
    readPED(fileIn = pedPath, header = TRUE, dataType = integer(), folderOut = path)
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

test_that("it loads BGData objects created by readPED.matrix", {

    # Create dummy BGData object
    path <- paste0(tmpPath, "test-", randomString(), "/", "BGData.RData")
    dir.create(dirname(path))
    BGData <- readPED.matrix(fileIn = pedPath, header = TRUE, dataType = integer())
    save(BGData, file = path)
    rm(BGData)
    expect_true(!("BGData" %in% ls()))

    # Load BGData object
    load.BGData(path)
    expect_true("BGData" %in% ls())
    expect_equal(dim(BGData@geno), c(nRows, nCols))

})

test_that("it loads BGData objects created by readPED.big.matrix", {

    # Create dummy BGData object
    path <- paste0(tmpPath, "test-", randomString())
    readPED.big.matrix(fileIn = pedPath, header = TRUE, dataType = integer(), folderOut = path)
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
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))
    bedDims <- dim(bedMatrix)
    bedDNames <- dimnames(bedMatrix)
    bedRow <- bedMatrix[1, ]
    BGData <- BGData(geno = bedMatrix)

    # Save BGData object
    path <- paste0(tmpPath, "test-", randomString(), "/", "BGData.RData")
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
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))
    bgData <- as.BGData(bedMatrix)
    expect_is(bgData, "BGData")
    expect_equal(dim(bgData@geno), dim(bedMatrix))
    expect_equal(nrow(bgData@pheno), nrow(bedMatrix))
    expect_equal(nrow(bgData@map), ncol(bedMatrix))
})

test_that("it throws an error if an alternate phenotype file does not exist when converting a BEDMatrix object to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))
    expect_error(as.BGData(bedMatrix, alternatePhenotypeFile = "NOT_FOUND"))
})

test_that("it reads an alternate phenotype file when converting a BEDMatrix object to a BGData object", {
    bedMatrix <- BEDMatrix::BEDMatrix(system.file("extdata", "example.bed", package = "BEDMatrix"))
    bgData <- as.BGData(bedMatrix, alternatePhenotypeFile = system.file("extdata", "pheno.txt", package = "BGData"))
    expect_is(bgData, "BGData")
    # Test if pheno has two extra columns
    expect_equal(ncol(bgData@pheno), 8)
    # Test merging and NA handling
    expect_equal(bgData@pheno[3, 8], 24.12)
    expect_equal(nrow(bgData@pheno), nrow(bgData@geno))
    expect_true(all(is.na(bgData@pheno[2, c(7, 8)])))
})

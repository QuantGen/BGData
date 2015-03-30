library(BGData)

context("BGData")

# Create temporary directory
tmpPath <- paste0("/tmp/BGData-", randomString(), "/")
dir.create(tmpPath)

# Create example PED file
pedPath <- paste0(tmpPath, 'ped-', randomString(), '.txt')
nRows <- 3
nCols <- 3
phenotypes <- data.frame(FID = c('0', '0', '0'),
                         IID = c('id_1', 'id_2', 'id_3'),
                         PAT = c('NA', 'NA', 'NA'),
                         MAT = c('NA', 'NA', 'NA'),
                         SEX = c('NA', 'NA', 'NA'),
                         PHENOTYPE = c('NA', 'NA', 'NA'),
                         stringsAsFactors = FALSE)
phenotypes[] <- lapply(phenotypes, type.convert, as.is=TRUE)
rownames(phenotypes) <- paste0('id_', 1:3)
genotypes <- matrix(c(4, 4, 4, 3, 2, 3, 1, 2, 1), nrow = nRows, ncol = nCols)
colnames(genotypes) <- paste0('mrk_', 1:3)
rownames(genotypes) <- paste0('id_', 1:3)
ped <- cbind(phenotypes, genotypes)
write.table(ped, file = pedPath, quote = FALSE, row.names = FALSE)


context("read.PED.BGData.mmMatrix")

test_that("it complains if folderOut already exists", {
    dirExistsPath <- paste0(tmpPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(
        read.PED.BGData.mmMatrix(fileIn = pedPath, header = TRUE, dataType = 'integer',
                                 n = nRows, folderOut = dirExistsPath)
    )
})

test_that("it reads a PED file into a BGData object", {
    
    # With minimum number of parameters (with exception of folderOut)
    BGData <- read.PED.BGData.mmMatrix(fileIn = pedPath, header = TRUE, dataType = integer(),
                                       folderOut = paste0(tmpPath, 'BGData-', randomString()))
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
    # With n
    BGData <- read.PED.BGData.mmMatrix(fileIn = pedPath, header = TRUE,
                                       dataType = integer(), n = nRows,
                                       folderOut = paste0(tmpPath, 'BGData-', randomString()))
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
    # With p
    BGData <- read.PED.BGData.mmMatrix(fileIn = pedPath, header = TRUE,
                                       dataType = integer(), p = nCols,
                                       folderOut = paste0(tmpPath, 'BGData-', randomString()))
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
    # With both n and p
    BGData <- read.PED.BGData.mmMatrix(fileIn = pedPath, header = TRUE,
                                       dataType = integer(), n = nRows, p = nCols,
                                       folderOut = paste0(tmpPath, 'BGData-', randomString()))
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
})


context("read.PED.BGData.matrix")

test_that("it reads a PED file into a matrix object", {
    
    # With minimum number of parameters (with exception of folderOut)
    BGData <- read.PED.BGData.matrix(fileIn = pedPath, header = TRUE,
                                     dataType = integer())
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
    # With n
    BGData <- read.PED.BGData.matrix(fileIn = pedPath, header = TRUE,
                                     dataType = integer(), n = nRows)
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
    # With p
    BGData <- read.PED.BGData.matrix(fileIn = pedPath, header = TRUE,
                                     dataType = integer(), p = nCols)
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
    # With both n and p
    BGData <- read.PED.BGData.matrix(fileIn = pedPath, header = TRUE,
                                     dataType = integer(), n = nRows, p = nCols)
    expect_equal(all.equal(BGData@pheno, phenotypes), TRUE)
    expect_equal(all.equal(BGData@geno[], genotypes), TRUE)
    
})

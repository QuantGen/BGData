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

context("read.PED.BGData")

test_that("read.PED.BGData complains if folderOut already exists", {
    dirExistsPath <- paste0(tmpPath, "dirExists")
    dir.create(dirExistsPath, showWarnings = FALSE)
    expect_error(
        read.PED.BGData(fileIn = pedPath, dataType = 'integer', header = TRUE,
                        n = nRows, folderOut = dirExistsPath)
    )
})

library(testthat)
library(BGData)

# Create temporary directory
testPath <- paste0("/tmp/BGData-", BGData:::randomString(), "/")
dir.create(testPath)
options("testPath") <- testPath

test_check("BGData")

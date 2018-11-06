library(parallel)

testDir <- function() {
    paste0(tempdir(), "/BGData-", BGData:::randomString(), "/")
}

hasCores <- function(numCores) {
    # For CRAN
    if (Sys.getenv("_R_CHECK_LIMIT_CORES_") == TRUE || numCores > parallel::detectCores()) {
        skip("Not enough cores or number of cores capped for CRAN submission checks.")
    }
    # For WinBuilder
    if (.Platform$OS.type == "windows" && numCores > 1) {
        skip("mc.cores > 1 is not supported on Windows.")
    }
}

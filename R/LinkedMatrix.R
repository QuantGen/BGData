#' @export
ffRowLinkedMatrix <- function(n, p, vmode, nNodes = NULL, folderOut = paste0("BGData_", randomString()), dimorder = 2:1) {
    ffLinkedMatrix(class = "RowLinkedMatrix", n = n, p = p, vmode = vmode, nNodes = nNodes, folderOut = folderOut, dimorder = dimorder)
}

#' @export
ffColumnLinkedMatrix <- function(n, p, vmode, nNodes = NULL, folderOut = paste0("BGData_", randomString()), dimorder = 1:2) {
    ffLinkedMatrix(class = "ColumnLinkedMatrix", n = n, p = p, vmode = vmode, nNodes = nNodes, folderOut = folderOut, dimorder = dimorder)
}

ffLinkedMatrix <- function(class, n, p, vmode, nNodes = NULL, folderOut = paste0("BGData_", randomString()), dimorder = if (class == "RowLinkedMatrix") 2:1 else 1:2) {

    # Create output directory
    if (file.exists(folderOut)) {
        stop(paste("Output folder", folderOut, "already exists. Please move it or pick a different one."))
    }
    dir.create(folderOut)

    # Determine chunk size and number of nodes
    if (is.null(nNodes)) {
        if (class == "RowLinkedMatrix") {
            chunkSize <- min(n, floor(.Machine$integer.max / p / 1.2))
            nNodes <- ceiling(n / chunkSize)
        } else {
            chunkSize <- min(p, floor(.Machine$integer.max / n / 1.2))
            nNodes <- ceiling(p / chunkSize)
        }
    } else {
        if (class == "RowLinkedMatrix") {
            chunkSize <- ceiling(n / nNodes)
            if (chunkSize * p >= .Machine$integer.max / 1.2) {
              stop("More nodes are needed")
            }
        } else {
            chunkSize <- ceiling(p / nNodes)
            if (chunkSize * n >= .Machine$integer.max / 1.2) {
              stop("More nodes are needed")
            }
        }
    }

    # Initialize list
    geno <- new(class)
    end <- 0
    if (class == "RowLinkedMatrix") {
        for (i in seq_len(nNodes)) {
            ini <- end + 1
            end <- min(n, ini + chunkSize - 1)
            filename <- paste0("geno_", i, ".bin")
            geno[[i]] <- ff(vmode = vmode, dim = c((end - ini + 1), p), dimorder = dimorder, filename = paste0(folderOut, .Platform$file.sep, filename))
            # Change ff path to a relative one
            physical(geno[[i]])$pattern <- "ff"
            physical(geno[[i]])$filename <- filename
        }
    } else {
        for (i in seq_len(nNodes)) {
            ini <- end + 1
            end <- min(p, ini + chunkSize - 1)
            filename <- paste0("geno_", i, ".bin")
            geno[[i]] <- ff(vmode = vmode, dim = c(n, (end - ini + 1)), dimorder = dimorder, filename = paste0(folderOut, .Platform$file.sep, filename))
            # Change ff path to a relative one
            physical(geno[[i]])$pattern <- "ff"
            physical(geno[[i]])$filename <- filename
        }
    }

    return(geno)
}

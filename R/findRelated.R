findRelated <- function(x, ...) {
    UseMethod("findRelated")
}

findRelated.matrix <- function(x, cutoff = 0.03, ...) {
    x[lower.tri(x, diag = TRUE)] <- 0
    pairs <- which(x > cutoff, arr.ind = TRUE, useNames = FALSE)
    samples <- unique(pairs[, 1L])
    rownames(x)[samples]
}

findRelated.symDMatrix <- function(x, cutoff = 0.03, verbose = FALSE, ...) {
    n <- nBlocks(x)
    pairs <- lapply(seq_len(n), function(i) {
        lapply(seq(i, n), function(j) {
            if (verbose) {
                message("Working on block ", i, " ", j)
            }
            block <- x[[i]][[j]][]
            # Remove lower triangle in blocks that contain the diagonal
            if (i == j) {
                block[lower.tri(block, diag = TRUE)] <- 0
            }
            pairs <- which(block > cutoff, arr.ind = TRUE, useNames = FALSE)
            # Remap local indices to sample names
            remap <- matrix(character(), nrow = nrow(pairs), ncol = ncol(pairs))
            remap[, 1L] <- rownames(block)[pairs[, 1L]]
            remap[, 2L] <- colnames(block)[pairs[, 2L]]
            return(remap)
        })
    })
    pairs <- do.call(rbind, lapply(pairs, function(x) do.call(rbind, x)))
    unique(pairs[, 1L])
}

summarize <- function(X, i = seq_len(nrow(X)), j = seq_len(ncol(X)), chunkSize = 5000L, nCores = getOption("mc.cores", 2L), verbose = FALSE) {
    res <- chunkedMap(X = X, FUN = function(chunk) {
        summaries <- .Call(C_summarize, chunk)
        rownames(summaries) <- colnames(chunk)
        colnames(summaries) <- c("freq_na", "allele_freq", "sd")
        return(summaries)
    }, i = i, j = j, chunkSize = chunkSize, nCores = nCores, verbose = verbose)
    res <- do.call(rbind, res)
    as.data.frame(res)
}

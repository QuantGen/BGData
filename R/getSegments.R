getSegments <- function(x, chr, bp, names, threshold, lag, trim = FALSE, verbose = FALSE) {
    uniqueChr <- unique(chr)
    out <- vector(mode = "list", length = length(uniqueChr))
    for (curChr in uniqueChr) {
        if (verbose) {
            message("Working on chromosome ", curChr)
        }
        # Extract chromosome data
        chrFilter <- which(chr == curChr)
        xChr <- x[chrFilter]
        bpChr <- bp[chrFilter]
        namesChr <- names[chrFilter]
        # Determine variants below threshold
        discoverySet <- which(xChr <= threshold[1])
        # Set discoveries and all variants within +/- lag to 1, leave rest as 0
        signal <- rep(0, length(chrFilter))
        for (discovery in discoverySet) {
            signal[abs(bpChr - bpChr[discovery]) <= lag] <- 1
        }
        # Determine the runs in the 0/1 signal
        runs <- rle(signal)
        # Determine at what positions within the chromosome the runs start and
        # end while removing 0-runs
        runStart <- c(1, cumsum(runs$length[-length(runs$length)]) + 1)
        withinSegment <- runs$values == 1
        runStart <- runStart[withinSegment]
        runEnd <- runStart + runs$length[withinSegment] - 1
        runLength <- runs$length[withinSegment]
        # Determine name and value of smallest variant within segment, and
        # optionally trim segment (i.e., remove variants that are not internal
        # to the segment containing GWAS-significant variants)
        # Would be nice to vectorize this like the other operations ...
        minValue <- vector(mode = "numeric", length = length(runStart))
        minName <- vector(mode = "numeric", length = length(runStart))
        for (curSeg in seq_along(runStart)) {
            segFilter <- seq(runStart[curSeg], runEnd[curSeg])
            xSeg <- xChr[segFilter]
            namesSeg <- namesChr[segFilter]
            minValuePos <- which.min(xSeg)
            minValue[curSeg] <- xSeg[minValuePos]
            minName[curSeg] <- namesSeg[minValuePos]
            if (trim) {
                # Determine which variants in the segment passed the threshold
                # (the second discovery set, so to speak)
                significantVariants <- which(xSeg <= threshold[2])
                # Set start of run to first significant variant and end of run
                # to last significant variant
                runStart[curSeg] <- segFilter[significantVariants[1]]
                runEnd[curSeg] <- segFilter[significantVariants[length(significantVariants)]]
                runLength[curSeg] <- runEnd[curSeg] - runStart[curSeg] + 1
            }
        }
        # Determine at what base-pair positions the runs start and end
        bpStart <- bpChr[runStart]
        bpEnd <- bpChr[runEnd]
        bpLength <- bpEnd - bpStart + 1
        # Determine at what positions within x the runs start and end (more
        # useful information than chromosome by chromosome because it is easier
        # to extract)
        xStart <- chrFilter[runStart]
        xEnd <- chrFilter[runEnd]
        # Prepare chromosome summary (there might be no segments, so we cannot
        # rely on recycling)
        outChr <- data.frame(
            chr = rep(curChr, times = length(runStart)),
            start = xStart,
            end = xEnd,
            length = runLength,
            bpStart = bpStart,
            bpEnd = bpEnd,
            bpLength = bpLength,
            minName = minName,
            minValue = minValue
        )
        out[[curChr]] <- outChr
    }
    # Combine chromosomes
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    return(out)
}

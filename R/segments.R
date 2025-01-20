segments <- function(statistic, chr, bp, threshold, gap, trim = FALSE, verbose = FALSE,snpid=NULL) {
    if (length(unique(c(length(statistic), length(chr), length(bp)))) != 1) {
        stop("statistic, chr, and bp need to match in length")
    }
    if (!is.numeric(statistic)) {
        stop("'statistic' needs to be a numeric vector")
    }
    if (!(is.numeric(chr) || is.character(chr))) {
        stop("'chr' needs to be a either a character or numeric vector")
    }
    if (!is.numeric(bp)) {
        stop("'bp' needs to be a numeric vector")
    }
    if (!is.numeric(threshold)) {
        stop("'threshold' needs to a number")
    }
    if (!is.numeric(gap)) {
        stop("'gap' needs to a number")
    }
    uniqueChr <- unique(chr)
    out <- vector(mode = "list", length = length(uniqueChr))
    for (curChr in uniqueChr) {
        if (verbose) {
            message("Working on chromosome ", curChr)
        }
        # Extract chromosome data
        chrFilter <- which(chr == curChr)
        statisticChr <- statistic[chrFilter]
        bpChr <- bp[chrFilter]
        # Determine variants below threshold
        discoverySet <- which(statisticChr <= threshold)
        # Set discoveries and all variants within +/- gap to 1, leave rest as 0
        signal <- rep(0, length(chrFilter))
        for (discovery in discoverySet) {
            signal[abs(bpChr - bpChr[discovery]) <= gap] <- 1
        }
        # Determine the runs in the 0/1 signal
        runs <- rle(signal)
        # Determine at what positions within the chromosome the runs start and
        # end while removing 0-runs
        runStart <- c(1, cumsum(runs[["lengths"]][-length(runs[["lengths"]])]) + 1)
        withinSegment <- runs[["values"]] == 1
        runStart <- runStart[withinSegment]
        runEnd <- runStart + runs[["lengths"]][withinSegment] - 1
        runLength <- runs[["lengths"]][withinSegment]
        # Determine value and position of smallest variant within segment, and
        # optionally trim segment (i.e., remove variants that are not internal
        # to the segment containing GWAS-significant variants)
        # Would be nice to vectorize this like the other operations ...
        minValue <- vector(mode = "numeric", length = length(runStart))
        minValuePos <- vector(mode = "integer", length = length(runStart))
        for (curSeg in seq_along(runStart)) {
            segFilter <- seq(runStart[curSeg], runEnd[curSeg])
            statisticSeq <- statisticChr[segFilter]
            minValuePosSeg <- which.min(statisticSeq)
            minValue[curSeg] <- statisticSeq[minValuePosSeg]
            minValuePos[curSeg] <- chrFilter[1] + segFilter[1] + minValuePosSeg - 2
            if (trim) {
                # Determine which variants in the segment passed the threshold
                significantVariants <- which(statisticSeq <= threshold)
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
        # Prepare chromosome summary (there might be no segments, so do not
        # rely on recycling)
        outChr <- data.frame(
            chr = rep(curChr, times = length(runStart)),
            start = xStart,
            end = xEnd,
            length = runLength,
            bpStart = bpStart,
            bpEnd = bpEnd,
            bpLength = bpLength,
            minValue = minValue,
            minValuePos = minValuePos,
            minValueBp=bp[minValuePos]
        )
        if(!is.null(snpid)){
            minValueID=snpid[minValuePos]
        }
        out[[curChr]] <- outChr
    }
    # Combine chromosomes
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    return(out)
}

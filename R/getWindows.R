#' Computes sliding LD-based windows.
#'
#' For each SNP present in a genotpe matrix, obtains nearby SNPs in LD based on genotype correlations.
#'
#' @param X A matrix-like object, typically `@@geno` of a [BGData-class]
#' object.
#' @param center_snp Columns of `X` to compute windows. By default, does not attempt to obtain LD windows surrounding
#' each element of `X`.
#' @param rSq R-squared threshold for determining genotype correlations.
#' @param maxGaps The number of consequtive SNPs that may fall below `rSq` until LD is considered lost
#' @param ids Indices of individuals used to compute LD, if not all individuals in X are desired.
#' @return A list the length of the number of columns of `X` containing sliding LD-based windows. The elements of
#' the list are named for each "center" SNP used in LD determination. Each element is a vector of SNP names--SNPs
#' in LD with the "center" SNP.
#' @export
getWindows <- function(X, center_snp = 1:100, rSq = 0.1, maxGaps = 2, ids = 1:nrow(X)) {
    win <- lapply(center_snp, function(i) getBlock(X, i, rSq = rSq, maxGaps = maxGaps, ids = ids))
    names(win) <- colnames(X)[center_snp]
    return(win)
}


# Hidden function used by getWindows().
getBlock=function(X, center, rSq, maxGaps, ids = 1:nrow(X)){
    p=ncol(X)
    # sanity check
    if(center<=0 | center>p){
        stop("center must be between 1 and ncol(X)")
    }
    block=list() # this list will contain the position of the SNPs in the blcok

    block$left=list()
    block$center=center
    block$right=list()

    xi=X[ids,center]

    ###########################################################
    ## while loop to the right
    ###########################################################
    ready= (center>=p) # no loop to the right if center=p
    lag=1
    nGaps=0
    while(!ready){
        R2=stats::cor(xi,X[ids,center+lag],use='complete.obs')^2
        if(R2>rSq){
            block$right[[lag]]=center+lag
            nGaps=0
        }else{
            nGaps=nGaps+1
            ready=nGaps>maxGaps
            if(!ready){
                block$right[[lag]]=center+lag
            }
        }
        if(!ready){lag=lag+1}
        if(center+lag>p){ ready=TRUE }
    }

    ###########################################################
    ## while loop to the left
    ###########################################################
    ready<-(center<=1) #no loop to the left if center=1
    lag=1
    nGaps=0
    while(!ready){
        R2=stats::cor(xi,X[ids,center-lag],use='complete.obs')^2
        if(R2>rSq){
            block$left[[lag]]=center-lag
            nGaps=0
        }else{
            nGaps=nGaps+1
            ready=nGaps>maxGaps
            if(!ready){
                block$left[[lag]]=center-lag
            }
        }
        if(!ready){lag=lag+1}
        if(center-lag<=0){ ready=TRUE }
    }

    # Trim excess SNPs on either end that did not meet rSq threshold.
    # If length(block$right) or length(block$left) - maxGaps == 0, then remove
    # block$right or block$left, respectively. If both need removed, the center
    # SNP is determined to be a singleton.
    if (length(block$right) > maxGaps) {
        block$right <- block$right[1:(length(block$right)-maxGaps)]
    } else {
        block$right <- NULL
    }

    if (length(block$left) > maxGaps) {
        block$left <- block$left[1:(length(block$left)-maxGaps)]
    } else {
        block$left <- NULL
    }

    return(sort(unlist(block)))
}

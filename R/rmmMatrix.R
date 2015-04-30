#' An S4 class to represent a row-distributed \code{mmMatrix}
#'
#' \code{rmmMatrix} inherits from \code{\link{list}}. Each element of the list is
#' an \code{ff_matrix} object.
#'
#' @export rmmMatrix
#' @exportClass rmmMatrix
rmmMatrix<-setClass('rmmMatrix',contains='list')

#' @export
setMethod('initialize','rmmMatrix',function(.Object,nrow=1,ncol=1,vmode='byte',folderOut=NULL,nChunks=NULL,dimorder=c(2,1)){
    if(is.null(folderOut)){
        folderOut<-paste0(tempdir(),'/mmMatrix-',randomString())
    }
    if(file.exists(folderOut)){
        stop(paste('Output folder',folderOut,'already exists. Please move it or pick a different one.'))
    }
    dir.create(folderOut)
    if(is.null(nChunks)){
        chunkSize<-min(nrow,floor(.Machine$integer.max/ncol/1.2))
        nChunks<-ceiling(nrow/chunkSize)
    }else{
        chunkSize<-ceiling(nrow/nChunks)
        if(chunkSize*ncol >= .Machine$integer.max/1.2){
            stop('More chunks are needed')
        }
    }
    ffList<-list()
    end<-0
    for(i in 1:nChunks){
        ini<-end+1
        end<-min(nrow,ini+chunkSize-1)
        filename=paste0('geno_',i,'.bin')
        ffList[[i]]<-ff(vmode=vmode,dim=c((end-ini+1),ncol),dimorder=dimorder,filename=paste0(folderOut,'/',filename))
        # Change ff path to a relative one
        physical(ffList[[i]])$pattern<-'ff'
        physical(ffList[[i]])$filename<-filename
    }
    .Object<-callNextMethod(.Object,ffList)
    return(.Object)
})


subset.rmmMatrix<-function(x,i,j,drop){
    if(class(i)=='logical'){
        i<-which(i)
    }else if(class(i)=='character'){
        i<-sapply(i,function(name){
            which(rownames(x)==name)
        },USE.NAMES=FALSE)
    }
    if(class(j)=='logical'){
        j<-which(j)
    }else if(class(j)=='character'){
        j<-sapply(j,function(name){
            which(colnames(x)==name)
        },USE.NAMES=FALSE)
    }
    n<-length(i)
    p<-length(j)
    originalOrder<-(1:n)[order(i)]
    sortedRows<-sort(i)

    dimX<-dim(x)
    if( p>dimX[2] | n>dimX[1] ){
        stop('Either the number of columns or number of rows requested exceed the number of rows or columns in x, try dim(x)...')
    }

    Z<-matrix(nrow=n,ncol=p,NA)
    colnames(Z)<-colnames(x)[j]
    rownames(Z)<-rownames(x)[i]

    INDEXES<-rowindexes(x,rows=sortedRows)

    whatChunks<-unique(INDEXES[,1])
    end<-0
    for(k in whatChunks){
        TMP<-matrix(data=INDEXES[INDEXES[,1]==k,],ncol=3)
        ini<-end+1
        end<-ini+nrow(TMP)-1
        Z[ini:end,]<-x[[k]][TMP[,3],j,drop=FALSE]
    }
    if(length(originalOrder)>1){
        Z[]<-Z[originalOrder,]
    }
    if(drop==TRUE&&(n==1||p==1)){
        # Revert drop.
        return(Z[,])
    }else{
        return(Z)
    }
}

#' @export
setMethod("[",signature(x="rmmMatrix",i="ANY",j="ANY",drop="ANY"),subset.rmmMatrix)

#' @export
setMethod("[",signature(x="rmmMatrix",i="ANY",j="missing",drop="ANY"),function(x,i,drop){
    j<-1:ncol(x)
    subset.rmmMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="rmmMatrix",i="missing",j="ANY",drop="ANY"),function(x,j,drop) {
    i<-1:nrow(x)
    subset.rmmMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="rmmMatrix",i="missing",j="missing",drop="ANY"),function(x,drop) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    subset.rmmMatrix(x,i,j,drop)
})


replace.rmmMatrix<-function(x,i,j,...,value){
    Z<-matrix(nrow=length(i),ncol=length(j),data=value)
    CHUNKS<-chunks(x)
    for(k in 1:nrow(CHUNKS) ){
        rows_z<-(i>=CHUNKS[k,2])&(i<=CHUNKS[k,3])
        rowLocal<-rowindexes(x,i[rows_z])[,3]
        x[[k]][rowLocal,j]<-Z[rows_z,]
    }
    return(x)
}

#' @export
setReplaceMethod("[",signature(x="rmmMatrix",i="ANY",j="ANY",value="ANY"),replace.rmmMatrix)

#' @export
setReplaceMethod("[",signature(x="rmmMatrix",i="ANY",j="missing",value="ANY"),function(x,i,value){
    j<-1:ncol(x)
    replace.rmmMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="rmmMatrix",i="missing",j="ANY",value="ANY"),function(x,j,value) {
    i<-1:nrow(x)
    replace.rmmMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="rmmMatrix",i="missing",j="missing",value="ANY"),function(x,value) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    replace.rmmMatrix(x,i,j,value=value)
})


dim.rmmMatrix<-function(x){
    p<-ncol(x[[1]])
    n<-0
    for(i in 1:length(x)){
        n<-n+nrow(x[[i]])
    }
    return(c(n,p))
}

#' @export
setMethod("dim",signature("rmmMatrix"),dim.rmmMatrix)


get.colnames.rmmMatrix<-function(x){
    out<-colnames(x[[1]])
    return(out)
}

set.colnames.rmmMatrix<-function(x,value){
    for(i in 1:length(x)){
        colnames(x[[i]])<-value
    }
    x
}

#' @export
setMethod("colnames",signature("rmmMatrix"),get.colnames.rmmMatrix)

#' @export
setMethod("colnames<-",signature("rmmMatrix"),set.colnames.rmmMatrix)


get.rownames.rmmMatrix<-function(x){
    out<-NULL
    if(!is.null(rownames(x[[1]]))){
        n<-dim(x)[1]
        out<-rep('',n)
        TMP<-chunks(x)
        for(i in 1:nrow(TMP)){
            out[(TMP[i,2]:TMP[i,3])]<-rownames(x[[i]])
        }
    }
    return(out)
}

set.rownames.rmmMatrix<-function(x,value){
    TMP<-chunks(x)
    for(i in 1:nrow(TMP)){
        rownames(x[[i]])<-value[(TMP[i,2]:TMP[i,3])]
    }
    x
}

#' @export
setMethod("rownames",signature("rmmMatrix"),get.rownames.rmmMatrix)

#' @export
setMethod("rownames<-",signature("rmmMatrix"),set.rownames.rmmMatrix)


#' @export
as.matrix.rmmMatrix<-function(x){
    x[,,drop=FALSE]
}


#' Finds the position of a set of rows in an rmmMatrix object.
rowindexes<-function(x,rows){
    CHUNKS<-chunks(x)
    nRow<-CHUNKS[nrow(CHUNKS),3]

    INDEX<-matrix(nrow=nRow,ncol=3)
    colnames(INDEX)<-c('chunk','row.global','row.local')
    INDEX[,2]<-1:nRow
    end<-0
    for(i in 1:length(x)){
        ini<-end+1
        end<-ini+CHUNKS[i,3]-CHUNKS[i,2]
        INDEX[ini:end,1]<-i
        INDEX[ini:end,3]<-1:nrow(x[[i]])
    }
    if(!is.null(rows)){ INDEX<-INDEX[rows,] }
    if(is.vector(INDEX)){
        tmp<-names(INDEX)
        INDEX<-matrix(INDEX,ncol=3)
        colnames(INDEX)<-tmp
    }
    return(INDEX)
}

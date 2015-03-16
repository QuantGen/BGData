#' An S4 class to represent a row-distributed \code{mmMatrix}
#'
#' \code{rDMatrix} inherits from \code{\link{list}}. Each element of the list is
#' an \code{ff_matrix} object.
#'
#' @export rDMatrix
#' @exportClass rDMatrix
rDMatrix<-setClass('rDMatrix',contains='list')

#' @export
setMethod('initialize','rDMatrix',function(.Object,nrow=1,ncol=1,vmode='byte',folderOut=NULL,nChunks=NULL,dimorder=c(2,1)){
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
})


subset.rDMatrix<-function(x,i,j,drop){
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
setMethod("[",signature(x="rDMatrix",i="ANY",j="ANY",drop="ANY"),subset.rDMatrix)

#' @export
setMethod("[",signature(x="rDMatrix",i="ANY",j="missing",drop="ANY"),function(x,i,drop){
    j<-1:ncol(x)
    subset.rDMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="rDMatrix",i="missing",j="ANY",drop="ANY"),function(x,j,drop) {
    i<-1:nrow(x)
    subset.rDMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="rDMatrix",i="missing",j="missing",drop="ANY"),function(x,drop) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    subset.rDMatrix(x,i,j,drop)
})


replace.rDMatrix<-function(x,i,j,...,value){
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
setReplaceMethod("[",signature(x="rDMatrix",i="ANY",j="ANY",value="ANY"),replace.rDMatrix)

#' @export
setReplaceMethod("[",signature(x="rDMatrix",i="ANY",j="missing",value="ANY"),function(x,i,value){
    j<-1:ncol(x)
    replace.rDMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="rDMatrix",i="missing",j="ANY",value="ANY"),function(x,j,value) {
    i<-1:nrow(x)
    replace.rDMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="rDMatrix",i="missing",j="missing",value="ANY"),function(x,value) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    replace.rDMatrix(x,i,j,value=value)
})


dim.rDMatrix<-function(x){
    p<-ncol(x[[1]])
    n<-0
    for(i in 1:length(x)){
        n<-n+nrow(x[[i]])
    }
    return(c(n,p))
}

#' @export
setMethod("dim",signature("rDMatrix"),dim.rDMatrix)


get.colnames.rDMatrix<-function(x){
    out<-colnames(x[[1]])
    return(out)
}

set.colnames.rDMatrix<-function(x,value){
    for(i in 1:length(x)){
        colnames(x[[i]])<-value
    }
    x
}

#' @export
setMethod("colnames",signature("rDMatrix"),get.colnames.rDMatrix)

#' @export
setMethod("colnames<-",signature("rDMatrix"),set.colnames.rDMatrix)


get.rownames.rDMatrix<-function(x){
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

set.rownames.rDMatrix<-function(x,value){
    TMP<-chunks(x)
    for(i in 1:nrow(TMP)){
        rownames(x[[i]])<-value[(TMP[i,2]:TMP[i,3])]
    }
    x
}

#' @export
setMethod("rownames",signature("rDMatrix"),get.rownames.rDMatrix)

#' @export
setMethod("rownames<-",signature("rDMatrix"),set.rownames.rDMatrix)


#' Finds the position of a set of rows in an rDMatrix object.
rowindexes<-function(x,rows){
    TMP<-chunks(x)
    nRow<-(TMP[nrow(TMP),ncol(TMP)])

    INDEX<-matrix(nrow=nRow,ncol=3)
    colnames(INDEX)<-c('chunk','row.global','row.local')
    INDEX[,2]<-1:nRow
    end<-0
    for(i in 1:length(x)){
        ini<-end+1
        end<-ini+TMP[i,3]-TMP[i,2]
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

#' An S4 class to represent a column-distributed \code{mmMatrix}.
#'
#' \code{cmmMatrix} inherits from \code{\link{list}}. Each element of the list is
#' an \code{ff_matrix} object.
#'
#' @export cmmMatrix
#' @exportClass cmmMatrix
cmmMatrix<-setClass('cmmMatrix',contains='list')

#' @export
setMethod('initialize','cmmMatrix',function(.Object,nrow=1,ncol=1,vmode='byte',folderOut=NULL,nChunks=NULL,dimorder=c(2,1)){
    if(is.null(folderOut)){
        folderOut<-paste0(tempdir(),'/mmMatrix-',randomString())
    }
    if(file.exists(folderOut)){
        stop(paste('Output folder',folderOut,'already exists. Please move it or pick a different one.'))
    }
    dir.create(folderOut)
    if(is.null(nChunks)){
        chunkSize<-min(ncol,floor(.Machine$integer.max/nrow/1.2))
        nChunks<-ceiling(ncol/chunkSize)
    }else{
        chunkSize<-ceiling(ncol/nChunks)
        if(chunkSize*nrow >= .Machine$integer.max/1.2){
            stop('More chunks are needed')
        }
    }
    ffList<-list()
    end<-0
    for(i in 1:nChunks){
        ini<-end+1
        end<-min(ncol,ini+chunkSize-1)
        filename=paste0('geno_',i,'.bin')
        ffList[[i]]<-ff(vmode=vmode,dim=c(nrow,(end-ini+1)),dimorder=dimorder,filename=paste0(folderOut,'/',filename))
        # Change ff path to a relative one
        physical(ffList[[i]])$pattern<-'ff'
        physical(ffList[[i]])$filename<-filename
    }
    .Object<-callNextMethod(.Object,ffList)
    return(.Object)
})


subset.cmmMatrix<-function(x,i,j,drop){
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
    originalOrder<-(1:p)[order(j)]
    sortedColumns<-sort(j)

    dimX<-dim(x)
    if( p>dimX[2] | n>dimX[1] ){
        stop('Either the number of columns or number of rows requested exceed the number of rows or columns in x, try dim(x)...')
    }

    Z<-matrix(nrow=n,ncol=p,NA)
    colnames(Z)<-colnames(x)[j]
    rownames(Z)<-rownames(x)[i]

    INDEXES<-colindexes(x,columns=sortedColumns)

    whatChunks<-unique(INDEXES[,1])
    end<-0
    for(k in whatChunks){
        TMP<-matrix(data=INDEXES[INDEXES[,1]==k,],ncol=3)
        ini<-end+1
        end<-ini+nrow(TMP)-1
        Z[,ini:end]<-x[[k]][i,TMP[,3],drop=FALSE]
    }
    if(length(originalOrder)>1){
        Z[]<-Z[,originalOrder]
    }
    if(drop==TRUE&&(n==1||p==1)){
        # Revert drop.
        return(Z[,])
    }else{
        return(Z)
    }
}

#' @export
setMethod("[",signature(x="cmmMatrix",i="ANY",j="ANY",drop="ANY"),subset.cmmMatrix)

#' @export
setMethod("[",signature(x="cmmMatrix",i="ANY",j="missing",drop="ANY"),function(x,i,drop){
    j<-1:ncol(x)
    subset.cmmMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="cmmMatrix",i="missing",j="ANY",drop="ANY"),function(x,j,drop) {
    i<-1:nrow(x)
    subset.cmmMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="cmmMatrix",i="missing",j="missing",drop="ANY"),function(x,drop) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    subset.cmmMatrix(x,i,j,drop)
})


replace.cmmMatrix<-function(x,i,j,...,value){
    Z<-matrix(nrow=length(i),ncol=length(j),data=value)
    CHUNKS<-chunks(x)
    for(k in 1:nrow(CHUNKS) ){
        col_z<-(j>=CHUNKS[k,2])&(j<=CHUNKS[k,3])
        colLocal<-colindexes(x,j[col_z])[,3]
        x[[k]][i,colLocal]<-Z[,col_z]
    }
    return(x)
}

#' @export
setReplaceMethod("[",signature(x="cmmMatrix",i="ANY",j="ANY",value="ANY"),replace.cmmMatrix)

#' @export
setReplaceMethod("[",signature(x="cmmMatrix",i="ANY",j="missing",value="ANY"),function(x,i,value){
    j<-1:ncol(x)
    replace.cmmMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="cmmMatrix",i="missing",j="ANY",value="ANY"),function(x,j,value) {
    i<-1:nrow(x)
    replace.cmmMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="cmmMatrix",i="missing",j="missing",value="ANY"),function(x,value) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    replace.cmmMatrix(x,i,j,value=value)
})


dim.cmmMatrix<-function(x){
    n<-nrow(x[[1]])
    p<-0
    for(i in 1:length(x)){
        p<-p+ncol(x[[i]])
    }
    return(c(n,p))
}

#' @export
setMethod("dim",signature("cmmMatrix"),dim.cmmMatrix)


get.colnames.cmmMatrix<-function(x){
    out<-NULL
    if(!is.null(colnames(x[[1]]))){
        p<-dim(x)[2]
        out<-rep('',p)
        TMP<-chunks(x)
        for(i in 1:nrow(TMP)){
            out[(TMP[i,2]:TMP[i,3])]<-colnames(x[[i]])
        }
    }
    return(out)
}

set.colnames.cmmMatrix<-function(x,value){
    TMP<-chunks(x)
    for(i in 1:nrow(TMP)){
        colnames(x[[i]])<-value[(TMP[i,2]:TMP[i,3])]
    }
    x
}

#' @export
setMethod("colnames",signature("cmmMatrix"),get.colnames.cmmMatrix)

#' @export
setMethod("colnames<-",signature("cmmMatrix"),set.colnames.cmmMatrix)


get.rownames.cmmMatrix<-function(x){
    out<-rownames(x[[1]])
    return(out)
}

set.rownames.cmmMatrix<-function(x,value){
    for(i in 1:length(x)){
        rownames(x[[i]])<-value
    }
    x
}

#' @export
setMethod("rownames",signature("cmmMatrix"),get.rownames.cmmMatrix)

#' @export
setMethod("rownames<-",signature("cmmMatrix"),set.rownames.cmmMatrix)


#' Finds the position of a set of columns in a cmmMatrix object.
colindexes<-function(x,columns){

    TMP<-chunks(x)
    nCol<-(TMP[nrow(TMP),ncol(TMP)])

    INDEX<-matrix(nrow=nCol,ncol=3)
    colnames(INDEX)<-c('chunk','col.global','col.local')
    INDEX[,2]<-1:nCol
    end<-0
    for(i in 1:length(x)){
        ini<-end+1
        end<-ini+TMP[i,3]-TMP[i,2]
        INDEX[ini:end,1]<-i
        INDEX[ini:end,3]<-1:ncol(x[[i]])
    }
    if(!is.null(columns)){ INDEX<-INDEX[columns,] }
    if(is.vector(INDEX)){
        tmp<-names(INDEX)
        INDEX<-matrix(INDEX,ncol=3)
        colnames(INDEX)<-tmp
    }
    return(INDEX)
}

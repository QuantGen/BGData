#' An S4 class to represent a column-distributed \code{dMatrix}.
#' 
#' \code{cDMatrix} inherits from \code{\link{list}}. Each element of the list is
#' an \code{ff_matrix} object.
#' 
#' @export cDMatrix
#' @exportClass cDMatrix
cDMatrix<-setClass('cDMatrix',contains='list')

#' @export
setMethod('initialize','cDMatrix',function(.Object,nrow=1,ncol=1,vmode='byte',folderOut=NULL,nChunks=NULL,dimorder=c(2,1)){
    if(is.null(folderOut)){
        folderOut<-paste0(tempdir(),'/dMatrix-',randomString())
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


#' An S4 class to represent a row-distributed \code{dMatrix}
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
        folderOut<-paste0(tempdir(),'/dMatrix-',randomString())
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


setClassUnion('dMatrix',c('cDMatrix','rDMatrix'))


## Defines method dim() (extract # of rows and number of columns) of an object getnosFF ##
dim.cDMatrix<-function(x){
   n<-nrow(x[[1]])
   p<-0
   for(i in 1:length(x)){
     p<-p+ncol(x[[i]])
   }
   return(c(n,p))
}

dim.rDMatrix<-function(x){
   p<-ncol(x[[1]])
   n<-0
   for(i in 1:length(x)){
     n<-n+nrow(x[[i]])
   }
   return(c(n,p))
}

 # sets method dim for cDMatrix and rDMatrix objects

#' @export
setMethod("dim",signature("cDMatrix"),dim.cDMatrix)

#' @export
setMethod("dim",signature("rDMatrix"),dim.rDMatrix)

## end of dim ############################################################################


#' Provides information about how data is distributed into binary files.
#' 
#' Row-distributed (\code{\linkS4class{rDMatrix}}) and column-distributed 
#' matrices (\code{\linkS4class{cDMatrix}}) have the content of the array mapped
#' to possibly multiple binary files. Each chunk is an \code{ff_matrix} object. 
#' \code{chunks} gives, for each chunk, the row or column indexes at which each 
#' chunk start and ends.
#' 
#' @param x Either an \code{\linkS4class{rDMatrix}} or a 
#'   \code{\linkS4class{cDMatrix}} object
#' @return A matrix with information per chunk in rows.
#' @export
chunks<-function(x){
    if(class(x)=='cDMatrix'){
        n<-length(x)
        OUT<-matrix(nrow=n,ncol=3,NA)
        colnames(OUT)<-c('chunk','col.ini','col.end')
        end<-0
        for(i in 1:n){
            ini<-end+1
            end<-ini+ncol(x[[i]])-1
            OUT[i,]<-c(i,ini,end)
        }
        return(OUT)
    }else if(class(x)=='rDMatrix'){
        n<-length(x)
        OUT<-matrix(nrow=n,ncol=3,NA)
        colnames(OUT)<-c('chunk','row.ini','row.end')
        end<-0
        for(i in 1:n){
            ini<-end+1
            end<-ini+nrow(x[[i]])-1
            OUT[i,]<-c(i,ini,end)
        }
        return(OUT)
    }else{
        stop("Unsupported type.")
    }
}

## colnames method for cDMatrix and rDMatrix ###########################################################
get.colnames.cDMatrix<-function(x){
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

get.colnames.rDMatrix<-function(x){
    out<-colnames(x[[1]])
    return(out)
}

set.colnames.cDMatrix<-function(x,value){
    TMP<-chunks(x)
    for(i in 1:nrow(TMP)){
        colnames(x[[i]])<-value[(TMP[i,2]:TMP[i,3])]
    }
    x
}

set.colnames.rDMatrix<-function(x,value){
    for(i in 1:length(x)){
        colnames(x[[i]])<-value
    }
    x
}

#' @export
setMethod("colnames",signature("cDMatrix"),get.colnames.cDMatrix)

#' @export
setMethod("colnames<-",signature("cDMatrix"),set.colnames.cDMatrix)

#' @export
setMethod("colnames",signature("rDMatrix"),get.colnames.rDMatrix)

#' @export
setMethod("colnames<-",signature("rDMatrix"),set.colnames.rDMatrix)

## end of colnames #######################################################################

## rownames method for cDMatrix and rDMatrix ##########################################################
get.rownames.cDMatrix<-function(x){
    out<-rownames(x[[1]])
    return(out)
}

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

set.rownames.cDMatrix<-function(x,value){
    for(i in 1:length(x)){
        rownames(x[[i]])<-value
    }
    x
}

set.rownames.rDMatrix<-function(x,value){
    TMP<-chunks(x)
    for(i in 1:nrow(TMP)){
        rownames(x[[i]])<-value[(TMP[i,2]:TMP[i,3])]
    }
    x
}

#' @export
setMethod("rownames",signature("cDMatrix"),get.rownames.cDMatrix)

#' @export
setMethod("rownames<-",signature("cDMatrix"),set.rownames.cDMatrix)

#' @export
setMethod("rownames",signature("rDMatrix"),get.rownames.rDMatrix)

#' @export
setMethod("rownames<-",signature("rDMatrix"),set.rownames.rDMatrix)

# end of rownames ########################################################################

## dimnames method for cDMatrix and rDMatrix ##########################################################
get.dimnames<-function(x){
    list(rownames(x), colnames(x))
}

#' @export
setMethod("dimnames",signature("dMatrix"),get.dimnames)

# end of dimnames ########################################################################

## finds the position of a set of columns in an object cDMatrix ###########################
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
########################################################################################

## finds the position of a set of rows in an object rDMatrix ###########################
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
#########################################################################################

## Indexing for cDMatrix objects ########################################################
subset.cDMatrix<-function(x,i,j,drop){
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
setMethod("[",signature(x="cDMatrix",i="ANY",j="ANY",drop="ANY"),subset.cDMatrix)

#' @export
setMethod("[",signature(x="cDMatrix",i="ANY",j="missing",drop="ANY"),function(x,i,drop){
    j<-1:ncol(x)
    subset.cDMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="cDMatrix",i="missing",j="ANY",drop="ANY"),function(x,j,drop) {
    i<-1:nrow(x)
    subset.cDMatrix(x,i,j,drop)
})

#' @export
setMethod("[",signature(x="cDMatrix",i="missing",j="missing",drop="ANY"),function(x,drop) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    subset.cDMatrix(x,i,j,drop)
})


replace.cDMatrix<-function(x,i,j,...,value){
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
setReplaceMethod("[",signature(x="cDMatrix",i="ANY",j="ANY",value="ANY"),replace.cDMatrix)

#' @export
setReplaceMethod("[",signature(x="cDMatrix",i="ANY",j="missing",value="ANY"),function(x,i,value){
    j<-1:ncol(x)
    replace.cDMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="cDMatrix",i="missing",j="ANY",value="ANY"),function(x,j,value) {
    i<-1:nrow(x)
    replace.cDMatrix(x,i,j,value=value)
})

#' @export
setReplaceMethod("[",signature(x="cDMatrix",i="missing",j="missing",value="ANY"),function(x,value) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    replace.cDMatrix(x,i,j,value=value)
})
 
## end of indexing cDMatrix #################################################################### 

## Indexing for rDMatrix objects ##########################################################
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

## end of indexing cDMatrix #################################################################### 


apply.DMatrix<-function(X,MARGIN,FUN,chunkSize=1e3,verbose=TRUE,...){
    FUN<-match.fun(FUN)
    if(!(class(X)%in%c('rDMatrix','cDMatrix'))){ stop('X must be either dMatrix or rMatrix') }
    
    n<-ifelse(MARGIN==1,nrow(X),ncol(X))
 
    if(MARGIN==1){  x<-X[1,] }else{ x<-X[,1] }

    tmp<-FUN(x,...)
    
    if(is.atomic(tmp)){
	ANS<-matrix(nrow=length(tmp),ncol=n,NA)
	rownames(ANS)<-names(tmp)
	if(MARGIN==1){
		colnames(ANS)<-rownames(X)
	}else{
           colnames(ANS)<-colnames(X)	    
	}
	nChunks<-ceiling(n/chunkSize)
    	end<-0
    	for(i in 1:nChunks){
        	if(verbose){cat(i,' out of ',nChunks,' \n')}
        	ini<-end+1
        	end<-min(ini+chunkSize-1,n)
        	if(MARGIN==1){
            	Z<-X[ini:end,]
        	}else{
            	Z<-X[,ini:end]
        	}
         	ANS[,ini:end]<-apply(FUN=FUN,MARGIN=MARGIN, X=Z,...)
    	}
    }else{
    	ANS<-vector('list',n)
    	names(ANS)<-ifelse(MARGIN==1,rownames(X),colnames(X))
    	end<-0
    	for(i in 1:n){
  		if(verbose){cat(i,' out of ',n,' \n')}
		if(MARGIN==1){
				ANS[[i]]<-FUN(X[i,],...)
		}else{
			ANS[[i]]<-FUN(X[,i],...)
		}    		
    	}
    }
    return(ANS[,,drop=TRUE])
}

#' Apply function for \code{\linkS4class{rDMatrix}} or 
#' \code{\linkS4class{cDMatrix}} objects.
#' 
#' This function brings chunks of data (of size \code{chunkSize}) from the 
#' distributed array into RAM as \code{matrix} objects and calls \code{apply} of
#' the base package to obtain the summaries for the chunk. Results from all the 
#' chunks are collected and returned.
#' 
#' @param X Either an \code{\linkS4class{rDMatrix}} or a 
#'   \code{\linkS4class{cDMatrix}} object.
#' @param MARGIN Use 1 to obtain row summaries or 2 to obtain column summaries.
#' @param chunkSize The number of columns or rows that are processed at a time 
#'   (see Details).
#' @return Returns a \code{matrix} or a \code{list} with results from FUN.
#' @export
setMethod("apply",signature("dMatrix"),apply.DMatrix)


colMeans.DMatrix<-function(x,na.rm=TRUE,chunkSize=1e3,...){
    if(na.rm){
        warning('Ignoring missing values')
    }
    ANS<-apply.DMatrix(X=x,MARGIN=2,FUN=mean,chunkSize=chunkSize,na.rm=na.rm,...)
    return(ANS)
}

#' @export
setMethod("colMeans",signature("dMatrix"),colMeans.DMatrix)


colSums.DMatrix<-function(x,na.rm=TRUE,chunkSize=1e3,...){
    if(na.rm){
        warning('Ignoring missing values')
    }
    ANS<-apply.DMatrix(X=x,MARGIN=2,FUN=sum,chunkSize=chunkSize,na.rm=na.rm,...)
    return(ANS)
}

#' @export
setMethod("colSums",signature("dMatrix"),colSums.DMatrix)


rowMeans.DMatrix<-function(x,na.rm=TRUE,chunkSize=1e3,...){
    if(na.rm){
        warning('Ignoring missing values')
    }
    ANS<-apply.DMatrix(X=x,MARGIN=1,FUN=mean,chunkSize=chunkSize,na.rm=na.rm,...)
    return(ANS)
}

#' @export
setMethod("rowMeans",signature("dMatrix"),rowMeans.DMatrix)


rowSums.DMatrix<-function(x,na.rm=TRUE,chunkSize=1e3,...){
    if(na.rm){
        warning('Ignoring missing values')
    }
    ANS<-apply.DMatrix(X=x,MARGIN=1,FUN=sum,chunkSize=chunkSize,na.rm=na.rm,...)
    return(ANS)
}

#' @export
setMethod("rowSums",signature("dMatrix"),rowSums.DMatrix)


summary.num<-function(x){
    out<-c(range(x,na.rm=T),mean(x,na.rm=T),sd(x,na.rm=T),mean(is.na(x)))
    names(out)<-c('min','max','mean','sd','prop NAs')
    return(out)
}

summary.char<-function(x){
    out<-table(x,useNA='always')
    out<-out/length(x)
    return(out)
}

summary.DMatrix<-function(object,MARGIN=2,chunkSize=1e3,...){
    # If MARGIN==1 summaries of columns are provided, this is the default, otherwise, row-summaries are returned.
    if(is.numeric(object[1,1])){
        ANS<-apply.DMatrix(X=object,MARGIN=MARGIN,FUN=summary.num,chunkSize=chunkSize,...)
    }else{
       if(is.character(object[1,1])|is.logical(object[1,1])){
           ANS<-apply.DMatrix(X=object,MARGIN=MARGIN,FUN=summary.char,chunkSize=chunkSize,...)
       }else{
           ANS<-apply.DMatrix(X=object,MARGIN=MARGIN,FUN=summary,chunkSize=chunkSize,...)
       }

    }
    return(ANS)
}

#' @export
setMethod("summary",signature("dMatrix"),summary.DMatrix)


randomString<-function(){
    paste(sample(c(0:9,letters,LETTERS),size=5,replace=TRUE),collapse="")
}

#' @include cDMatrix.R rDMatrix.R
NULL


setClassUnion('dMatrix',c('cDMatrix','rDMatrix'))


get.dimnames<-function(x){
    list(rownames(x), colnames(x))
}

#' @export
setMethod("dimnames",signature("dMatrix"),get.dimnames)


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

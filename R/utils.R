
crossprods.chunk<-function(chunk,x,y=NULL,nChunks,use_tcrossprod=FALSE){
  ## Performs crossprod() or tcrossprod()
  #  for a chunk (set of columns or sets of rows) of x
 
  n<-ifelse(use_tcrossprod,ncol(x),nrow(x))
 
  if(!is.null(y)){ y=as.matrix(y) }
  chunkID=rep(1:nChunks,each=ceiling(n/nChunks))[1:n]
  if(!is.null(y)){ y=y[chunkID==chunk,,drop=FALSE]}

  if(use_tcrossprod){
      X=x[,chunkID==chunk,drop=FALSE]
      Xy=tcrossprod(X,y)
  }else{
      X=x[chunkID==chunk,,drop=FALSE]
      Xy=crossprod(X,y)
  }
  return(Xy)
}

#' Computes crossprod (x'y or x'x) or tcrossprod (xy' or xx', used when use_tcrossprod=TRUE) in parallel.
#' @param x matrix, ff_matrix, rmmMatrix or cmmMatrix
#' @param y, vector, matrix, ff_matrix, rmmMatrix or cmmMatrix. By default
#' @param nChunks the number of 'chunks' used when X and y are partitioned.
#' @param mc.cores  the number of cores (passed to mclapply)
#' @param use_tcrossprod  if FALSE crossprods computes x'y or x'x (if y=NULL), otherwise crossprods computes xy' or xx' (if y=NULL).
#' @return xx', x'x, x'y, or xy' depending on whether y is provided and on whether use_tcrossprod=TRUE/FALSE
#' @export
crossprods<-function(x,y=NULL,nChunks=detectCores(),mc.cores=detectCores(),use_tcrossprod=FALSE){
  # Computes crossprod(x,y) or tcrossprod(x,y)
  if(nChunks==1){
    if(use_tcrossprod){
     Xy=tcrossprod(x,y)
    }else{
      Xy=crossprod(x,y)
    }
  }else{ 
    tmpIndex=1:nChunks
    TMP=mclapply(X=tmpIndex,FUN=crossprods.chunk,x=x,y=y,nChunks=nChunks,mc.cores=mc.cores,use_tcrossprod=use_tcrossprod)
     ## We now need to add up chunks sequentially
     Xy=TMP[[1]]
     if(length(TMP)>1){
        for(i in 2:length(TMP)){
          Xy=Xy+TMP[[i]]
        }
     }
   }
   return(Xy)
}

#' Computes a genomic relationship matrix G=xx'.
#' 
#' Offers options for centering and scaling the columns of x before computing
#' xx'. If \code{centerCol=FALSE}, \code{scaleCol=FALSE} and
#' \code{scaleG=FALSE}, \code{getG} produces the same outcome than
#' \code{tcrossprod}.
#' 
#' @param x matrix, ff_matrix, rmmMatrix or cmmMatrix
#' @param nChunks The number of columns that are processed at a time.
#' @param scaleCol TRUE/FALSE whether columns must be scaled before computing
#'   xx'.
#' @param scaleG TRUE/FALSE whether columns must be scaled before computing xx'.
#' @param i (integer, boolean or character) Indicates which rows should be used.
#'   By default, all rows are used.
#' @param j (integer, boolean or character) Indicates which columns should be 
#'   used. By default, all columns are used.
#' @return A positive semi-definite symmetric numeric matrix.
#' @export
getG<-function(x,nChunks=ceiling(ncol(x)/1e3),scaleCol=TRUE,scaleG=TRUE,verbose=TRUE,i=1:nrow(x),j=1:ncol(x),minVar=1e-5,
               nChunks2=detectCores(),mc.cores=detectCores()){
    nX<-nrow(x); pX<-ncol(x); centerCol=TRUE # if this is made a parameter the imputation od NAs need to be modified.
    
    # converting boolean to integer index (it leads to a more efficient subsetting than booleans)
    if(is.logical(i)){ i<-which(i) }
    if(is.logical(j)){ j<-which(j) }
    
    n<-length(i); 	p<-length(j)
    
    if(n>nX|p>pX){ stop('Index out of bounds')}
    
    if(is.numeric(i)){ if( (min(i)<1)|(max(i)>nX)){ stop('Index out of bounds') }}
    if(is.numeric(j)){ if( (min(j)<1)|(max(j)>pX)){ stop('Index out of bounds') }}
    
    tmp<-x[i,1:2]
    n<-nrow(tmp)
    
    G<-matrix(0,nrow=n,ncol=n)
    rownames(G)<-rownames(tmp)
    colnames(G)<-rownames(G)
    
    end<-0;
    delta<-ceiling(p/nChunks);
    
    for(k in 1:nChunks){
        ini<-end+1;
        if(ini<=p){
            end<-min(p,ini+delta-1)
            if(verbose){
        	    cat("Chunk: ",k," (markers ", ini,":",end," ~",round(100*end/p,1),"% done)\n",sep="");
                cat("  =>Acquiring genotypes...\n")
            }
        
            # subset
            tmp<-j[ini:end]
            X=x[i,tmp,drop=FALSE];
        
            if(scaleCol){
                VAR<-apply(X=X,FUN=var,MARGIN=2,na.rm=TRUE)
                tmp<-which(VAR<minVar)
                if(length(tmp)>0){
                    X<-X[,-tmp]
                    VAR<-VAR[-tmp]
                }
            }
        
            if(ncol(X)>0){
                if(verbose){ cat("  =>Computing...\n") }
                if(centerCol|scaleCol){
                    X<-scale(X,center=centerCol,scale=scaleCol)
                }
                TMP<-is.na(X)
                if(any(TMP)){    X<-ifelse(TMP,0,X) }

              if(nChunks2>1){
                TMP<-crossprods(x=X,use_tcrossprod=TRUE,nChunks=nChunks2,mc.cores=mc.cores)
              }else{
                TMP<-tcrossprod(X)
              }
              G<-G+TMP
            }
        }
    }
    if(scaleG){
        tmp<-mean(diag(G))
        G<-G/tmp
    }
    
    return(G)
}

#' Generate and store a simulated plaintext raw PED file (see \code{--recodeA}
#' in PLINK) or PED-like file for testing purposes.
#' 
#' @param filename The path where to save the generated file.
#' @param n The number of observations to generate.
#' @param p The number of markers to generate.
#' @param genoChars The alphabet used to generate the genotypes.
#' @param na.string The symbol used to denote missing values.
#' @param propNA The probability of generating NAs.
#' @param returnGenos Whether to return the genotypes from the function.
#' @export
simPED<-function(filename,n,p,genoChars=0:2,na.string=NA,propNA=.02,returnGenos=FALSE){
    if(file.exists(filename)){
        stop(paste('File',filename,'already exists. Please move it or pick a different name.'))
    }
    markerNames<-paste0('mrk_',1:p)
    subjectNames<-paste0('id_',1:n)
    if(returnGenos){
        OUT<-matrix(nrow=n,ncol=p,NA)
        colnames(OUT)<-markerNames
        rownames(OUT)<-subjectNames
    }
    fileOut<-file(filename,open='w')
    pedP<-6+p
    header<-c(c('FID','IID','PAT','MAT','SEX','PHENOTYPE'),markerNames)
    write(header,ncolumns=pedP,append=TRUE,file=fileOut)
    for(i in 1:n){
        geno<-sample(genoChars,size=p,replace=TRUE)
        geno[runif(p)<propNA]<-na.string
        pheno<-c(0,subjectNames[i],rep(NA,4))
        x<-c(pheno,geno)
        write(x,ncolumns=pedP,append=TRUE,file=fileOut)
        if(returnGenos){
            OUT[i,]<-geno
        }
    }
    close(fileOut)
    if(returnGenos){
        return(OUT)
    }
}


randomString<-function(){
    paste(sample(c(0:9,letters,LETTERS),size=5,replace=TRUE),collapse="")
}


normalizeType<-function(val){
    type<-typeof(val)
    # detect strings
    if(type=='character'&&length(val)>0){
        # convert to type if type and value match
        convert<-try(vector(mode=val),silent=TRUE)
        if(class(convert)=='try-error'){
            # return a character type if conversion failed
            warning('could no convert type, using character instead')
            character()
        }else{
            # return conversion result otherwise
            convert
        }
        # value doesn't contain type information and can be handled by typeof
    }else{
        val
    }
}

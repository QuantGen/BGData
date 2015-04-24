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
getG<-function(x,nChunks=ceiling(ncol(x)/1e3),scaleCol=TRUE,scaleG=TRUE,verbose=TRUE,i=1:nrow(x),j=1:ncol(x),minVar=1e-5){
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
                G<-G+tcrossprod(X)
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
    write(header,ncol=pedP,append=TRUE,file=fileOut)
    for(i in 1:n){
        geno<-sample(genoChars,size=p,replace=TRUE)
        geno[runif(p)<propNA]<-na.string
        pheno<-c(0,subjectNames[i],rep(NA,4))
        x<-c(pheno,geno)
        write(x,ncol=pedP,append=TRUE,file=fileOut)
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


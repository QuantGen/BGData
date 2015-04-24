#' Computes a genomic relationship matrix G=xx'.
#' 
#' Offers options for centering and scaling the columns of x before computing
#' xx'. If \code{centerCol=FALSE}, \code{scaleCol=FALSE} and
#' \code{scaleG=FALSE}, \code{getG} produces the same outcome than
#' \code{tcrossprod}.
#' 
#' @param x matrix, ff_matrix, rmmMatrix or cmmMatrix
#' @param nChunks The number of columns that are processed at a time.
#' @param centerCol TRUE/FALSE whether columns must be centered before computing
#'   xx'.
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


### 
#' Creates a BGData object (reads and recode) from GBS data in the format provided by Jesse Poland's lab
#' @return A BGData object.
#' @export

readGBS<-function( fileIn, rsCol,allelesCol,nColSkip,header=TRUE,
                   class='matrix',n=NULL,p=NULL,na.strings='N',returnData=TRUE,
                   verbose=TRUE,nChunks=NULL,vmode='byte',dimorder=1:2,
                   folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",
                                   basename(fileIn)),sep='') 
                  ){
    ########
    # Reads and recode GBS data
    # Value: a BGData object with slots @map and @geno
    ########
    if(is.null(n)){
        tmp<-scan(gzfile(fileIn),nlines=1,what=character(),quiet=T)
        n<-length(tmp)-nColSkip
    }
    if(is.null(p)){
       tmp<-gzfile(fileIn,open='r')
       p <- 0
       while(length(scan(tmp,what=character(),nlines=1,quiet=TRUE))>0){
            p<-p+1
            if(verbose){
              cat('Determining number of markers (reading line ',p, ')\n')
            }
       }
       if(header){  p<-p-1}   
       close(tmp)
     }    
     if(header){
        pedFile<-gzfile(fileIn,open='r')
        tmp<-scan(pedFile,nlines=1,what=character(),quiet=TRUE)
        GID<-tmp[-(1:nColSkip)]
    }else{ GID<-paste0('ID_',1:n)}
    if(class=='matrix'){
        geno<-matrix(nrow=n,ncol=p)
    }else{
        geno<-new(class,nrow=n,ncol=p,vmode=vmode,folderOut=folderOut,
                  nChunks=nChunks,dimorder=dimorder)
    }
    rownames(geno)<-GID
    fileIn<-gzfile(fileIn,open='r')
    map<-matrix(nrow=p,ncol=nColSkip,NA)
    if(header){ tmp<-scan(fileIn,nlines=1,quiet=TRUE,what=character())
                colnames(map)<-tmp[1:nColSkip] 
              }

    for(i in 1:p){
        time<-proc.time()        
        x<-scan(fileIn,nlines=1,what=character(),quiet=TRUE)
        alleleOne<-unlist(strsplit(x[allelesCol],split='/'))[1]
        map[i,]<-x[(1:nColSkip)]
        x<-x[-(1:nColSkip)]
        z<-rep(0,n)
        z[x==alleleOne]<-2
        z[x=='H']<-1
        z[x==na.strings]<-NA
        geno[,i]<-z
        if(verbose){
            cat('Marker ',i,' (of ',p, '),',round(proc.time()[3]-time[3],3),
                 ' sec./marker.\n',sep='')
        }
    }
    close(fileIn)

    # Adding names
    colnames(geno)<-map[,rsCol]
  
    pheno<-data.frame(GID=GID,stringsAsFactors=FALSE)
    MAP<-data.frame(type.convert(map[,1],as.is=TRUE),stringsAsFactors=FALSE)
    if(ncol(map)>2){
    	for(i in 2:ncol(map)){
    		MAP<-cbind(MAP,type.convert(map[,i],as.is=TRUE))
    	}
    }
    colnames(MAP)<-colnames(map)
 
    BGData<-new('BGData',geno=geno,pheno=pheno,map=MAP)
    if(class!='matrix'){
        attr(BGData,'origFile')<-list(path=fileIn,dataType=typeof(dataType))
        attr(BGData,'dateCreated')<-date()
        save(BGData,file=paste(folderOut,'/BGData.RData',sep=''))
    }
    if(returnData){
        return(BGData)
    }
}

setOldClass('ff_matrix') # Convert ff_matrix into an S4 class
setClassUnion('geno',c('dMatrix','matrix','ff_matrix'))


#' An S4 class to represent GWAS data.
#' 
#' @slot pheno A \code{\link{data.frame}} that contains phenotypes.
#' @slot map A \code{\link{data.frame}} that contains a genetic map.
#' @slot geno A \code{geno} object (\code{dMatrix}, \code{ff_matrix}, or
#'   \code{\link{matrix}}) that contains genotypes.
#' @export genData
#' @exportClass genData
genData<-setClass('genData',slots=c(pheno='data.frame',map='data.frame',geno='geno'))

#' @export
setMethod('initialize','genData',function(.Object,geno,pheno,map){
    if(!is(geno,'geno')){
        stop("Only dMatrix, ff_matrix, or regular matrix objects are allowed for geno.")
    }
    if(is.null(colnames(geno))){
        colnames(geno)<-paste0('mrk_',1:ncol(geno))
    }
    if(is.null(rownames(geno))){
        rownames(geno)<-paste0('id_',1:nrow(geno))
    }
    if(missing(pheno)){
        pheno<-data.frame(IID=rownames(geno))
        rownames(pheno)<-rownames(geno)
    }
    if(missing(map)){
        map<-data.frame(mrk=colnames(geno))
        rownames(map)<-colnames(geno)
    }
    .Object@geno<-geno
    .Object@pheno<-pheno
    .Object@map<-map
    return(.Object)
})


#' Creates a \code{\linkS4class{genData}} object from a plaintext file.
#' 
#' \code{setGenData} assumes that the plaintext file (\code{fileIn}) contains 
#' records of individuals in rows, and phenotypes, covariates and markers in 
#' columns. The columns included in columns \code{1:nColSkip} are used to 
#' populate the slot \code{\code{@@pheno}} of a \code{\linkS4class{genData}}
#' object, and the remaining columns are used to fill the slot
#' \code{\code{@@geno}}. If the first row contains a header
#' (\code{header=TRUE}), data in this row is used to determine variables names
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}.
#' Genotypes are stored in a distributed matrix (\code{dMatrix}). By default a
#' column-distributed (\code{\linkS4class{cDMatrix}}) is used for \code{@@geno},
#' but the user can modify this using the \code{distributed.by} argument. The
#' number of chunks is either specified by the user (use \code{nChunks} when
#' calling \code{setGenData}) or determined internally so that each
#' \code{ff_matrix} object has a number of cells that is smaller than 
#' \code{.Machine$integer.max/1.2}. \code{setGenData} creates a folder 
#' (\code{folderOut}) that contains the binary flat files (\code{geno_*.bin}) 
#' and the \code{\linkS4class{genData}} object (typically named
#' \code{genData.RData}. Optionally (if \code{returnData} is TRUE) it returns
#' the \code{\linkS4class{genData}} object to the environment. The filename of
#' the \code{ff_matrix} objects are saved as relative names. Therefore, to be
#' able to access the content of the data included in \code{@@geno} the working
#' directory must either be the folder where these files are saved
#' (\code{folderOut}) or the object must be loaded using the \code{loadGenData}
#' function included in the package.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use 'character' for A/C/G/T or 
#'   'integer' for numeric coding.
#' @param distributed.by If columns a column-distributed matrix 
#'   (\code{\linkS4class{cDMatrix}}) is created, if rows a row-distributed 
#'   matrix (\code{\linkS4class{rDMatrix}}).
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param folderOut The path to the folder where to save the binary files.
#' @param returnData If TRUE, the function returns a 
#'   \code{\linkS4class{genData}} object.
#' @param na.strings The character string use to denote missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype 
#'   information in the file.
#' @param idCol The index of the ID column.
#' @param verbose If TRUE, progress updates will be posted.
#' @param nChunks The number of chunks to create.
#' @param dimorder The physical layout of the chunks.
#' @return If \code{returnData} is TRUE, a \code{\linkS4class{genData}} object 
#'   is returned.
#' @export
setGenData<-function(fileIn,header,dataType,distributed.by='columns',n=NULL,p=NULL,
                     folderOut=paste('genData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep=''),
                     returnData=TRUE,na.strings='NA',nColSkip=6,idCol=2,verbose=FALSE,nChunks=NULL,
                     dimorder=if(distributed.by=='rows') 2:1 else 1:2){
    
    if(file.exists(folderOut)){
        stop(paste('Output folder',folderOut,'already exists. Please move it or pick a different one.'))
    }
    if(!dataType%in%c('character','integer','numeric')){
        stop('dataType must be either character, integer or numeric')
    }
    if(!distributed.by%in%c('columns','rows')){
        stop('distributed.by must be either columns or rows')
    }
    
    vMode<-ifelse(dataType%in%c('character','integer'),'byte','double')
    
    if(is.null(n)){
        # gzfile and readLines throw some warnings, but since it works, let's
        # disable warnings for this block
        warnLevel<-unlist(options('warn'))
        options(warn=-1)
        detN<-gzfile(fileIn,open='r')
        n <- 0
        while(length(readLines(detN,n=1))>0){
            n<-n+1
        }
        if(header){
            n<-n-1
        }
        close(detN)
        # restore previous warning level
        options(warn=warnLevel)
    }
    if(header){
        pedFile<-gzfile(fileIn,open='r')
        tmp<-scan(pedFile,nlines=1,what=character(),quiet=TRUE)
        p<-length(tmp)-nColSkip
        phtNames<-tmp[1:nColSkip]
        mrkNames<-tmp[-(1:nColSkip)]
    }else{
        if(is.null(p)){
            detP<-gzfile(fileIn,open='r')
            tmp<-scan(detP,nlines=1,what=character(),quiet=TRUE)
            p<-length(tmp)-nColSkip
            close(detP)
        }
        pedFile<-gzfile(fileIn,open='r')
        phtNames<-paste('v_',1:nColSkip,sep='')
        mrkNames<-paste('mrk_',1:p,sep='')
    }
    
    IDs<-rep(NA,n)
    
    pheno<-matrix(nrow=n,ncol=nColSkip)
    colnames(pheno)<-phtNames
    
    geno<-new(ifelse(distributed.by=='columns','cDMatrix','rDMatrix'),nrow=n,ncol=p,vmode=vMode,folderOut=folderOut,nChunks=nChunks,dimorder=dimorder)
    colnames(geno)<-mrkNames
    
    for(i in 1:n){
        time<-proc.time()
        xSkip<-scan(pedFile,n=nColSkip,what=character(),na.strings=na.strings,quiet=TRUE)
        x<-scan(pedFile,n=p,what=dataType,na.strings=na.strings,quiet=TRUE)
        pheno[i,]<-xSkip
        IDs[i]<-xSkip[idCol]
        geno[i,]<-x
        if(verbose){
            cat('Subject',i,' ',round(proc.time()[3]-time[3],3),'sec/subject.','\n')
        }
    }
    close(pedFile)
    
    # Adding names
    rownames(geno)<-IDs
    rownames(pheno)<-IDs
    
    pheno<-as.data.frame(pheno,stringsAsFactors=FALSE)
    pheno[]<-lapply(pheno,type.convert,as.is=TRUE)
    
    genData<-new('genData',geno=geno,pheno=pheno)
    
    attr(genData,'origFile')<-list(path=fileIn,dataType=dataType)
    attr(genData,'dateCreated')<-date()
    
    save(genData,file=paste(folderOut,'/genData.RData',sep=''))
    
    if(returnData){ return(genData) }
}


#' @export
loadGenData<-function(path,envir=.GlobalEnv){
    ##
    # Use: to load a genData object using the name of the folder where the meta-data and data are stored.
    # path: the name of the folder where the data and meta data are stored.
    # envir: the name of the environment where the object is returned.
    # See also: load2() and setGenData()
    ##
    if('genData'%in%ls(envir=envir)){
        stop('There is already an object called genData in the environment. Please move it.')
    }
    if(!file.exists(paste0(path,'/genData.RData'))){
        stop(paste('Could not find a genData object in path',path))
    }
    cwd<-getwd()
    setwd(path)
    load('genData.RData',envir)
    cat('Loaded genData object into environment under name genData')
    # Open all chunks for reading (we do not store absolute paths to ff files,
    # so this has to happen in the same working directory)
    chunks<-chunks(genData@geno)
    for(i in 1:nrow(chunks)){
        open(genData@geno[[i]])
    }
    # Restore working directory
    setwd(cwd)
}


#' @export
load2<-function(file,envir=parent.frame(),verbose=TRUE){
    ##
    # Function to load genData or dMatrix objects
    # file: the name of the .RData file to be loaded (and possibly a path)
    # envir: the environment where to load the data
    # verbose: TRUE/FALSE
    # See also: loadGenData()
    ##
    
    # determining the object name
    lsOLD<-ls(); 
    load(file=file)
    lsNEW<-ls();
    objectName<-lsNEW[(!lsNEW%in%lsOLD)&(lsNEW!='lsOLD')];  
    
    # determining path and filename
    path<-dirname(file)
    fname<-basename(file)
    
    # stores current working directiory and sets working directory to path
    cwd<-getwd()
    setwd(path)
    
    # determining object class
    objectClass<-class(eval(parse(text=objectName)))
    
    if(verbose){ 
        cat(' Meta data (',fname,') and its data were stored at folder ',path,'.\n',sep='')
        cat(' Object Name: ',objectName,'\n',sep='')
        cat(' Object Class: ',objectClass,'\n',sep='')
    }
    if(!(objectClass%in%c('genData','rDMatrix','cDMatrix'))){ stop( ' Object class must be either genData, cDMatrix or rDMatrix')}
    
    # Determining number of chunks
    if(objectClass=='genData'){
        tmpChunks<-chunks(eval(parse(text=paste0(objectName,'@geno'))))
    }else{
        tmpChunks<-chunks(eval(parse(text=objectName)))
    }
    
    # opening files
    for(i in 1:nrow(tmpChunks)){
        if(verbose){ cat(' Opening flat file ', i,'\n')  }
        if(objectClass=='genData'){
            open(eval(parse(text=paste0(objectName,'@geno[[',i,']]'))))
        }else{
            open(eval(parse(text=paste0(objectName,'[[',i,']]'))))
        }
    }  
    # sending the object to envir
    assign(objectName,get(objectName), envir=envir)
    
    # restoring the working directory
    setwd(cwd)
    if(verbose){ cat(' Original directory (',getwd(),') restored \n',sep='')}
}


#' Conducts an association study (GWAS) using a \code{\linkS4class{genData}}
#' object.
#' 
#' This function conducts an association test using the formula provided by the 
#' user (\code{formula}) plus one column of \code{@@geno}, one column at a time.
#' The data from the association tests is obtained from a 
#' \code{\linkS4class{genData}} object.
#' 
#' @param formula A formula (e.g. weight~sex+age) with the response on the 
#'   left-hand side and predictors (all the covariates except the markers) on 
#'   the right-hand side. The variables included in the formula must be in the 
#'   \code{@@pheno} object of the \code{\linkS4class{genData}}.
#' @param data A \code{\linkS4class{genData}} object.
#' @param method The regression method to be used. Currently, the following 
#'   methods are implemented: \code{\link{lm}}, \code{\link{lm.fit}}, 
#'   \code{\link{lsfit}}, \code{\link{glm}} and \code{\link[lme4]{lmer}}.
#' @param plot If TRUE a Manhattan plot is produced and filled with points as 
#'   the association tests are run.
#' @param verbose If TRUE more messages are printed.
#' @param min.pValue Numeric, the minimum p-value expected, used to determine 
#'   the limits of the vertical axis of the Manhattan plot.
#' @param chunkSize Represents the number of columns of \code{@@geno} that are 
#'   brought into RAM for processing (10 by default).
#' @return Returns a matrix with estimates, SE, p-value, etc.
#' @export
GWAS<-function(formula,data,method,plot=FALSE,verbose=FALSE,min.pValue=1e-10,chunkSize=10,...){
    if(class(data)!='genData'){ stop('data must genData')}
    
    if(!method%in%c('lm','lm.fit','lsfit','glm','lmer','SKAT')){
        stop('Only lm, glm, lmer and SKAT have been implemented so far.')
    }
    ## We can have 'specialized methods, for instance for OLS it is better to use lsfit that is what GWAS.ols do
    if(method%in%c('lm','lm.fit','lsfit','SKAT')){
        if(method%in%c('lm','lm.fit','lsfit')){
            OUT<-GWAS.ols(formula=formula,data=data,plot=plot,verbose=verbose,min.pValue=min.pValue,chunkSize=,chunkSize,...)	
        }
        if(method%in%c('SKAT')){
            OUT<-GWAS.SKAT(formula=formula,data=data,plot=plot,verbose=verbose,min.pValue=min.pValue,...)	
        }           
    }else{
        FUN<-match.fun(method)
        # could subset based on NAs so that subsetting does not take place in each iteration of the GWAS loop
        pheno<-data@pheno
        fm<-FUN(formula,data=pheno,...)
        tmp<-getCoefficients(fm)
        p<-ncol(data@geno)
        OUT<-matrix(nrow=p,ncol=length(tmp),NA)
        rownames(OUT)<-colnames(data@geno)
        colnames(OUT)<-colnames(tmp)
        GWAS.model<-update(as.formula(formula),'.~z+.')
        if(plot){
            tmp<-paste(as.character(GWAS.model[2]),as.character(GWAS.model[3]),sep='~')
            plot(numeric()~numeric(),xlim=c(0,p),ylim=c(0,-log(min.pValue,base=10)),ylab='-log(p-value)',xlab='Marker',main=tmp)
        }
        nChunks<-ceiling(p/chunkSize)
        end<-0
        tmpRow<-0
        
        for(i in 1:nChunks){
            time.in<-proc.time()[3]
            ini<-end+1
            end<-min(ini+chunkSize-1,p)
            Z<-data@geno[,ini:end,drop=FALSE]
            
            for(j in 1:(end-ini+1)){
                pheno$z<-Z[,j]
                fm<-FUN(GWAS.model,data=pheno,...)
                tmp<-getCoefficients(fm)
                tmpRow<-tmpRow+1
                OUT[tmpRow,]<-tmp
                if(plot){
                    x=c(tmpRow-1,tmpRow)
                    y=-log(OUT[c(tmpRow-1,tmpRow),4],base=10)
                    if(tmpRow>1){ lines(x=x,y=y,col=8,lwd=.5) }
                    points(y=-log(tmp[4],base=10),col=2,cex=.5,x=tmpRow)
                }
            }
            if(verbose){ cat(sep='','Chunk ',i,' of ', nChunks,' (',round(proc.time()[3]-time.in,2),' seconds/chunk, ',round(i/nChunks*100,3),'% done )\n') }
        }
    }
    return(OUT)
}

getCoefficients<-function(x){
    UseMethod('getCoefficients')
}
getCoefficients.lm<-function(x){
    summary(x)$coef[2,]
}
getCoefficients.glm<-function(x){
    summary(x)$coef[2,]
}
getCoefficients.lmerMod<-function(x){
    ans<-summary(x)$coef[2,]
    ans<-c(ans,c(1-pnorm(ans[3])))
    return(ans)
}

## GWAS 'Ordinary least squares' (e.g., lsfit lm.fit lm)
GWAS.ols<-function(formula,data,plot=FALSE,verbose=FALSE,min.pValue=1e-10,chunkSize=10,...){
    ##
    # formula: the formula for the GWAS model without including the marker, e.g., y~1 or y~factor(sex)+age
    # all the variables in the formula must be in data@pheno
    # data (genData) containing slots @pheno and @geno
    ##
    
    X <- model.matrix(formula,data@pheno)
    X <- X[match(rownames(data@pheno),rownames(X)),]
    y<-data@pheno[,as.character(terms(formula)[[2]])]
    p<-ncol(data@geno)
    tmp<-ls.print(lsfit(x=X,y=y,intercept=FALSE),print=FALSE)$coef.table[[1]]
    OUT<-matrix(nrow=p,ncol=ncol(tmp),NA)
    colnames(OUT)<-colnames(tmp)
    rownames(OUT)<-colnames(data@geno)
    X<-cbind(0,X)
    
    if(plot){
        tmp<-paste(as.character(formula[2]),as.character(formula[3]),sep='~')
        plot(numeric()~numeric(),xlim=c(0,p),ylim=c(0,-log(min.pValue,base=10)),ylab='-log(p-value)',xlab='Marker',main=tmp)
    }
    nChunks<-ceiling(p/chunkSize)
    end<-0
    tmpRow<-0
    
    for(i in 1:nChunks){
        time.in<-proc.time()[3]
        ini<-end+1
        end<-min(ini+chunkSize-1,p)
        Z<-data@geno[,ini:end,drop=FALSE]
        
        for(j in 1:(end-ini+1)){
            X[,1]<-Z[,j]
            tmpRow<-tmpRow+1
            fm<-lsfit(x=X,y=y,intercept=FALSE)
            tmp<-ls.print(fm,print=FALSE)$coef.table[[1]][1,]
            OUT[tmpRow,]<-tmp
            
            if(plot){
                tmp.x=c(tmpRow-1,tmpRow)
                tmp.y=-log(OUT[c(tmpRow-1,tmpRow),4],base=10)
                if(tmpRow>1){ lines(x=tmp.x,y=tmp.y,col=8,lwd=.5) }
                points(y=-log(tmp[4],base=10),col=2,cex=.5,x=tmpRow)
            }
        }
        if(verbose){ cat(sep='','Chunk ',i,' of ', nChunks,' (',round(proc.time()[3]-time.in,2),' seconds/chunk, ',round(i/nChunks*100,3),'% done )\n') }
    }
    return(OUT)
}

GWAS.SKAT<-function(formula,data,groups,plot=FALSE,verbose=FALSE,min.pValue=1e-10,...){
    ##
    # formula: the formula for the GWAS model without including the markers, e.g., y~1 or y~factor(sex)+age
    # all the variables in the formula must be in data@pheno
    # data (genData) containing slots @pheno and @geno
    # groups: a vector mapping markers into groups (can be integer, character or factor).
    ##
    
    library(SKAT)
    
    p<-length(unique(groups))
    
    OUT<-matrix(nrow=p,ncol=2,NA)
    colnames(OUT)<-c('nMrk','p-value')
    levels<-unique(groups)
    rownames(OUT)<-levels
    
    H0<-SKAT_Null_Model(formula,data=data@pheno,...)
    
    if(plot){
        tmp<-paste(as.character(formula[2]),as.character(formula[3]),sep='~')
        plot(numeric()~numeric(),xlim=c(0,p),ylim=c(0,-log(min.pValue,base=10)),ylab='-log(p-value)',xlab='Marker',main=tmp)
    }
    
    for(i in 1:p){
        Z<-genData@geno[,groups==levels[i],drop=FALSE]
        fm<-SKAT(Z=Z,obj=H0,...)
        OUT[i,]<-c(ncol(Z),fm$p.value)
        
        if(plot){
            tmp.x=c(i-1,i)
            tmp.y=-log(OUT[tmp.x,2],base=10)
            if(i>1){ lines(x=tmp.x,y=tmp.y,col=8,lwd=.5) }
            points(y=tmp.y[2],col=2,cex=.5,x=i)
        }
    }
    if(verbose){
        cat(sep='','Group ',i,' of ', p,' (',round(proc.time()[3]-time.in,2),' seconds/chunk, ',round(i/p*100,3),'% done )\n')
    }
    return(OUT)
}


#' @export
getG<-function(x,nChunks=3,scaleCol=TRUE,scaleG=TRUE,verbose=TRUE,i=1:nrow(x),j=1:ncol(x),minVar=1e-5){
    ###
    # Computes a genomic relationship matrix G=XX'
    # Offers options for centering and scaling G=WW' where W=scale(X,center=centerCol,scale=scaleCol)
    # And of scaling the final output so that the average diagonal value is equal to one (scaleG=TRUE)
    # If scaleCol=centerCol=scaleG=FALSE it behaves as tcrossprod(X)
    # Arguments:
    #   x: matrix, ff_matrix, rDMatrix or cDMatrix
    #   nChunks: the number of columns that are processed at a time.
    #   scaleCol, centerCol: TRUE/FALSE whether columns must be centered and scaled before computing XX'
    #   i,j: (integer, boolean or character) indicating which columsn and which rows should be used.
    #        By default all columns and rows are used.
    #  Genomic relationship matrix
    # Value: G=XX'
    ###
    
    nX<-nrow(x);   	pX<-ncol(x)
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
        end<-min(p,end+delta-1)
        if(verbose){
            cat("Submatrix: ",k," (out of",nChunks,")\n");
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
            if(scaleCol){
                X<-scale(X,center=TRUE,scale=scaleCol)
            }
            TMP<-is.na(X)
            if(any(TMP)){    X<-ifelse(TMP,0,X) }
            G<-G+tcrossprod(X)
        }
    }
    if(scaleG){
        tmp<-mean(diag(G))
        G<-G/tmp
    }
    return(G)
}


##  Utils

#' @export
simPED<-function(filename,n,p,genoChars=1:4,na.string=0,propNA=.02,returnGenos=FALSE){
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

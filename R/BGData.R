#' @include mmMatrix.R
NULL


setOldClass('ff_matrix') # Convert ff_matrix into an S4 class
setClassUnion('geno',c('mmMatrix','matrix','ff_matrix'))


#' An S4 class to represent GWAS data.
#' 
#' @slot pheno A \code{\link{data.frame}} that contains phenotypes.
#' @slot map A \code{\link{data.frame}} that contains a genetic map.
#' @slot geno A \code{geno} object (\code{mmMatrix}, \code{ff_matrix}, or 
#'   \code{\link{matrix}}) that contains genotypes.
#' @export BGData
#' @exportClass BGData
BGData<-setClass('BGData',slots=c(pheno='data.frame',map='data.frame',geno='geno'))

#' @export
setMethod('initialize','BGData',function(.Object,geno,pheno,map){
    if(!is(geno,'geno')){
        stop("Only mmMatrix, ff_matrix, or regular matrix objects are allowed for geno.")
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


#' Creates a memory-mapped \code{\linkS4class{BGData}} object from a plaintext
#' raw PED file (generated with \code{--recodeA} in PLINK) or PED-like file.
#' 
#' \code{readPED} assumes that the plaintext file (\code{fileIn}) contains 
#' records of individuals in rows, and phenotypes, covariates and markers in 
#' columns. The columns included in columns \code{1:nColSkip} are used to 
#' populate the slot \code{\code{@@pheno}} of a \code{\linkS4class{BGData}} 
#' object, and the remaining columns are used to fill the slot 
#' \code{\code{@@geno}}. If the first row contains a header 
#' (\code{header=TRUE}), data in this row is used to determine variables names 
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}. 
#' Genotypes are stored in a distributed matrix (\code{mmMatrix}). By default a 
#' column-distributed (\code{\linkS4class{cmmMatrix}}) is used for 
#' \code{@@geno}, but the user can modify this using the \code{distributed.by} 
#' argument. The number of chunks is either specified by the user (use 
#' \code{nChunks} when calling \code{readPED}) or determined internally so that 
#' each \code{ff_matrix} object has a number of cells that is smaller than 
#' \code{.Machine$integer.max/1.2}. \code{readPED} creates a folder 
#' (\code{folderOut}) that contains the binary flat files (\code{geno_*.bin}) 
#' and the \code{\linkS4class{BGData}} object (typically named 
#' \code{BGData.RData}. Optionally (if \code{returnData} is TRUE) it returns the
#' \code{\linkS4class{BGData}} object to the environment. The filename of the 
#' \code{ff_matrix} objects are saved as relative names. Therefore, to be able 
#' to access the content of the data included in \code{@@geno} the working 
#' directory must either be the folder where these files are saved 
#' (\code{folderOut}) or the object must be loaded either using 
#' \code{loadBGData} or \code{load2}.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use 'character()' for A/C/G/T or 
#'   'integer()' for numeric coding.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param na.strings The character string used in the plaintext file to denote 
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype 
#'   information in the file.
#' @param idCol The index of the ID column.
#' @param returnData If TRUE, the function returns a \code{\linkS4class{BGData}}
#'   object.
#' @param verbose If TRUE, progress updates will be posted.
#' @param nChunks The number of chunks to create.
#' @param distributed.by If columns a column-distributed matrix 
#'   (\code{\linkS4class{cmmMatrix}}) is created, if rows a row-distributed 
#'   matrix (\code{\linkS4class{rmmMatrix}}).
#' @param folderOut The path to the folder where to save the binary files.
#' @param dimorder The physical layout of the chunks.
#' @return If \code{returnData} is TRUE, a \code{\linkS4class{BGData}} object is
#'   returned.
#' @seealso \code{\linkS4class{BGData}}, \code{mmMatrix}, 
#'   \code{\linkS4class{rmmMatrix}}, \code{\linkS4class{cmmMatrix}}, 
#'   \code{\link[ff]{ff}}
#' @export
readPED<-function(fileIn,header,dataType,n=NULL,p=NULL,na.strings='NA',
                  nColSkip=6,idCol=2,returnData=TRUE,verbose=FALSE,
                  nChunks=NULL,distributed.by='columns',
                  folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep=''),
                  dimorder=if(distributed.by=='rows') 2:1 else 1:2){

    if(file.exists(folderOut)){
        stop(paste('Output folder',folderOut,'already exists. Please move it or pick a different one.'))
    }
    if(!typeof(dataType)%in%c('character','integer','numeric')){
        stop('dataType must be either character(), integer() or numeric()')
    }
    if(!distributed.by%in%c('columns','rows')){
        stop('distributed.by must be either columns or rows')
    }

    class<-ifelse(distributed.by=='columns','cmmMatrix','rmmMatrix')
    vmode<-ifelse(typeof(dataType)%in%c('character','integer'),'byte','double')

    readPED.default(fileIn=fileIn,header=header,dataType=dataType,class=class,
                    n=n,p=p,na.strings=na.strings,nColSkip=nColSkip,idCol=idCol,
                    returnData=returnData,verbose=verbose,nChunks=nChunks,
                    vmode=vmode,folderOut=folderOut,dimorder=dimorder)
}

#' Creates a \code{\linkS4class{BGData}} object from a plaintext PED-like file.
#' 
#' \code{readPED.matrix} assumes that the plaintext file (\code{fileIn}) 
#' contains records of individuals in rows, and phenotypes, covariates and 
#' markers in columns. The columns included in columns \code{1:nColSkip} are 
#' used to populate the slot \code{\code{@@pheno}} of a 
#' \code{\linkS4class{BGData}} object, and the remaining columns are used to 
#' fill the slot \code{\code{@@geno}}. If the first row contains a header 
#' (\code{header=TRUE}), data in this row is used to determine variables names 
#' for \code{@@pheno} and marker names for \code{@@map} and \code{@@geno}.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use 'character()' for A/C/G/T or 
#'   'integer()' for numeric coding.
#' @param n The number of individuals.
#' @param p The number of markers.
#' @param na.strings The character string used in the plaintext file to denote 
#'   missing value.
#' @param nColSkip The number of columns to be skipped to reach the genotype 
#'   information in the file.
#' @param idCol The index of the ID column.
#' @param verbose If TRUE, progress updates will be posted.
#' @return Returns a \code{\linkS4class{BGData}} object.
#' @seealso \code{\linkS4class{BGData}}
#' @export
readPED.matrix<-function(fileIn,header,dataType,n=NULL,p=NULL,
                         na.strings='NA',nColSkip=6,idCol=2,
                         verbose=FALSE){

    readPED.default(fileIn=fileIn,header=header,dataType=dataType,class="matrix",
                    n=n,p=p,na.strings=na.strings,nColSkip=nColSkip,idCol=idCol,
                    returnData=TRUE,verbose=verbose)
}

readPED.default<-function(fileIn,header,dataType,class,n=NULL,p=NULL,na.strings='NA',
                          nColSkip=6,idCol=2,returnData=TRUE,verbose=FALSE,nChunks=NULL,
                          vmode=NULL,folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep=''),
                          dimorder=NULL){

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

    if(class=='matrix'){
        geno<-matrix(nrow=n,ncol=p)
    }else{
        if(is.null(dimorder)){
            if(class=='rmmMatrix'){
                dimorder<-2:1
            }else{
                dimorder<-1:2
            }
        }
        geno<-new(class,nrow=n,ncol=p,vmode=vmode,folderOut=folderOut,nChunks=nChunks,dimorder=dimorder)
    }

    colnames(geno)<-mrkNames

    for(i in 1:n){
        time<-proc.time()
        xSkip<-scan(pedFile,n=nColSkip,what=character(),quiet=TRUE)
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

    BGData<-new('BGData',geno=geno,pheno=pheno)

    if(class!='matrix'){
        attr(BGData,'origFile')<-list(path=fileIn,dataType=typeof(dataType))
        attr(BGData,'dateCreated')<-date()
        save(BGData,file=paste(folderOut,'/BGData.RData',sep=''))
    }

    if(returnData){
        return(BGData)
    }
}


#' Loads a BGData object using the name of the folder where the metadata and
#' data are stored.
#' 
#' @param path The name of the folder where the data and meta data are stored.
#' @param envir The name of the environment where the object is returned.
#' @seealso \code{\link{load2}} and \code{\link{readPED}}
#' @export
loadBGData<-function(path,envir=.GlobalEnv){
    if('BGData'%in%ls(envir=envir)){
        stop('There is already an object called BGData in the environment. Please move it.')
    }
    if(!file.exists(paste0(path,'/BGData.RData'))){
        stop(paste('Could not find a BGData object in path',path))
    }
    cwd<-getwd()
    setwd(path)
    load('BGData.RData',envir)
    cat('Loaded BGData object into environment under name BGData')
    # Open all chunks for reading (we do not store absolute paths to ff files,
    # so this has to happen in the same working directory)
    chunks<-chunks(get('BGData', envir=envir)@geno)
    for(i in 1:nrow(chunks)){
        open(get('BGData', envir=envir)@geno[[i]])
    }
    # Restore working directory
    setwd(cwd)
}


#' Loads BGData or mmMatrix objects.
#' 
#' @param file The name of the .RData file to be loaded (and possibly a path).
#' @param envir The environment where to load the data.
#' @param verbose TRUE/FALSE
#' @seealso \code{\link{loadBGData}}
#' @export
load2<-function(file,envir=parent.frame(),verbose=TRUE){
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
    if(!(objectClass%in%c('BGData','rmmMatrix','cmmMatrix'))){ stop( ' Object class must be either BGData, cmmMatrix or rmmMatrix')}

    # Determining number of chunks
    if(objectClass=='BGData'){
        tmpChunks<-chunks(eval(parse(text=paste0(objectName,'@geno'))))
    }else{
        tmpChunks<-chunks(eval(parse(text=objectName)))
    }

    # opening files
    for(i in 1:nrow(tmpChunks)){
        if(verbose){ cat(' Opening flat file ', i,'\n')  }
        if(objectClass=='BGData'){
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


#' Performs single marker regressions using a \code{\linkS4class{BGData}} 
#' object.
#' 
#' Implements single marker regressions. The regression model includes all the
#' covariates specified in the right-hand-side of the \code{formula} plus one
#' column of \code{@@geno}, one column at a time. The data from the association
#' tests is obtained from a \code{\linkS4class{BGData}} object.
#' 
#' @param formula A formula (e.g. weight~sex+age) with the response on the 
#'   left-hand side and predictors (all the covariates except the markers) on 
#'   the right-hand side. The variables included in the formula must be in the 
#'   \code{@@pheno} object of the \code{\linkS4class{BGData}}.
#' @param data A \code{\linkS4class{BGData}} object.
#' @param method The regression method to be used. Currently, the following 
#'   methods are implemented: \code{\link{lm}}, \code{\link{lm.fit}}, 
#'   \code{\link{lsfit}}, \code{\link{glm}} and \code{\link[lme4]{lmer}}.
#' @param plot If TRUE a Manhattan plot is produced and filled with points as 
#'   the association tests are run.
#' @param verbose If TRUE more messages are printed.
#' @param min.pValue Numeric, the minimum p-value expected, used to determine 
#'   the limits of the vertical axis of the Manhattan plot.
#' @param chunkSize Represents the number of columns of \code{@@geno} that are 
#'   brought into RAM for processing (5000 by default).
#' @return Returns a matrix with estimates, SE, p-value, etc.
#' @export
GWAS<-function(formula,data,method,plot=FALSE,verbose=FALSE,min.pValue=1e-10,chunkSize=5000,...){
    if(class(data)!='BGData'){ stop('data must BGData')}

    if(!method%in%c('lm','lm.fit','lsfit','glm','lmer','SKAT')){
        stop('Only lm, glm, lmer and SKAT have been implemented so far.')
    }
    ## We can have 'specialized methods, for instance for OLS it is better to use lsfit that is what GWAS.ols do
    if(method%in%c('lm','lm.fit','lsfit','SKAT')){
        if(method%in%c('lm','lm.fit','lsfit')){
            OUT<-GWAS.ols(formula=formula,data=data,plot=plot,verbose=verbose,min.pValue=min.pValue,chunkSize=chunkSize,...)
        }
        if(method=='SKAT'){
            OUT<-GWAS.SKAT(formula=formula,data=data,plot=plot,verbose=verbose,min.pValue=min.pValue,...)
        }
    }else{
        if(method=='lmer'&&!requireNamespace("lme4",quietly=TRUE)){
            stop("lme4 needed for this function to work. Please install it.",call.=FALSE)
        }
        if(method=='lmer'){
            FUN<-lme4::lmer
        }else{
            FUN<-match.fun(method)
        }
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
    # data (BGData) containing slots @pheno and @geno
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
    # data (BGData) containing slots @pheno and @geno
    # groups: a vector mapping markers into groups (can be integer, character or factor).
    ##

    if(!requireNamespace("SKAT",quietly=TRUE)){
        stop("SKAT needed for this function to work. Please install it.",call.=FALSE)
    }

    p<-length(unique(groups))

    OUT<-matrix(nrow=p,ncol=2,NA)
    colnames(OUT)<-c('nMrk','p-value')
    levels<-unique(groups)
    rownames(OUT)<-levels

    H0<-SKAT::SKAT_Null_Model(formula,data=data@pheno,...)

    if(plot){
        tmp<-paste(as.character(formula[2]),as.character(formula[3]),sep='~')
        plot(numeric()~numeric(),xlim=c(0,p),ylim=c(0,-log(min.pValue,base=10)),ylab='-log(p-value)',xlab='Marker',main=tmp)
    }

    for(i in 1:p){
        Z<-data@geno[,groups==levels[i],drop=FALSE]
        fm<-SKAT::SKAT(Z=Z,obj=H0,...)
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

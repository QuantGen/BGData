#' @include package.R
NULL

setOldClass('ff_matrix') # Convert ff_matrix into an S4 class
setOldClass('BEDMatrix') # Convert BEDMatrix into an S4 class

setClassUnion('geno',c('LinkedMatrix','BEDMatrix','matrix','ff_matrix'))


#' An S4 class to represent GWAS data.
#' 
#' @slot pheno A \code{\link{data.frame}} that contains phenotypes.
#' @slot map A \code{\link{data.frame}} that contains a genetic map.
#' @slot geno A \code{geno} object (\code{LinkedMatrix}, \code{ff_matrix}, or
#'   \code{\link{matrix}}) that contains genotypes.
#' @export BGData
#' @exportClass BGData
BGData<-setClass('BGData',slots=c(pheno='data.frame',map='data.frame',geno='geno'))

#' @export
setMethod('initialize','BGData',function(.Object,geno,pheno,map){
    if(!is(geno,'geno')){
        stop("Only LinkedMatrix, BEDMatrix, ff_matrix, or regular matrix objects are allowed for geno.")
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
#' populate the slot \code{@@pheno} of a \code{\linkS4class{BGData}} object, and
#' the remaining columns are used to fill the slot \code{@@geno}. If the first
#' row contains a header (\code{header=TRUE}), data in this row is used to
#' determine variables names for \code{@@pheno} and marker names for
#' \code{@@map} and \code{@@geno}. Genotypes are stored in a distributed matrix
#' (\code{LinkedMatrix}). By default a column-distributed
#' (\code{\linkS4class{ColumnLinkedMatrix}}) is used for \code{@@geno}, but the user can
#' modify this using the \code{distributed.by} argument. The number of chunks is
#' either specified by the user (use \code{nChunks} when calling \code{readPED})
#' or determined internally so that each \code{ff_matrix} object has a number of
#' cells that is smaller than \code{.Machine$integer.max/1.2}. \code{readPED}
#' creates a folder (\code{folderOut}) that contains the binary flat files
#' (\code{geno_*.bin}) and the \code{\linkS4class{BGData}} object (typically
#' named \code{BGData.RData}. Optionally (if \code{returnData} is TRUE) it
#' returns the \code{\linkS4class{BGData}} object to the environment. The
#' filename of the \code{ff_matrix} objects are saved as relative names.
#' Therefore, to be able to access the content of the data included in
#' \code{@@geno} the working directory must either be the folder where these
#' files are saved (\code{folderOut}) or the object must be loaded using
#' \code{load.BGData}.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{integer()} or 
#'   \code{double()} for numeric coding. Character coding is currently not 
#'   supported: use the \code{--recodeA} option of PLINK to convert the PED file
#'   into a raw file.
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
#'   (\code{\linkS4class{ColumnLinkedMatrix}}) is created, if rows a row-distributed
#'   matrix (\code{\linkS4class{RowLinkedMatrix}}).
#' @param folderOut The path to the folder where to save the binary files.
#' @param dimorder The physical layout of the chunks.
#' @return If \code{returnData} is TRUE, a \code{\linkS4class{BGData}} object is
#'   returned.
#' @seealso \code{\linkS4class{BGData}}, \code{LinkedMatrix},
#'   \code{\linkS4class{ColumnLinkedMatrix}}, \code{\linkS4class{RowLinkedMatrix}},
#'   \code{\link[ff]{ff}}
#' @export
readPED<-function(fileIn,header,dataType,n=NULL,p=NULL,na.strings='NA',
                  nColSkip=6,idCol=2,returnData=TRUE,verbose=FALSE,
                  nChunks=NULL,distributed.by='rows',
                  folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep=''),
                  dimorder=if(distributed.by=='rows') 2:1 else 1:2){

    if(file.exists(folderOut)){
        stop(paste('Output folder',folderOut,'already exists. Please move it or pick a different one.'))
    }

    dataType<-normalizeType(dataType)
    if(!typeof(dataType)%in%c('integer','double')){
        stop('dataType must be either integer() or double()')
    }

    if(!distributed.by%in%c('columns','rows')){
        stop('distributed.by must be either columns or rows')
    }

    class<-ifelse(distributed.by=='columns','ColumnLinkedMatrix','RowLinkedMatrix')
    vmode<-ifelse(typeof(dataType)=='integer','byte','double')

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
#' used to populate the slot \code{@@pheno} of a \code{\linkS4class{BGData}}
#' object, and the remaining columns are used to fill the slot \code{@@geno}. If
#' the first row contains a header (\code{header=TRUE}), data in this row is
#' used to determine variables names for \code{@@pheno} and marker names for
#' \code{@@map} and \code{@@geno}.
#' 
#' @param fileIn The path to the plaintext file.
#' @param header If TRUE, the file contains a header.
#' @param dataType The coding of genotypes. Use \code{character()} for A/C/G/T 
#'   or \code{integer()} for numeric coding.
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

    readPED.default(fileIn=fileIn,header=header,dataType=normalizeType(dataType),
                    class="matrix",n=n,p=p,na.strings=na.strings,nColSkip=nColSkip,
                    idCol=idCol,returnData=TRUE,verbose=verbose)
}

readPED.default<-function(fileIn,header,dataType,class,n=NULL,p=NULL,na.strings='NA',
                          nColSkip=6,idCol=2,returnData=TRUE,verbose=FALSE,nChunks=NULL,
                          vmode=NULL,folderOut=paste('BGData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep=''),
                          dimorder=NULL){

    if(is.null(n)){
        n<-getLineCount(fileIn,header)
    }
    if(header){
        headerLine<-getFileHeader(fileIn)
        p<-length(headerLine)-nColSkip
        phtNames<-headerLine[1:nColSkip]
        mrkNames<-headerLine[-(1:nColSkip)]
    }else{
        if(is.null(p)){
            p<-getColumnCount(fileIn)-nColSkip
        }
        phtNames<-paste('v_',1:nColSkip,sep='')
        mrkNames<-paste('mrk_',1:p,sep='')
    }

    # Prepare pheno and add colnames
    pheno<-matrix(nrow=n,ncol=nColSkip)
    colnames(pheno)<-phtNames

    # Prepare geno and add colnames
    if(class=='matrix'){
        geno<-matrix(nrow=n,ncol=p)
    }else{
        if(is.null(dimorder)){
            if(class=='RowLinkedMatrix'){
                dimorder<-2:1
            }else{
                dimorder<-1:2
            }
        }

        # Handle chunking logic
        if(is.null(folderOut)){
            folderOut<-paste0(tempdir(),'/BGData-',randomString())
        }
        if(file.exists(folderOut)){
            stop(paste('Output folder',folderOut,'already exists. Please move it or pick a different one.'))
        }
        dir.create(folderOut)

        # Determine chunk size and number of chunks
        if(is.null(nChunks)){
            if(class=='RowLinkedMatrix'){
                chunkSize<-min(n,floor(.Machine$integer.max/p/1.2))
                nChunks<-ceiling(n/chunkSize)
            }else{
                chunkSize<-min(p,floor(.Machine$integer.max/n/1.2))
                nChunks<-ceiling(p/chunkSize)
            }
        }else{
            if(class=='RowLinkedMatrix'){
                chunkSize<-ceiling(n/nChunks)
                if(chunkSize*p >= .Machine$integer.max/1.2){
                    stop('More chunks are needed')
                }
            }else{
                chunkSize<-ceiling(p/nChunks)
                if(chunkSize*n >= .Machine$integer.max/1.2){
                    stop('More chunks are needed')
                }
            }
        }

        # Initialize list
        geno<-new(class)
        end<-0
        if(class=='RowLinkedMatrix'){
            for(i in 1:nChunks){
                ini<-end+1
                end<-min(n,ini+chunkSize-1)
                filename=paste0('geno_',i,'.bin')
                geno[[i]]<-ff(vmode=vmode,dim=c((end-ini+1),p),dimorder=dimorder,filename=paste0(folderOut,.Platform$file.sep,filename))
                # Change ff path to a relative one
                physical(geno[[i]])$pattern<-'ff'
                physical(geno[[i]])$filename<-filename
            }
        }else{
            for(i in 1:nChunks){
                ini<-end+1
                end<-min(p,ini+chunkSize-1)
                filename=paste0('geno_',i,'.bin')
                geno[[i]]<-ff(vmode=vmode,dim=c(n,(end-ini+1)),dimorder=dimorder,filename=paste0(folderOut,.Platform$file.sep,filename))
                # Change ff path to a relative one
                physical(geno[[i]])$pattern<-'ff'
                physical(geno[[i]])$filename<-filename
            }
        }

        # Generate index
        index<-LinkedMatrix::index(geno)
    }
    colnames(geno)<-mrkNames

    # Parse PED file
    pedFile<-gzfile(fileIn,open='r')
    if(header){
        scan(pedFile,nlines=1,what=character(),quiet=TRUE)
    }
    for(i in 1:n){
        time<-proc.time()
        xSkip<-scan(pedFile,n=nColSkip,what=character(),quiet=TRUE)
        x<-scan(pedFile,n=p,what=dataType,na.strings=na.strings,quiet=TRUE)
        pheno[i,]<-xSkip
        if(class=='matrix'){
            geno[i,]<-x
        }else{
            geno<-`[<-`(geno,i,1:ncol(geno),index=index,value=x)
        }
        if(verbose){
            cat('Subject',i,' ',round(proc.time()[3]-time[3],3),'sec/subject.','\n')
        }
    }
    close(pedFile)

    # Add rownames
    rownames(geno)<-pheno[, idCol]
    rownames(pheno)<-pheno[, idCol]

    # Convert pheno to a data.frame
    pheno<-as.data.frame(pheno,stringsAsFactors=FALSE)
    pheno[]<-lapply(pheno,type.convert,as.is=TRUE)

    # Construct BGData object
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


#' Loads BGData objects.
#' 
#' @param file The name of the .RData file to be loaded.
#' @param envir The environment where to load the data.
#' @param verbose TRUE/FALSE
#' @export
load.BGData<-function(file,envir=parent.frame(),verbose=TRUE){

    # Determine object name
    lsOLD<-ls()
    load(file=file)
    lsNEW<-ls()
    objectName<-lsNEW[(!lsNEW%in%lsOLD)&(lsNEW!='lsOLD')]

    # Determine path and filename
    path<-dirname(file)
    fname<-basename(file)

    # Store current working directory and set working directory to path
    cwd<-getwd()
    setwd(path)

    # Determine object class
    objectClass<-class(get(objectName))
    if(objectClass!='BGData'){
        stop('Object class must be BGData')
    }

    if(verbose){
        cat('Meta data (',fname,') and its data were stored at folder ',path,'.\n',sep='')
        cat('Object Name: ',objectName,'\n',sep='')
        cat('Object Class: ',objectClass,'\n',sep='')
    }

    # Determine number of chunks
    chunks<-LinkedMatrix::chunks(get(objectName)@geno)

    # Open all chunks for reading (we do not store absolute paths to ff files,
    # so this has to happen in the same working directory)
    for(i in 1:nrow(chunks)){
        if(verbose){
            cat('Opening flat file ',i,'\n')
        }
        open(get(objectName)@geno[[i]])
    }

    # Send the object to envir
    assign(objectName,get(objectName),envir=envir)

    # Restore the working directory
    setwd(cwd)
    if(verbose){
        cat('Original directory (',getwd(),') restored \n',sep='')
    }
}

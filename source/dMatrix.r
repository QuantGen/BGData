# Load dependencies
library(ff)

## Defines a classes cDMatrix & rDMatrix ###############################################################
# The class inherits from list, each element of the list is an FF object
# cDMatrix splits the matrix by columns, rDMatrix by rows
setClass('cDMatrix',contains='list')
setClass('rDMatrix',contains='list')

## Idea we can define in the future a class for a collection of rDMatrices or cDMatrices (dDatabase)
# Defiens the class genData the object has three slots:
setClass('genData',slots=c(pheno='data.frame',map='data.frame',geno='list') )


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
 setMethod("dim",signature("cDMatrix"),dim.cDMatrix)
 setMethod("dim",signature("rDMatrix"),dim.rDMatrix)
## end of dim ############################################################################

## chunks: returns which columns are contained in each chunk of a cDMatrix or rDMatrix object ##########
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

setMethod("colnames",signature("cDMatrix"),get.colnames.cDMatrix)
setMethod("colnames<-",signature("cDMatrix"),set.colnames.cDMatrix)
setMethod("colnames",signature("rDMatrix"),get.colnames.rDMatrix)
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

setMethod("rownames",signature("cDMatrix"),get.rownames.cDMatrix)
setMethod("rownames<-",signature("cDMatrix"),set.rownames.cDMatrix)
setMethod("rownames",signature("rDMatrix"),get.rownames.rDMatrix)
setMethod("rownames<-",signature("rDMatrix"),set.rownames.rDMatrix)
# end of rownames ########################################################################

## dimnames method for cDMatrix and rDMatrix ##########################################################
get.dimnames<-function(x){
    list(rownames(x), colnames(x))
}

setMethod("dimnames",signature("cDMatrix"),get.dimnames)
setMethod("dimnames",signature("rDMatrix"),get.dimnames)
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
subset.cDMatrix<-function(x,i,j){
        rows<-i
        columns<-j
        n<-length(rows)
        p<-length(columns)
        originalOrder<-(1:p)[order(columns)]
        columns<-sort(columns)

        dimX<-dim(x)
        if( p>dimX[2] | n>dimX[1] ){
                stop('Either the number of columns or number of rows requested exceed the number of rows or columns in x, try dim(x)...')
        }

        Z<-matrix(nrow=n,ncol=p,NA)
        colnames(Z)<-colnames(x)[columns]
        rownames(Z)<-rownames(x)[rows]

        INDEXES<-colindexes(x,columns=columns)

        whatChunks<-unique(INDEXES[,1])
        end<-0
        for(i in whatChunks){
                TMP<-matrix(data=INDEXES[INDEXES[,1]==i,],ncol=3)
                ini<-end+1; end<-ini+nrow(TMP)-1
                Z[,ini:end]<-x[[i]][rows,TMP[,3],drop=FALSE]
        }
        if(length(originalOrder)>1){
            Z[]<-Z[,originalOrder]
        }
        if(n==1||p==1){
            return(as.vector(Z))
        }else{
            return(Z)
        }
 }

replace.cDMatrix<-function(x,i=(1:nrow(x)),j=(1:ncol(x)),...,value){
    indices<-colindexes(x,j)
    for(k in 1:nrow(indices)){
        index<-indices[k,]
        chunk<-index[1]
        col<-index[3]
        x[[chunk]][i,col]<-value
    }
    x
}

setMethod("[",signature(x="cDMatrix",i="numeric",j="numeric"),subset.cDMatrix)
setMethod("[",signature(x="cDMatrix",i="numeric",j="missing"),function(x,i){
    j<-1:ncol(x)
    subset.cDMatrix(x,i,j)
})
setMethod("[",signature(x="cDMatrix",i="missing",j="numeric"),function(x,j) {
    i<-1:nrow(x)
    subset.cDMatrix(x,i,j)
})
setMethod("[",signature(x="cDMatrix",i="missing",j="missing"),function(x) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    subset.cDMatrix(x,i,j)
})
setReplaceMethod("[",signature("cDMatrix"),replace.cDMatrix)
 
## end of indexing cDMatrix #################################################################### 

#*# check this one
## Indexing for rDMatrix objects ##########################################################
subset.rDMatrix<-function(x,i,j){
        rows<-i
        columns<-j
        n<-length(rows)
        p<-length(columns)
        originalOrder<-(1:n)[order(rows)]
        rows<-sort(rows)

        dimX<-dim(x)
        if( p>dimX[2] | n>dimX[1] ){
                stop('Either the number of columns or number of rows requested exceed the number of rows or columns in x, try dim(x)...')
        }

        Z<-matrix(nrow=n,ncol=p,NA)
        colnames(Z)<-colnames(x)[columns]
        rownames(Z)<-rownames(x)[rows]

        INDEXES<-rowindexes(x,rows=rows)

        whatChunks<-unique(INDEXES[,1])
        end<-0
        for(i in whatChunks){
                TMP<-matrix(data=INDEXES[INDEXES[,1]==i,],ncol=3)
                ini<-end+1; end<-ini+nrow(TMP)-1
                Z[ini:end,]<-x[[i]][TMP[,3],columns,drop=FALSE]
        }
        if(length(originalOrder)>1){
            Z[]<-Z[originalOrder,]
        }
        if(n==1||p==1){
            return(as.vector(Z))
        }else{
            return(Z)
        }
 }

replace.rDMatrix<-function(x,i=(1:nrow(x)),j=(1:ncol(x)),...,value){
    indices<-rowindexes(x,i)
    for(k in 1:nrow(indices)){
        index<-indices[k,]
        chunk<-index[1]
        row<-index[3]
        x[[chunk]][row,j]<-value
    }
    x
}

setMethod("[",signature(x="rDMatrix",i="numeric",j="numeric"),subset.rDMatrix)
setMethod("[",signature(x="rDMatrix",i="numeric",j="missing"),function(x,i){
    j<-1:ncol(x)
    subset.rDMatrix(x,i,j)
})
setMethod("[",signature(x="rDMatrix",i="missing",j="numeric"),function(x,j) {
    i<-1:nrow(x)
    subset.rDMatrix(x,i,j)
})
setMethod("[",signature(x="rDMatrix",i="missing",j="missing"),function(x) {
    i<-1:nrow(x)
    j<-1:ncol(x)
    subset.rDMatrix(x,i,j)
})
setReplaceMethod("[",signature("rDMatrix"),replace.rDMatrix)

## end of indexing cDMatrix #################################################################### 


## Creates and rDMatrix or cDMatrix from a ped file

setGenData<-function(fileIn,n,header,dataType,distributed.by='rows',p=NULL,
                    folderOut=paste('genData_',sub("\\.[[:alnum:]]+$","",basename(fileIn)),sep=''),
                    returnData=TRUE,na.strings='NA',nColSkip=6,idCol=2,verbose=FALSE,nChunks=NULL,
                    dimorder=if(distributed.by=='rows') 2:1 else 1:2){
        ###
        # Use: creates (returns, saves or both) a genData object from an ASCII file
        # fileIn (character): the name of the ped file.
        # n (integer): the number of individuals.
        # dataType : the coding of genotypes, use 'character' for A/C/G/T or 'integer' for numeric coding.
        #            if type='numeric' vmode needs must be 'double'.
        # folderOut (charater): the name of the folder where to save the binary files.
        #                     NOTE: associated to this file there will be files containing the acutal data genos_*.bin
        # returnData (logical): if TRUE the function returns a list with genotypes (X), MAP file and subject information (phenos)
        # header (logical): TRUE if the 1st line of the ped file is a header
        # na.strings (character): the character string use to denote missing value.
        # nColSkip (integer): the number of columsn to be skipped.
        # idCol (integer): the column that contains the subject ID
        # Requires: package ff
        ###

    if(file.exists(folderOut)){
        stop(paste('Output folder',folderOut,'already exists. Please remove it or pick a different one.'))
    }
    dir.create(folderOut)

    vMode<-ifelse( dataType%in%c('character','integer'),'byte','double')

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

	if(!distributed.by%in%c('columns','rows')){stop('distributed.by must be either columns or rows') }
	
    if(is.null(nChunks)){
		if(distributed.by=='columns'){
			chunkSize<-min(p,floor(.Machine$integer.max/n/1.2))
			nChunks<-ceiling(p/chunkSize)
		}else{
			chunkSize<-min(n,floor(.Machine$integer.max/p/1.2))
			nChunks<-ceiling(n/chunkSize)		
		}

	}else{
		if(distributed.by=='columns'){
			chunkSize<-ceiling(p/nChunks)
			if(chunkSize*n >= .Machine$integer.max/1.2){ stop(' More chunks are needed')}
		}else{
			chunkSize<-ceiling(n/nChunks)
			if(chunkSize*p >=  .Machine$integer.max/1.2){ stop(' More chunks are needed')}
		}
	}
	
    genosList<-list()

    end<-0
    for(i in 1:nChunks){
        if(distributed.by=='columns'){
			ini<-end+1
			end<-min(p,ini+chunkSize-1)
			genosList[[i]]<-ff(vmode=vMode,dim=c(n,(end-ini+1)),dimorder=dimorder,filename=paste(folderOut,'/geno_',i,'.bin',sep=''))
			colnames(genosList[[i]])<-mrkNames[ini:end]
		}else{
			ini<-end+1
			end<-min(n,ini+chunkSize-1)
			genosList[[i]]<-ff(vmode=vMode,dim=c((end-ini+1),p),dimorder=dimorder,filename=paste(folderOut,'/geno_',i,'.bin',sep=''))
			colnames(genosList[[i]])<-mrkNames
		}
    }

    for(i in 1:n){
		#if(verbose){ cat(' Subject ',i,'\n')}
		time1<-proc.time()
		xSkip<-scan(pedFile,n=nColSkip,what=character(),na.strings=na.strings,quiet=TRUE)
		x<-scan(pedFile,n=p,what=dataType,na.strings=na.strings,quiet=TRUE)
		
		pheno[i,]<-xSkip
		
		time2<-proc.time()
        IDs[i]<-xSkip[idCol]
        ## now we split x into its chunks
		end<-0
		time3<-proc.time()
		
		if(distributed.by=='columns'){
			for(j in 1:nChunks){
				ini<-end+1
				end<-ini+ncol(genosList[[j]])-1
				genosList[[j]][i,]<-x[ini:end]
			}
		}else{
			tmpChunk<-ceiling(i/chunkSize)
			tmpRow<-i-(tmpChunk-1)*chunkSize
			genosList[[tmpChunk]][tmpRow,]<-x
		}
		
		time4<-proc.time()
        if(verbose){ cat(' Subject ',i,'  ',round(time4[3]-time1[3],3),' sec/subject.','\n')}
    }
    close(pedFile)

	# Adding names
	if(distributed.by=='rows'){
	    end<-0
		for(i in 1:nChunks){
			ini<-end+1
			end<-min(ini+nrow(genosList[[i]])-1,n)
			rownames(genosList[[i]])<-IDs[ini:end]
		}
	}else{
		for(i in 1:nChunks){
			rownames(genosList[[i]])<-IDs
		}
	}
	
	pheno<-as.data.frame(pheno,stringsAsFactors=FALSE)
	pheno[]<-lapply(pheno,type.convert,as.is=TRUE)
	map<-data.frame(mrk=mrkNames,maf=as.numeric(NA),freqNA=as.numeric(NA),stringsAsFactors=FALSE)
	
	geno<-new(ifelse(distributed.by=='columns','cDMatrix','rDMatrix'),genosList)
	genData<-new('genData',geno=geno,map=map,pheno=pheno)

    for(i in 1:nChunks){
        attr(attributes(genData@geno[[i]])$physical,"pattern")<-'ff'
        attr(attributes(genData@geno[[i]])$physical,"filename")<-paste('geno_',i,'.bin',sep='')
    }
    save(genData,file=paste(folderOut,'/genData.RData',sep=''))

    if(returnData){ return(genData) }
}
 
## END OF MAKE makeGenosFF ###############################################################


apply.DMatrix<-function(X,MARGIN,FUN,chunkSize=1e3,...){
    FUN<-match.fun(FUN)
    if(!(class(X)%in%c('rDMatrix','cDMatrix'))){ stop('X must be either dMatrix or rMatrix') }

    nCol<-ifelse(MARGIN==1,nrow(X),ncol(X))

    if(MARGIN==1){
        x<-X[1,]
    }else{
        x<-X[,1]
    }
    tmp<-FUN(x,...)

    ANS<-matrix(nrow=length(tmp),ncol=nCol,NA)
    if(MARGIN==1){
        rownames(ANS)<-names(tmp)
        colnames(ANS)<-rownames(X)
    }else{
        rownames(ANS)<-names(tmp)
        colnames(ANS)<-colnames(X)
    }

    nChunks<-floor(nCol/chunkSize)
    end<-0

    for(i in 1:nChunks){
        cat(i,' out of ',nChunks,' \n')
        ini<-end+1
        end<-ini+chunkSize-1
        if(MARGIN==1){
            Z<-X[ini:end,]
        }else{
            Z<-X[,ini:end]
        }
        ANS[,ini:end]<-apply(FUN=FUN,MARGIN=MARGIN, X=Z,...)
    }
    if(end<nCol){
        ini<-end+1
        end<-nCol
        if(MARGIN==1){
            Z<-X[ini:end,]
        }else{
            Z<-X[,ini:end]
        }
        ANS[,ini:end]<-apply(FUN=FUN,MARGIN=MARGIN, X=Z,...)
    }
    return(ANS)
}

setMethod("apply",signature("rDMatrix"),apply.DMatrix)
setMethod("apply",signature("cDMatrix"),apply.DMatrix)


colMeans.DMatrix<-function(x,chunkSize=1e3,...){
    ANS<-apply.DMatrix(X=x,MARGIN=2,FUN=mean,chunkSize=chunkSize,...)
    return(ANS)
}

setMethod("colMeans",signature("rDMatrix"),colMeans.DMatrix)
setMethod("colMeans",signature("cDMatrix"),colMeans.DMatrix)


colSums.DMatrix<-function(x,chunkSize=1e3,...){
    ANS<-apply.DMatrix(X=x,MARGIN=2,FUN=sum,chunkSize=chunkSize,...)
    return(ANS)
}

setMethod("colSums",signature("rDMatrix"),colSums.DMatrix)
setMethod("colSums",signature("cDMatrix"),colSums.DMatrix)


rowMeans.DMatrix<-function(x,chunkSize=1e3,...){
    ANS<-apply.DMatrix(X=x,MARGIN=1,FUN=mean,chunkSize=chunkSize,...)
    return(ANS)
}

setMethod("rowMeans",signature("rDMatrix"),rowMeans.DMatrix)
setMethod("rowMeans",signature("cDMatrix"),rowMeans.DMatrix)


rowSums.DMatrix<-function(x,chunkSize=1e3,...){
    ANS<-apply.DMatrix(X=x,MARGIN=1,FUN=sum,chunkSize=chunkSize,...)
    return(ANS)
}

setMethod("rowSums",signature("rDMatrix"),rowSums.DMatrix)
setMethod("rowSums",signature("cDMatrix"),rowSums.DMatrix)


summary.DMatrix<-function(object,MARGIN=2,chunkSize=1e3,...){
    # If MARGIN==1 summaries of columns are provided, this is the default, otherwise, row-summaries are returned.
    ANS<-apply.DMatrix(X=object,MARGIN=MARGIN,FUN=summary,chunkSize=chunkSize,...)
    return(ANS)
}

setMethod("summary",signature("rDMatrix"),summary.DMatrix)
setMethod("summary",signature("cDMatrix"),summary.DMatrix)


## Example: GWAS using function lm

GWAS<-function(formula,data,method,plot=FALSE,verbose=FALSE,min.pValue=1e-10,chunkSize=10,...){
        ##
        # formula: the formula for the GWAS model without including the marker, e.g., y~1 or y~factor(sex)+age
        # all the variables in the formula must be in data@pheno
        # data (genData) containing slots @pheno and @geno
        # method: a descritpion of the regression method (e.g.,lm, glm...)
        ##
        if(class(data)!='genData'){ stop('data must genData')}
        
        if(!method%in%c('lm','lm.fit','lsfit','glm','lmer')){
                stop('Only lm, glm and lmer have been implemented so far.')
        }
        ## We can have 'specialized methods, for instance for OLS it is better to use lsfit that is what GWAS.ols do
        if(method%in%c('lm','lm.fit','lsfit')){
        	OUT<-GWAS.ols(formula=formula,data=data,plot=plot,verbose=verbose,min.pValue=min.pValue,chunkSize=,chunkSize,...)	
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
        		Z<-as.matrix(data@geno[,ini:end])

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
        Z<-as.matrix(data@geno[,ini:end])

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


## Computes a Genomic Relationship Matrix
getG<-function(x,n_submatrix=3,scaleCol=TRUE,verbose=TRUE,minMAF=1/100){
	#This function takes as input an object of the class cFF or rFF and computes the
	#Genomic relationship matrix
	#Arguments:
	#the cFF object that holds the matrix with Markers (e.g. SNPs or )

	#Obtain the number of rows and columns for the matrix with molecular markers
	rows_cols=dim(x)
	n=rows_cols[1]
	p=rows_cols[2]
	
	to_column=0;
	delta=ceiling(p/n_submatrix);
	G=matrix(0,nrow=n,ncol=n)
	rownames(G)<-rownames(x)
	colnames(G)<-rownames(x)
	
	K<-0
	from_column<-0
    for(k in 1:n_submatrix){
		from_column=to_column+1;
		to_column=min(p,from_column+delta-1)
		if(verbose){
			cat("Submatrix: ",k," (",from_column,":",to_column,")\n");
			cat("  =>Aquiring genotypes...\n")
		}
		X=x[,from_column:to_column];
		tmp<-colMeans(X,na.rm=TRUE)/2
		maf<-ifelse(tmp>.5,1-tmp,tmp)
		VAR<-apply(X=X,FUN=var,MARGIN=2,na.rm=TRUE)
		
		tmp<-which((maf<minMAF)|(VAR<1e-8))
		if(length(tmp)>0){
			X<-X[,-tmp]
			VAR<-VAR[-tmp]
		}

		if(ncol(X)>0){
			cat("  =>Computing...\n")
			if(scaleCol){
				X<-scale(X,center=TRUE,scale=scaleCol)
			}
			X<-ifelse(is.na(X),0,X)
			K<-K+ifelse(scaleCol,ncol(X)*(n-1)/n,sum(VAR))
			G<-G+tcrossprod(X)
		}
		
   }
   G<-G/K
   return(G)
}


##  Utils

simPED<-function(filename,n,p,propNA=.02){
   fileOut<-file(filename,open='w')
   pedP<-6+p
   for(i in 1:n){
        geno<-sample(1:4,size=p,replace=TRUE)
        geno[runif(p)<propNA]<-0
        pheno<-c(0,paste0('id_',i),rep(NA,4))
        x<-c(pheno,geno)
        write(x,ncol=pedP,append=TRUE,file=fileOut)
   }
   close(fileOut)
}



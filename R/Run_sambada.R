#' @title Run sambada on parallel cores
#' @description Read sambadas input file to retrieve necessary information (num indiv etc...), split the dataset using SamBada's Supervision tool, run sambada on the splitted dataset and merge all using Supervision.
#' @author Solange Gaillard, Sylvie Stucki
#' @param genoFile The name of the file in the current directory of genetic information, complient with samBada's format (use prepare_geno to transform it)
#' @param envFile  The name of the file in the current directory of environmental information (use create_env to create it and prepare_env to reduce the correlated dataset and check order)
#' @param idGeno Name of the column in the \ref{genoFile} corresponding to the id of the animals
#' @param idEnv Name of the column in the \ref{envFile} corresponding to the id of the animals
#' @param dimMax Maximum number of environmental variables included in the logistic models. Use 1 for univariate models, 2 for univariate and bivariates models
#' @param cores Number of cores to use. If NULL, the #cores-1 will be used where #cores corresponds to all cores available on your computer.
#' @param All additional parameters in samBada: see documentation. In case you have to specify several words, you can either specify them in one string and separate them with a space or add a vector string
#' @param keepAllFiles logical If TRUE, all parameter files and splitted genofile and log-files are not removed. Default FALSE
sambada_parallel = function(genoFile, envFile, idGeno, idEnv, dimMax=1, cores=NULL, wordDelim=' ', saveType='END BEST 0.05', spatial=NULL, autoCorr=NULL, ShapeFile=NULL, outputFile=NULL, colSupEnv=NULL, colSupMark=NULL, subsetVarEnv=NULL, subsetVarMark=NULL, headers=TRUE, keepAllFiles=FALSE){
  #a faire
  #documenter inputs spéciaux + check saveType and autoCorr
  
  ### Check inputs
  required_string=c('genoFile', 'envFile', 'idGeno', 'idEnv')
  for(i in 1:length(required_string)){
    if(typeof(eval(parse(text = required_string[i])))!='character') stop(paste(required_string[i],"argument supposed to be character"))
  }
  if (!file.exists(genoFile)) stop("genoFile not found.")
  if (!file.exists(envFile)) stop("envFile not found.")
  if(typeof(dimMax)!='double') stop("dimMax supposed to be integer")
  if(dimMax%%1!=0) stop("dimMax supposed to be integer")
  if(dimMax>2) warning('dimMax is >2. Are you sure this is what you want?')
  if(typeof(wordDelim)!='character') stop("wordDelim supposed to be a single string")
  if(nchar(wordDelim)!=1) stop("wordDelim supposed to be a single string")
  add_opt=c('saveType', 'spatial', 'autoCorr', 'shapeFile','outputFile', 'colSupEnv', 'colSupMark', 'subsetVarEnv', 'subsetVarMark')
  for(i in 1:length(add_opt)){
    if(!is.null(eval(parse(text = add_opt[i])))){
      if(typeof(eval(parse(text = add_opt[i])))!='character') stop(paste(add_opt[i],"argument supposed to be character or NULL"))
      eval(parse(text = add_opt[i]))=unlist(strsplit(eval(parse(text = add_opt[i])), split=" ")) #to build character vector if all in one string (if already vector, won't change)
    }
  }
  #check saveType
  saveType=toupper(saveType)
  if(length(saveType)!=3) stop("saveType should be composed of 3 words (either in vector or delimited by a space)")
  if(!(saveType[1] %in% c('END','REAL'))) stop("The first word of saveType should be either END or REAL")
  if(!(saveType[2] %in% c('ALL','SIGNIF','BEST'))) stop("The second word of saveType should be either ALL, SIGNIF or BEST")
  if(!(grepl("[[:digit:]\\.-]",saveType[3]))) stop("The third word of saveType should be a number between 0 and 1")
  if(as.numeric(saveType[3])>1 | as.numeric(saveType[3])<0) stop("The third word of saveType should be a number between 0 and 1")
  
  #check spatial (first two columns checked later)
  spatial[3:5]=toupper(spatial[3:5])
  if(length(spatial)!=5) stop("spatial should be composed of 5 words (either in vector or delimited by a space)")
  if(!(spatial[3] %in% c('SPHERICAL','CARTESIAN'))) stop("The third word of spatial should be either SPHERICAL or CARTESIAN")
  if(!(spatial[4] %in% c('DISTANCE','GAUSSIAN','BISQUARE','NEAREST'))) stop("The fourth word of spatial should be either DISTANCE,GAUSSIAN,BISQUARE or NEAREST")
  if(!(grepl("[[:digit:]\\.-]",spatial[5]))) stop("The fifth word of spatial should be a number")
  if(spatial[4]=='NEAREST' & as.numeric(spatial[5])%%1!=0 ) stop("If the fourth word of spatial is NEAREST, then the fifth word an integer representing the number of neighbours")
  
  #Check autoCorr
  if(length(autoCorr)!=3) stop("autoCorr should be composed of 3 words (either in vector or delimited by a space)")
  if(!(autoCorr[1] %in% c('GLOBAL','LOCAL','BOTH'))) stop("The first word of autoCorr should be either GLOBAL, LOCAL or BOTH")
  if(!(autoCorr[2] %in% c('ENV','MARK','BOTH'))) stop("The second word of autoCorr should be either ENV, MARK or BOTH")
  if(!(grepl("[[:digit:]\\.-]",autoCorr[3]))) stop("The third word of autoCorr should be an integer representing the number of permutation")
  if(as.numeric(autoCorr[3])%%1!=0) stop("The third word of autoCorr should be an integer representing the number of permutation")
  
  #Check shapeFile
  shapeFile=toupper(shapeFile)
  if(!(shapeFile %in% c('YES','NO'))) stop("shapeFile should be either YES or NO")
  
  ### Load required libraries if cores=NULL or cores>1
  if(is.null(cores)){
    loadParallel=TRUE
  } else if (cores>1) {
    loadParallel=TRUE
  } else {
    loadParallel=FALSE
  }
  if(loadParallel==TRUE){
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package \"doParallel\" needed for this function to work. Please install it.", call. = FALSE)
    }
    if (!requireNamespace("foreach", quietly = TRUE)) {
      stop("Package \"foreach\" needed for this function to work. Please install it.", call. = FALSE)
    }
  }

  ### Open genofile to count number of mark and indiv
  print("Setting up parameter file for Supervision") 
  
  # Count num mark and check that idGeno, colSupMark and subsetVarMark present in genoFile
  con <- file(genoFile,"r")
  first_line <- readLines(con,n=1)
  #Check column present
  if(regexpr(idGeno, first_line)[1]==-1) stop("idGeno not present in header line of genoFile")
  if(!is.null(colSupMark)){ #Check that all colSupMark in header of envFile
    if(sum(colSupMark %in% first_line)!=length(colSupMark)) stop('not all colSupMark in header line of genoFile')
  }
  if(!is.null(subsetVarMark)){ #Check that all subsetVarMark in header of envFile
    if(sum(subsetVarMark %in% first_line)!=length(subsetVarMark)) stop('not all subsetVarMark in header line of genoFile')
  }
  #Count num mark
  nummark=lengths(gregexpr(wordDelim,first_line))-1 
  rm(first_line)
  close(con)
  
  # numVarEnv
  env=read.csv(envFile, sep=wordDelim) 
  numVarEnv=ncol(env)
  if(!(idEnv %in% colnames(env))) stop('idEnv not in envFile')
  if(!is.null(colSupEnv)){ #Check that all colSupEnv in header of envFile
    if(sum(colSupEnv %in% colnames(env))!=length(colSupEnv)) stop('not all colSupEnv in header line of envFile')
  }
  if(!is.null(subsetVarEnv)){ #Check that all subsetVarEnv in header of envFile
    if(sum(subsetVarEnv %in% colnames(env))!=length(subsetVarEnv)) stop('not all subsetVarEnv in header line of envFile')
  }
  if(!is.null(spatial)){ #Check that all spatial in header of envFile
    if(sum(spatial[1:2] %in% colnames(env))!=length(spatial)) stop('The first two words of spatial parameter should be in header line of envFile')
  }
  
  if(nrow(env)!=numIndiv){
    stop('length of envfile and genofile do not match!')
  }
  
  # Count numIndiv from genoFile
  con <- file(genoFile,"rb")
  numIndiv = 0L
  while (length(chunk <- readBin(con, "raw", 20000)) > 0) {
    numIndiv = numIndiv + sum(chunk == as.raw(10L))
  }
  numIndiv=numIndiv - 1
  close(con)
  
  genoFileShort=sub('.csv','',sub('-recode','',genoFile)) #without .csv and without -recode if from prepare_geno
  

  
  ### Check order ID. If not same order use prepare_env. Could also use the -env from supervision
  
  con <- file(genoFile,"r")
  #Read file line by line. Since big file, more rapid than read.table
  rm(IDgeno2)
  j=1
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0){
    if(j==1){ #header line
      col_pos=regexpr(idGeno, line) #position beginning of Id column in header line (already checked that present)
      if(col_pos==1){
        col_pos2=0
      } else{
        col_pos2=lengths(gregexpr(wordDelim, substr(line, 1, col_pos))) #number of columns before ID
      }
    } else{ # line after header
      if(col_pos2==0){ #Id in first column
        length_col=regexpr(wordDelim,line)[1]-1
        ID=substr(line, 1, length_col)       
      } else{ #id in subsequent column
        col_delim=gregexpr(wordDelim,line)[[1]] #stores position of column delimiter
        ID=substr(line, col_delim[col_pos2]+1, col_delim[col_pos2+1]-1) 
      }

      if(j==2){ #First line after header, not use rbind
        IDgeno2=ID
      } else{
        IDgeno2=rbind(IDgeno2, ID) #If already exists, rbind IDgeno2
      } 
    } 
    j=j+1
  }
  rm(line)
  close(con)
  
  #If Id from geno and env file do not match => error
  if(!identical(as.vector(IDgeno2), as.vector(env[,idEnv]))){
    stop('IDs or ID order between the envfile and genofile do not match. Please double check your idGeno and idEnv input or use prepare_env function to reorder your envfile')
  }
  
  ### Preparing samBada's parameter file
  
  #Cleaning supervisions mess on.exit
  if(keepAllFiles==FALSE){
    #file.remove(paste0(genoFileShort,'_paramSupervision.txt'), showWarnings = FALSE)
    on.exit(unlink(paste0(genoFileShort,'_paramSupervision.txt')))
    on.exit(unlink(paste0(genoFileShort,'_param',0:(cores-1),'.txt')))
    on.exit(unlink(paste0(genoFileShort,'_param.txt')))
    on.exit(unlink(paste0(genoFileShort,'-mark-',0:(cores-1),'-',(0:(cores-1))*sizeBlock,'.csv'))) 
    on.exit(unlink(paste0(genoFileShort,'-mark-',0:(cores-1),'-',(0:(cores-1))*sizeBlock,'-log.csv'))) 
    for(dimI in 0:dimMax){
      on.exit(unlink(paste0(genoFileShort,'-mark-',0:(cores-1),'-',(0:(cores-1))*sizeBlock,'-Out-',dimI,'.csv'))) 
    }
  }
  
  # Sambada base parameter (nummark will be changed if supervision is used)
  params=c()
  params["HEADERS"]="Yes"
  if(wordDelim!=' '){
    params["WORDDELIM"]=wordDelim
  }
  params["NUMVARENV"]=numVarEnv
  params["NUMMARK"]=nummark
  params["NUMINDIV"]=numIndiv
  params["IDINDIV"]=paste(idEnv, idGeno)
  #write(paste("COLSUPENV ID ID_gen Code Breed Latitude Longitude Country",sep=" "), file = fileParam, append = TRUE, sep = "\n")
  add_opt=c('dimMax', 'saveType', 'spatial', 'autoCorr', 'shapeFile','outputFile', 'colSupEnv', 'colSupMark', 'subsetVarEnv', 'subsetVarMark')
  for(i in 1:length(add_opt)){
    if(!is.null(eval(parse(text = add_opt[i])))){
      params[toupper(add_opt[i])]=paste(eval(parse(text = add_opt[i])), collapse=' ')
    }
  }
  
  ### Detect number of clusters for parallel computing
  
  #Define the number of cores to be used and load libraries for parallel computing
  if(is.null(cores)){
    cores=parallel::detectCores()
    cores=cores[1]-1
  }
  
  # If cores=1 run directly sambada
  if(cores==1 | !is.null(autoCorr)){
    #write sambada parameterfile
    print("Number of cores used: 1, running sambada without splitting marker file")
    paramFile=paste0(genoFileShort,'-param.txt') #name of paramFile: genoFile (without -recode and without .csv) + -param.txt
    write_sambada_parameter(paramFile, params, 'sambada')
    sambada(paramFile, envFile, genoFile)
    exit()
  }
  
  ### Split datafile with supervision
  
  #Write parameter file for supervision
  fileSupervision=paste0(genoFileShort,'_paramSupervision.txt')
  paramsSuper=c()
  paramsSuper["genofile"]=genoFile
  paramsSuper["paramFile"]="null"
  paramsSuper["numEnv"]=1
  paramsSuper["numMark"]=nummark
  paramsSuper["numLigne"]=numIndiv+1 #For the header line
  sizeBlock=ceiling(nummark/cores)
  paramsSuper["sizeBlock"]=sizeBlock
  write_sambada_parameter(fileSupervision, paramsSuper, 'supervision')

  #Run supervision
  print("Running Supervision to divide the parameter file")
  a=supervision(nomFichier=fileSupervision)
  a=supervision(nomFichier=fileSupervision, numBlock = 0, blockSize = 0, maxDimension = 0, selScore = "", scoreThreshold = 0, sortScore = "", wordDelim = ' ')
  
  ### Run sambada
  
  #Prepare Sambada's parameter file for each supervision split (file 1:n-1 are the same)
  params["NUMMARK"]=paste(sizeBlock, nummark,sep=" ") #Sabada need both the number of marker of the block + total number of markers
  params=params[names(params)!="COLSUPENV" & names(params)!="IDINDIV"] #Supervision has deleted the id from molecular file, and we can't have ID for env file only
  params["COLSUPENV"]=paste(append("Code", colSupEnv), collapse=" ")
    
  for(i in 0:(cores-1)){
    fileParam=paste0(genoFileShort,'_param',i,'.txt')
    if(i==(cores-1)){
      SizeLastBlock=nummark-((cores-1)*sizeBlock) #Last bock, number of markers might differ
      params["NUMMARK"]=paste(SizeLastBlock, nummark,sep=" ")
    }
    write_sambada_parameter(fileParam, params, 'sambada')
  }
  
  #Build cluster
  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  #Close cluster on exit
  on.exit(tryCatch({ stopCluster(cl)}, error=function(e){})) #If already closed, do nothing
  
  ### Run sambada on parallel 
  #sambada('ADAPTmap2_paramSupervision_0.txt',envFile,'ADAPTmap2_recode-mark-0-0.csv')
  # use foreach function
  finalMatrix = foreach(i=0:(cores-1), .combine=cbind, .packages='test2') %dopar% {tempMatrix = sambada(paste0(genoFileShort,'_param',i,'.txt'),envFile,paste0(genoFileShort,'-mark-',i,'-',i*sizeBlock,'.csv'))}

  #Close cluser
  stopCluster(cl)
  
  #Supervision merge
  print("Running supervision to merge the files")
  #supervision( base-name.txt, numBlock, blockSize, maxDimension, selScore, scoreThreshold, sortScore, wordDelim)
  supervision(genoFile, cores, sizeBlock, dimMax, selScore = "Both", scoreThreshold = 0.001, sortScore = "Wald", wordDelim = ' ')

  for(dim in 0:dimMax){
    file.rename(paste0(genoFileShort,'-res-',dim,'.csv'), paste0(genoFileShort,'-Out-',dim,'.csv'))
  }
  print(paste("Output in ",paste0(genoFileShort,'-Out-',0:dimMax,'.csv', collapse=' and ')))
}

write_sambada_parameter = function(paramFile, paramMatrix, programType){
  #For sambada, write the corresponding category at beginning of lines
  if(programType=='sambada'){
    write(paste(names(paramMatrix[1]),paramMatrix[1]), paramFile, append=FALSE)
    for(i in 2:length(paramMatrix)){
      write(paste(names(paramMatrix[i]),paramMatrix[i]), paramFile, append=TRUE, sep="\n")
    }    
  # For supervision, just write the value of the parameter (without category)
  } else if(programType=='supervision'){
    write(paramMatrix[1], paramFile, append=FALSE)
    for(i in 2:length(paramMatrix)){
      write(paramMatrix[i], paramFile, append=TRUE, sep="\n")
    }  
  }
}

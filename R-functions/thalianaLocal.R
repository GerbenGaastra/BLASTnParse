# Written by GT Gaastra (r)
# created last 24-02-2012
# modified last 24-02-2012





## Parse raw gene list and return gene synonym, start&end gene
geneRange <- function(geneFile) {
  rawMatrix <- as.matrix(read.csv(geneFile,skip=2,header=TRUE,sep="\t"))
  output <- matrix(NA,nrow(rawMatrix),3)
  colnames(output) <- c("gene","start","stop")
  for( i in 1:nrow(rawMatrix)) {
    range <- unlist(strsplit(rawMatrix[i,"Location"],"[.]{2}"))
    output[i,] <- c(rawMatrix[i,"Synonym"],range)
  }
  output
}

check <- function(geneName,entries,geneStart,geneEnd,output) {
  blastRes <- as.matrix(read.table(blastFile,sep="\t",header=FALSE))
  for(j in j:nrow(output) ) {
    # probes should be unique
    currID <-output[i,"ID"]
    BlastsFromID <- which(blastRes[,1] == currID)
    if( length(BlastsFromID) == 1) {
      # probes should fall into gene range
      if( blastRes[BlastsFromID,9] > geneStart & blastRes[BlastsFromID,10] < geneEnd ) {
        output[i,"unique"] = "yes"
      } else {
        output[i,"unique"] = "outOfRange"
      }
    } else {
      output[i,"unique"] = "non unique"
    }
  }
  output
}

## function will break when dirInput contains dirctories
blastThaliana <- function(dirInput,dirBlast,dirOutput,genome,geneFile,eValue) {
  #setting up parameters
  if( substring(dirOutput,nchar(dirOutput),nchar(dirOutput)) != "/" ) {
    dirOutput <- paste( dirOutput, "/", sep="")
  }
  if( substring(dirBlast,nchar(dirBlast),nchar(dirBlast)) != "/" ) {
    dirBlast <- paste( dirBlast, "/", sep="")
  }
  if( substring(dirInput,nchar(dirInput),nchar(dirInput)) != "/" ) {
    dirInput <- paste( dirInput, "/", sep="")
  }
  ## preparing matrix with information about start&stop bp of genes on chromosome
  geneInfo <- geneRange(geneFile)
  ## lists of files with or without directory
  lFilesnPath <- list.files(path=dirInput,full.names=TRUE)
  lFiles <- list.files()
  for(i in 1:length(lFiles)) {
    ## reading in input file and rewriting for blast
    inputFile <- as.matrix(read.table(lFilesnPath[i],skip=1,header=FALSE))
    temp <- fastaRewrite(inputFile[,c(1,3)])
    ## setting up parameters for blast
    fastaFile <- paste(dirBlast,"tempfile.fasta",sep="")
    blastFile <- paste(dirBlast,lFiles[i],".out",sep="")
    write(temp, fastaFile)
    # performing blast
    blastLocal(genome,fastaFile,blastFile,evalue=eValue)
    #prepping parameters for check function
    GeneName <- substring(lFiles[i],17,nchar(lFiles[i])-4) # parsing gene name
    entries <- which(geneInfo[,"gene"] == GeneName) # getting lines containing current gene
    geneStart <- as.numeric(min(geneInfo[entries,"start"])) - 400 #flanking correction
    geneEnd <- as.numeric(max(geneInfo[entries,"stop"])) + 400 #flanking correction
    output<- cbind(inputFile,rep(NA,nrow(inputFile)))
    colnames(output) <- c("ID","direction","seq","unique")
    res <- check(geneName,entries,geneStart,geneEnd,output) # performing check
    # writing result to file
    write(res,paste(dirOutput,lFiles[i],".out",sep=""))
  }
}



inputDir <- "D:/CompMolBio/BlastLocal/input"
blastDir <- "D:/CompMolBio/BlastLocal/blastFiles"
outputDir <- "D:/CompMolBio/Blastlocal/output"
genome <- "D:/CompMolBio/BlastLocal/TAIR10_chr_all.fas"
geneFile <- "D:/CompMolBio/BlastLocal/NC_003070.ptt"
blastThaliana(inputDir,blastDir,outputDir,genome,geneFile,0.001)







ids <- unique(blastRes[,1])
nr <- lapply(ids,function(x){length(which(blastRes[,1]==x))})
combined <- matrix(c(ids,nr,rep("NA",length(ids))),,3)
combined[,3] 

  
  
  blastLocal(genome_name="C:/Documents and Settings/gbic/My Documents/CHR_I.fna",temp_filename="C:/Documents and Settings/gbic/My Documents/query.faa")

  
  mInput <- as.matrix(read.table("AGI_WUR_Probes.txt",sep="\t",na.strings = c("NA","","-"),header=TRUE))
# Written by GT Gaastra (r)
# created last 24-02-2012
# modified last 24-02-2012



write.table(output,"fasta.input",row.names=FALSE,col.names=FALSE,quote=FALSE)

listUnique <- function(dirInput,dirOutput,genome,eValue) {
  ## function will break when dirInput contains dirctories
  if( substring(dirOutput,nchar(dirOutput),nchar(dirOutput)) != "/" ) {
    dirOutput <- paste( dirOutput, "/", sep="")
  }
  if( substring(dirInput,nchar(dirOutput),nchar(dirInput)) != "/" ) {
    dirInput <- paste( dirInput, "/", sep="")
  }
  lFilesnPath <- list.files(path=dirInput,full.names=TRUE)
  lFiles <- list.files()
  for(i in 1:length(lFiles)) {
    inputFile <- as.matrix(read.table(lFilesnPath[i],skip=1,header=FALSE))
    temp <- fastaRewrite(inputFile[,c(1,3)])
    fastaFile <- paste(dirOutput,"tempfile.fasta",sep="")
    blastFile <- paste(dirOutput,lFiles[i],".out",sep="")
    cat("fastaFile: ",fastaFile,", blastFile ",blastFile,"\n",sep=" ")
    write(temp, fastaFile)
    blastLocal(genome,fastaFile,blastFile,evalue=eValue)
  }

}


dirInput <- "D:/CompMolBio/Centrotype/Sequentie_probes"
dirOutput <- "D:/CompMolBio/Centrotype/output"
genome <- "D:/CompMolBio/BlastLocal/TAIR10_chr_all.fas"
listUnique(dirInput,dirOutput,genome,0.001)






ids <- unique(blastRes[,1])
nr <- lapply(ids,function(x){length(which(blastRes[,1]==x))})
combined <- matrix(c(ids,nr,rep("NA",length(ids))),,3)
combined[,3] 

  
  
  blastLocal(genome_name="C:/Documents and Settings/gbic/My Documents/CHR_I.fna",temp_filename="C:/Documents and Settings/gbic/My Documents/query.faa")

  
  mInput <- as.matrix(read.table("AGI_WUR_Probes.txt",sep="\t",na.strings = c("NA","","-"),header=TRUE))
# Written by GT Gaastra (r)
# created last 24-02-2012
# modified last 24-02-2012
# Wrapper for Wormblast

blastLocal <- function(genome_name="mygenome.fasta",temp_filename = "seq.aligned",task = "blastn", evalue = 0.001){
  command <- paste("blastn -subject=" , genome_name , " -query=", temp_filename, " -out=", temp_filename,".out -outfmt=6 -task=",task," -evalue=",evalue,sep="")
  cat(command,"\n")
  shell(command)
}

formatFASTA <- function(input) {
  if( class(input) != "matrix") { stop("'input' is not an matrix") }
  if( ncol(input) != 2) { stop("'input' has not 2 collumns") }
  output <- NULL
  temp <- apply(mInput,1,function(x){
    temp <- as.vector(c(paste(">",x[1],sep=""),as.character(x[2]) ))
    #output <- as.vector(c(output,temp))
    #output
  })
  for(i in 1:ncol(temp)  ){
    output <- c(output,temp[,i])
  }
}
write.table(output,"fasta.input",row.names=FALSE,col.names=FALSE,quote=FALSE)

  # for(j in 1:nrow(matrixFile)) {
    # temp1 <- paste(">",matrixFile[j,1],sep="") #row name
    # temp2 <- paste(temp1,matrixFile[j,2],sep="\n") # sequence
    # query <- paste(query,temp2,sep="\n")
  # }
  
 blastLocal <- function(genome_name="mygenome.fasta",temp_filename = "seq.aligned",task = "blastn", evalue = 0.0001){
  command <- paste("blastn  -subject=" , genome_name ,
    " -query=", temp_filename,
    " -out=", temp_filename,".out",
    " -outfmt=\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"",
    " -task=",task,
    " -evalue=",evalue,
    sep="")
  shell(command)
}

ids <- unique(blastRes[,1])
nr <- lapply(ids,function(x){length(which(blastRes[,1]==x))})
combined <- matrix(c(ids,nr,rep("NA",length(ids))),,3)
combined[,3] 

  
  
  blastLocal(genome_name="C:/Documents and Settings/gbic/My Documents/CHR_I.fna",temp_filename="C:/Documents and Settings/gbic/My Documents/query.faa")

  
  mInput <- as.matrix(read.table("AGI_WUR_Probes.txt",sep="\t",na.strings = c("NA","","-"),header=TRUE))
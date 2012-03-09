# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Blasts nucluotide seqeunce vs nucleotide database
# on http://www.wormbase.org/db/searches/blast_blat
# To-Do: Retrieves first n sequences

getPosThaliana <- function(query,eValue="1E+0",daba="elegans_genome") {
  # setting up parameters
  uriToPost <- "http://www.arabidopsis.org/cgi-bin/Blast/TAIRblast.pl"
  # names list containing all fields and their values
  postValues <- new("list",
    Algorithm="blastn",
    BlastTargetSet="ATH1_seq",
    textbox="seq",
    QueryText=query,
    Matrix="Blosum62",
    Expectation="0.001",
    ReplyVia="BROWSER",
    ReplyFormat="TABULAR")
    #submit="Run BLAST")
  # Post and download Form
  Sys.sleep(runif(1))
  formData <- postForm(uriToPost, .params = postValues,style="POST")
  
  formData
}


blastLocal <- function(genome_name="mygenome.fasta",temp_filename = "seq.aligned",task = "blastn", evalue = 0.001){
  command <- paste("blastn -subject=" , genome_name , " -query=", temp_filename, " -out=", temp_filename,".out -outfmt=6 -task=",task," -evalue=",evalue,sep="")
  cat(command,"\n")
  shell(command)
}
blastLocal(genome_name="C:/Documents and Settings/gbic/My Documents/CHR_I.fna",temp_filename="C:/Documents and Settings/gbic/My Documents/query.faa")

## read in files
allFiles <- list.files()

for(i in 1:length(allFiles)) {
  ### reading in data file
  ##################input <- as.matrix(read.table(allFiles[[i]],header=TRUE,row.names=NULL))
  input <- as.matrix(read.table(allFiles[[i]],header=TRUE,row.names=NULL))
  
  ############################3
  matrixFile <- cbind(input,rep(NA,nrow(input)))
  colnames(matrixFile) <- c(colnames(input),"unique")
  query <- NULL
  #rewriting it to semi FASTA
  for(j in 1:nrow(matrixFile)) {
    temp1 <- paste(">",matrixFile[j,1],sep="") #row name
    temp2 <- paste(temp1,matrixFile[j,3],sep="\n") # sequence
    query <- paste(query,temp2,sep="\n")
  }
  # submitting BLAST
  output <- getPosThaliana(query)
  # rewriting BLAST output
  output <- substr(output,23,nchar(output))
  output <- strsplit(output,"\n")
  blastFile <- strsplit(output[[1]],"\t"))
  # parsing output
  
  for(k in 1:nrow(matrixFile)) {
    bUnique <- 1 #full match (tested below)
    rowTemp = which(blastFile[,1] == matrixFile[k,1]) 
    if( length(rowTemp) != 0) {
      fullMatched <- which(blastFile[rowTemp,4] == "25" & blastFile[rowTemp,5] == "0" & blastFile[rowTemp,6] == "0") 
      if( length(unique(blastFile[rowTemp[fullMatched],9])) != 1 & length(unique(blastFile[rowTemp[fullMatched],10])) != 1) {
        bUnique <- 0 # non unique match
      }
    } else {
      bUnique <- -1 # incomplete or non-match
    }
    matrixFile[k,4] <- bUnique
  }
  outName <- paste("./output/c_",allFiles[[i]],sep="")
  write.table(matrixFile,file=outName,sep="\t")
}

# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Blasts nucluotide seqeunce vs nucleotide database
# on http://www.wormbase.org/db/searches/blast_blat
# To-Do: Retrieves first n sequences

query <- ">486029\nCTTCGTTTCCCTCTTCTGCGATTTC\n>486032\nGATTGCACCTTCGATGGCCCTGAAA\n>486033\nGCCACATAGTGACGGCGATGGAAAC\n>486042\nTCGCGCTTCACGAAACGCCTTTCGT\n>486049\nGCTTCCAGAGATAGCCCGACGGCGG\n>486050\nACTCAGCAGCAGTCGTCACGACGTG"


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
  formData <- postForm(uriToPost, .params = postValues,style="POST")
  formData
}

txt <- getPosThaliana()
write(txt,"outputT.txt")


## read in files
setwd("D:/CompMolBio/Centrotype/Sequentie_probes")
allFiles <- list.files()

#for(i in 1:length(allFiles)) {
  ### rewrite base files
  matrixFile <- as.matrix(read.table("sequentie_probesAT1G01010.txt",header=TRUE))
  query <- NULL
  for(j in 1:nrow(matrixFile)) {
    temp1 <- paste(">",rownames(matrixFile)[j],sep="")
    temp2 <- paste(temp1,matrixFile[j,2],sep="\n")
    query <- paste(query,temp2,sep="\n")
  }
  oneFile <- getPosThaliana(query)
  oneFile <- substr(oneFile,23,nchar(oneFile))
  tempB <- strsplit(oneFile,"\n")
  tempC <- strsplit(tempB[[1]],"\t")
  tempC
#}

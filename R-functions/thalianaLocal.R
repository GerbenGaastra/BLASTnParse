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

check <- function(geneName,geneInfo,blastFile,output) {
  blastRes <- as.matrix(read.table("blastFile.txt",sep="\t",header=FALSE))
  for(j in 1:nrow(output) ) {
    
    entries <- which(geneInfo[,"Synonym"] == output[j,"gene") # getting lines containing current gene
    geneStart <- as.numeric(min(geneInfo[entries,"start"])) - 400 #flanking correction
    geneEnd <- as.numeric(max(geneInfo[entries,"stop"])) + 400 #flanking correction
    # probes should be unique
    currID <-output[j,"row.names"]
    BlastsFromID <- which(blastRes[,1] == currID)
    if( length(BlastsFromID) == 1) {
      # probes should fall into gene range
      if( blastRes[BlastsFromID,9] > geneStart & blastRes[BlastsFromID,10] < geneEnd ) {
        output[j,"unique"] = "yes"
      } else {
        output[j,"unique"] = "outOfRange"
      }
    } else {
      output[j,"unique"] = "non unique"
    }
  }
  output
}

###############################################################
input <- "D:/CompMolBio/Centrotype/Sequentie_probes/combined.txt"
evalue <- 0.001
genome <- "D:/CompMolBio/BlastLocal/TAIR10_chr_all.fas"
geneFile <- "D:/CompMolBio/BlastLocal/NC_003070.ptt"
###########################################################

cat("reading input file:", input,"\n",sep=" ")
inputFile <- read.csv(input,sep="\t")
cat("rewriting input to fasta format\n")
## Two collumns needed, one with fasta ID, one containing the sequences
fastaRewrite(inputFile[,c("ID","seq")])

cat("starting BLASTs\n")
blastLocal(genome,"input2.fasta","blastFile.txt",evalue=evalue)
cat("all BLASTs performed\n")
## preparing matrix with information about start&stop bp of genes on chromosome
geneInfo <- geneRange(geneFile)
cat("extracted starting & ending bp of all genes")
output<- cbind(inputFile,rep(NA,nrow(inputFile)))
colnames(output) <- c("gene","row.names","direction","seq","unique")
res <- check(geneName,entries,geneStart,geneEnd,blastFile,output)
write.table(res,"output.txt",sep="\t")
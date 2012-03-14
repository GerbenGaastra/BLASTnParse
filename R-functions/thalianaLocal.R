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

check <- function(output,blastRes,geneInfo) {
  blastID <- matrix(NA,nrow(blastRes),2)
  colnames(blastID) <- c("gene","ID")
  # seperating genename and probe id
  for(i in 1:nrow(blastRes) ) {
    temp <- unlist(strsplit(blastRes[i,1],"-"))
    blastID[i,1] <- temp[1]
    blastID[i,2] <- temp[2]
  }
  for(j in 1:nrow(output) ) {
    if(j %% 1000 == 0) {
      cat("computing row", j, "\n", sep=" ")
    }
    # probes should be unique
    currID <-sub("^ +","",output[j,"row.names"]) #probe ID of current row
    #cat("currID",currID,"\n",sep=" ")
    BlastsFromID <- which(blastID[,"ID"] == currID) #list of times ID gave blast result
    #cat("BlastsFromID",BlastsFromID,"\n",sep=" ")
    if( length(BlastsFromID) == 1) {
      # probes should fall into gene range
      entries <- which(geneInfo[,"gene"] == output[j,"gene"]) # getting lines containing current gene
      #cat("entries",entries,"\n",sep=" ")
      if( length(entries) != 0 ) {
        geneStart <- as.numeric(min(geneInfo[entries,"start"])) - 2400 #flanking correction
        #cat("geneStart",geneStart,"\n",sep=" ")
        geneEnd <- as.numeric(max(geneInfo[entries,"stop"])) + 2400 #flanking correction
        #cat("geneEnd",geneEnd,"\n",sep=" ")
        if( as.numeric(blastRes[BlastsFromID,9]) > geneStart & as.numeric(blastRes[BlastsFromID,10]) < geneEnd ) {
          output[j,"unique"] = "y"
        } else {
          output[j,"unique"] = "y-outOfRange"
        }
      } else {
        output[j,"unique"] = "y-rangeUnknown"
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
inputFile <- as.matrix(read.csv(input,sep="\t"))
cat("rewriting input to fasta format\n")
## Two collumns needed, one with fasta ID, one containing the sequences
fastaID <<- apply(inputFile[,c("gene","row.names")],1,function(x) {
  paste(">",paste(x["gene"],sub("^ +","",x["row.names"]),sep="-"),sep="")
})
forFasta <- matrix( c(fastaID,inputFile[,"seq"]),,2)
fastaRewrite(forFasta)
cat("starting BLASTs\n")
blastLocal(genome,"output.fasta","blastFile.txt",evalue=evalue)
cat("all BLASTs performed\n")
## preparing matrix with information about start&stop bp of genes on chromosome
geneInfo <- geneRange(geneFile)
cat("extracted starting & ending bp of all genes\n")
## Prepping output matrix
output<- cbind(inputFile,rep(NA,nrow(inputFile)))
blastRes <- as.matrix(read.table("blastFile.txt",sep="\t",header=FALSE))
colnames(output) <- c("gene","row.names","direction","seq","unique")
res <- check(output,blastRes,geneInfo)
write.table(res,"output.txt",sep="\t")
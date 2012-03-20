# modified by GerbenGaastra
# created 20-02-2012
# modified last 22-02-2012
# Listing gene ranges from ncbi .ptt file

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

# Written by GT Gaastra (r)
# created last 12-03-2012
# modified last 12-03-2012
# rewrites input into fasta format

## function will rebuild its input (form 'matrix[,2]') returns 
## a vector which can be written to file as a valid (ncbi) fasta file
## input collumn 1 must contain the fasta ID
## input collumn 2 must contain the sequence
fastaRewrite<- function(input) {
  #if( class(input) != "matrix") { stop("'input' is not an matrix") }
  if( ncol(input) != 2 ) { stop("'input' has does not have 2 collumns") }
  output <- NULL
  fp <- file("output.fasta","w")
  for( i in 1:nrow(input) ) {
    #if( (i %% 1000) == 0 ) {
    #  cat("parsing row",i,"\n",sep=" ")
    #  write(output,"input.fasta",append=TRUE,sep="\t")
    #  output <- NULL
    #}
    writeLines(as.character(input[i,1]), con=fp)
    writeLines(as.character(input[i,2]), con=fp)
    
    
  }
  close(fp)
}
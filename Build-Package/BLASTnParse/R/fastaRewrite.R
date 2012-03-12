# Written by GT Gaastra (r)
# created last 12-03-2012
# modified last 12-03-2012
# rewrites input into fasta format

## function will rebuild its input (form 'matrix[,2]') returns 
## a vector which can be written to file as a valid (ncbi) fasta file
## input collumn 1 must contain the fasta description
## input collumn 2 must contain the sequence
fastaRewrite<- function(input) {
  if( class(input) != "matrix") { stop("'input' is not an matrix") }
  if( ncol(input) != 2) { stop("'input' has not 2 collumns") }
  output <- NULL
  for( i in 1:nrow(input) ) {
    output <- c(output,paste(">",input[i,1],sep=""),as.character(input[i,2]) )
  }
  output
}
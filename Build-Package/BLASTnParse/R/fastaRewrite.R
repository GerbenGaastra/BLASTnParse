# Written by GT Gaastra (r)
# created last 12-03-2012
# modified last 12-03-2012
# rewrites input into fasta format

## function will rebuild its input (form 'matrix[,2]') returns 
## a vector which can be written to file as a valid (ncbi) fasta file
## input collumn 1 must contain the fasta ID
## input collumn 2 must contain the sequence
fastaRewrite<- function(input,outFile,checkID=FALSE) {
  if( missing(input) ) { stop("argument 'input' is missing") }
  if( missing(outFile) ) { stop("argument 'outFile' is missing") }
  if( ncol(input) != 2 ) { stop("'input' has does not have 2 collumns") }
  if( mode(outFile) != "character" ) { stop("mode 'outFile' is not a char string") }
  fp <- file(outFile,"w")
  apply( input,1,function(x) {
    if(checkID) {
      if( substr(as.character(x[1]),1,1) != ">") {
        defline <- paste(">",as.character(x[1]),collapse="",sep="")
        writeLines(defline, con=fp)
      } else {
        writeLines(as.character(x[1]), con=fp)
      }
    } else {
      writeLines(as.character(x[1]), con=fp)
    }
    writeLines(as.character(x[2]), con=fp)
  })
  close(fp)
}

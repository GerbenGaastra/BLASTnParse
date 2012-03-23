\name{fastaRewrite}
\alias{fastaRewrite}
\docType{package}
\title{
  Rewrite a matrix into vector in fasta format
}
\description{
  The function \emph{fastaRewrite} converts its input into a fasta format file,
  The function expects a 2 column dataframe of matrix as input. Each row should 
  contain on column 1 the fasta identifier, and on column 2 the corresponding
  sequence. Some very basic error checking on the fasta defline is performed is checkID is true.
}
\usage{
  fastaRewrite(input,outFile,checkID=FALSE) 
}
\arguments{
  \item{input}{Matrix containing on column 1 the fasta identifiers and on collumn 2 the corresponding sequence.}
  \item{outFile}{String containing (path and) filename of the desired output file}
  \item{checkID}{A logical value, indicating whether the fasta identifier should be checked and rewritten}
}
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml}, 
\url{www.rug.nl/gbic}
}
\examples{
## filling example data
input <- c("WBGene00007709", "WBGene00017675", "tcgtttattatttgtcaccgggt", "ctccatcccacacacacacactt")
dim(input) <- c(2,2)
outFile <- "output.fasta"

## Formatting input and writing output to file: 'input.fasta' 
fastaRewrite(input,outFile)
}
\keyword{ manip }

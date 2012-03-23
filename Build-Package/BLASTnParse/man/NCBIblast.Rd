\name{NCBIblast}
\alias{NCBIblast}
\title{
  Run batch BLAST NCBI via qblast
}
\description{
  Performs batch BLAST on sequences provided in input matrix, using the ncbi qblast.
  NCBIblast returns matrix containing the BLAST results
}
\usage{
  NCBIblast(input,eValue,program="blastn",db="nt",verbose=FALSE)
}
\arguments{
  \item{input}{Matrix with 2 collumns containing in collumn 1 the IDs and in column 2 the sequences}
  \item{eValue}{cuttoff value to filter BLAST results}
  \item{program}{BLAST program to run, blastn is standard}
  \item{db}{string, specifying BLAST database}
  \item{verbose}{logical, states whether progres should be printed to terminal}
}
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html}, 
\url{www.rug.nl/gbic}
}
\examples{
  ## building input
  input <- matrix(NA,2,2)
  input <- c("ID1","ID2","CTTCGTTTCCCTCTTCTGCGATTTC","GATTGCACCTTCGATGGCCCTGAAA")
  dim(input) <- c(2,2)
  colnames(input) <- c("ID","query")
  eValue <- "1e-2"
  
  ## running request
  NCBIblast(input,eValue)
}
\seealso{
  \code{\link{NCBIsubmit}} - Submit single BLAST request to the NCBI servers, 
  \code{\link{NCBIcheck}} - Check status of request on NCBI server, 
  \code{\link{NCBIdownload}} - Download ready BLAST request
}
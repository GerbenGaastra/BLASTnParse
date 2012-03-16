\name{NCBIblast}
\alias{NCBIblast}
\title{
  Run batch BLAST NCBI
}
\description{
  Performs batch blast on sequences provided in input matrix, 
  NCBIblast returns matrix containing BLAST results
}
\usage{
  NCBIblast(input,eValue,program="blastn",db="nt")
}
\arguments{
  \item{input}{Matrix with 2 collumns containing in collumn 1 the IDs and in column 2 the sequences}
  \item{eValue}{cuttoff value to filter BLAST results}
  \item{program}{BLAST program to run, blastn is standard}
  \item{db}{string, specifying BLAST database}
}
\details{
\tabular{ll}{
Package: \tab firstPackage\cr
Type: \tab Package\cr
Version: \tab 0.0.1-0\cr
Date: \tab 2012-02-21\cr
License: \tab GPL (>= 3)\cr
}
}
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{www.rug.nl/gbic}
}
\keyword{ package }
\examples{
  input <- matrix(NA,2,2)
  input <- c("ID1","ID2","CTTCGTTTCCCTCTTCTGCGATTTC","GATTGCACCTTCGATGGCCCTGAAA")
  dim(input) <- c(2,2)
  colnames(input) <- c("ID","query")
  eValue <- "1e-2"
  NCBIblast(input,eValue)
}
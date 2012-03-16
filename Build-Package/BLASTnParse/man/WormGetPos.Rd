\name{wormGetPos}
\alias{wormGetPos}
\docType{package}
\title{
  Get position 
}
\description{
  Post BLAST on 'wormblast.org' and returns the best match. 
  Return contains chromosome, bp position, number of identities.
}
\usage{
 wormGetPos(sequence,eValue,db,handle) 
}
\arguments{
  \item{sequence}{string containing bp-sequence to BLAST}
  \item{eValue}{String containing minimal expectancy of result}
  \item{db}{database which wormblast should use}
  \item{handle}{Specifying curl handle}
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
## String containing sequence to allign
sequenceOfInterest <- "tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"    

## Sending sequence to wormblast, returning relevant data into list.
list <- wormGetPos(sequenceOfInterest)
}
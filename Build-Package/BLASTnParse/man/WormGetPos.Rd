\name{wormGetPos}
\alias{wormGetPos}
\title{
  Perform BLAST request to wormbase.org
}
\description{
  Post BLAST on 'wormblast.org' and returns the best match. 
  Return contains chromosome, bp position, number of identities.
}
\usage{
 wormGetPos(query,eValue,db,handle) 
}
\arguments{
  \item{query}{string containing  query sequence to be BLASTed.}
  \item{eValue}{String containing the minimal expectancy value, used to filter BLAST results.}
  \item{db}{catabase against which should be BLASTed.}
  \item{handle}{the curl handle to be used, if no handle is provided, one is created at runtime.}
}   
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{http://www.wormbase.org/db/searches/blast_blat}, 
\url{www.rug.nl/gbic}
}

\examples{
## String containing sequence to allign
sequenceOfInterest <- "tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"    

## Sending sequence to wormblast, returning relevant data into list.
list <- wormGetPos(sequenceOfInterest)
}

\seealso{
  \code{\link{wormParse}} - Parse chromosome, bp matched and allignmet lenght from wormBLAST output.R
  \code{\link{wormDownload}} - Post and download one BLAST request from wormbase.org.
}


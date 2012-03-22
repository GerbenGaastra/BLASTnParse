\name{wormGetPos}
\alias{wormGetPos}
\title{
  Send BLAST request to wormbase.org
}
\description{
  Post BLAST on 'wormblast.org' and returns the best match. 
  Return contains chromosome, bp position, number of identities.
}
\usage{
 wormGetPos(query,eValue,db,handle) 
}
\arguments{
  \item{query}{string containing  query sequence to be BLASTed}
  \item{eValue}{String containing the minimal expectancy value to  result}
  \item{db}{Database against which should be BLASTed}
  \item{handle}{The curl handle to bu used, if no handle is provided, one is created at runtime}
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

\examples{
## String containing sequence to allign
sequenceOfInterest <- "tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"    

## Sending sequence to wormblast, returning relevant data into list.
list <- wormGetPos(sequenceOfInterest)
}

\seealso{
  \code{\link{wormParse}} - Parse chromosome, bp matched and allignmet lenght from wormBLAST output.
  \code{\link{wormDownload}} - Post and download one BLAST request from wormbase.org.
}


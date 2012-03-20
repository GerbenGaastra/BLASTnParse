\name{NCBIcheck}
\alias{NCBIcheck}
\title{
  Check status of request on NCBI server
}
\description{
  Checks the status of an earlyer sended BLAST request on the NCBI server using 
  the RID as the identifier.
}
\usage{
  NCBIcheck(reqID) 
}
\arguments{
  \item{reqID}{String containing the RID of which the status will be returned}
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
  vectorResult <- NCBIsubmit("GATTGCACCTTCGATGGCCCTGAAA","0.01")
  #Sys.sleep(vectorResult$RTOE)
  NCBIcheck(vectorResult$RID)
}
\seealso{
  \code{\link{NCBIblast}} - Run batch BLAST NCBI
  \code{\link{NCBIsubmit}} - Submit single BLAST request to the NCBI servers
  \code{\link{NCBIdownload}} - Download ready BLAST request
}

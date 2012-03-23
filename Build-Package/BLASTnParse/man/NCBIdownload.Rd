\name{NCBIdownload}
\alias{NCBIdownload}
\title{
  Download ready BLAST request
}
\description{
  Downloads completed BLAST result, and returns it in matrix format.
}
\usage{
  NCBIdownload( reqID ) 
}
\arguments{
  \item{reqID}{string containing the RID value of the BLAST request}
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
    if( NCBIcheck(vectorResult$RID) == 1 ) {
      NCBIdownload(vectorResult$RTOE)
    }
}
\seealso{
  \code{\link{NCBIblast}} - Run batch BLAST NCBI
  \code{\link{NCBIsubmit}} - Submit single BLAST request to the NCBI servers
  \code{\link{NCBIcheck}} - Check status of request on NCBI server
}
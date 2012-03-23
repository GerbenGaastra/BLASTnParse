\name{NCBIcheck}
\alias{NCBIcheck}
\title{
  Check status of request on NCBI server
}
\description{
  Checks the status of an earlier sended BLAST request on the NCBI server using 
  the RID as the identifier.
}
\usage{
  NCBIcheck(reqID) 
}
\arguments{
  \item{reqID}{String containing the RID of which the status will be returned}
}
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
value{
  A '0' is returned when the status of the RID is unfinished, a '2' is returned for a finished job which has no results, a '3' is returned for a failure on the NCBI server (please contact blast-help [at] ncbi.nlm.nih.gov, a returned '4' indicates an unknown RID, this means a wrom or an expired RID. 
  
  Finally an returned 1 indicates an finished job with results.
\references{
\url{http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html}, 
\url{www.rug.nl/gbic}
}
\examples{
  vectorResult <- NCBIsubmit("GATTGCACCTTCGATGGCCCTGAAA","0.01")
  #Sys.sleep(vectorResult$RTOE)
  NCBIcheck(vectorResult$RID)
}
\seealso{
  \code{\link{NCBIblast}} - Run batch BLAST NCBI, 
  \code{\link{NCBIsubmit}} - Submit single BLAST request to the NCBI servers, 
  \code{\link{NCBIdownload}} - Download ready BLAST request
}

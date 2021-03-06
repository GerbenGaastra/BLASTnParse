\name{NCBIsubmit}
\alias{NCBIsubmit}
\title{
  Submit single BLAST request to the NCBI servers
}
\description{
  NCBIsubmit stores a single BLAST request on the NCBI servers, and returns a vector with the RID (Request ID) and RTOE (request time of execution)
}
\usage{
  NCBIsubmit(query,eValue,program="blastn",db="nt") 
}
\arguments{
  \item{query}{string containing the query sequence}
  \item{eValue}{cuttoff value to filter BLAST results}
  \item{program}{BLAST program to be run by the NCBI servers}
  \item{db}{string, specifying BLAST database}
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
  ### Send a BLAST request and retrieve the RID and RTOe
  vectorResult <- NCBIsubmit("GATTGCACCTTCGATGGCCCTGAAA","0.01","blastn","nt")
}
\seealso{
  \code{\link{NCBIblast}} - Run batch BLAST NCBI, 
  \code{\link{NCBIcheck}} - Check status of request on NCBI server, 
  \code{\link{NCBIdownload}} - Download ready BLAST request
}
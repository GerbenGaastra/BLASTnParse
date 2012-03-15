\name{NCBIsubmit}
\alias{NCBIsubmit}
\title{
  Submit BLAST request to the NCBI servers
}
\description{
  Store an BLAST request on the NCBI servers, and returns a vector with the RID (Request ID) and RTOE (request time of execution)
}
\usage{
  NCBIsubmit(query,eValue) 
}
\arguments{
  \item{query}{string containing the query sequence}
  \item{eValue}{cuttoff value to filter BLAST results}
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
  ### Send a BLAST request and retrieve the RID and RTOe
  vectorResult <- NCBIsubmit("GATTGCACCTTCGATGGCCCTGAAA","0.01")
}
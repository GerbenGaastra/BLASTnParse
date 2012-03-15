\name{NCBIdownload}
\alias{NCBIdownload}
\title{
  Download BLAST request
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
    vectorResult <- NCBIsubmit("GATTGCACCTTCGATGGCCCTGAAA","0.01")
    #Sys.sleep(vectorResult$RTOE)
    if( NCBIcheck(vectorResult$RID) == 1 ) {
      NCBIdownload(vectorResult$RTOE)
    }
}
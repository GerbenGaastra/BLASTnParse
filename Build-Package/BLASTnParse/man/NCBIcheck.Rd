\name{NCBIcheck}
\alias{NCBIcheck}
\title{
  Check status request
}
\description{
  Checks the status of an earlyer sended BLAST request.
}
\usage{
  NCBIcheck(reqID) 
}
\arguments{
  \item{reqID}{String containing the RID of which the status should be reported}
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
  NCBIcheck(vectorResult$RID)
}

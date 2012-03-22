\name{BLASTthaliana}
\alias{BLASTthaliana}
\title{
  Send a blast request to wormbase en parse out position
}
\description{
  Post BLAST on 'wormblast.org' and returns the best match. 
  Return contains chromosome, bp position, number of identities.
}
\usage{
  BLASTthaliana(ID,query,eValue="1E+0",db="ATH1_bacs_con")
}
\arguments{
  \item{ID}{string, query id}
  \item{query}{String, containing query sequence}
  \item{eValue}{String, containing the minimal expectancy value to  result, in the form "1E+0"}
  \item{db}{Database against which should be BLASTed}
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
## sample data
ID <- "query1"
query <- "GACCCGAGAAAATCCAAGACCTATG"

## Running blast
BLASTthaliana(ID,query)
}



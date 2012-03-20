\name{wormParse}
\alias{wormParse}
\title{
  Parse chromosome, bp matched and allignmet lenght from wormBLAST output
}
\description{
  Parses relevant information from returned by wormDownload. Returns some of
  the information in a named list. The names of each entry is the name of the
  corresponding chromosome. At this stage not all available information in the 
  blast result is parsed.
}
\usage{
  wormParse(BLASTresult) 
}
\arguments{
  \item{BLASTresult}{HTML data from downloadForm, containing output
    from wormblast form.}
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
## List of fields and values
query <-"tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"
eValue <- "1E-16"
db <- "elegans_genome"
handle <- getCurlHandle()

## Post form and retrieve result
BLASTresult <- wormDownload(query, eValue, db, handle)

## Post form and retrieve result
list <- wormParse(BLASTresult)
}
\seealso{
  \code{\link{wormGetPos}} - Get data of the first BLAST result from wormbase.org.
  \code{\link{wormDownload}} - Post and download one BLAST request from wormbase.org.
}
\keyword{ package }
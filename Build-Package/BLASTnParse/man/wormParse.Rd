\name{wormParse}
\alias{wormParse}
\title{
  Parse chromosome, bp matched and allignmet lenght from wormBLAST output
}
\description{
  Parses relevant information from BLAST result
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
## Much simplified return from downloadForm
BLASTresult <- "<HTML><BODY><H1>CHROMOSOME_I</H1><P>Sbjct:&nbsp;10184580.</P><BODY></HTML>"    

## Post form and retrieve result
list <- wormParse(BLASTresult)
}
\seealso{
  \code{\link{wormGetPos}} - Get data of the first BLAST result from wormbase.org.
  \code{\link{wormDownload}} - Post and download one BLAST request from wormbase.org.
}
\keyword{ package }
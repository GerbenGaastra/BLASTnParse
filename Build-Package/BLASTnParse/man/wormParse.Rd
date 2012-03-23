\name{wormParse}
\alias{wormParse}
\title{
  Parse chromosome, bp matched and allignment lenght from wormBLAST output
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
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{http://www.wormbase.org/db/searches/blast_blat}, 
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
  \code{\link{wormGetPos}} - Get data of the first BLAST result from wormbase.org, 
  \code{\link{wormDownload}} - Post and download one BLAST request from wormbase.org.
}
\keyword{ package }
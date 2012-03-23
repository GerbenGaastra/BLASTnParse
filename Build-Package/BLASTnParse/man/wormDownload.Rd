\name{wormDownload}
\alias{wormDownload}
\title{
  Post and download one BLAST request from wormbase.org
}
\description{
  The function \emph{wormDownlaod} posts the form on arabidopsis.org and 
  returns the downloaded BLAST result as a string.
}
\usage{
  wormDownload(query, eValue="1E+0", db="elegans_genome", handle=getCurlHandle()) 
}
\arguments{
  \item{query}{the query sequence, as a character string}
  \item{eValue}{expectancy cutoff value to filter BLAST results}
  \item{db}{database to BLAST against, default is the C. elegans genome}
  \item{handle}{the curl handle to bu used, if no handle is provided, one is created at runtime}
}    
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{http://www.wormbase.org/db/searches/blast_blat}
\url{www.rug.nl/gbic}
}
\examples{
  ## List of fields and values
  query <-"tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"
  eValue <- "1E-16"
  db <- "elegans_genome"
  handle <- getCurlHandle()
  
  ## Post form and retrieve result
  htmlReturn <- wormDownload(query, eValue, db, handle)
}
\seealso{
  \code{\link{wormGetPos}} - Get data of the first BLAST result from wormbase.org.
  \code{\link{wormParse}} - Parse chromosome, bp matched and allignmet lenght from wormBLAST output.
}
\keyword{ package }

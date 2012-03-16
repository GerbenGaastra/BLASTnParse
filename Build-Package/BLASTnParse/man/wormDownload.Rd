\name{wormDownload}
\alias{wormDownload}
\docType{package}
\title{
  download Form
}
\description{
  Wrapper with headercheck around htmlForm()
}
\usage{
  wormDownload(sequence, eValue="1E+0", db="elegans_genome", handle= getCurlHandle()) 
}
\arguments{
  \item{sequence}{string containing query sequence}
  \item{eValue}{expectancy cutoff value to filter BLAST results}
  \item{db}{database to BLAST against, default is the C. elegans genome}
  \item{handle}{curl handle}
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
  ## List of fields and values
  sequence <-"tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"
  eValue <- "1E-16"
  db <- "elegans_genome"
  handle <- getCurlHandle()
  
  ## Post form and retrieve result
  htmlReturn <- wormDownload(sequence, eValue, db, handle)
}

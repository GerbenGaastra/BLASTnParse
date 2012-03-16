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
  wormDownload(uri, postValues, handle) 
}
\arguments{
  \item{uri}{string containing url of form.}
  \item{postValues}{names list with imput fields of the form 
    and their values}
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
postValues <- new("list",
    query_sequence="tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact",
    query_type="nucl",
    blast_app="blastn",
    db="nucl",
    blast_e_value="1E-16",
    database="elegans_genome",
    search_type="blast",
    submit="Submit")
    
## URL of webform
uriToPost <- "http://www.wormbase.org/db/searches/blast_blat"    

## Post form and retrieve result
handle = getCurlHandle()
htmlReturn <- wormDownload(uriToPost, postValues, handle)
}

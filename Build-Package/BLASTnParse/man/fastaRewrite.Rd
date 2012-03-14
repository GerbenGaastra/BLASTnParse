\name{fastaRewrite}
\alias{fastaRewrite}
\docType{package}
\title{
  Rewrite a matrix into vector in fasta format
}
\description{
  rewrites 2 collumn matrix into a vector which can be written to
  file in (NCBI) fasta format
}
\usage{
  fastaRewrite(input) 
}
\arguments{
  \item{input}{Matrix containing per row on collumn 1 the fasta 
  identifiers and on collumn 2 the corresponding sequence.}
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
## filling example data
input <- c("WBGene00007709", "WBGene00017675", "tcgtttattatttgtcaccgggt", "ctccatcccacacacacacactt")
dim(input) <- c(2,2)

## Formatting input and writing output to file: 'input.fasta' 
fastaRewrite(input)
}

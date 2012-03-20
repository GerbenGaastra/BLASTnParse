\name{geneRange}
\alias{geneRange}
\title{
  Parse raw gene list and return gene synonym, start&end gene
}
\description{
  Using a NCBI *.ptt or *.rnt file as input, 
}
\usage{
  geneRange( geneFile )
}
\arguments{
  \item{geneFile}{text file in NCBI .ptt or .rnt format}
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
## Example .ptt location
pttUri <- "ftp://ftp.ncbi.nih.gov/genomes/Caenorhabditis_elegans/CHR_I/NC_003279.ptt"
##
matrixGenes <- geneRange( pttUri )
}

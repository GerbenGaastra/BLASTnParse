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

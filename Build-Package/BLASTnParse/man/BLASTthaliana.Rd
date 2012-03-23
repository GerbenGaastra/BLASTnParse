\name{BLASTthaliana}
\alias{BLASTthaliana}
\title{
  Send a blast request to Arabidopsis.org
}
\description{
  Post an BLAST request on 'arabidopsis.org', downloads the results, and parses the downloaded data into a matrix.
}
\usage{
  BLASTthaliana(ID, query, eValue="1E+0", db="ATH1_bacs_con")
}
\arguments{
  \item{ID}{string, containing query id}
  \item{query}{String, containing query sequence}
  \item{eValue}{String, containing the minimal expectancy value which is used to filter the BLAST hits. Must provided in this format: "1E+0"}
  \item{db}{Database against which should be BLASTed}
}   
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{http://www.arabidopsis.org/Blast/},
\url{www.rug.nl/gbic}
}

\examples{
## sample data
ID <- "query1"
query <- "GACCCGAGAAAATCCAAGACCTATG"

## Running blast
BLASTthaliana(ID,query)
}



\name{BLASTnParse-package}
\alias{BLASTnParse-package}
\alias{BLASTnParse}
\docType{package}
\title{
  Collection of functions to automate BLAST requests
}
\description{
  This package aimes to implement some functions to make specific BLAST requests more
  scalable. Most functions take an matrix with probe identifiers and sequences as input
  and return the BLAST results as a matrix.
}
\details{
\tabular{ll}{
Package: \tab BLASTnParse\cr
Type: \tab Package\cr
Version: \tab 0.0.1-0\cr
Date: \tab 2012-02-21\cr
Authors@R: \tab c(person("Gerben", "Gaastra", role=c("aut", "cre"), email =
"gaastra@gtga.nl"))\cr
Depends: \tab R (>= 1.14.1), xtable, RCurl\cr
License: \tab GPL (>= 3)\cr
URL: \tab http://www.r-project.org, http://www.rug.nl/gbic\cr
Packaged: \tab 2012-03-06 11:34:32 UTC; gbic\cr
Built: \tab R 2.14.1; ; 2012-03-23 12:16:26 UTC; windows\cr
}

Index:
\preformatted{
Webblast:
BLASTthaliana           Send a blast request to Arabidopsis.org and
                        recieve a named vector
wormGetPos              Send BLAST request to wormbase.org
-wormDownload           Post and download one BLAST request from
                        wormbase.org
-wormParse              Parse chromosome, bp matched and allignmet
                        lenght from wormBLAST output

Webblast using qblast:
NCBIblast               Run batch BLAST NCBI
-NCBIcheck              Check status of request on NCBI server
-NCBIdownload           Download ready BLAST request
-NCBIsubmit             Submit single BLAST request to the NCBI servers

BLAST on local system:
BLASTlocal              Run BLAST locally using NCBI blast+ executable

Miscellaneous:
fastaRewrite            Rewrite a matrix into vector in fasta format
geneRange               Parse raw gene list and return gene synonym,
                        start&end gene
}
}
\author{
Author: GT Gaastra
Maintainer: GT Gaastra <gaastra@gtga.nl>

}
\keyword{ package }



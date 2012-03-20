\name{BLASTlocal}
\alias{BLASTlocal}
\title{
  Run BLAST locally using NCBI blast+ executable
}
\description{
  BLASTs input fasta file against a subject (for example a genome file), other blast options still have to be implemented
}
\usage{
  BLASTlocal(subject.fasta,input.fasta,outputPath,task = "blastn", evalue = 0.0001)
}
\arguments{
  \item{subject.fasta}{Provide a fasta file with which to serve as subject, for example a genome file}
  \item{input.fasta}{Provide a fasta file with contains one or more queries}
  \item{outputPath}{Provide a string containing a valid output path}
  \item{task}{specify which BLAST program to use, defaults to blastn}
  \item{evalue}{specify the cutoff value for filtering BLAST results}
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
  ## some example data
  subject <- c(">subject1","TCGAGATACCGTATGAAGAGAGCGAGATAGGAGATCAGAT")
  eValue <- "0.1"
  input <- c(">probe1","AGATACCGTATGAAGGAGC",">probe2",
    "CGTATGAAGAGAGCGAGATAGGAGATCA",">probe3","GATACCGTATGAAGAGACCGAGATAGGAGA")
  ## BLASTlocal only accepts files from drive
  write(subject,"subject.fasta")
  write(input,"input.fasta")
  ## running blasts
  BLASTlocal("subject.fasta","input.fasta",output="output.txt",evalue=eValue)
}
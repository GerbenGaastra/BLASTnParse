\name{BLASTlocal}
\alias{BLASTlocal}
\title{
  Run BLAST locally using NCBI blast+ executable
}
\description{
  BLASTlocal call the 'blastn' executable, which is part of the blast+ install, from within the R environment. The function takes an input file (with for example probes and subject file (for example a genome file). Both have to be valid fasta files conform the NCBI guidelines. Other types of BLAST requests still have to be implemented.
}
\usage{
  BLASTlocal(subject.fasta, input.fasta, outputPath, task = "blastn", evalue = 0.0001)
}
\arguments{
  \item{subject.fasta}{provide a fasta file with which to serve as subject, for example a genome file}
  \item{input.fasta}{provide a fasta file with contains one or more queries}
  \item{outputPath}{provide a string containing a valid output path}
  \item{task}{specify which BLAST program to use, defaults to blastn}
  \item{evalue}{specify the minimal expectancy value used as a cutoff value for filtering BLAST results}
}
\author{
Author: Gerben Gaastra
Maintainer: Gerben Gaastra <gaastra@gtga.nl>
}
\references{
\url{www.rug.nl/gbic}\\
\url{http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download}
}
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
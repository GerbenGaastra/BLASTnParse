# Adapted from Danny Arends
# modified by GerbenGaastra
# adopted 16-02-2012
# modified last 16-02-2012
# Wrapper around shell('blastn [args]')

## function will use the installed blastn program to BLAST the input along the provided genome
BLASTlocal <- function(subject.fasta,input.fasta,outputPath,task = "blastn", evalue = 0.0001){
  if( missing(subject.fasta) ) stop ("Please provide a valid genome filename")
  if( missing(input.fasta) ) stop ("Please provide a valid input filename")
  if( missing(outputPath) ) stop ("Please provide a valid output filename")
  if( !file.exists(subject.fasta) ) stop("Genome file not found")
  if( !file.exists(input.fasta) ) stop("fasta input file not found")
  installed <- system("blastn -help",ignore.stderr=TRUE,ignore.stdout=TRUE,show.output.on.console=FALSE)
  if( installed != 0 ) { stop("Please install the NCBI blast+ program") }
  command <- paste("blastn ",
    " -subject=" , subject.fasta , ## see '("blastn -help")' for all options
    " -query=", input.fasta,
    " -out=", outputPath,
    " -outfmt=\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"",
    " -task=",task,
    " -evalue=",evalue,
    sep="")
  shell(command)
}

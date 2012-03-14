# Written by GT Gaastra (r)
# created last 24-02-2012
# modified last 24-02-2012
# Wrapper around shell('blastn [args]')

## function will use the installed blastn program to BLAST the input along the provided genome
blastLocal <- function(genome_name,fasta_input,outputPath,task = "blastn", evalue = 0.0001){
  if( missing(genome_name) ) stop ("Please provide a valid genome filename")
  if( missing(fasta_input) ) stop ("Please provide a valid input filename")
  if( !file.exists(genome_name) ) stop("Genome file not found")
  if( !file.exists(fasta_input) ) stop("fasta input file not found")
  command <- paste("blastn ",
    " -subject=" , genome_name , ## see '("blastn -help")' for all options
    " -query=", fasta_input,
    " -out=", outputPath,
    " -outfmt=\"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue\"",
    " -task=",task,
    " -evalue=",evalue,
    sep="")
  shell(command)
}

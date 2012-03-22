# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Blasts nucluotide seqeunce vs nucleotide database
# on http://www.wormbase.org/db/searches/blast_blat
# To-Do: Retrieves first n sequences

BLASTthaliana <- function(ID,query,eValue="1E+0",db="ATH1_bacs_con") {
  temp1 <- paste(">",ID,sep="") #row name
  fasta <- paste(temp1,query,sep="\n")
  # setting up parameters
  uriToPost <- "http://www.arabidopsis.org/cgi-bin/Blast/TAIRblast.pl"
  # names list containing all fields and their values
  postValues <- new("list",
    Algorithm="blastn",
    BlastTargetSet=db,
    textbox="seq",
    QueryText=fasta,
    Matrix="Blosum62",
    ReplyVia="BROWSER",
    Expectation=eValue,
    ReplyFormat="TABULAR")
  #submit="Run BLAST")
  # Post and download Form
  Sys.sleep(runif(1))
  formData <- postForm(uriToPost, .params = postValues,style="POST")
  
  # rebuild output data to table
  comments <- c("Query id", "Subject id", "identity", "alignment length", "mismatches", "gap openings", "q. start", "q. end", "s. start", "s. end", "e-value", "bit score")
  formData <- substr(formData,23,nchar(formData)-3)
  formData <- strsplit(formData,"\n")
  blastFile <- unlist(strsplit(formData[[1]],"\t"))
  names(blastFile) <- comments
  blastFile
}

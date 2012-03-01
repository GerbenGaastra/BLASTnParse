# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Retrieve Position and Chromosome from post-data

parseBLAST <- function(BLASTresult) {
  out <- vector("list",3)
  resChr <- regexpr("CHROMOSOME_([IVX]+)?",BLASTresult) #Chromosome
  out[[1]] <- as.character(substr(BLASTresult, resChr+11, resChr + attr(resChr, "match.length") -1))
  resPos <- regexpr("(Sbjct+)(.{7})([0123456789]{2,})",BLASTresult) #bp position
  out[[2]] <- as.character(substr(BLASTresult, resPos+12, resPos + attr(resPos, "match.length") -1))
  iFullMatch  <- gregexpr("(Identities&nbsp;=&nbsp;60/60+).+?",BLASTresult)[[1]]
  if( length(iFullMatch) > 1) {
    out[[2]] == -2
  }
  resLen <- regexpr("(Identities&nbsp;=&nbsp;)+([0-9]{1,})",BLASTresult) #length matched 
  out[[3]] <- as.character(substr(BLASTresult, resLen + 23, resLen + attr(resLen, "match.length")- 1))
  
  names(out) <- c("chr","pos","nrOfMatches")
  out
}
parseBLAST(txt)
#####################################
#####  meuk #########################

uriToPost <- "http://www.wormbase.org/db/searches/blast_blat"
q_seq = "tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact"
postValues <- new("list",
    query_sequence=q_seq,
    query_type="nucl",
    blast_app="blastn",
    db="nucl",
    blast_e_value="1E+0",
    database="elegans_genome",
    search_type="blast",
    submit="Submit")



txt <- downloadForm(uriToPost,postValues) 

parseBLAST(txt)
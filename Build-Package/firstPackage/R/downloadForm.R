# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Wrapper for Wormblast

# adding some control to postForm
downloadForm <- function(uri, postValues) {
  if(url.exists(uri)) {
    HTMLreturn = postForm(uri, .params = postValues)
  } else {
    error <- cat("url: ",uri,"not available.\n",sep=" ")
    stop(error)
  }
  HTMLreturn
}

# Retrieving position, match length en chromosome from the best hit, 
# '""' returned if no results present
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

# Performing 1 blast and returns parsed results
getPosition <- function(q_seq,eValue="1E+0",daba="elegans_genome") {
  # setting up parameters
  uriToPost <- "http://www.wormbase.org/db/searches/blast_blat"
  # names list containing all fields and their values
  postValues <- new("list",
    query_sequence=q_seq,
    query_type="nucl",
    blast_app="blastn",
    db="nucl",
    blast_e_value=eValue,
    database=daba,
    search_type="blast",
    submit="Submit")
  # Post and download Form
  formData <-downloadForm(uriToPost,postValues)
  # Parsing Post-data
  result <- parseBLAST(formData)
  result
}

###

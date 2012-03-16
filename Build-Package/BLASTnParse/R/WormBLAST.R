# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Wrapper for Wormblast


# adding some control to postForm
wormDownload <- function(uri, postValues, handle) {
  #if(url.exists(uri)) {
    HTMLreturn = postForm(uri, .params = postValues, curl = handle)
  #} else {
  #  error <- cat("url: ",uri,"not available.\n",sep=" ")
  #  stop(error)
  #}
  HTMLreturn
}

# Retrieving position, match length en chromosome from the best hit, 
# '""' returned if no results present
wormParse <- function(BLASTresult) {
  out <- vector("list",3)
  resChr <- regexpr("CHROMOSOME_([IVX]+)?",BLASTresult) #Chromosome
  out[[1]] <- as.character(substr(BLASTresult, resChr+11, resChr + attr(resChr, "match.length") -1))
  resPos <- regexpr("(Sbjct+)(.{7})([0123456789]{2,})",BLASTresult) #bp position
  out[[2]] <- as.character(substr(BLASTresult, resPos+12, resPos + attr(resPos, "match.length") -1))
  iFullMatch  <- gregexpr("(Identities&nbsp;=&nbsp;60/60+).+?",BLASTresult)[[1]]
  if( length(iFullMatch) > 1) {
    out[[2]] == -2 ## error code, more than 1 full match
  }
  resLen <- regexpr("(Identities&nbsp;=&nbsp;)+([0-9]{1,})",BLASTresult) #length matched 
  out[[3]] <- as.character(substr(BLASTresult, resLen + 23, resLen + attr(resLen, "match.length")- 1))
  
  names(out) <- c("chr","pos","nrOfMatches")
  out
}

# Performing 1 blast and returns parsed results
wormGetPos <- function(sequence, eValue="1E+0",db="elegans_genome", handle = getCurlHandle()) {
  # setting up parameters
  if(missing(sequence)) stop("No sequence to query for, please provide a sequence")
  if("RCurl" %in% rownames(installed.packages())){
    require("RCurl")
  }else{
    stop("Please install the RCurl package (install.packages(\"RCurl\")")
  }
  uriToPost <- "http://www.wormbase.org/db/searches/blast_blat"
  # names list containing all fields and their values
  postValues <- new("list",
    query_sequence=sequence,
    query_type="nucl",
    blast_app="blastn",
    db="nucl",
    blast_e_value=eValue,
    database=db,
    search_type="blast",
    submit="Submit")
  # Post and download Form
  formData <-wormDownload(uriToPost,postValues,handle = handle)
  # Parsing Post-data
  result <- wormParse(formData)
  result
}

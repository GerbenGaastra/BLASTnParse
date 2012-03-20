# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Wrapper for Wormblast

# adding some control to postForm
wormDownload <- function(query, eValue="1E+0", db="elegans_genome", handle = getCurlHandle()) {
  if( missing(query) ) { stop("Please provide a query sequence") }
  if( mode(eValue) != "character" ) { stop("'eValue' must be a character string") }
  
  #if( mode(db) != "character" ) { stop("'db' must be a character string") }
  uri <- "http://www.wormbase.org/db/searches/blast_blat"
  # names list containing all fields and their values
  postValues <- new("list",
    query_sequence=query,
    query_type="nucl",
    blast_app="blastn",
    db="nucl",
    blast_e_value=eValue,
    database=db,
    search_type="blast",
    submit="Submit")
  if(url.exists(uri)) {
    HTMLreturn = postForm(uri, .params = postValues, curl = handle)
  } else {
    error <- cat("url: ",uri,"not available.\n",sep=" ")
    stop(error)
  }
  if( regexpr("(An error occurred:)+",HTMLreturn) != -1) {
    res <- gregexpr("(<div class=\"spacer\">\n    &nbsp;\n</div>    \n\n)+",HTMLreturn)
    error <- substr(HTMLreturn,(res[[1]][1]+ attr(res[[1]],"match.length")[1]),res[[1]][2]-7)
    stop(paste("wormbase error:",error,sep=" "))
  }
  HTMLreturn
}

# Retrieving position, match length en chromosome from BLAST result, 
# If there are no hits, an empty list is returned
wormParse <- function(BLASTresult) {
  ## seperating results per chromosome and discarding header
  report_hit <- strsplit(BLASTresult,"(report_hit)+")[[1]][-1]
  namesChr <-unlist(lapply(report_hit,function(x) { 
    resChr <-regexpr("CHROMOSOME_([IVX]+)?",x)
    as.character(substr(x, resChr+11, resChr + attr(resChr, "match.length") -1))
  }))
  names(report_hit) <- namesChr
  ## seperating each hit per chromosomere
  hits <- lapply(report_hit,function(x) {strsplit(x,"(Score)+")[[1]][-1] })
  lapply(hits,function(x) {
    do.call( rbind, lapply(x,
      function(y) {
        out <- vector("list",2)
        resPos <- regexpr("(Sbjct+)(.{7})([0123456789]{2,})",y) #bp position
        out[1] <- as.character(substr(y, resPos+12, resPos + attr(resPos, "match.length") -1))
        resLen <- regexpr("(Identities&nbsp;=&nbsp;)+([0-9]{1,})",y) #length matched 
        out[[2]] <- as.character(substr(y, resLen + 23, resLen + attr(resLen, "match.length")- 1))
        names(out) <- c("bp-mapped","length matched")
        out
      })
    )
  })
}

# Performing 1 blast and returns parsed results of the best hit
wormGetPos <- function(query, eValue="1E+0",db="elegans_genome", handle = getCurlHandle()) {
  # setting up parameters
  if(missing(sequence)) stop("No sequence to query for, please provide a sequence")
  if("RCurl" %in% rownames(installed.packages())){
    require("RCurl")
  }else{
    stop("Please install the RCurl package (install.packages(\"RCurl\")")
  }
  # Post and download Form
  formData <-wormDownload(query,eValue,db,handle = handle)
  # Parsing Post-data
  result <- wormParse(formData)
  result
}

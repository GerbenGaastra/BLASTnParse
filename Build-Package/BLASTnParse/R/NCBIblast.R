# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# NCBI blast, Qblast

# send one query and return the RID and RTOE
NCBIsubmit <- function(query,eValue,program="blastn",db="nt") {
  # building and sending request
  h = basicTextGatherer()
  inputFields <- paste("CMD=Put&PROGRAM=",program,"&DATABASE=",db,"&EXPECT=",eValue,"&QUERY=",query,collapse="",sep="")
  #cat(inputFields,"\n")
  curlPerform(url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
    httpheader=c(Accept="text/xml", Accept="multipart/*", 'Content-Type' = "application/x-www-form-urlencoded"),
    postfields=inputFields,
    writefunction = h$update,
    verbose = TRUE
  )
  #parse h$value for RID (request ID) and RTOE (Request time of Execution)
  out <- vector("list",2)
  res1 <- regexpr( "(RID%3D+)([A-Za-z0-9-]{2,})",h$value() )
  res2 <- regexpr( "(RTOE%3D+)([0-9]{1,})",h$value() )
  out[[1]] <- substr(h$value(),res1+6, res1 + attr(res1, "match.length") -1 )
  out[[2]] <- substr(h$value(),res2+7, res2 + attr(res2, "match.length") -1 )
  names(out) <- c("RID","RTOE")
  out
}

#checking status of RID
NCBIcheck <- function( reqID ) {
  tempres <- getURL(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=",reqID,sep=""))
  res <- regexpr( "(Status=+)([A-Za-z]{1,})",tempres )
  status <- substr(tempres,res+7,res+attr(res,"match.length")-1)
  #cat("status",status,"\n")
  out <- 0
  if( status == "READY") {
    if ( regexpr("(ThereAreHits=yes+)",tempres ) != -1 ) {
    out <- 1 #BLAST with results
    names(out) <- "READY+hits"
    } else { 
      out <- 2 #BLAST without results
      names(out) <- "READY-0Hits"
    } 
  }
  if( status == "FAILED") {
    cat("Search", reqID, "failed; please report to blast-help [at] ncbi.nlm.nih.gov.\n", sep=" ")
    out <- 3 #Unknown error
    names(out) <- "FAILED"
  }
  if( status == "UNKNOWN" ) {
    cat("Search", reqID, "not found on server on server", sep=" ")
    out <- 4 #RID already deleted on NCBI server (results lasts 24 hours, of 30minutes for large requests)names
    names(out) <- "UNKNOWN"
  }
  out
}

#download results of RID, return as table 
NCBIdownload <- function(reqID) {
   uri <- paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID=",reqID,sep="")
   res <- getURL(uri) #retrieving from NCBI server
   res <- substr(res,regexpr("(hits found)+",res)+11,regexpr("(</PRE>)+",res)-2) # removing data before/after table
   res <- strsplit(res,"\\n") # reformatting table
   res <- do.call(rbind, strsplit(res[[1]],"\\t"))
   colnames(res) <- c("query id", "subject ids", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
   res
}

############################################
# input <- matrix(NA,2,2)
# input <- c("ID1","ID2","CTTCGTTTCCCTCTTCTGCGATTTC","GATTGCACCTTCGATGGCCCTGAAA")
# dim(input) <- c(2,2)
# colnames(input) <- c("ID","query")
# eValue <- "1e-2"
##########################################

NCBIblast <- function(input,eValue,program="blastn",db="nt",verbose=FALSE) {
  mRID <- matrix(NA,nrow(input),6)
  rownames(mRID) <- rownames(input)
  colnames(mRID) <- c("ID","query","RID","RTOE","ready","result")
  mRID[,c("ID","query")] = input
  mRID[,c("ready","result")] = 0
  cnt <- 1
  #send ALL the queries!!
  lapply(mRID[,"query"],function(jj){
    resSub <- NCBIsubmit(jj,eValue,program,db)
    mRID[cnt,"RID"] <<- resSub$RID
    mRID[cnt,"RTOE"] <<- resSub$RTOE
    #cat(cnt,jj,resSub$RID,resSub$RTOE,"\n",sep=" - ")
    cnt <<- cnt + 1
  })
  #starting check, and if ready, download loop
  output <- NULL # defining output column
  while ( any(mRID[,"ready"] == 0)) {
    Sys.sleep(runif(1)*5)
    ## make list of unfinished BLASTs
    toDo <- which( mRID[,"ready"] == 0)
    if(verbose) { cat(length(toDo), "left to do\n", sep=" ") }
    ## apply 'NCBIcheck' to unfinished BLASTs
    mRID[toDo,"ready"] <- unlist(lapply(mRID[toDo,"RID"],NCBIcheck))
    j <- 1
    apply(mRID,1,function(x) {
      if(x["ready"] == 1 & x["result"] != 1) {
        res <- NCBIdownload(x["RID"])
        sequenc <- matrix( c( rep( mRID[j,"ID"],nrow(res)), rep(mRID[j,"query"],nrow(res)) ), nrow(res), 2)
        colnames(sequenc) <- c("ID","sequence")
        res1 <- cbind(sequenc,res) #combining blast result and queries 
        output <<- rbind(output,res1)
        mRID[j,"result"] <<- 1 
        
      }
      j <<- j + 1
    })
  }
  output
}
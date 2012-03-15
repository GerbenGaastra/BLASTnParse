# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# NCBI blast, Qblast

NCBIsubmit <- function(query,eValue) {
  # building and sending request
  h = basicTextGatherer()
  hardCoded <- "CMD=Put&PROGRAM=blastn&DATABASE=nt"
  inputFields <- paste(hardCoded,"&EXPECT=",eValue,"&QUERY=",query,collapse="",sep="")
  cat(inputFields)
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
    if ( grep("(ThereAreHits=yes+)",tempres ) ) {
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
    (out) <- "UNKNOWN"
  }
  out
}

# download results of RID, return as table 
# NCBIdownload <- function(reqID) {
   # uri <- paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID=",reqID,sep="")
   # res <- getURL(uri) #retrieving from NCBI server
   # res <- substr(res,regexpr("(hits found)+",res)+11,regexpr("(</PRE>)+",res)-2) # removing data before/after table
   # res <- strsplit(res,"\\n") # reformatting table
   # res <- do.call(rbind, strsplit(res[[1]],"\\t"))
   # colnames(res) <- c("query id", "subject ids", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
   # res
# }

############################################
# input <- matrix(NA,2,1)
# rownames(input) <- c("ID1","ID2")
# colnames(input) <- c("query")
# input[1,1] = "CTTCGTTTCCCTCTTCTGCGATTTC"
# input[2,1] = "GATTGCACCTTCGATGGCCCTGAAA"
# eValue <- "1e-2"
##########################################

# NCBIblast <- function(input,eValue) {
  # mRID <- matrix(NA,nrow(input),5)
  # rownames(mRID) <- rownames(input)
  # colnames(mRID) <- c("query","RID","RTOE","ready","result")
  # mRID[,"query"] = input[,"query"]
  # mRID[,c("ready","result")] = 0
  # cnt <- 1
  # lapply(input[,"query"],function(jj,cnt){
    # resSub <- submit(jj)
    # mRID[cnt,"RID"] <- resSub$RID
    # mRID[cnt,"RTOE"] <- resSub$RTOE
    # cat(cnt,jj,resSub$RID,resSub$RTOE,"\n",sep=" - ")
    # cnt <<- cnt + 1
  # },cnt)
  # send ALL the queries!!
  # for( i in 1:nrow(input)) {
    # resSub <- submit(input[i,"query"],eValue)
    #cat(res[[1]],res[[2]],"\n",sep=" - ")
    # mRID[i,"RID"] <- resSub$RID
    # mRID[i,"RTOE"] <- resSub$RTOE
  # }
  # output <- NULL # defining output column
  # while ( any(mRID[,"ready"] == 0)) {
    # cat( length(which(mRID[,"ready"] == 0)), "left to do\n", sep=" ")
    # Sys.sleep(runif(1)*5)
    # toDo <- which( mRID[,"ready"] == 0)
    # mRID[toDo,"ready"] <- unlist(lapply(mRID[toDo,"RID"],check))
    # cat(unlist(lapply(mRID[toDo,"RID"],check)))
    # print(mRID)
    # for( j in 1:nrow(mRID) ) {
      # if(mRID[j,"ready"] == 1 & mRID[j,"result"] != 1) {
        # res <- download(mRID[j,"RID"])
        # sequenc <- matrix(c(rep(rownames(mRID)[j],nrow(res)), rep(mRID[j,"query"], nrow(res)) ), nrow(res),2)
        # colnames(sequenc) <- c("ID","sequence")
        # res1 <- cbind(sequenc,res[,-1]) #combining blast result and queries (removing non informational collumn with -1)
        # output <- rbind(output,res1)
        # mRID[j,"result"] <- 1 
      # }
    # }
  # }
  # output
# }

# input <- matrix(NA,2,1)
# rownames(input) <- c("AGI17","AGI18")
# colnames(input) <- c("query")
# input[1,1] = "CCTAAATTTCTGATTTTCAGAGTTTGAGACCGTTTCGATTCAAACCCCCACCGAACCCAA"
# input[2,1] = "ACAGCTAGGGGAAATGGATCATCAGTAGCCGAGGAGCTCAATTCAAATTCAAAGGAAAAA"
# eValue <- "1e-2"
# runBatchBlast(input,eValue)
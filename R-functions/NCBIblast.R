# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# NCBI blast, Qblast




submit <- function(query) {
  
  # send ALL the queries!!
  h = basicTextGatherer()
  curlPerform(url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
    httpheader=c(Accept="text/xml", Accept="multipart/*", 'Content-Type' = "application/x-www-form-urlencoded"),
    postfields=paste("CMD=Put&PROGRAM=blastn&DATABASE_PREFIX=genomes/c_elegans_chr&DATABASE=nt&QUERY=",query,sep=""),
    writefunction = h$update,
    verbose = TRUE
  )
  ## ara_chr        arabidopsis chr
  ## 3702           Arabidopsis thaliana est
  ## c_elegans_chr  c_elegans chr
  
  #parse h$value for RID
  out <- vector("list",2)
  res1 <- regexpr( "(RID%3D+)([A-Za-z0-9-]{2,})",h$value() )
  res2 <- regexpr( "(RTOE%3D+)([0-9]{1,})",h$value() )
  out[[1]] <- substr(h$value(),res1+6, res1 + attr(res1, "match.length") -1 )
  out[[2]] <- substr(h$value(),res2+7, res2 + attr(res2, "match.length") -1 )
  names(out) <- c("RID","RTOE")
  out
}

check <- function( reqID ) {
  tempres <- getURL(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=",reqID,sep=""))
  res <- regexpr( "(Status=+)([A-Za-z]{1,})",tempres )
  status <- substr(tempres,res+7,res+attr(res,"match.length")-1)
  if( status == "READY") {
    out <- 1
  } else {
    out <- 0
  }
  out
}

download <- function(reqID) {
   uri <- paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID=",reqID,sep="")
   res <- getURL(uri)
   res <- substr(res,regexpr("(hits found)+",res)+11,regexpr("(</PRE>)+",res)-2)
   res <- strsplit(res,"\\n")
   res <- do.call(rbind, strsplit(res[[1]],"\\t"))
   colnames(res) <- c("query id", "subject ids", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score")
   res
}




# input <- matrix(NA,2,1)
# rownames(input) <- c("ID1","ID2")
# colnames(input) <- c("query")
# input[1,1] = "CTTCGTTTCCCTCTTCTGCGATTTC"
# input[2,1] = "GATTGCACCTTCGATGGCCCTGAAA"

output <- matrix(NA,1,13)
runBatchBlast <- function(input) {
  mRID <- matrix(NA,nrow(input),5)
  rownames(mRID) <- rownames(input)
  colnames(mRID) <- c("query","RID","RTOE","ready","result")
  mRID[,"query"] = input[,"query"]
  mRID[,c("ready","result")] = 0
  cnt <- 1
  apply(input,1,function(jj){
    resSub <- submit(jj["query"])
    mRID[cnt,"RID"] <- resSub$RID
    mRID[cnt,"RTOE"] <- resSub$RTOE
    cnt <<- cnt + 1
  })
  for( i in 1:nrow(input)) {
    resSub <- submit(input[i,"query"])
    #cat(res[[1]],res[[2]],"\n",sep=" - ")
    mRID[i,"RID"] <- resSub$RID
    mRID[i,"RTOE"] <- resSub$RTOE
  }
  while ( any(mRID[,"ready"] == 0)) {
    cat( length(which(mRID[,"ready"] == 0)), "left to do\n", sep=" ")
    Sys.sleep(runif(1)*10)
    
    mRID[,"ready"] <- lapply(mRID[,"RID"],check)
    
    mRID[,"ready"] <- apply(mRID,1,function(x){
      check(x["RID"])
    })
    output <- NULL
    for( j in 1:nrow(mRID) ) {
      if(mRID[j,"ready"] == 1 & mRID[j,"result"] != 1) {
        res <- download(mRID[j,"RID"])
        sequenc <- matrix(c(rep(rownames(mRID)[j],nrow(res)), rep(mRID[j,"query"], nrow(res)) ), nrow(res),2)
        colnames(sequenc) <- c("ID","sequence")
        res1 <- cbind(sequenc,res[,-1]) #combining blast result and queries (removing non informational collumn with -1)
        output <- rbind(output,res1)
        mRID[j,"result"] <- 1 
      }
    }
  }
  output
}
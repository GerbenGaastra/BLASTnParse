# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# NCBI blast, Qblast


# temp <- httpPOST("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi","CMD=Put&QUERY=MKN&DATABASE=nr&PROGRAM=blastp&FILTER=L&ViewReport=View%20report&HITLIST_SZE=500",handle)
# temp <- httpPUT("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&QUERY=tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcactTCTATCTATTATATCCTATACGTGTGTGATAGTACACACAA")
# temp2 <- httpGET("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?\
  # CMD=Get&RID=954517013-7639-11119",handle)
  
  # write(temp2,"temp.html")
  
  # h = basicTextGatherer()
  # curlPerform(url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
              # httpheader=c(Accept="text/xml", Accept="multipart/*", 'Content-Type' = "application/x-www-form-urlencoded"),
              # postfields="CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&QUERY=tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact",
              # writefunction = h$update,
              # verbose = TRUE
             # )
     # write(h$value(),"temp.html")          
  
  # reqID <- "N1JT4RUU016" 
  
  
  # tempres <- getURL(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=",reqID,sep=""))
   # write(tempres,"temp.html")          
   
   # http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=N1JT4RUU016
   
   ###    tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact
   
###################################################
submit <- function(query) {
  
  # send ALL the queries!!
  h = basicTextGatherer()
  curlPerform(url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
    httpheader=c(Accept="text/xml", Accept="multipart/*", 'Content-Type' = "application/x-www-form-urlencoded"),
    postfields=paste("CMD=Put&PROGRAM=blastn&DATABASE_PREFIX=genomes/ara_chr&DATABASE=nt&QUERY=",query,sep=""),
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

res <- submit(query)
check(res$RID)
txt <- download(res$RID)
write(txt,"temp.html")

download <- function(reqID) {
   uri <- paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID=",reqID,sep="")
   res <- getURL(uri)
   res
}




input <- matrix(NA,2,1)
rownames(input) <- c("ID1","ID2")
colnames(input) <- c("query")
input[1,1] = "CTTCGTTTCCCTCTTCTGCGATTTC"
input[2,1] = "GATTGCACCTTCGATGGCCCTGAAA"


runBatchBlast <- function(input) {
  output <- matrix(NA,nrow(input),5)
  rownames(output) <- rownames(input)
  colnames(output) <- c("query","RID","RTOE","ready","result")
  output[,"query"] = input[,"query"]
  output[,"ready"] = 0
  for( i in 1:nrow(input)) {
    resSub <- submit(input[i,"query"])
    #cat(res[[1]],res[[2]],"\n",sep=" - ")
    output[i,"RID"] <- resSub$RID
    output[i,"RTOE"] <- resSub$RTOE
  }
  
   while ( any(output[,"ready"] == 0)) {
    Sys.sleep(runif(1)*10)
    for( i in 1:nrow(output)) {
      iStatus <- check(output[i,"RID"])
      output[i,"ready"] <- iStatus
    }
    
  }
}
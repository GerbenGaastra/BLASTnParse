# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# NCBI blast, Qblast


handle <- getCurlHandle()
temp <- httpPOST("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi","CMD=Put&QUERY=MKN&DATABASE=nr&PROGRAM=blastp&FILTER=L&ViewReport=View%20report&HITLIST_SZE=500",handle)
temp <- httpPUT("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&QUERY=tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcactTCTATCTATTATATCCTATACGTGTGTGATAGTACACACAA")
temp2 <- httpGET("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?\
  CMD=Get&RID=954517013-7639-11119",handle)
  
  write(temp2,"temp.html")
  
  h = basicTextGatherer()
  curlPerform(url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
              httpheader=c(Accept="text/xml", Accept="multipart/*", 'Content-Type' = "application/x-www-form-urlencoded"),
              postfields="CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&QUERY=tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact",
              writefunction = h$update,
              verbose = TRUE
             )
     write(h$value(),"temp.html")          
  
  reqID <- "N1JT4RUU016" 
  
  
  tempres <- getURL(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=",reqID,sep=""))
   write(tempres,"temp.html")          
   
   http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=N1JT4RUU016
   
   
   
###################################################
submit <- function( list of queries ) {
  # send ALL the queries!!
  
  
  # basic idea
  h = basicTextGatherer()
  curlPerform(url = "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi",
              httpheader=c(Accept="text/xml", Accept="multipart/*", 'Content-Type' = "application/x-www-form-urlencoded"),
              postfields="CMD=Put&PROGRAM=blastn&MEGABLAST=on&DATABASE=nt&QUERY=tcgtttattatttgtcaccgggttccatcccccttacgtttgacaatcattgcactcact",
              writefunction = h$update,
              verbose = TRUE
             )
  parse h$value() to get reqID
}

check <- function( list of reqID ) {
  tempres <- getURL(paste("http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=",reqID,sep=""))
  
  parse tempres to get status
}

download <- function( checked list of reqID) {
   get >>> http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=N1JT4RUU016
}




runBatchBlast <- function() {
  submit()
  
  while ( any(list reqID[,2]) {
    check()
    download()
  }
}
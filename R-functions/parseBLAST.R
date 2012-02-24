# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Retrieve Position and Chromosome from post-data

parseBLAST <- function(BLAST) {
  out <- vector("list",2)
  resChr <- regexpr("CHROMOSOME_([IVX]+)?",BLAST)
  out[[1]] <- as.character(substr(BLAST, resChr+11, resChr + attr(resChr, "match.length") -1))
  resPos <- regexpr("(Sbjct+)(.{7})([0123456789]{2,})",BLAST)
  out[[2]] <- as.character(substr(BLAST, resPos+12, resPos + attr(resPos, "match.length") -1))
  names(out) <- c("chr","pos")
  out
}
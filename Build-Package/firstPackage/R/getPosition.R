# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Blasts nucluotide seqeunce vs nucleotide database
# on http://www.wormbase.org/db/searches/blast_blat
# To-Do: Retrieves first n sequences

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




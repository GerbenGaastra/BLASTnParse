# (c) by GT Gaastra
# created last 24-02-2012
# modified last 24-02-2012
# Wrapper for postForm

downloadForm <- function(uri, postValues) {
  if(url.exists(uri)) {
    HTMLreturn = postForm(uri, .params = postValues)
  } else {
    error <- cat("url: ",uri,"not available.\n",sep=" ")
    stop(error)
  }
  HTMLreturn
}
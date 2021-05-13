upload.file <- function(filename, featurecol, samplecol) {
  fh <- read.csv(filename, header = TRUE)
  fileparsed <- fh[,c(featurecol, samplecol)]
  fileparsed <- as.data.frame(t(fileparsed)); 
  colnames(fileparsed) <- fileparsed[1,]; fileparsed <- fileparsed[-1,]
  return(fileparsed)
}


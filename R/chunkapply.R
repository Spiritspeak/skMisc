
# Use package abind: asub for subsetting, abind for binding
chunkapply <- function(x,MARGIN,FUN,...,chunksize=1000){
  chunks <- subdivide(seq_len(dim(x)[MARGIN]), divlen=chunksize)
  
}


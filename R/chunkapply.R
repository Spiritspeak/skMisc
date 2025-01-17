

#' Iterate over chunks of a vector
#' 
#' This returns an iterator that counts from one, returning integer sequences 
#' of a pre-specified length with each iteration up until an optional maximum value. 
#'
#' @param chunksize How long should each sequence be?
#' @param count What should be the maximum value returned by the iterator? 
#' Above this value, the iterator will return \code{NULL}.
#'
#' @return An iterator function.
#' @note This does not return an iterator in the manner used by the iterators package,
#' since those objects currently misbehave in an Rstudio live session.
#' 
#' @export
#'
#' @examples
#' myit <- ichunk(chunksize=10,count=15)
#' myit()
#' myit()
#' myit()
#' 
ichunk <- function(chunksize,count=NULL){
  if (!is.numeric(chunksize) || length(chunksize) != 1){
    stop("chunksize must be a single numeric value")
  }else if(length(count)>1 || (!is.numeric(count) && !is.null(count))){
    stop("count must be a single numeric value or NULL")
  }
  chunksize <- as.integer(chunksize)
  i <- 0L
  nextEl <- function(){
    if (is.null(count) || i*chunksize < count){
      i <<- i + 1L
      iterseq <- seq(from=1+(i-1)*chunksize,
                     to=min(count,i*chunksize))
      return(iterseq)
    }
  }
  nextEl
}


chunkapply.vector <- function(x, FUN, ..., chunksize=1000){
  stopifnot(is.vector(x), is.function(FUN), is.numeric(chunksize))
  out <- vector(length=length(x))
  counter <- ichunk(chunksize=chunksize, count=length(x))
  warnlength <- FALSE
  while(T){
    currchunk <- counter()
    if(length(currchunk) == 0){ break; }
    currout <- forceAndCall(1, FUN, currx, ...)
    if(length(currout) != length()){ warnlength <- TRUE }
    out[currchunk] <- currout
  }
  if(warnlength){ warning("Output length did not match input length on at least one function call") }
  return(out)
}




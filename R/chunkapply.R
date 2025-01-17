

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

arrsubset <- function(x,dim,idx,drop=FALSE){
  newx <- Quote(x[ ,drop=drop])
  newx <- newx[c(1,2,rep(3,length(dim(x))),4)]
  newx[2+dim] <- list(idx)
  return(eval(newx))
}

chunkapply <- function(x, FUN, ..., chunksize=1000, MARGIN=1){
  stopifnot(is.vector(x), is.function(FUN), is.numeric(chunksize))
  out <- vector(length=length(x))
  counter <- ichunk(chunksize=chunksize, count=length(x))
  grabber <- NULL
  if(is.null(dim(x))){
    grabber <- function(idx){ x[idx] }
  }else{
    stopifnot(length(MARGIN)==1,any(dim(x)==MARGIN))
    grabber <- function(idx){ arrsubset(x=x,dim=MARGIN,idx=idx) }
  }
  while(T){
    curridx <- counter()
    if(length(curridx) == 0){ break; }
    currx <- grabber(curridx)
    currout <- forceAndCall(1, FUN, currx, ...)
    if(length(currout) != length(curridx)){ 
      stop("Function output length (",length(currout),
           ") was not the same as input length (",length(currx),")")
    }else if(!is.vector(currout)){
      stop("Function output not in form of a vector")
    }
    out[curridx] <- currout
  }
  return(out)
}
#chunkapply.vector(x=1:1000,FUN=as.character,chunksize=5)


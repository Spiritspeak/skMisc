
dropLeadingZero_old <- function(x){
  xnew <- c()
  for(i in x){
    if(isTRUE(i==0)){
      xnew <- c(xnew,"0")
    } else if (isTRUE(i>=1) | isTRUE(i<=-1)){
      xnew <- c(xnew, as.character(i))
    } else
      xnew <- c(xnew, gsub("(?<![0-9])0+(?=\\.)", "", i, perl = TRUE))
  }
  return(xnew)
}


#' Compound tokens without overflowing memory and crashing R
#' @description A wrapper around \link[quanteda]{tokens_compound} that processes your tokens in chunks, 
#' set by argument \code{stepsize}. See \link[quanteda]{tokens_compound} for more info.
#'
#' @export
#' @examples 
#' toks<-tokens(data_corpus_inaugural)
#' compounded<-tokens_compound_stepwise(x=toks,pattern="I am",stepsize=10)
#' 
#' #note: does not work?
#' 
tokens_compound_stepwise<-function(x, pattern, stepsize=100, concatenator = "_", 
                                   valuetype = c("glob", "regex", "fixed"), 
                                   case_insensitive = TRUE, join = TRUE){
  valuetype<-match.arg(valuetype)
  output<-list()
  lx<-length(x)
  mxsteps<-ceiling(lx/stepsize)
  for(i in seq_len(mxsteps)){
    cat("\rCompounding tokens... ",round(100*(i-1)/mxsteps,digits=2),"%",sep="")
    range<-(1+(i-1)*stepsize):min(i*stepsize, lx)
    currcomp<-quanteda::tokens_compound(x[range],pattern,concatenator,valuetype,case_insensitive,join)
    output<-c(output,currcomp)
  }
  cat("\rCompounding tokens... finished.")
  return(output)
}



#' The logistic function
#'
#' @param x A numeric value
#'
#' @return \code{x} with the logistic function applied
#' @export
#'
#' @examples
#' logistic(-10:10)
#' 
logistic<-function(x){
  1/(1+exp(-x))
}
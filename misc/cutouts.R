
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


#wtd.median
#' Weighted Median
#'
#' @param x an input vector
#' @param wts a vector of weights
#' @param na.rm Logical indicating whether NA values in the input and weight vectors should be stripped. 
#'
#' @return A weighted median of the input values and weights.
#' @export
#'
#' @examples
#' wtd.median(1:5,c(.5,4,1,2,1))
#' 
wtd.median<-function(x,wts,na.rm=T){
  #clean
  if(na.rm){
    to.include<-which(!(is.na(x) | is.na(wts)))
    x<-x[to.include]
    wts<-wts[to.include]
  }
  
  #sort
  xord<-order(x)
  x<-x[xord]
  wts<-wts[xord]
  
  #standardize
  wts<-wts/sum(wts)
  
  #find middle
  cumwts<-cumsum(wts)-.5
  sig<-sign(cumwts)
  if(!any(sig==0)){
    midx<-which(sig==1)[1]
    out<-x[midx]
  }else{
    midx<-c(which(sig==1)[1],which(sig==0)[1])
    out<-mean(x[midx])
  }
  return(out)
}


#' Downweight outliers
#' @description Computes weights; trials within certain bounds of the mean receive the maximum weight while trials
#' outside these bounds are downweighted to 0 or an optional minimum.
#'
#' @param x A numeric vector
#' @param mean An optional mean of the vector
#' @param s An optional standard deviation of the vector
#' @param sdist The number of standard deviations beyond which values should be downweighted
#' @param taper A number indicating how strongly values exceeding the standard deviation should taper off
#' @param scale How the weight vector should be scaled: "norm" sets the sum to 1, "max" sets the maximum to 1.
#' @param min A minimum weight. 
#'
#' @return A numeric vector of weights
#' @export
#'
#' @examples
logit.weightfun<-function(x,mean=mean(x),s=sd(x),
                          sdist=3,taper=10,
                          scale=c("max","norm"),min=0){
  scale <- match.arg(scale)
  zx <- (x-m)/s * taper
  out <- inv.logit((zx-sdist*taper)) * inv.logit((-zx-sdist*taper))
  out <- (out/max(out)) * (1-min) + min
  if(scale== "norm"){
    out < -out / sum(out)
  }
  return(out)
}


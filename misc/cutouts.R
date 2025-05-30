
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



#' Influence function of the Pearson correlation coefficient
#'
#' @param x,y Numeric vectors
#'
#' @return Influence values of all observations.
#' @export
#'
#' @examples
#' outlier<-numeric(100)
#' outlier[1]<-1000
#' cor.influence(rnorm(100)+outlier,rnorm(100)+outlier)
#' 
cor.influence<-function(x,y){
  x<-x-mean(x)
  y<-y-mean(y)
  x*y-(x^2+y^2)/2*cor(x,y)
}


#' Read and merge all .csv files in a folder
#'
#' @param folder path to a folder
#' @param readfunc list of functions that will be used to read the files; 
#' if the first function fails, the second function will be used, etc.
#'
#' @return A data.frame containing all merged .csv files 
#' @export
#' 
#' @examples 
#' 
read.csv.folder <- function(folder="./", readfunc=list(read.csv,read.csv2,read.table)){
  flist<-list.files(folder)
  flist<-flist[grepl(".csv",flist)]
  datlist<-list()
  ct<-0
  for(file in flist){
    ct<-ct+1
    datlist[[ct]]<-"fail"
    for(i in seq_len(length(readfunc))){
      #try(
      datlist[[ct]]<-do.call(readfunc[[i]],list(file=paste0(folder,file),stringsAsFactors=F))
      #,silent=T)
      if(any(datlist[[ct]]!="fail")){ break; }
    }
    if(any(datlist[[ct]]=="fail")){
      warning("Failed to read ",file)
    }
  }
  combodat<-datlist[[1]]
  if(length(datlist)>1){
    for(i in 2:length(datlist)){
      combodat<-rbind(combodat,datlist[[i]])
    }
  }
  return(combodat)
}

#' Check convergence of a \code{brmsfit} object
#'
#' @param x A \code{brmsfit} object
#' @param min.ess The minimum effective sample size for any parameter. 
#' Values below will trigger a warning.
#' @param max.rhat The maximum R-hat value for any parameter. 
#' Values above will trigger a warning.
#'
#' @return Silently returns a \code{data.frame} with parameter names, 
#' ESS, and rhat values.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- brm(time ~ age * sex, data = kidney)
#' brms.check(fit)
#' }
#' 
brms.check <- function(x,min.ess=400,max.rhat=1.05){
  info <- data.frame(par=dimnames(x$fit)$parameters,
                     ess=posterior::ess_basic(x),
                     rhat=posterior::rhat(x))
  
  toolow.ess <- info$ess < min.ess
  if(any(toolow.ess)){
    warning("The following parameter(s) have an excessively low ",
            "effective sample size, consider raising the iteration count: ",
            paste0(info$par[toolow.ess],collapse=", "))
  }
  
  toohigh.rhat <- info$rhat > max.rhat
  if(any(toohigh.rhat)){
    warning("The following parameter(s) have an excessively high ",
            "Rhat, consider changing the model or adding more warmup samples: ",
            paste0(info$par[toohigh.rhat],collapse=", "))
  }
  
  return(invisible(info))
}


# Add option to define positive, negative, and zero color; 
# add support for simple vector input
colorEdges<-function(x,maxedge=NULL){
  edgevec<-as.vector(x)
  if(is.null(maxedge)){ 
    edgevec <- edgevec/max(abs(edgevec)) 
  }else{
    edgevec <- edgevec/abs(maxedge)
  }
  cols<-hsv(h=ifelse(sign(edgevec)==1,2/3,0),
            s=abs(edgevec)*.9+ifelse(edgevec==0,0,.1),
            v=1)
  matrix(cols,ncol=ncol(x),nrow=nrow(x))
}

angular.mean <- function(x,period=2*pi){
  period/(2*pi)*atan2(mean(sin(x/period*2*pi)),mean(cos(x/period*2*pi)))
}

# angular.mean(c(12,11,7),period=12)
# angular.mean(c(12,11,8),period=12)



#' Generate and concatenate multiple integer sequences
#' 
#' This generates multiple integer sequences and concatenates them.
#'
#' @param from,to the starting and (maximal) end values of the sequences. 
#' Multiple can be given.
#'
#' @returns An integer vector of multiple concatenated integer sequences.
#' @export
#' @author Sercan Kahveci
#'
#' @examples
#' seq_composite(from=c(1,7),to=c(3,8))
#' 
seq_composite <- function(from,to){
  unlist(mapply(\(x,y){x:y},from,to))
}

#' t-distribution fitter
#' This is a wrapper around [MASS::fitdistr()] specifically intended to
#' fit the t-distribution.
#' 
#' @param x Vector of values to fit the t-distribution to
#' @param df Starting df value
#'
#' @return An object of class \code{"fitdistr"}.
#' @export
#'
#' @examples
#' h<-rt(1000,df=3)*3+10
#' tpars(h)
#' 
tpars <- function(x, df=30){
  MASS::fitdistr(x=x,
                 densfun="t",
                 start=list(m=mean(x), s=sd(x), df=df),
                 lower=c(m=-Inf, s=0, df=1))
}
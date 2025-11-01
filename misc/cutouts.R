
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


#' Carry non-NA values forward into NA values
#' 
#' This replaces every \code{NA} value with their last preceding non-\code{NA} value.
#'
#' @param x A vector.
#'
#' @returns \code{x} with all \code{NA} values replaced with their 
#' last preceding non-\code{NA} value.
#' @export
#' @author Sercan Kahveci
#'
#' @examples
#' carryforward(c(NA,1,NA,NA,2))
#' 
carryforward <- function(x){
  runend <- which(!is.na(x) & is.na(quicklag(x))) -1L
  runstart <- which(is.na(x) & !is.na(quicklag(x)))
  if(length(runstart)==0){ return(x) }
  if(runstart[1] > runend[1]){
    runend <- runend[-1]
  }
  if(length(runend)==0){
    runend <- c(runend, length(x))
  }else if(runstart[length(runstart)] > runend[length(runend)]){
    runend <- c(runend, length(x))
  }
  
  x[seq_composite(runstart, runend)] <- rep(x[runstart-1], times=runend-runstart+1)
  return(x)
}


# TODO: change format to: numeric = c("age","height")
# Add support for tidyselectors

#' Change classes of columns in a data.frame
#' 
#' @description \code{retype()} changes the class of specific columns; 
#' \code{retype_all()} changes the class of all columns of a given class.
#'
#' @param df a data frame
#' @param ... Unquoted column names, paired with the desired class, e.g. 
#' 
#' \code{age = numeric(), language = character()}
#'
#' @export
#' @author Sercan Kahveci
#'
#' @examples 
#' sapply(ToothGrowth,class)
#' NewToothGrowth <- retype(ToothGrowth, supp = character(), dose = factor())
#' sapply(NewToothGrowth,class)
#' 
retype <- function(df, ...){
  args <- list(...)
  
  varnames <- names(args)
  vartypes <- sapply(args, class)
  
  effcols <- names(df)[names(df) %in% varnames]
  
  for(effcol in effcols){
    df[,effcol] <- as(df[,effcol], vartypes[which(varnames == effcol)])
  }
  return(df)
}

#' @rdname retype
#' @param df A \code{data.frame}.
#' @param from An empty vector of the class to convert from, or a string. 
#' Columns sharing the class of argument \code{from} will be converted 
#' to the class of argument \code{to}.
#' @param to An empty vector of the class to convert to, or a string. 
#' Columns sharing the class of argument \code{from} will be converted 
#' to the class of argument \code{to}.
#'
#' @export
#'
#' @examples 
#'
#' sapply(mtcars,class)
#' newmtcars <- retype_all(mtcars,from="numeric",to="character")
#' sapply(newmtcars,class)
#' 
retype_all <- function(df, from, to){
  for(i in which(sapply(df, class) == from)){
    df[[i]] <- as(df[[i]], to)
  }
  df
}

#' Verify variable types in bulk
#'
#' @param ... Named arguments, where the argument is the object to be checked and 
#' the name of the argument is the mode (numeric, list, character, etc).
#'
#' @return Returns true on success, causes error if not.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' try(verify_types(character="test",numeric=0000,character=12345))
#' 
verify_types <- function(...){
  args <- list(...)
  call <- as.list(match.call()[-1])
  types <- unique(names(args))
  for(type in types){
    ids <- which(type == names(args))
    for(id in ids){
      if(!do.call(paste0("is.",type),list(args[[id]]))){
        stop("Variable ",as.character(call[[id]])," is not of type ",type)
      }
    }
  }
  return(T)
}



ComputeLowerModels <- function(form, data, group="", ...){
  args <- list(...)
  cluster <- makeCluster(detectCores()-1)
  registerDoParallel(cluster)
  testforms <- RemoveTopTerms(form,group)
  
  results <-
    foreach(currform=testforms, .packages="lmerTest") %dopar% {
      do.call(lmer, c(list(formula=currform, data=data), args))
    }
  
  stopCluster(cluster)
  return(results)
}

ComputeLowerModels2 <- function(model, data, group="", ...){
  form <- formula(model)
  data <- model.frame(model)
  #data <- model@call$data
  args <- list(...)
  testforms <- RemoveTopTerms(form, group)
  
  cluster <- makeCluster(detectCores())
  registerDoParallel(cluster)
  results<-
    foreach(currform=testforms,.packages="lmerTest") %dopar% {
      do.call(lmer,c(list(formula=currform, data=data, REML=F), args))
    }
  stopCluster(cluster)
  names(results) <- names(testforms)
  
  warns <- sapply(results, function(x){ paste(x@optinfo$warnings, collapse="\n") })
  message(paste0("Warning in model ", names(results)[warns!=""], ": ", warns[warns != ""]))
  
  anovatable <- AnovaTable(model, results)
  
  print(anovatable)
  return(invisible(list(anovatable=anovatable, models=results)))
}


arrsubset <- function(x, dim, idx, drop=FALSE){
  newx <- Quote(x[ ,drop=drop])
  newx <- newx[c(1,2,rep(3,length(dim(x))),4)]
  newx[2+dim] <- list(idx)
  return(eval(newx))
}

chunkapply <- function(x, FUN, ..., chunksize=1000, MARGIN=1){
  stopifnot(is.function(FUN), is.numeric(chunksize))
  out <- vector(length=length(x))
  counter <- ichunk(chunksize=chunksize, count=length(x))
  
  chunklength <- prod(dim(x)[-MARGIN])
  
  grabber <- NULL
  if(is.null(dim(x))){
    grabber <- function(idx){ x[idx] }
  }else{
    stopifnot(length(MARGIN)==1,any(seq_along(dim(x))==MARGIN))
    grabber <- function(idx){ arrsubset(x=x,dim=MARGIN,idx=idx) }
  }
  while(T){
    curridx <- counter()
    if(length(curridx) == 0){ break; }
    currx <- grabber(curridx)
    #currout <- forceAndCall(1, FUN, currx, ...)
    currout <- forceAndCall(1, FUN, currx)
    if(dim(currout) != length(curridx)){
      stop("Function output length (",length(currout),
           ") was not the same as input length (",length(currx),")")
    }else if(!is.vector(currout)){
      stop("Function output not in form of a vector")
    }
    out[curridx] <- currout
  }
  return(out)
}
#chunkapply(x=1:1000,FUN=as.character,chunksize=5)
#chunkapply(x=matrix(1:1000,nrow=100),FUN=as.character,chunksize=5)


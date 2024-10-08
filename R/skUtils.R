# Files that are "up to standard": holdout.R and formatting.R


.onLoad<-function(libname, pkgname){
  packageStartupMessage("Thank you for loading skMisc v0.01")
}

#I don't want to import rlang, so it will be done this way instead.
args2strings <- function(...) sapply(substitute({ ... })[-1], deparse)

#' clamp
#'
#' @param val The vector/matrix to clamp
#' @param minval Minimum value; all lower values are clamped to this value
#' @param maxval Maximum value; all higher values are clamped to this value
#'
#' @return Clamped vector.
#' @export 
#'
#' @examples 
#' clamp(0:10,2,8)
#' clamp0(rnorm(10))
clamp <- function(val,minval=-Inf,maxval=Inf){
  val[val<minval]<-minval
  val[val>maxval]<-maxval
  val
}

#' @rdname clamp
#' @export
clamp0 <- function(val,minval=0,maxval=1){
  val[val<minval]<-minval
  val[val>maxval]<-maxval
  val
}

# duplicated() determines which elements of a vector or data frame are 
# duplicates of elements with smaller subscripts, and returns a logical vector 
# indicating which elements (rows) are duplicates.

#' Count duplicate values in a vector
#' 
#' which.duplidate() determines for each element of a vector 
#' how many times it has occurred so far.
#' It works similarly to [base::duplicated()] which only determines 
#' whether a value has occurred before and not how many times.
#' 
#' @param x A vector
#'
#' @return A vector of the same length as \code{x}, where each value has been replaced
#' with the number of times its value has been repeated so far.
#' @export
#'
#' @examples
#' which.duplicate(c(1,6,5,2,1,1,8,6,5))
which.duplicate<-function(x){
  vals<-unique(x)
  repvec<-numeric(length(vals))
  for(v in vals){
    repvec[x==v]<-seq_len(sum(x==v))
  }
  return(repvec)
}

#' Scale a vector
#' 
#' Like scale() but returns a vector and is faster
#' 
#' @param x Numeric vector to standardize
#'
#' @return Scaled numeric vector with mean of 0 and sd of 1
#' @export
#'
#' @examples
#' vec.scale(1:10)
#' 
vec.scale<-function(x){
  xt<-na.omit(x)
  m<-mean.default(xt)
  (x-m)/sqrt((sum((xt-m)^2)/(length(xt)-1)))
}

#' Smooth a numeric vector using a moving window algorithm
#'
#' @param vect Numeric vector to be smoothened
#' @param width Over how many values should the vector be averaged?  
#' @param both.sides If TRUE (default), takes the mean of \code{width} values before and after the current index. If FALSE, only takes values ahead of the current index.
#' @param alg Method by which to smooth the vector. 'mean' or 'gauss' are supported.
#'
#' @return Smoothed numeric vector
#' @export
#'
#' @examples temp<- smoothvect(beaver1$temp)
#' plot(temp,type="l")
smoothvect<-function(vect,width=2,both.sides=T,alg=c("mean","gauss")){
  output<-numeric()
  
  normsum<-function(x){x/sum(x)}
  
  if(alg=="mean"){
    for(i in seq_len(length(vect))){
      output[i]<-mean(vect[ max(i-width*both.sides,1):min(i+width,length(vect))],na.rm=T)
    }
  }
  if(alg=="gauss"){
    winvect<- (-width*both.sides):(width)
    window<-dnorm(winvect/length(vect)*3)
    for(i in seq_len(length(vect))){
      output[i]<-sum(
        vect[max(i-width*both.sides,1):min(i+width,length(vect))] * 
          normsum(window[ (max(i-width,1)-i):(min(i+width,length(vect))-i) +width+1])
        ,na.rm=T)
    }
  }
  output
}

#' Initiate an empty data frame
#'
#' @param namelist A character vector of column names.
#'
#' @return A data.frame with 0 rows.
#' @export
df.init<-function(namelist){
  setNames(data.frame(matrix(ncol = length(namelist), nrow = 0)), namelist)
}

#' Change classes of columns in a data.frame
#' @description \code{retype()} changes the class of specific columns; \code{retype_all()} changes the class of all columns of a given class.
#'
#' @param df a data frame
#' @param ... Unquoted column names, paired with the desired class, e.g. 
#' 
#' \code{age = numeric(), language = character()}
#'
#' @export
#'
#' @examples 
#' sapply(ToothGrowth,class)
#' NewToothGrowth <- retype(ToothGrowth, supp = character(), dose = factor())
#' sapply(NewToothGrowth,class)
retype<-function(df, ...){
  args<-list(...)

  varnames<-names(args)
  vartypes<-sapply(args, class)

  effcols<-names(df)[names(df) %in% varnames]

  for(effcol in effcols){
    df[,effcol]<-as(df[,effcol], vartypes[which(varnames==effcol)])
  }
  return(df)
}

#' @rdname retype
#' @param df A data.frame
#' @param from An empty vector of the class to convert from, or a string. Columns sharing the class of argument \code{from} will be converted to the class of argument \code{to}.
#' @param to An empty vector of the class to convert to, or a string. Columns sharing the class of argument \code{from} will be converted to the class of argument \code{to}.
#'
#' @export
#'
#' @examples 
#'
#' sapply(mtcars,class)
#' newmtcars <- retype_all(mtcars,from="numeric",to="character")
#' sapply(newmtcars,class)
retype_all<-function(df,from,to){
  for(i in which(sapply(df,class)==from)){
    df[[i]]<-as(df[[i]],to)
  }
  df
}

#' Verify variable types in bulk
#'
#' @param ... Named arguments, where the argument is the object to be checked and the name of the argument is the mode (numeric, list, character, etc)
#'
#' @return Returns true on success, causes error if not.
#' @export
#'
#' @examples
#' try(verify_types(character="test",numeric=0000,character=12345))
verify_types<-function(...){
  args<-list(...)
  call<-as.list(match.call()[-1])
  types<-unique(names(args))
  for(type in types){
    ids<-which(type == names(args))
    for(id in ids){
      if(!do.call(paste0("is.",type),list(args[[id]]))){
        stop("Variable ",as.character(call[[id]])," is not of type ",type)
      }
    }
  }
  return(T)
}


#' Read and merge all .csv files in a folder
#'
#' @param folder path to a folder
#' @param readfunc list of functions that will be used to read the files; if the first function fails, the second function will be used, etc.
#'
#' @return A data.frame containing all merged .csv files 
#' @export
read.csv.folder<-function(folder="./", readfunc=list(read.csv,read.csv2,read.table)){
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

#' Install packages if neccesary, then load them.
#' @param ... Unquoted names of packages to try loading, and if unable, install and load.
#'
#' @examples trypackages(stats,utils,compiler)
#' @export
trypackages<-function(...){
  packs<-args2strings(...)
  for(pack in packs){
    if(!require(pack,character.only=T)){
      install.packages(pack)
      require(pack,character.only=T)
    }
  }
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


#' Set column and row names of an object
#' These are convenience functions that return an object with its column or row names changed.
#' Use it in pipes.
#' 
#' @param x an object
#' @param names column or row names to be assigned to the object
#' 
#' @export
#' @examples 
#' setColNames(ToothGrowth,c("length","supplement","dosage"))
#' setRowNames(BOD,BOD$Time)
setColNames<-function(x,names){ colnames(x)<-names; return(x) }
#' @export
#' @rdname setColNames
setRowNames<-function(x,names){ rownames(x)<-names; return(x) }

#' Levenshtein distance
#' 
#' Counts the number of single character deletions, insertions, and substitutions 
#' that need to be performed to turn the source string into the target string.
#'
#' @param source,target Strings to be compared.
#'
#' @return The Levenshtein distance between the two strings.
#' @export
#'
#' @examples LevenshteinDistance("Yoghurt","Youtube")
LevenshteinDistance<-function(source,target){
  source<-strsplit(source,"")[[1]]
  target<-strsplit(target,"")[[1]]
  sl<-length(source)
  tl<-length(target)
  d<- matrix(nrow=sl+1,ncol=tl+1)

  d[,1]<-seq_len(sl+1)-1
  d[1,]<-seq_len(tl+1)-1
  
  for(i in seq_len(sl+1)[-1]){
    for(j in seq_len(tl+1)[-1]){
      d[i, j]<- min(
        d[i-1, j] + 1,
        d[i, j-1] + 1,
        d[i-1, j-1] + (source[i-1] == target[j-1])
      )
    }
  }
  return(d[sl+1,tl+1])
}

#' Split a character column into multiple values
#'
#' @param x a character vector to split into columns
#' @param sep a caracter separating the different values
#'
#' @return a \code{data.frame} of boolean values, with rows representing the unpacked
#' vector entries and columns indicating whether the specific value
#' @export
#'
#' @examples
#' unsplit<-c("flour;salt;baking soda;steak;sugar;water;sauce;vinegar",
#' "flour;sauce;mustard;salt;pepper;vinegar;baking soda;water;tomatoes;onion;steak")
#' splitColumn(unsplit)
splitColumn<-function(x,sep=";"){
  vals<-lapply(x,function(y){strsplit(y,sep)[[1]]})
  uniques<-unique(unlist(vals))
  idx<-t(sapply(vals,function(y){uniques %in% y}))
  colnames(idx)<-ifelse(is.na(uniques),"NA",uniques)
  return(as.data.frame(idx))
}


#' Merge Multiple Data Frames
#' 
#' This function makes calls to \code{merge()} to merge every other dataset 
#' with the one next to it, repeating until only one dataset remains. 
#'
#' @param x a list of data frames
#' @param ... all other arguments for \code{merge} can be provided here
#'
#' @return A single, merged \code{data.frame}
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' #generate test data
#' testlist<-list()
#' lsize<-50
#' for(i in 1:lsize){
#'   testlist[[i]]<-data.frame(key=sample(1:500,100),
#'                             junk=letters[sample(1:26,100,replace=T)])
#'   colnames(testlist[[i]])[2]<-paste0("info",i)
#' }
#' multimerge(testlist,by="key",all=T)
multimerge<-function(x,...){
  while(length(x)>1){
    out<-list()
    while(length(x)>0){
      if(length(x)>=2){
        out[[length(out)+1]]<-merge(x[[1]],x[[2]],...)
        x[[1]]<-NULL
        x[[1]]<-NULL
      }else{
        out[[length(out)+1]]<-x[[1]]
        x[[1]]<-NULL
      }
    }
    x<-out
  }
  return(x)
}


#' Divide a vector or list
#' 
#' Divide a vector or list into parts of (preferably) equal length.
#' Either the length or the number of the parts can be set.
#'
#' @param x the to-be-divided object
#' @param divs,divlen The number of divisions and the preferred length of divisions. 
#' One and only one of \code{divs} and \code{divlen} must be given.
#'
#' @return A list consisting of \code{x}, divided in parts.
#' @export
#'
#' @examples
#' DivideSeries(letters,divs=5)
#' DivideSeries(1:10,divlen=3)
DivideSeries<-function(x,divs,divlen){
  if(missing(divlen)){
    divlen<-(length(x)+1)/divs
    stopifnot(length(x)/divs>=1)
    mapply(a=ceiling(cumsum(c(0,rep(divlen,divs-1)))),
           b=ceiling(cumsum(rep(divlen,divs))-1),
           FUN=function(a,b){x[a:b]},SIMPLIFY=F)
  }else if(missing(divs)){
    reps<-ceiling(length(x)/divlen)
    mapply(a=cumsum(c(1,rep(divlen,reps -1))),
           b=replace(cumsum(rep(divlen,reps)),reps,length(x)),
           FUN=function(a,b){x[a:b]},SIMPLIFY=F)
  }
}

divide2<-function(x,chunksize=NULL,chunks=NULL){
  lx<-length(x)
  if(!is.null(chunksize)){
    chunks<-ceiling(lx/chunksize)
  }else if(!is.null(chunks)){
    chunksize<-ceiling(lx/chunks)
  }
  chvec<-rep(seq_len(chunks),each=chunksize,length.out=lx)
  newx<-unname(split(x,chvec))
  return(newx)
}

# divide(1:16,length=5)
# 
# DivideSeries(letters[1:16],divs=7)
# DivideSeries(1:16,divlen=7)



# Remove rows with OLs from data frame
removeOLs<-function(.tbl,olvars,groups=NULL){
  newtbl<-.tbl %>% group_by(across(all_of(groups))) %>% 
    filter(if_all(.cols=(olvars),.fns=function(x){abs(scale(x))<3})) %>% ungroup()
  message("Filtered ",nrow(.tbl)-nrow(newtbl)," rows")
  return(newtbl)
}

# Replace OLs with NA
maskOLs<-function(.tbl,olvars,groups=NULL){
  if(!is.null(groups)){
    groupvar<-interaction(.tbl[groups])
  }else{
    groupvar<-rep(1,nrow(.tbl))
  }
  for(olvar in olvars){
    key<-tapply(.tbl[[olvar]],groupvar,function(x) abs(x-mean(x,na.rm=T))/sd(x,na.rm=T) >3) %>% unlist() %>% which()
    .tbl[[olvar]][key]<-NA
    message("Masked ",length(key), " outliers from variable ",olvar)
  }
  return(.tbl)
}

# Remove OLs in a vector
vec.removeOLs<-function(x){
  excl<-which(abs(scale(x))>3)
  message("Excluding ",length(excl)," observations from vector")
  if(length(excl)>0)  return(x[-excl])
  else return(x)
}

# Formatted t-test comparing to zero
zerodiff<-function(x){
  x<-na.omit(x)
  tt<-t.test(x)
  cat("t (",tt$parameter,") = ",tt$statistic, ", p = ",tt$p.value,
      ", g = ",format(CohenD(x,correct=T),digits=5),", M = ",format(mean(x),digits=5),"\n",sep="")
}

# Formatted t-test
twodiff<-function(form,data,paired=F){
  term<-as.character(form)[3]
  outcome<-as.character(form)[2]
  data<-data[,c(term,outcome)]
  data<-na.omit(data)
  
  tt<-t.test(form,data,paired=paired)
  x<-data[[outcome]][data[[term]]==unique(data[[term]])[1]]
  y<-data[[outcome]][data[[term]]==unique(data[[term]])[2]]
  cat("t (",tt$parameter,") = ",tt$statistic, ", p = ",tt$p.value,
      ", g = ",format(CohenD(x=x,y=y,correct=T),digits=5),", Mdiff = ",format(mean(x)-mean(y),digits=5),"\n",sep="")
}

# Formatted t-test for two inputs
twodiff2<-function(x,y,paired=F){
  tt<-t.test(x=x,y=y,paired=paired)
  cat("t (",tt$parameter,") = ",tt$statistic, ", p = ",tt$p.value,
      ", g = ",format(CohenD(x=x,y=y,correct=T),digits=5),", Mdiff = ",format(mean(x)-mean(y),digits=5),"\n",sep="")
}

npr.zerodiff<-function(x){
  x<-na.omit(x)
  test<-coin::wilcoxsign_test(rep(0,length(x))~x,exact=T)
  cat("Z = ",format(test@statistic@teststatistic,digits=5),
      ", p = ",format(test@distribution@pvalue(test@statistic@teststatistic),digits=5),
      ", Mdiff = ",format(mean(x),digits=5),"\n",sep="")
}

# Formatted Wilcoxon test
npr.twodiff<-function(form,data){
  term<-as.character(form)[3]
  outcome<-as.character(form)[2]
  data<-data[,c(term,outcome)]
  data[[term]]<-factor(data[[term]])
  data<-na.omit(data)
  test<-coin::wilcox_test(form,data, distribution="exact",conf.int=T)
  
  x<-data[[outcome]][data[[term]]==unique(data[[term]])[1]]
  y<-data[[outcome]][data[[term]]==unique(data[[term]])[2]]
  
  cat("Z = ",format(test@statistic@teststatistic,digits=5),
      ", p = ",format(test@distribution@pvalue(test@statistic@teststatistic),digits=5),
      ", Mdiff = ",format(mean(x)-mean(y),digits=5),"\n",sep="")
}

# formatted correlation
print.cor<-function(x,y,method="pearson"){
  h<-cor.test(x,y,method=method)
  cat(method," r (",h$parameter,") = ",dropLeadingZero(format(h$estimate,digits=2)),
      ", p = ",dropLeadingZero(format(h$p.value,digits=3)),"\n",sep="")
}



getMetrics<-function(origdata,preddata){
  cm<-table(origdata,preddata)
  return(c(acc=acc<-sum(diag(cm))/sum(cm),
           chance_acc=chance_acc<-max(rowSums(cm))/sum(cm),
           kappa=(acc-chance_acc)/(1-chance_acc),
           ppv=ppv<-cm[2,2]/sum(cm[,2]),
           chance_ppv=chance_ppv<-sum(cm[2,])/sum(cm),
           kappa_ppv=(ppv-chance_ppv)/(1-chance_ppv),
           npv=npv<-cm[1,1]/sum(cm[,1]),
           chance_npv=chance_npv<-sum(cm[1,])/sum(cm),
           kappa_npv=(npv-chance_npv)/(1-chance_npv),
           sens=cm[2,2]/sum(cm[2,]),
           spec=cm[1,1]/sum(cm[1,])))
}


#' Multiple correlation
#' Computes the \href{https://en.wikipedia.org/wiki/Multiple_correlation}{multiple correlation coefficient}
#' of variables in \code{ymat} with the variable \code{x}
#' @param x Either a matrix of variables whose multiple correlation with each other is to be estimated; or a vector of which the multiple correlation with variables in \code{ymat} is to be estimated
#' @param ymat a matrix or data.frame of variables of which the multiple correlation with \code{x} is to be estimated
#' @param use optional character indicating how to handle missing values (see \link{cor})
#'
#' @return The multiple correlation coefficient
#' @export
#' @seealso https://www.personality-project.org/r/book/chapter5.pdf
#'
#' @examples
#' multiple.cor(mtcars[,1],mtcars[,2:4])
multiple.cor<-function(x,ymat,use="everything"){
  if(missing(ymat)){
    cv<-cor(x,use=use)
    corvec<-numeric(ncol(x))
    for(i in seq_along(corvec)){
      gfvec<-cv[(1:nrow(cv))[-i],i]
      dcm<-cv[(1:nrow(cv))[-i],(1:ncol(cv))[-i]]
      rsq<-t(gfvec) %*% solve(dcm) %*% gfvec
      corvec[i]<-sqrt(as.vector(rsq))
    }
    names(corvec)<-colnames(cv)
    return(corvec)
  }else{
    cv<-cor(cbind(x,ymat),use=use)
    gfvec<-cv[2:nrow(cv),1]
    dcm<-cv[2:nrow(cv),2:ncol(cv)]
    rsq<-t(gfvec) %*% solve(dcm) %*% gfvec
    return(sqrt(as.vector(rsq)))
  }
}


# Asymmetric CorMat with CorrCrunch ####
#' Create a Correlation Table
#'
#' @param df A data.frame. 
#' @param rowids,columnids character vectors containing column names from \code{df} that need to be correlated.  
#' @param rowdf,columndf data.frames whose columns need to be correlated. 
#' Either \code{df, rowids, & columnids} or \code{rowdf & columndf} are required.
#'
#' @return A formatted markdown table containing correlation coefficients, p-values, and 
#' the number and percentage of cases that need to be removed to flip the sign of each correlation coefficient.
#' @export
#'
#' @examples CorTable(mtcars,rowids=c("mpg","disp","hp"),columnids=c("drat","wt","qsec"))
#' 
#' CorTable(rowdf=mtcars[,c(1,3,4)],columndf=mtcars[,5:7])
CorTable<-function(df,rowids,columnids,rowdf,columndf){
  if(missing(df) | missing(rowids) | missing(columnids)){
    df<-cbind(rowdf,columndf)
    rowids<-colnames(rowdf)
    columnids<-colnames(columndf)
  }
  
  cormat<-matrix(NA,nrow=length(rowids),ncol=length(columnids))
  colnames(cormat)<-columnids
  rownames(cormat)<-rowids
  dfmat<-pmat<-hmat<-cormat
  
  for(i in rowids){
    for(j in columnids){
      corobj<-cor.test(df[,i],df[,j])
      cormat[i,j]<-corobj$estimate
      pmat[i,j]<-corobj$p.value
      hmat[i,j]<-CorrCrunch(df[,i],df[,j])$h
      dfmat[i,j]<-corobj$parameter
    }
  }
  
  outmat<-matrix("",ncol=length(columnids),nrow=length(rowids)*5-1)
  colnames(outmat)<-abbreviate(columnids)
  outrows<-rep("",length(rowids)*5-1)
  outrows[5*(0:(length(rowids)-1))+1]<-rowids
  rownames(outmat)<-outrows
  
  outmat[5*(0:(length(rowids)-1))+1,]<-gsub("0\\.","\\.",paste0("r=  ",format(cormat,nsmall=2)))
  outmat[5*(0:(length(rowids)-1))+2,]<-gsub("0\\.","\\.",paste0("p=  ",format(pmat,nsmall=3)))
  outmat[5*(0:(length(rowids)-1))+3,]<-gsub("0\\.","\\.",paste0("h=    ",format(hmat)))
  outmat[5*(0:(length(rowids)-1))+4,]<-gsub("0\\.","\\.",paste0("h/df=",format(hmat/dfmat,nsmall=2)))
  knitr::kable(outmat,digits=2,align="r")
}


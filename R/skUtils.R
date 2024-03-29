.onLoad<-function(libname, pkgname){
  packageStartupMessage("Thank you for loading skMisc v0.01")
}

#I don't want to import rlang, so it will be done this way instead.
args2strings <- function(...) sapply(substitute({ ... })[-1], deparse)

#' Create unique pairs
#' @description Combines vectors such that unique unordered sets are derived from the vectors' cross sections. 
#' @param ... two or more vectors of equal length
#'
#' @return a character vector consisting of all input vectors concatenated term-by-term and in alphabetic order.
#' @export
#'
#' @examples 
#' pair(1:4,4:1)
#' #[1] "1-4" "2-3" "2-3" "1-4"
pair<-function(...){
  args<-list(...)
  mat<-matrix(unlist(args),ncol=length(args))
  apply(mat,MARGIN=1,FUN=function(x){paste(sort(x),collapse="-")})
}

#' Crunch Outliers
#'
#' @param x Numeric vector to remove outliers from
#' @param DS A positive numeric value. If value exceeds this many standard deviations, it is counted as an outlier
#' @param hardlimit A numeric vector with two values. If set, values below the first value and above the second will be counted as outliers, and means/standard deviations will be computed from values within these bounds only.
#'
#' @return Vector with outlying values set to NA
#' @export
#' 
#' 
OLcrunch<-function(x,DS=3,rerun=F,hardlimit=NULL){
  if(!missing(hardlimit)){
    x[x<hardlimit[1] | x>hardlimit[2]]<-NA
  }
  m<-mean(x,na.rm=T)
  s<-sd(x,na.rm=T)
  x[(x>m+s*DS) | (x<m-s*DS)]<-NA
  if(rerun){
    while(any( abs(scale(x))>DS, na.omit=T )){
      x[which(abs(scale(x))>DS)]<-NA
    }
  }
  return(x)
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
#' #      len      supp      dose 
#' #"numeric"  "factor" "numeric" 
#' NewToothGrowth <- retype(ToothGrowth, supp = character(), dose = factor())
#' sapply(NewToothGrowth,class)
#' #      len        supp        dose 
#' #"numeric" "character"    "factor" 
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
#' #       mpg       cyl      disp        hp      drat        wt
#' # "numeric" "numeric" "numeric" "numeric" "numeric" "numeric"
#' #      qsec        vs        am      gear      carb 
#' # "numeric" "numeric" "numeric" "numeric" "numeric" 
#' 
#' newmtcars <- retype_all(mtcars,from="numeric",to="character")
#' sapply(newmtcars,class)
#' #         mpg         cyl        disp          hp        drat
#' # "character" "character" "character" "character" "character"
#' #          wt        qsec          vs          am        gear        carb 
#' # "character" "character" "character" "character" "character" "character" 
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


# CorrCrunch ####

#' Analyse the robustness of a correlation
#' @description \code{CorrCrunch()} computes the minimum number of cases that need to be removed from a dataset to flip the sign of a correlation coefficient.
#' This can be useful in distinguishing genuine correlations from spurious findings that hinge on one or two outliers.
#' Cases are removed iteratively; in each iteration the case that maximally shrinks the correlation coefficient is removed.
#' @param x,y Numeric vectors to correlate.
#' @param verbose if TRUE, prints verbose output.
#'
#' @return A list containing the number of cases that need to be removed to flip the sign of the correlation coefficient; 
#' the proportion removed cases in the data; and a data.frame without these cases.
#' @export
#'
#' @examples CorrCrunch(mtcars$mpg,mtcars$wt)
#' #Holdout needed to flip the sign: 19 (63.33%)
#' #Final r: 0.01181141
CorrCrunch<-function(x,y,verbose=F){
  iter<-0
  delidx<-0
  iterdf<-data.frame(x=x,y=y)
  iterdf<-iterdf[!is.na(rowSums(iterdf)),]
  origrval<-cor(iterdf$x,iterdf$y)
  rval<-origrval
  corsign<-sign(origrval)
  if(verbose){ cat(sep="","Holdout: ",0,", r: ",rval,"\n") }
  while(sign(rval)==corsign & iter<1000){
    iter<-iter+1
    for(i in 1:nrow(iterdf)){
      loopdf<-iterdf[-i,]
      if(corsign*cor(loopdf$x,loopdf$y)<corsign*rval){
        delidx<-i
        rval<-cor(loopdf$x,loopdf$y)
      }
    }
    if(verbose){ cat(sep="","Holdout: ",iter,", r: ",rval,"\n") }
    iterdf<-iterdf[-delidx,]
  }
  structure(.Data=list(h=iter,h.prop=iter/(length(x)-2),lastdf=iterdf),class="CorrCrunch")
}

print.CorrCrunch<-function(x){
  cat("Holdout needed to flip the sign: ",x$h,
      " (",round(x$h.prop*100,digits=2),"%)\n",sep="")
  cat("Final r: ",cor(x$lastdf$x,x$lastdf$y),"\n",sep="")
}
registerS3method("print", "CorrCrunch", print.CorrCrunch)


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
  
  outmat[5*(0:(length(rowids)-1))+1,]<-gsub("0\\.","\\.",paste0("r=  ",format(cormat,digits=0,nsmall=2)))
  outmat[5*(0:(length(rowids)-1))+2,]<-gsub("0\\.","\\.",paste0("p=  ",format(pmat,digits=0,nsmall=3)))
  outmat[5*(0:(length(rowids)-1))+3,]<-gsub("0\\.","\\.",paste0("h=    ",format(hmat,digits=0)))
  outmat[5*(0:(length(rowids)-1))+4,]<-gsub("0\\.","\\.",paste0("h/df=",format(hmat/dfmat,digits=0,nsmall=2)))
  knitr::kable(outmat,digits=2,align="r")
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

#' Compute column and row variances
#' @param x an input matrix of data.frame
#' @param na.rm Logical indicating whether NA values should be omitted before variance computation
#'
#' @export
#' @examples 
#' colVars(WorldPhones)
#' rowVars(WorldPhones)
colVars<-function(x,na.rm=T){
  apply(x,MARGIN=2,FUN=var,na.rm=T)
}
#' @export
#' @rdname colVars
rowVars<-function(x,na.rm=T){
  apply(x,MARGIN=1,FUN=var,na.rm=T)
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

#' Generate a matrix of combinations of values
#'
#' @param ... Character vectors, named or unnamed, or unquoted names of named arguments. 
#' Character vectors will be used to generate a matrix where each row represents a unique combination
#' of all values, akin to \code{expand.grid()}. Arguments which are unquoted names of named arguments
#' will become copies of the column generated by the eponymous named character vector. 
#'
#' @return A matrix.
#' @export
#'
#' @examples
#' hh<-c("a","b")
#' comboTable(a=letters[1:3], b=2,a,b,c=c("e","f"),d,c,d=hh,"huh",a,hh)
comboTable<-function(...){
  cl<-match.call(expand.dots = FALSE)[[2]]
  
  #deconstruct arguments
  named<- names(cl)!=""
  repeated<- as.character(cl) %in% names(cl)[named]
  unnamed<- names(cl)=="" & !repeated

  #ensure named args occur only once
  if(any(duped<-named & duplicated(names(cl)))){
    stop("When there are named arguments with identical names, every such argument",
         "after the first one should have no value associated with it. Violating args: ",
         paste(which(duped),collapse=", "))
  }
  
  #build args
  e.args<-as.list(vector(length=length(cl)))
  e.args[named]<-lapply(cl[named],eval)
  e.args[repeated]<-NA
  e.args[unnamed]<-lapply(cl[unnamed],eval)
  
  #get table
  mat<-unname(as.matrix(do.call(expand.grid,e.args)))
  
  #enter duplicates
  for(i in which(repeated)){
    col<-which(names(cl) %in% as.character(cl[i]) & named)
    mat[,i]<-mat[,col]
  }
  
  #return
  return(mat)
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


suppressor<-function(x,soft,hard, strength){
  hard<-hard-soft
  y<-x
  y[x< -soft]<- hard*(y[x< -soft]+soft)/(strength*abs(y[x< -soft]+soft)+hard)-soft
  y[x> soft]<- hard*(y[x> soft]-soft)/(strength*abs(y[x> soft]-soft)+hard)+soft
  return(y)
}


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

# ggplot theme appropriate for manuscripts
apatheme<-function(){
  theme_bw() + 
    theme(legend.position="bottom",panel.grid=element_blank(),
          panel.border = element_blank(),
          axis.line=element_line(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold",size=unit(14,"pt")),
          axis.title=element_text(face="bold",size=unit(14,"pt"),margin=unit(rep(0,4),"pt")),
          legend.text = element_text(size=unit(14,"pt")),
          axis.ticks.length = unit(-5,"pt"),
          axis.text=element_text(color="black",size=unit(13,"pt")),
          axis.text.x=element_text(margin=unit(c(8,0,0,0),"pt")),
          axis.text.y=element_text(margin=unit(c(0,8,0,0),"pt")))
}

# Remove leading zero from formatted numbers (taken from stackoverflow, and edited)
dropLeadingZero <- function(l){
  lnew <- c()
  for(i in l){
    if(isTRUE(i==0)){ #zeros stay zero
      lnew <- c(lnew,"0")
    } else if (isTRUE(i>=1) | isTRUE(i<=-1)){ #above one stays the same
      lnew <- c(lnew, as.character(i))
    } else
      lnew <- c(lnew, gsub("(?<![0-9])0+(?=\\.)", "", i, perl = TRUE))
  }
  as.character(lnew)
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
cor.print<-function(x,y,method="pearson"){
  h<-cor.test(x,y,method=method)
  cat(method," r (",h$parameter,") = ",dropLeadingZero(format(h$estimate,digits=2)),
      ", p = ",dropLeadingZero(format(h$p.value,digits=3)),"\n",sep="")
}



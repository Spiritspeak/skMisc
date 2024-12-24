# Files that are "up to standard": holdout.R and formatting.R
# Add a match.merge like function that just gives the matching indices
# Consider changing where.duplicated to which.unique() and add an option to limit it to only duplicated values


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
#' 
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

#' Negative %in%
#' Returns which values in the left vector are not in the right vector.
#' @param x Values whose presence will be checked for in \code{table}
#' @param table Values that will yield a \code{FALSE} if they exist in \code{x}
#'
#' @return A logical vector indicating whether each value of \code{x} lacks a match in \code{table}
#' @export
#'
#' @examples
#' (1:5) %nin% (1:3)
#' 
`%nin%` <- function(x,table){
  !(x %in% table)
}

#' Count duplicate values in a vector
#' 
#' which.duplidate() determines for each element of a vector 
#' how many times it has occurred so far.
#' It works similarly to [base::duplicated()] which only determines 
#' whether a value has occurred before and not how many times.
#' 
#' @param x A vector
#'
#' @return A vector of the same length as \code{x}, where each element represents
#' the number of times its value in \code{x} has been repeated so far.
#' @export
#'
#' @examples
#' which.duplicate(c(1,6,5,2,1,1,8,6,5))
#' 
which.duplicate<-function(x){
  vals<-unique(x)
  repvec<-numeric(length(vals))
  for(v in vals){
    repvec[x==v]<-seq_len(sum(x==v))
  }
  return(repvec)
}


#' Find the index of the first occurrence of each value
#' This replaces all values with the index of their first occurrence
#' 
#' @param x A vector
#' @param na.first Logical value indicating whether to replace the first
#' occurrence of each value with its index (if \code{FALSE}) or with NA (if \code{TRUE}).
#'
#' @return A vector of the length of \code{x} with each value representing either 
#' \code{NA} if it is the first occurrence of a unique value, or 
#' the index of the first instance of each value in \code{x} if it is a duplicated value
#' @export
#'
#' @examples
#' which.first(c("a","b","a","k"))
#' 
which.first <- function(x, na.first=FALSE){
  ux<-unique(x)
  newx<-rep(NA,length(x))
  if(na.first){
    for(u in ux){
      idvec<-which(x==u)
      newx[idvec[-1]]<-idvec[1]
    }
  }else{
    for(u in ux){
      idvec<-which(x==u)
      newx[idvec]<-idvec[1]
    }
  }
  return(newx)
}

unique.set<-function(x){
  x<-t(apply(x,1,sort))
  x<-apply(x,1,paste,collapse="-")
  duplicated(x)
}


#' Generate unique pairs or N-tuplets
#'
#' @param nval Number of values to arrange into unique tuplets
#' @param ntuplet N-tuplets to arrange the values uniquely into
#'
#' @return A matrix where each row is a unique N-tuplet
#' @export
#'
#' @examples
#' allpairs(nval=20,ntuplet=3)
#' 
allpairs<-function(nval,ntuplet=2){
  currmat<-matrix(seq_len(nval),ncol=1)
  for(tuple in seq_len(ntuplet)[-1]){
    newmats<-list()
    for(i in seq_len(nrow(currmat))){
      startval<-currmat[i,tuple-1]+1
      if(startval<=nval){
        itervec<-startval:nval
        newmats[[i]]<-do.call(cbind,c(as.list(currmat[i,]),list(itervec)))
      }
    }
    currmat<-do.call(rbind,newmats)
  }
  return(currmat)
}

#' Convert dates to Nth weekday of the month values
#' Computes which weekday of the month each date represents. 
#' E.g., for each "second monday of the month", it gives 2.
#' @param x Date(s) to convert
#'
#' @return A vector of Nth weekday of the month values
#' @export
#'
#' @examples
#' nthWeekdayOfMonth("2000-04-08")
#' 
nthWeekdayOfMonth<-function(x){
  x <- as.Date(x)
  out <- numeric(length(x))
  for(i in seq_along(x)){
    monthday <- as.numeric(format(x[i],'%d'))
    wdayvec <- weekdays(x[i]-0:monthday)
    out[i] <- sum(wdayvec[length(wdayvec)] == wdayvec)
  }
  return(out)
}

#' Quantize vector
#' Will replace vector values with their quantile.
#' 
#' @param x A numeric vector to quantize
#' @param n The number of quantiles to divide it into
#'
#' @return A vector with values replaced by their quantile membership
#' @export
#'
#' @examples
#' a<-rnorm(100)
#' quantize(a,5)
#' 
quantize<-function(x,n){
  quants<-quantile(x,seq_len(n-1)/n,na.rm=T)
  quants<-c(-Inf,quants,Inf)
  cut(x,breaks=quants,labels=seq_len(n))
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
#' subdivide(letters,divs=5)
#' subdivide(1:10,divlen=3)
#' 
subdivide<-function(x,divs,divlen){
  xl<-length(x)
  if(missing(divlen)){
    divlen<-(xl+1)/divs
    stopifnot(xl/divs>=1)
    cs<-ceiling(cumsum(rep(divlen,divs-1)))
    a<-c(0,cs)
    b<-c(cs-1,xl)
  }else if(missing(divs)){
    divs<-ceiling(xl/divlen)
    cs<-cumsum(c(1,rep(divlen,divs-1)))
    a<-cs
    b<-c(cs[-1],xl)
  }
  exli<-vector(mode="list",length=divs)
  for(i in seq_len(divs)){
    exli[[i]]<-x[a[i]:b[i]]
  }
  return(exli)
}

#' Turn a matrix into a long-format data.frame
#'
#' @param x A matrix
#'
#' @return a \code{data.frame} with three columns: \code{row} and \code{col} indicating the
#' row and column names, and \code{value} indicating the respective value in the matrix. 
#' If no row or column names are available, the row or column number is used instead.
#' @export
#'
#' @examples
#' mymatrix<-matrix(1:80,ncol=8,nrow=10)
#' unwrap.matrix(mymatrix)
#' 
#' carmatrix<-as.matrix(mtcars)
#' unwrap.matrix(carmatrix)
#' 
unwrap.matrix<-function(x){
  dn<-dimnames(x)
  unwrap<-expand.grid(row=if(is.null(dn[[1]])){ seq_len(nrow(x)) }else{ dn[[1]] },
                      col=if(is.null(dn[[2]])){ seq_len(ncol(x)) }else{ dn[[2]] },
                      stringsAsFactors=F)
  unwrap[["value"]]<-as.vector(x)
  return(unwrap)
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

#' Initiate an empty data frame
#'
#' @param x A character vector of column names.
#'
#' @return A data.frame with 0 rows.
#' @export
df.init<-function(x){
  setNames(data.frame(matrix(ncol = length(x), nrow = 0)), x)
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
#' 
setColNames<-function(x,names){ colnames(x)<-names; return(x) }
#' @export
#' @rdname setColNames
setRowNames<-function(x,names){ rownames(x)<-names; return(x) }

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


#' Change classes of columns in a data.frame
#' @description \code{retype()} changes the class of specific columns; 
#' \code{retype_all()} changes the class of all columns of a given class.
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
#' 
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
#' 
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

#' Create substrings with a maximal length by splitting at specific characters
#' 
#' This function splits a string into substrings of length \code{width} or shorter.
#' The splitting is done at the characters specified in \code{split}, in order of preference.
#' 
#' This combines the functionality of [base::strwrap()] and [base::strsplit()]; 
#' instead of a string wrapped with newlines, the result is multiple substrings.
#'
#' @param x A character vector of length 1.
#' @param width The maximum character length to break the vector at.
#' @param split A vector of regular expressions to match a character to break the string at.
#' The function will try to break the string at the first value specified in this argument;
#' if that fails, it will move on to the second, then the third, etc.
#'
#' @return A character vector consisting of strings of length \code{width} or shorter, 
#' and split at the characters specified in \code{split}.
#' 
#' @export
#' @md
#'
#' @examples
#' thanks <- paste(readLines(file.path(R.home("doc"), "THANKS")), collapse = "\n")
#' strsplit.wrap(thanks,width=80)
#' 
#' alphabet <- paste0(letters,collapse="")
#' strsplit.wrap(alphabet,width=3)
strsplit.wrap <- function(x, width=2000, split=c("\n"," ",",","")){
  output <- character()
  if(!any(split=="")){ split <- c(split,"") }
  while(nchar(x) > 0){
    cstr <- substr(x,1,width)
    if(nchar(cstr) < width){
      output[length(output)+1] <- x
      x <- ""
    }else{
      for(splitchar in split){
        
        if(nzchar(splitchar)){
          nls <- gregexpr(splitchar,cstr)[[1]]
          end <- nls[length(nls)]
          if(end!=-1){
            output[length(output)+1] <- 
              trimws(substr(cstr,1,end),whitespace=splitchar)
            x <- substr(x, end+1, nchar(x))
            break
          }
        }else{
          output[length(output)+1] <- substr(cstr,1,width)
          x <- substr(x, width+1, nchar(x))
        }
        
      }
    }
  }
  return(output)
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
#' 
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
#' @param sep a character separating the different values
#'
#' @return a \code{data.frame} of boolean values, with each row representing 
#' a value of x and each column representing a unique value 
#' in \code{x} following splitting. A column is marked TRUE in a specific row if
#' the value representing that column was present in that row.
#' 
#' @export
#'
#' @examples
#' unsplit<-c("flour;salt;baking soda;steak;sugar;water;sauce;vinegar",
#' "flour;sauce;mustard;salt;pepper;vinegar;baking soda;water;tomatoes;onion;steak")
#' vector2columns(unsplit)
#' 
vector2columns<-function(x, sep=";"){
  vals <- strsplit(x, sep)
  uniques <- unique(unlist(vals))
  idx <- t(sapply(vals, function(y){ uniques %in% y }))
  colnames(idx) <- ifelse(is.na(uniques), "NA", uniques)
  out <- as.data.frame(idx)
  return(out)
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
#' 
multimerge<-function(x, ...){
  while(length(x) > 1){
    out <- list()
    while(length(x) > 0){
      if(length(x) >= 2){
        out[[length(out) + 1]] <- merge(x[[1]], x[[2]], ...)
        x[[1]] <- NULL
        x[[1]] <- NULL
      }else{
        out[[length(out) + 1]] <- x[[1]]
        x[[1]] <- NULL
      }
    }
    x <- out
  }
  return(x)
}



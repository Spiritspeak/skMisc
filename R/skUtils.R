

#I don't want to import rlang, so it will be done this way instead.
args2strings <- function(...) sapply(substitute({ ... })[-1], deparse)

#' @name clamp
#' @title Clamp
#' @description Clamp a numeric vector between a minimum and maximum value.
#' 
#' @param x The vector/matrix to clamp.
#' @param minval Minimum value; all lower values are clamped to this value.
#' @param maxval Maximum value; all higher values are clamped to this value.
#' @param na.rm for \code{scale2range}, whether \code{NA}s should be removed before
#' applying the function (default \code{TRUE}. 
#' If \code{FALSE} and \code{x} contains \code{NA}s, then the
#' complete vector will be \code{NA}.
#'
#' @return Clamped vector.
#' @author Sercan Kahveci
#' @export 
#'
#' @examples 
#' clamp(0:10,2,8)
#' clamp0(rnorm(10))
#' scale2range(rnorm(10))
#' 
NULL

#' @describeIn clamp Set all values exceeding the minimum and maximum to 
#' the minimum and maximum value, respectively.
#' @export
clamp <- function(x, minval=-Inf, maxval=Inf){
  x[x < minval] <- minval
  x[x > maxval] <- maxval
  x
}

#' @describeIn clamp Same as \code{clamp())} but with a default range of 0 to 1.
#' @export
clamp0 <- function(x, minval=0, maxval=1){
  x[x < minval] <- minval
  x[x > maxval] <- maxval
  x
}

#' @describeIn clamp Rescales the vector such that it fits neatly between 
#' the given minimum and maximum values.
#' @export
scale2range <- function(x, minval=0, maxval=1, na.rm=TRUE){
  x <- x - min(x, na.rm=na.rm)
  x <- x / max(x, na.rm=na.rm) * maxval + minval
  x
}

#' Transform all data.frame columns into ordered IDs
#' 
#' Transforms every column into a factor and then into an integer.
#'
#' @param x A \code{data.frame}.
#'
#' @return \code{x} with each column replaced by integer IDs.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' cols2ids(mtcars)
#' 
cols2ids <- function(x){
  for(col in seq_along(x)){
    x[,col] <- as.numeric(as.factor(x[[col]]))
  }
  x
}

# TODO: add counting direction as argument

#' Count duplicate values in a vector
#' 
#' which.duplicate() determines for each element of a vector 
#' how many times its value has occurred so far.
#' 
#' It works similarly to [base::duplicated()] which only determines 
#' whether a value has occurred before and not how many times.
#' 
#' @param x A vector.
#'
#' @return A vector of the same length as \code{x}, where each element represents
#' the number of times the value in the same position in \code{x} has been repeated so far.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' which.duplicate(c(1,6,5,2,1,1,8,6,5))
#' 
which.duplicate <- function(x){
  vals <- unique(x)
  repvec <- numeric(length(vals))
  for(v in vals){
    repvec[x == v] <- seq_len(sum(x == v))
  }
  return(repvec)
}


#' Find the index of the first occurrence of each value
#' 
#' This replaces all values with the index of their first occurrence.
#' 
#' @param x A vector.
#' @param na.first Logical value indicating whether to replace the first
#' occurrence of each value with its index (if \code{FALSE}) or with NA (if \code{TRUE}).
#'
#' @return A vector of the length of \code{x} with each value representing either 
#' \code{NA} if it is the first occurrence of a unique value, or 
#' the index of the first instance of each value in \code{x} if it is a duplicated value.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' duplicateof(c("a","b","a","k"))
#' 
duplicateof <- function(x, na.first=FALSE){
  ux <- unique(x)
  newx <- rep(NA, length(x))
  if(na.first){
    for(u in ux){
      idvec <- which(x == u)
      newx[idvec[-1]] <- idvec[1]
    }
  }else{
    for(u in ux){
      idvec <- which(x == u)
      newx[idvec] <- idvec[1]
    }
  }
  return(newx)
}

# TODO: add option to ignore duplicate values

#' Set-based unique and duplicate detection
#' 
#' These functions are like [base::unique()] and [base::duplicated()] except they only look at
#' whether two list elements contain the same values - the order does not matter.
#'
#' @param x A list of vectors
#'
#' @return For \code{setunique()}, a list of unique sets. 
#' For \code{setduplicated()}, a logical vector indicating whether
#' a set occurred previously.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' mysets <- list(a=1:3,b=2,c=3:1,d=c(1,3))
#' setunique(mysets)
#' setduplicated(mysets)
#' 
setunique <- function(x){
  x[!duplicated(lapply(x, sort))]
}

#' @rdname setunique
#' @export
#' 
setduplicated <- function(x){
  duplicated(lapply(x, sort))
}


#' Produce consecutive ranks
#' 
#' Ranks a vector such that ranks are always 1 apart if they succeed each other
#' or 0 apart if they are tied. 
#'
#' @param x A numeric vector.
#'
#' @return An integer vector of the same length as \code{x}.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' testvec <- c(0,2,3,3,7,9,9,0)
#' cons.rank(testvec)
#' 
cons.rank <- function(x){
  ux <- unique(x)
  ur <- as.integer(rank(ux))
  out <- integer(length(x))
  names(out) <- names(x)
  for(i in seq_along(ux)){
    out[x == ux[i]] <- ur[i]
  }
  return(out)
}


#' Identify runs of identical values in a vector
#'
#' @param x A vector.
#'
#' @return An integer vector of the same length as \code{x}, 
#' where each value identifies the run to which
#' the same-position value in \code{x} belongs.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' myvec <- sample(letters[1:3],size=20,replace=TRUE)
#' cbind(myvec,runs(myvec))
#' 
runs <- function(x){
  ln <- length(x)
  runbound <- which(x[-1] != x[-ln])
  runend <- c(runbound,ln)
  runstart <- c(1,runbound+1)
  out <- integer(ln)
  for(i in seq_along(runstart)){
    out[runstart[i]:runend[i]] <- i
  }
  return(out)
}

#' Generate unique pairs or N-tuplets
#'
#' @param nval Number of values to arrange into unique tuplets.
#' @param ntuplet N-tuplets to arrange the values uniquely into.
#' @param incl.self Determines whether a value can be paired with itself.
#'
#' @return A matrix where each row is a unique N-tuplet.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' allpairs(nval=5,ntuplet=3)
#' 
allpairs <- function(nval, ntuplet=2, incl.self=FALSE){
  currmat <- matrix(seq_len(nval), ncol=1)
  offset <- ifelse(incl.self, 0, 1)
  for(tuple in seq_len(ntuplet)[-1]){
    newmats <- list()
    for(i in seq_len(nrow(currmat))){
      startval <- currmat[i,tuple-1] + offset
      if(startval <= nval){
        itervec <- startval:nval
        newmats[[i]] <- do.call(cbind, c(as.list(currmat[i,]), list(itervec)))
      }
    }
    currmat <- do.call(rbind, newmats)
  }
  return(currmat)
}

#' Interweave vectors
#' 
#' Interweave vectors, such that the result features the first value from vector A, 
#' then the first value from vector B, then the second value from vector A, etc.
#'
#' @param ... Any number of vectors of equal length or length 1.
#'
#' @return The input vectors interwoven into one.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' # Interweave 3 equal-length vectors 
#' interweave(1:3, 4:6, 7:9)
#' 
#' # Interweave vector with single value
#' interweave(1:5, 0)
#' 
interweave <- function(...){
  nvecs <- ...length()
  veclens <- sapply(seq_len(nvecs), \(x){length(...elt(x))})
  maxlen <- max(veclens)
  if(!all(veclens == 1 | veclens == maxlen)){
    stop("Inputs not of equal length")
  }
  if(nvecs==0 | maxlen==0){
    return(NULL)
  }
  out <- vector(length=nvecs * maxlen)
  for(i in seq_len(nvecs)){
    out[i+nvecs*(0:(maxlen-1))] <- ...elt(i)
  }
  return(out)
}

quicklag <- function(x){ c(NA, x[-length(x)]) }


#' Convert dates to Nth weekday of the month values
#' 
#' Computes which weekday of the month each date represents. 
#' E.g., for each "second monday of the month", it gives 2.
#' 
#' @param x Date(s) to convert.
#'
#' @return A vector of Nth weekday of the month values.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' nthWeekdayOfMonth("2000-04-08")
#' 
nthWeekdayOfMonth <- function(x){
  x <- as.Date(x)
  out <- numeric(length(x))
  for(i in seq_along(x)){
    monthday <- as.numeric(format(x[i],'%d'))
    wdayvec <- weekdays(x[i]-0:monthday)
    out[i] <- sum(wdayvec[length(wdayvec)] == wdayvec)
  }
  return(out)
}

#' Find quantiles of vector
#' 
#' Will replace vector values with their quantile.
#' 
#' @param x A numeric vector to replace with the quantile of its values.
#' @param n The number of quantiles to divide it into.
#'
#' @return A vector with values replaced by their quantile membership.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' a<-rnorm(100)
#' findQuantile(a,5)
#' 
findQuantile <- function(x, n){
  quants <- c(-Inf,quantile(x, seq_len(n-1)/n, na.rm=T), Inf)
  .bincode(x, breaks=quants)
}

#' Divide a vector or list
#' 
#' Divide a vector or list into parts of (preferably) equal length.
#' Either the length or the number of the parts can be set.
#'
#' @param x the to-be-divided object.
#' @param divs,divlen The number of divisions and the preferred length of divisions. 
#' One and only one of \code{divs} and \code{divlen} must be given.
#'
#' @return A list consisting of \code{x}, divided in parts.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' subdivide(letters, divs=5)
#' subdivide(1:10, divlen=3)
#' 
subdivide <- function(x, divs, divlen){
  xl <- length(x)
  if(missing(divlen)){
    divlen <- (xl+1)/divs
    stopifnot(xl/divs >= 1)
    cs <- ceiling(cumsum(rep(divlen, divs-1)))
    a <- c(0, cs)
    b <- c(cs-1, xl)
  }else if(missing(divs)){
    divs <- ceiling(xl/divlen)
    cs <- cumsum(c(1,rep(divlen, divs-1)))
    a <- cs
    b <- c(cs[-1], xl)
  }
  exli <- vector(mode="list",length=divs)
  for(i in seq_len(divs)){
    exli[[i]] <- x[a[i]:b[i]]
  }
  return(exli)
}

#' Assign numeric values to bundles summing up to a set value
#'
#' @param x A numeric vector to bundle up.
#' @param maxval The maximum value a bundle can have.
#' @param group.high What should be done if an individual value exceeds \code{maxval}?
#' \code{"own"} assigns it to its own bundle, \code{"na"} assigns it to \code{NA}, and 
#' \code{"error"} gives an error.
#' @param group.na What should be done if a value is \code{NA}? 
#' \code{"own"} assigns it to its own bundle, \code{"na"} assigns it to \code{NA}, 
#' \code{"same"} assigns it to the same bundle as the previous value, and 
#' \code{"error"} gives an error.
#'
#' @return An integer vector of the same length as \code{x}, where each value indicates
#' to which bundle the value in \code{x} has been assigned.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' testvec <- sample(c(NA, NA, 1:14))
#' t(cbind(testvec, bundle_up(testvec, maxval=12, group.high="own", group.na="own")))
#' 
#' bundle(testvec, maxval=10, group.high="na", group.na="na")
#' 
bundle <- function(x, 
                   maxval, 
                   group.high=c("own", "na", "error"),
                   group.na=c("own", "na", "same", "error")){
  group.high <- match.arg(group.high)
  group.na <- match.arg(group.na)
  cs <- 0
  grp <- 1
  out <- numeric(length(x))
  
  naval <- is.na(x)
  na.addna <- 0
  na.grpbump <- 2
  if(any(naval) & group.na=="error"){
    stop("NA values not permitted")
  }else if(group.na=="na"){
    na.addna <- NA
    na.grpbump <- 1
  }else if(group.na=="same"){
    x[naval] <- 0
    naval[naval] <- F
  }
  
  excess <- x > maxval
  ex.addna <- 0
  ex.grpbump <- 2
  if(any(excess) & group.high=="error"){
    stop("Value exceeds maximum")
  }else if(group.high=="na"){
    ex.addna <- NA
    ex.grpbump <- 1
  }
  
  for(i in seq_along(x)){
    if(naval[i]){
      grp <- grp + na.grpbump
      cs <- 0
      out[i] <- grp -1 + na.addna
    }else if(excess[i]){
      grp <- grp + ex.grpbump
      cs <- 0
      out[i] <- grp -1 + ex.addna
    }else{
      cs <- cs + x[i]
      if(cs > maxval){
        cs <- x[i]
        grp <- grp + 1
      }
      out[i] <- grp
    }
  }
  nonna <- which(!is.na(out))
  if(out[nonna[1]] != 1){ out <- out - out[nonna[1]] +1 }
  out
}

#' Initiate an empty data frame
#'
#' @param x A character vector of column names.
#' @param nrow Number of rows (defaults to 0).
#'
#' @return A data.frame filled with \code{NA}.
#' @author Sercan Kahveci
#' @export
#' @examples
#' test <- df.init(c("A","B","C"))
#' 
df.init <- function(x, nrow=0){
  data.frame(matrix(ncol = length(x), nrow = nrow, dimnames=list(NULL, x)))
}

#' Scale a vector
#' 
#' Like scale() but returns a vector and is faster.
#' 
#' @param x Numeric vector to standardize.
#'
#' @return Scaled numeric vector with mean of 0 and standard deviation of 1.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' vec.scale(1:10)
#' 
vec.scale <- function(x){
  xt <- na.omit(x)
  m <- mean.default(xt)
  (x-m) / sqrt( sum((xt-m)^2)/(length(xt)-1) )
}

#' Split a character column into multiple values
#'
#' @param x a character vector to split into columns.
#' @param sep a character separating the different values.
#'
#' @return a \code{data.frame} of boolean values, with each row representing 
#' a value of x and each column representing a unique value 
#' in \code{x} following splitting. A column is marked TRUE in a specific row if
#' the value representing that column was present in that row.
#' @author Sercan Kahveci 
#' @export
#'
#' @examples
#' unsplit<-c("flour;salt;baking soda;steak;sugar;water;sauce;vinegar",
#' "flour;sauce;mustard;salt;pepper;vinegar;baking soda;water;tomatoes;onion;steak")
#' vec2columns(unsplit)
#' 
vec2columns <- function(x, sep=";"){
  vals <- strsplit(x, sep)
  uniques <- unique(unlist(vals))
  idx <- t(sapply(vals, function(y){ uniques %in% y }))
  colnames(idx) <- ifelse(is.na(uniques), "NA", uniques)
  out <- as.data.frame(idx)
  return(out)
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

#' Levenshtein distance
#' 
#' Counts the number of single character deletions, insertions, and substitutions 
#' that need to be performed to turn the source string into the target string.
#'
#' @param x,y Strings to be compared.
#'
#' @return The Levenshtein distance between the two strings.
#' @author Sercan Kahveci
#' @export
#'
#' @examples 
#' LevenshteinDistance("sitting","kitten")
#' 
#' LevenshteinDistance(c("cheese","fish"),c("child","flan"))
#' 
LevenshteinDistance <- function(x,y){
  x <- strsplit(x,"")
  y <- strsplit(y,"")
  xl <- sapply(x,length)
  yl <- sapply(y,length)
  diffs <- integer(length(x))
  
  for(k in seq_along(diffs)){
    d <- matrix(NA,nrow=xl[k] + 1, ncol=yl[k] + 1)
  
    d[,1] <- seq_len(xl[k] + 1) - 1
    d[1,] <- seq_len(yl[k] + 1) - 1
    
    for(j in seq_len(yl[k] + 1)[-1]){
      for(i in seq_len(xl[k] + 1)[-1]){
        d[i, j] <- min(
          d[i-1, j] + 1,
          d[i, j-1] + 1,
          d[i-1, j-1] + (x[[k]][i-1] != y[[k]][j-1])
        )
      }
    }
    diffs[k] <- d[xl[k] + 1, yl[k] + 1]
  }
  return(diffs)
}

##########################################
# Ops unrelated to objects or statistics #
##########################################

#' Install packages if neccesary, then load them.
#' 
#' @param ... Unquoted names of packages to try loading, 
#' and if unable, install and load.
#'
#' @examples trypackages(stats,utils,compiler)
#' @export
#' @author Sercan Kahveci
#' 
trypackages <- function(...){
  packs <- args2strings(...)
  for(pack in packs){
    if(!do.call("require", list(pack, character.only=T))){
      install.packages(pack)
      do.call("require", list(pack, character.only=T))
    }
  }
}

#' Retry running a function until it succeeds
#' 
#' An expression is executed using [base::try()] and re-run until
#' it raises no more errors or until a maximum number of evaluations is reached.
#' In the latter case, it raises an error of its own.
#'
#' @param expr An expression to try to execute.
#' @param n Maximum number of times to try to run the expression.
#'
#' @return The output of the expression.
#' @details 
#' This was primarily developed to robustly query websites. Sometimes
#' a query fails or produces output which raises an error; 
#' this function issues a retry in such a case.
#' @author Sercan Kahveci 
#' @export
#'
#' @examples
#' k <- 0
#' retry({k <<- k+1; if(k<4){ stop()}})
#' 
#' 
retry <- function(expr, n=5){
  expr <- substitute(expr)
  outcome <- NULL
  for(k in seq_len(n)){
    if(k > 1){ message("Attempt: ",k) }
    outcome <- try(eval(expr))
    if(!any("try-error" == attributes(outcome)$class)){ break }
  }
  if(!any("try-error" == attributes(outcome)$class)){
    outcome
  }else{
    stop("Expression failed maximum number of attempts")
  }
}

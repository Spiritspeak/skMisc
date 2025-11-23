

#I don't want to import rlang, so it will be done this way instead.
args2strings <- function(...) sapply(substitute({ ... })[-1], deparse)


#' Semi-randomize a sequence
#' 
#' This algorithm first uses [base::sample()] to randomize the sequence and then 
#' repetitively eliminates excessively long subsequences of the same value 
#' by swapping individual values.
#'
#' @param x A vector to randomize.
#' @param maxrep Maximum permissible number of sequential occurrences of the same value.
#'
#' @returns Vector \code{x} in semirandomized order, with no more than \code{maxrep}
#' sequential occurrences of the same value.
#' @details
#' This gives an error if semirandomization is not possible, 
#' for example when one value in \code{x} is too numerous.
#' 
#' @export
#'
#' @examples
#' Semirandomize(rep(0:1,each=96),maxrep=2)
#' 
#' Semirandomize(rep(0:2,c(1,4,1)),maxrep=2)
#' 
Semirandomize <- function(x, maxrep){
  # Check if semirandom sequence is possible
  tx <- table(x)
  maxprop<-max(tx/sum(tx))
  if(maxprop > maxrep/(1+maxrep)){
    stop("Semirandom sequence with given parameters not possible. ",
         "One of the values is too abundant.")
  }
  
  # Do the randomization
  x <- sample(x)
  counter <- 0L
  while(TRUE){
    counter<-counter+1L
    # Detect runs with excessive length
    rlex <- rle(x)
    repsize <- rep(rlex$lengths,rlex$lengths)
    repids <- rep(seq_along(rlex$lengths),rlex$lengths)
    excessreps <- unique(repids[repsize>maxrep])
    
    # Quit if no runs with excessive length
    if(length(excessreps)<1L){
      break
    # Reshuffle if too many swaps occurred
    }else if(counter > length(x)*maxprop){
      counter <- 0L
      x <- sample(x)
    }
    # From each run, pick one value to swap
    swapids <- sapply(excessreps,\(x){
      sample(which(repids==x),size=1L)
    })
    swapvals <- unique(x[swapids])
    # Per unique value, swap each picked element with that value with 
    # a picked element of another value; if that is not possible, swap with
    # an element of value that was not picked
    for(i in swapvals){
      currswapids <- swapids[x[swapids]==i]
      currswapvals <- swapids[x[swapids]!=i][sample(seq_along(currswapids))]
      swapids<-swapids[x[swapids]!=i & !(swapids %in% currswapvals[!is.na(currswapvals)])]
      currswapvals[is.na(currswapvals)]<-sample(which(x!=i),sum(is.na(currswapvals)))
      x[currswapids] <- x[currswapvals]
      x[currswapvals] <- i
    }
  }
  return(x)
}


#' Quota-fair Sampler
#' 
#' This samples \code{size} times from \code{x} such that the proportions given by
#' \code{prob} are met _exactly_. When they cannot be met exactly, sampling occurs as
#' fairly as possible; in that case, only the observations that cannot be allocated fairly
#' are randomly sampled with probabilities such that many repetitions of the sampling
#' algorithm would yield the desired proportions. If no proportions are given,
#' all observations are sampled equally often.
#'
#' @param x A vector to sample from.
#' @param size The desired length of the output vector.
#' @param prob (Optional) With what proportion should each value 
#' in \code{x} occur in the output vector?
#'
#' @returns A vector of length \code{size} with elements fairly sampled from \code{x}.
#' @author Sercan Kahveci
#' @export
#' @md
#'
#' @examples
#' allocate(1:2,size=50)
#' 
#' test <- replicate(10000,allocate(1:5,size=23,prob=c(.1,.3,.2,.1,.3)))
#' table(test) / length(test)
#' 
allocate <- function(x, size, prob=1){
  if(length(prob)==length(x)){
    prob <- prob/sum(prob)
  }else if(length(prob)==1L){
    prob <- rep(1/length(x),length(x))
  }else{
    stop("Prob must be of length 1 or the same length as x")
  }
  xsizes <- size * prob
  floorsizes <- floor(xsizes)
  miscsizes <- xsizes-floorsizes
  newx <- rep(x,floorsizes)
  if(sum(miscsizes)>0L){
    addedx <- sample(x,size=round(sum(miscsizes)),prob=miscsizes/sum(miscsizes),replace=F)
    newx <- c(newx,addedx)
  }
  sample(newx)
}

#' Return most common element(s) of vector
#'
#' @param x A vector.
#' @param n Number of most common elements to return.
#'
#' @returns A vector of the most common elements of \code{x}, with length \code{n}.
#' @export
#'
#' @examples
#' mostcommon(c(1,2,2,3,7,6,4,4,4))
#' mostcommon(c(1,2,2,3,7,6,4,4,4), n=2)
#' 
mostcommon <- function(x,n=1){
  xnames <- as.character(x)
  xfreq <- table(x)
  xord <- order(xfreq,decreasing=TRUE)
  x[match(names(xfreq[xord][seq_len(n)]),xnames)]
}

#' Convert empty vectors to NA
#' 
#' This converts vectors of length 0 to \code{NA}.
#'
#' @param x A vector.
#'
#' @returns \code{NA} if \code{x} is length 0, else \code{x} itself.
#' @export
#'
#' @examples
#' NONE2NA(NULL)
#' NONE2NA(logical(0))
#' NONE2NA(1:5)
#' 
NONE2NA <- function(x){
  if(length(x) == 0){
    NA
  }else{
    x
  }
}

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

#' Cumulative occurrence count of each value in a vector
#' 
#' occurrence() determines for each element of a vector 
#' how many times its value has occurred so far.
#' 
#' It works similarly to [base::duplicated()] which only determines 
#' whether a value has occurred before and not how many times.
#' 
#' @param x A vector.
#' @param fromLast if \code{TRUE}, the last occurrence of a value will be counted as the first,
#' and any earlier occurrences will be counted as second, third, etc.
#'
#' @return A vector of the same length as \code{x}, where each element represents
#' the number of times the value in the same position in \code{x} has been repeated so far.
#' @author Sercan Kahveci
#' @seealso [firstoccurrence()]
#' @export
#'
#' @examples
#' occurrence(c(1,6,5,2,1,1,8,6,5))
#' 
occurrence <- function(x, fromLast=FALSE){
  vals <- unique(x)
  repvec <- numeric(length(vals))
  if(!fromLast){
    for(v in vals){
      repvec[x == v] <- seq_len(sum(x == v))
    }
  }else{
    for(v in vals){
      repvec[x == v] <- rev(seq_len(sum(x == v)))
    }
  }
  return(repvec)
}


#' Replace elements with the index of the first occurrence of their value
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
#' @seealso [occurrence()]
#' @export
#'
#' @examples
#' firstoccurrence(c("a","b","a","k"))
#' 
firstoccurrence <- function(x, na.first=FALSE){
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


#' Set-based unique and duplicate detection
#' 
#' These functions are like [base::unique()] and [base::duplicated()] except they only look at
#' whether two list elements contain the same values - the order does not matter.
#'
#' @param x A list of vectors.
#' @param unique.only Should duplicated values within sets be ignored? Defaults to \code{FALSE}.
#'
#' @return For \code{setunique()}, a list of unique sets. 
#' For \code{setduplicated()}, a logical vector indicating whether
#' a set occurred previously.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' mysets <- list(a=c(1,2,3),b=2,c=c(3,2,1),d=c(1,3),e=c(1,1,2,3,3))
#' setunique(mysets)
#' setunique(mysets, unique.only=TRUE)
#' 
#' setduplicated(mysets)
#' setduplicated(mysets, unique.only=TRUE)
#' 
setunique <- function(x, unique.only=FALSE){
  x[!duplicated(lapply(x, function(y){
    if(unique.only){ y <- unique(y) }
    sort(y)
  }))]
}

#' @rdname setunique
#' @export
#' 
setduplicated <- function(x, unique.only=FALSE){
  duplicated(lapply(x, function(y){
    if(unique.only){ y <- unique(y) }
    sort(y)
  }))
}


#' Rotate a specific value to the front of vectors
#' 
#' This function rotates vectors in a list such that the first value in each vector 
#' is the first value in \code{first} that occurs in that vector.
#' 
#' This can be used to identify cyclic vectors that are duplicated 
#' with \code{duplicated(alignchains(x))}.
#' 
#' @param x List of vectors or a single vector to align through rotation.
#' @param first Optional vector indicating which values would be preferred as first value, 
#' in order of preference.
#' If left blank, this is set to the unique values in \code{x} in order of occurrence.
#' @param matches.only Should vectors only be rotated if at least one value 
#' has a match in \code{first}? Default is \code{FALSE}. If set to \code{FALSE},
#' then vectors without a match will be rotated such that the first value is the first
#' that occurs in \code{x}.
#'
#' @returns If \code{x} was a list, this returns 
#' a list with all constituent vectors rotated in accordance with the algorithm.
#' If \code{x} was a vector, this returns 
#' a vector rotated in accordance with the algorithm.
#' @md 
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' vecs <- list(a=c(5,3,6,7,3,2),
#'              b=c(7,8,4,2,3),
#'              c=c(3,4,5,6,1,2),
#'              d=c(1,5,3),
#'              e=c(5,3,1),
#'              f=c(3,1,5))
#' alignvecs <- rotate2front(vecs,first=1:8)
#' duplicated(alignvecs)
#' 
#' rotate2front(list(1:3,letters[3:1],letters[1:3],3:1))
#' rotate2front(c(4,2,3,1),first=1)
#' 
rotate2front <- function(x, first=NULL, matches.only=FALSE){
  if(!matches.only || is.null(first)){
    first <- unique(c(first, unlist(x)))
  }
  if(!is.list(x) && is.vector(x)){
    x <- list(x)
    vec.out <- T
  }else{
    vec.out <- F
  }
  out <- vector(mode="list", length=length(x))
  for(i in seq_along(x)){
    preferredinitial <- first[which(first %in% x[[i]])[1]]
    if(!is.na(preferredinitial)){
      newinitial <- which(x[[i]] == preferredinitial)[1]
      out[[i]] <- x[[i]][c(newinitial:length(x[[i]]),
                           (1:newinitial)[-newinitial])]
    }else{
      out[[i]] <- x[[i]]
    }
  }
  if(vec.out){
    out <- out[[1]]
  }
  out
}


#' Produce consecutive ranks
#' 
#' Ranks values of a vector such that ranks are always 1 apart if they succeed each other
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
  match(x,sort(unique(x)))
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
  if(!all(veclens == 1L | veclens == maxlen)){
    stop("Inputs not of equal length")
  }
  if(nvecs==0L || maxlen==0L){
    return(NULL)
  }
  out <- vector(length=nvecs * maxlen)
  for(i in seq_len(nvecs)){
    out[i+nvecs*(0:(maxlen-1L))] <- ...elt(i)
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
#' @param maxval The maximum sum a bundle can have.
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
#' t(cbind(testvec, bundle(testvec, maxval=12, group.high="own", group.na="own")))
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
  if(any(naval) && group.na=="error"){
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
  if(any(excess) && group.high=="error"){
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

#' Standardize a vector
#' 
#' A method for scale() that takes and returns a vector.
#' 
#' @param x Numeric vector to standardize.
#' @param center whether centering should occur (default \code{TRUE}).
#' @param scale Whether scaling should occur (default \code{TRUE}).
#'
#' @return Numeric vector with mean of 0 (if centered) and standard deviation of 1 (if scaled).
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' scale(1:10)
#' 
scale.vector <- function(x, center=TRUE, scale=TRUE){
  xt <- na.omit(x)
  if(center){
    m <- mean.default(xt)
    x <- x-m
  }
  if(scale){
    x / sqrt( sum((xt-m)^2)/(length(xt)-1L) )
  }else{
    x
  }
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
LevenshteinDistance <- function(x, y){
  x <- strsplit(x,"")
  y <- strsplit(y,"")
  xl <- lengths(x)
  yl <- lengths(y)
  diffs <- integer(length(x))
  
  for(k in seq_along(diffs)){
    d <- matrix(NA,nrow=xl[k] + 1L, ncol=yl[k] + 1L)
  
    d[,1] <- seq_len(xl[k] + 1L) - 1L
    d[1,] <- seq_len(yl[k] + 1L) - 1L
    
    for(j in seq_len(yl[k] + 1L)[-1L]){
      for(i in seq_len(xl[k] + 1L)[-1L]){
        d[i, j] <- min(
          d[i-1L, j] + 1L,
          d[i, j-1L] + 1L,
          d[i-1L, j-1L] + (x[[k]][i-1L] != y[[k]][j-1L])
        )
      }
    }
    diffs[k] <- d[xl[k] + 1L, yl[k] + 1L]
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

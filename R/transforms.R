
# testdat<-rgamma(10,shape=2)*50
# skewness(testdat)
# e1071::skewness(testdat,type=2) # mathematically identical to this
skewness <- function(x){
  n <- length(x)
  x <- x-mean(x)
  sqrt(n*(n-1))/(n-2)*mean(x^3)/(mean(x^2)^(3/2))
}

# testdat<-rt(500,100)
# kurtosis(testdat)
# e1071::kurtosis(testdat,type=2) # mathematically identical to this
kurtosis <- function(x){
  n <- length(x)
  x <- x-mean(x)
  (n-1)/((n-2)*(n-3))*((n+1)*(mean(x^4)/mean(x^2)^2-3)+6)
}

#' Vector transformations
#' 
#' Normalizing transformations for skewed and leptokurtic variables
#'
#' @param x Vector to be transformed.
#' @param type Can be "none" (no transform), "log" (\code{log(x)}), "log1p" (\code{log(x+1)}), 
#' "sqrt" (\code{sqrt(x)}), "inv" (\code{1/x}), "atan" (\code{atan(x)}),
#' "asinh" (\code{asinh(x-mean(x)) + mean(x)}),
#' "logshift" (\code{log(x+shift)-log(shift)} if \code{shift>0}, else \code{log(x+shift)}),
#' "invshift" (\code{(1+shift)/(x+1+shift)}),
#' "asinhrate" (\code{asinh((x-mean(x)) * rate) / rate + mean(x)}),
#' "bcp" (Box-Cox transform from [car::bcPower()]), 
#' "bcnp" (Shifted Box-Cox transform from [car::bcnPower()]),
#' or "yj" (Yeo-Johnson transformation from [car::yjPower()]). Non-positive values in \code{x} are
#' only allowed in "none", "asinh", "asinhrate", "bcnp" and "yj".
#' @param par A parameter that determines the shape of the transformation. If no parameter is given, then the function finds
#' * for "logshift" and "invshift", the parameter that minimizes the skew
#' * for "asinhrate", the parameter that minimizes the excess kurtosis
#' * for "bcp", "bcnp", and "yj", the parameter that gives the most normally distributed values, using [[car::powerTransform]]
#' For \code{"logshift"} and \code{"invshift"}, this is the shift; 
#' For \code{"asinhrate"}, this is the rate, with higher values representing a stronger reduction of kurtosis;
#' for \code{"bcp"} and \code{"yj"}, this is the lambda parameter; 
#' for \code{"bcnp"}, this is two values: the lambda and gamma parameters.
#' For all other transformations, \code{par} is ignored.
#' 
#' @details
#' * For skewed distributions, consider "log", "sqrt", "inv", "atan", "logshift", "invshift", "bcp", "bcnp", or "yjp".
#' 
#' * For distributions with excess kurtosis, consider "asinh" or "asinhrate".
#' 
#' @md
#' @author Sercan Kahveci
#' @returns A transformed vector. 
#' Attribute \code{"transform.method"} indicates the used method, and 
#' \code{"transform.par"} indicates the transformation parameter(s).
#' @seealso [removeOLs()], [vec.removeOLs()]
#' @export
#'
#' @examples
#' # Normality comparison
#' set.seed(1)
#' u<-rgamma(1000,2)
#' logtrans<-trans(u,"log")
#' shapiro.test(logtrans) # Not normal
#' logshifttrans<-trans(u,"logshift")
#' shapiro.test(logshifttrans) # Not normal but less bad
#' bctrans<-trans(u,"bcp")
#' shapiro.test(bctrans) # Normal
#' 
#' # Shifted distribution
#' u<-rgamma(1000,2)+2
#' logshifttrans<-trans(u,"logshift")
#' shapiro.test(logshifttrans) # Not normal 
#' invshifttrans<-trans(u,"invshift")
#' shapiro.test(invshifttrans) # Not normal
#' 
#'
#' u<-rt(1000,2)+10 
#' yjtrans<-trans(u,"yjp")
#' shapiro.test(yjtrans) # Not normal; still leptokurtic
#' asinhtrans<-trans(u,"asinh")
#' shapiro.test(asinhtrans) # Normal!
#' 
trans <- function(x,type=c("none","log","log1p","sqrt","inv","atan","asinh",
                           "logshift","invshift","asinhrate","bcp","bcnp","yjp"),par=NULL){
  type<-match.arg(type)
  if(type=="none"){
    out<-x
  }else if(type=="log"){
    out<-log(x)
  }else if(type=="log1p"){
    out<-log1p(x)
  }else if(type=="sqrt"){
    out<-sqrt(x)
  }else if(type=="inv"){
    out<-1/x
  }else if(type=="atan"){
    out<-atan(x)
  }else if(type=="asinh"){
    out<-asinh(x-mean(x))+mean(x)
  }else if(type=="logshift"){
    logshift<-function(x,p){
      log(x+p)-if(p<=0){0}else{log(p)}
    }
    if(is.null(par)){
      minpar<-min(x)
      scorer<-function(p){
        trval<-logshift(na.omit(x),exp(p)-minpar)
        skewness(trval)^2
      }
      par<-exp(optim(par=0,fn=scorer,method="BFGS")$par)-minpar
    }
    out<-logshift(x,par[1])
  }else if(type=="invshift"){
    invshift<-function(x,p){
      (p+1)/(x+p+1)
    }
    if(is.null(par)){
      minpar<-min(x)
      scorer<-function(p){
        trval<-invshift(na.omit(x),exp(p)-minpar)
        skewness(trval)^2
      }
      par<-exp(optim(par=0,fn=scorer,method="BFGS")$par)-minpar
    }
    out<-invshift(x,par[1])
  }else if(type=="asinhrate"){
    asinhrate<-function(x,p){
      asinh((x-mean(x))*p)/p+mean(x)
    }
    if(is.null(par)){
      scorer<-function(p){
        trval<-asinhrate(na.omit(x),exp(p))
        kurtosis(trval)^2
      }
      par<-exp(optim(par=0,fn=scorer,method="BFGS")$par)
    }
    out<-asinhrate(x,par[1])
  }else if(type=="bcp"){
    if(is.null(par)){
      par<-car::powerTransform(x,family="bcPower")$lambda
    }
    out<-car::bcPower(U=x,lambda=par[1])
  }else if(type=="bcnp"){
    if(is.null(par)){
      res<-car::powerTransform(x,family="bcnPower")
      par<-c(res$lambda,res$gamma)
    }
    out<-car::bcnPower(U=x,lambda=par[1],gamma=par[2])
  }else if(type=="yjp"){
    if(is.null(par)){
      par<-car::powerTransform(x,family="yjPower")$lambda
    }
    out<-car::yjPower(U=x,lambda=par[1])
  }
  attr(out,"transform.method")<-type
  attr(out,"transform.par")<-par
  return(out)
}


#' Remove outlying observations from a data.frame or matrix
#' 
#' Outliers are defined as values deviating more than X standard deviations (SDs) from the mean.
#' 
#' @param .tbl A \code{data.frame} or matrix to exclude outliers from
#' @param olvars Names or numeric index of the variables to detect outliers in. 
#' If \code{NULL}, all variables will be checked for outliers.
#' @param groups (optional) name or numeric index of the variable identifying groups of observations; 
#' outlier detection will be performed separately per group.
#' @param s If a value deviates more SDs from the mean than this value, it is marked as an outlier
#' @param make.na If \code{FALSE} (default), excludes all rows that have an outlier in
#' at least one variable in \code{olvars} (listwise).
#' If \code{TRUE}, the function instead turns the individual outlying values into \code{NA},
#' and does not exclude any rows.
#' @param verbose Report in console how many values have been removed per variable?
#' 
#' @details This does not detect any outliers in groups with less than 3 non-NA observations.
#'
#' @return The input \code{data.frame} or matrix with outliers excluded.
#' @author Sercan Kahveci
#' @export
#' @seealso [vec.removeOLs()] for the same outlier exclusion applied to a single vector. 
#' [trans()] for normalizing skewed or leptokurtic distributions.
#'
#' @examples
#' # Standard deviation limits can be set with argument s
#' removeOLs(mtcars, olvars=c("mpg", "disp", "hp"))
#' removeOLs(mtcars, olvars=c("mpg", "disp", "hp"), s=1)
#' 
#' # Replace OLs with NA with argument make.na
#' testdata <- mtcars
#' testdata$mpg[1] <- 40
#' testdata$hp[2] <- 500
#' removeOLs(testdata, olvars=c("mpg", "disp", "hp"), groups="vs", make.na=TRUE)
#' 
#' # Also works on matrices
#' testmat <- matrix(rnorm(1000), ncol=5)
#' testmat[cbind(sample(1:200,5),1:5)]<-1000
#' removeOLs(testmat)
#' 
removeOLs <- function(.tbl, olvars=NULL, groups=NULL, s=3, make.na=FALSE, verbose=TRUE){
  if(is.null(olvars)){ # If no olvars are given, evaluate all columns except those given by the group argument
    olvars <- setdiff(colnames(.tbl),groups)
    if(is.null(olvars)){ # same procedure but for when no colnames exist
      olvars <- setdiff(seq_len(NCOL(.tbl)),groups)
    }
  }
  # Error if not all olvars are numeric
  stopifnot(all(sapply(olvars,function(x)is.numeric(.tbl[,x]))))
  # determine olvars colname type
  if(is.numeric(olvars)){ 
    stopifnot(all(olvars <= NCOL(.tbl)))
    coltype <- "col. "
  }else{
    stopifnot(all(olvars %in% colnames(.tbl)))
    coltype <- ""
  }
  # Create single crossed groupvar
  if(!is.null(groups)){
    groupvar <- interaction(.tbl[,groups])
  }else{
    groupvar <- rep(1, NROW(.tbl))
  }
  # For each of olvars, determine which values are outlying, and NA-ify if applicable
  keylist <- list()
  for(olvar in olvars){
    key <- ave(x=.tbl[,olvar,drop=T], 
               groupvar,
               FUN=function(x){ 
                 if(sum(!is.na(x)) > 2){
                   abs(scale.vector(x)) > s
                 }else{
                   rep_len(NA, length(x))
                 }
               }) |> as.logical() |> which()
    if(make.na){
      .tbl[key,olvar] <- NA
      if(verbose){ message("Masked ", length(key), " outliers from ",coltype, olvar) }
    }else{
      keylist[[olvar]] <- key
    }
  }
  # If entire rows are to be removed, then do that
  if(!make.na){
    # Remove outlying rows
    keys <- unique(unlist(keylist))
    if(length(keys) > 0){
      .tbl <- .tbl[-keys, ]
    }
    
    # Report what you did
    if(verbose){ 
      nols<-sapply(keylist,length)
      endstr <- paste0(nols[1]," outliers in ",coltype, olvars[1])
      if(length(olvars) > 1){
        laststr <- paste0("and ",nols[length(nols)]," in ", coltype, olvars[length(olvars)])
        if(length(olvars) == 2){
          endstr <- paste(endstr, laststr)
        }else{
          endstr <- paste(endstr,
                          paste0(nols[-c(1,length(nols))]," in ", coltype, 
                                 olvars[-c(1,length(olvars))], collapse=", "),
                          laststr,
                          sep=", ")
        }
      }
      message("Excluded ", length(keys), " rows, due to ", endstr)
    }
  }
  return(.tbl)
}


#' Remove outlying observations from a vector
#'
#' @param x Vector to remove outliers from
#' @param s If a value deviates more SDs from the mean than this value, it is marked as an outlier
#' @param make.na If \code{FALSE}, excludes the outliers. 
#' If \code{TRUE}, replaces them with \code{NA}.
#' @param verbose Report in console how many values have been removed?
#'
#' @return A vector with outliers removed or replaced with \code{NA}.
#' @author Sercan Kahveci
#' @export
#' @seealso [removeOLs()]
#'
#' @examples
#' testvec <- c(1,3,5,7,9,11,13,15,17,19,100000)
#' vec.removeOLs(testvec)
#' vec.removeOLs(testvec,make.na=TRUE)
#' 
vec.removeOLs <- function(x, s=3, make.na=FALSE, verbose=TRUE){
  excl <- which(abs(scale.vector(x)) > s)
  if(verbose){ message("Excluded ", length(excl), " observations from vector") }
  if(length(excl) > 0){
    if(make.na){
      x[excl]<-NA
      return(x)
    }else{
      return(x[-excl])
    }
  }else{ 
    return(x)
  }
}





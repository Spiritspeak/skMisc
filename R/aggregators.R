

hedgesg <- function(x,y=NULL,paired=FALSE){
  if(is.null(y) | (paired & !is.null(y))){
    if(paired){
      x <- x-y
    }
    mean(x)/sd(x)
  }else if(!is.null(y) & !paired){
    dfx <- length(x)
    dfy <- length(y)
    (mean(x)-mean(y)) / 
      sqrt((dfx*var(x)+dfy*var(y))/(dfx+dfy))
  }
}


#' Standard error
#' 
#' Compute the standard error of the mean.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed?
#'
#' @return The standard error of \code{x}. Returns \code{NA} for
#' vectors with non-removed missing values, 
#' as well as those with 1 or less non-\code{NA} values.
#' @author Sercan Kahveci#' 
#' @export
#'
#' @examples
#' serr(rnorm(100))
#' 
serr <- function(x, na.rm=FALSE){
  sd(x, na.rm=na.rm) / sqrt(if(na.rm) sum(!is.na(x)) else length(x))
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

#' Tukey's Trimean
#' 
#' A robust mean estimator more efficient than the median, first proposed by Tukey.
#'
#' @param x A numeric vector
#' @param na.rm Whether to remove missing values. 
#' The function will error if \code{NA}s are present while \code{FALSE}.
#'
#' @return The trimean of \code{x}.
#' @export
#'
#' @examples
#' a <- rnorm(100)
#' a[1] <- 1000
#' trimean(a)
#' mean(a)
#' 
trimean <- function(x, na.rm=FALSE){
  if(na.rm){ x <- x[!is.na(x)] }
  sum(c(1,2,1) * quantile(x,c(.25,.5,.75)), na.rm=na.rm)/4
}

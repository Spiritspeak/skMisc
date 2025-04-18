

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

#' Title
#' 
#' @param x A numeric vector of modular quantities to be averaged.
#' @param mod The modulus of \code{x}, i.e. where the numbers "loop back" to zero. 
#' On a clock this is 12, in degrees it is 360, etc.
#' @param check Should a check be performed to detect whether there is no identifiable mean?
#' This should be set to \code{FALSE} when such a situation is not expected in the data, since
#' this check can be computationally costly.
#' @param na.rm Should \code{NA} values be removed? Defaults to \code{TRUE}.
#' 
#' @details
#' This function finds the value that has the smallest squared distance 
#' to all values in \code{x}.
#' 
#' @returns
#' @export
#'
#' @examples
#' 
#' 
#' # Situations where no modular mean exists (points are evenly divided around the circle)
#' modular.mean(c(0,3,6,9),mod=12)
#' modular.mean(c(0,1,3,4,6,7,9,10),mod=12)
#' 
modular.mean <- function(x, mod, check=TRUE, na.rm=TRUE){
  nas <- is.na(x)
  if(any(nas)){
    if(na.rm){
      x <- x[!nas]
    }else{
      return(NA)
    }
  }
  x <- x %% mod
  
  # Detect if the mean location on a circle is in the middle
  if(length(x) > 1 && check && length(unique(x))>1){
    if(abs(mean(cos(x/mod*2*pi)))<.Machine$double.eps && 
       abs(mean(sin(x/mod*2*pi)))<.Machine$double.eps){
      return(NA)
    }
  }
  scorefun <- function(par){
    scv <- abs(par - x)
    scv <- pmin(scv, mod - scv)
    sum(scv ^ 2)
  }
  est <- optimize(f=scorefun, 
                  interval=c(0, mod - .Machine$double.eps),
                  tol=sqrt(.Machine$double.eps))$minimum
  return(est)
}

#' Get classification accuracy metrics
#'
#' @param x The true class memberships
#' @param y The predicted class memberships
#'
#' @md
#' @author Sercan Kahveci
#' @return A vector with the following accuracy metrics:
#' * Accuracy: \code{acc}
#' * Chance accuracy: \code{chance_acc}
#' * Cohen's Kappa: \code{kappa}
#' * Positive predictive value: \code{ppv}
#' * Chance positive predictive value: \code{chance_ppv}
#' * Kappa of the positive predictive value: \code{kappa_ppv}
#' * Negative predictive value: \code{npv}
#' * Chance negative predictive value: \code{chance_npv}
#' * Kappa of the negative predictive value: \code{kappa_npv}
#' * Sensitivity: \code{sens}
#' * Specificity: \code{spec}
#' 
#' @export
#'
#' @author Sercan Kahveci
#' @examples
#' mymod <- glm(am ~ cyl + disp + hp,
#'              family="binomial", 
#'              data = mtcars)
#' 
#' newpred <- as.numeric(predict(mymod, 
#'                               data=mtcars,
#'                               type="response") > 0.5)
#' 
#' ClassMetrics(x=mtcars$am, y=newpred)
#' 
ClassMetrics <- function(x, y){
  cm <- table(x, y)
  return(c(acc=acc <- sum(diag(cm)) / sum(cm),
           chance_acc=chance_acc <- max(rowSums(cm)) / sum(cm),
           kappa=(acc-chance_acc) / (1 - chance_acc),
           ppv=ppv <- cm[2,2] / sum(cm[,2]),
           chance_ppv=chance_ppv <- sum(cm[2,]) / sum(cm),
           kappa_ppv=(ppv - chance_ppv) / (1 - chance_ppv),
           npv=npv <- cm[1,1] / sum(cm[,1]),
           chance_npv=chance_npv <- sum(cm[1,]) / sum(cm),
           kappa_npv=(npv - chance_npv) / (1 - chance_npv),
           sens=cm[2,2] / sum(cm[2,]),
           spec=cm[1,1] / sum(cm[1,]) ))
}

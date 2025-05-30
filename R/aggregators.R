

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


#' @name std.errors
#' @title Standard error estimators
#' 
#' @description Compute the standard errors (SE) of various estimators, 
#' algebraically or through bootstrapping.
#'
#' @param x A numeric vector.
#' @param boot If 0, SE will be determined algebraically. Otherwise,
#' this specifies the number of bootstraps used to compute the SE.
#' @param na.rm Logical. Should missing values be removed?
#'
#' @return The standard error of \code{x}. Returns \code{NA} for
#' vectors with non-removed missing values, 
#' as well as those with 1 or less non-\code{NA} values.
#' @note 
#' The SEs of the skewness and kurtosis are computed for 
#' the second skewness and kurtosis formulas reported by 
#' Joannes and Gill (1998).
#' 
#' These estimators all assume a normal distribution. 
#' Results may be inaccurate if the data is distributed differently. 
#' In such occasions, bootstrapping should be preferred.
#' 

#' @references Harding, B., Tremblay, C., & Cousineau, D. (2014). 
#' Standard errors: A review and evaluation of standard error estimators 
#' using Monte Carlo simulations. 
#' The Quantitative Methods for Psychology, 10(2), 107-123. 
#' doi: 10.20982/tqmp.10.2.p107
#' 
#' D. N. Joanes and C. A. Gill (1998). 
#' Comparing measures of sample skewness and kurtosis. 
#' The Statistician, 47, 183–189.
#' 
#' @export
#'
#' @examples
#' serr(rnorm(100))
#' 
#' serr.sd(rnorm(100))
#' 
#' serr.var(rnorm(100))
#' 
NULL

#' @describeIn std.errors SE of the mean
#' @export
#' 
serr <- function(x, boot=0, na.rm=FALSE){
  if(boot==0){
    sd(x, na.rm=na.rm) / sqrt(if(na.rm) sum(!is.na(x)) else length(x))
  }else{
    if(na.rm){ x <- x[!is.na(x)]}
    sd(replicate(boot, mean(sample(x, replace=TRUE))))
  }
}

#' @describeIn std.errors SE of the standard deviation
#' @export
#' 
serr.sd <- function(x, boot=0, na.rm=FALSE){
  if(boot==0){
    sd(x, na.rm=na.rm) / sqrt(2 * ifelse(na.rm, sum(!is.na(x)),length(x)) - 2)
  }else{
    if(na.rm){ x <- x[!is.na(x)]}
    sd(replicate(boot, sd(sample(x, replace=TRUE))))
  }
}

#' @describeIn std.errors SE of the variance
#' @export
#' 
serr.var <- function(x, boot=0, na.rm=FALSE){
  if(boot==0){
    var(x, na.rm=na.rm) * sqrt(2 / (ifelse(na.rm, sum(!is.na(x)), length(x)) - 1))
  }else{
    if(na.rm){ x <- x[!is.na(x)]}
    sd(replicate(boot, var(sample(x, replace=TRUE))))
  }
}

#' @describeIn std.errors SE of the skewness
#' @export
#'
serr.skewness <- function(x, na.rm=FALSE){
  n <- ifelse(na.rm, sum(is.na(x)), length(x))
  sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)))
}

#' @describeIn std.errors SE of the kurtosis
#' @export
#'
serr.kurtosis <- function(x, na.rm=FALSE){
  n <- ifelse(na.rm, sum(is.na(x)), length(x))
  2 * sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)))*
    sqrt((n^2-1)/((n-3)*(n+5)))
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
  if(check && length(unique(x))>1){
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


circular.mean <- function(x, mod=2*pi, check=TRUE, na.rm=TRUE){
  nas <- is.na(x)
  if(any(nas)){
    if(na.rm){
      x <- x[!nas]
    }else{
      return(NA)
    }
  }
  
  xtrans <- x/mod*(2*pi)
  xsin <- mean.default(sin(xtrans))
  xcos <- mean.default(cos(xtrans))
  if(check && length(x)>1){
    if(abs(xsin)<.Machine$double.eps && abs(xcos)<.Machine$double.eps){
      return(NA)
    }
  }
  return(atan2(xsin,xcos))
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

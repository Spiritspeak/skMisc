

#' Hedges' g
#' 
#' A measure of effect size similar to Cohen's d with a more accurate 
#' pooled standard deviation of the contrasted groups.
#'
#' @param x,y Numeric vectors. \code{y} can be omitted if 
#' the effect size should be computed for the difference from zero.
#' @param paired Are the values in \code{x} and \code{y} paired? 
#' Then Hedges' g will be computed for \code{x-y}.
#'
#' @return Hedges' g, positive if x > y.
#' @export
#'
#' @examples
#' # Difference from zero
#' hedgesg(x=rnorm(100,m=3))
#' 
#' # Difference from another group
#' hedgesg(x=rnorm(100,m=0),y=rnorm(100,m=3))
#' 
#' # Difference from paired value
#' a <- rnorm(100,m=3)
#' b <- a+rnorm(100,m=3)
#' hedgesg(x=a,y=b)
#' 
hedgesg <- function(x, y=NULL, paired=FALSE){
  if(paired){
    if(!is.null(y)){
      x <- x-y
      mean(x)/sd(x)
    }else{
      stop("No y given to compute a paired effect size for.")
    }
  }else{
    if(is.null(y)){
      mean(x)/sd(x)
    }else{
      dfx <- length(x)
      dfy <- length(y)
      (mean(x)-mean(y)) / 
        sqrt((dfx*var(x)+dfy*var(y))/(dfx+dfy))
    }
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
#' The algebraic estimators all assume a normal distribution. 
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


# TODO: add formula in latex

#' Mean or median of a modular quantity
#' 
#' @param x A numeric vector of modular quantities to be aggregated
#' @param mod The modulus of \code{x}, i.e. where the numbers "loop back" to zero. 
#' On a clock this is 12, in degrees it is 360, etc.
#' @param check Should a check be performed to detect whether there is 
#' no identifiable mean or median?
#' Set to \code{FALSE} when this is not expected to occur, since this check can be costly.
#' @param na.rm Should \code{NA} values be removed? Defaults to \code{TRUE}.
#' 
#' @details
#' These functions should not be confused with the [circular.mean()], which is computed
#' in a different way. 
#' 
#' @returns \code{modular.mean()} and \code{modular.median()} respectively find 
#' the value with the smallest squared or absolute distance 
#' to all values in \code{x}, like the regular mean and median do.
#' @export
#'
#' @examples
#' modular.mean(c(0,1,2,10),mod=12)
#' 
#' # Situations where no modular mean exists (points are evenly divided around the circle)
#' modular.mean(c(0,3,6,9),mod=12)
#' modular.mean(c(0,1,3,4,6,7,9,10),mod=12)
#' 
modular.mean <- function(x, mod, check=TRUE, na.rm=FALSE){
  scorefun <- function(par){
    scv <- abs(par - x)
    scv <- pmin(scv, mod - scv)
    sum(scv ^ 2)
  }
  compute.modular.aggregate(x, mod, check, na.rm, scorefun=scorefun)
}

#' @rdname modular.mean
#' @export
#' 
modular.median <- function(x, mod, check=TRUE, na.rm=FALSE){
  scorefun <- function(par){
    scv <- abs(par - x)
    scv <- pmin(scv, mod - scv)
    sum(scv)
  }
  compute.modular.aggregate(x, mod, check, na.rm, scorefun=scorefun)
}

compute.modular.aggregate <- function(x, mod, check, na.rm, scorefun){
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
  if(check){
    if(length(unique(x)) > 1){
      if(abs(mean(cos(x/mod*2*pi))) < .Machine$double.eps && 
         abs(mean(sin(x/mod*2*pi))) < .Machine$double.eps){
        return(NA)
      }
    }
  }
  est <- optimize(f=scorefun, 
                  interval=c(0, mod - .Machine$double.eps),
                  tol=sqrt(.Machine$double.eps))$minimum
  return(est)
}

# TODO: add formula in latex

#' Circular mean
#'
#' @param x A numeric vector.
#' @param mod The modulus of \code{x}, i.e. where the numbers "loop back" to zero. 
#' On a clock this is 12, in degrees it is 360, etc.
#' @param check Should a check be performed to detect whether there is 
#' no identifiable mean or median?
#' @param na.rm 
#'
#' @note This should not be confused with the [modular.mean()], which is computed differently.
#' @return
#' @export
#'
#' @examples
#' # Bedtimes for the past week
#' bedtimes <- c(23,23,00,01,23,02,03)
#' circular.mean(bedtimes,mod=24)
#' 
#' # Intuitively, these may require different outcomes, but they have the same circular mean
#' circular.mean(c(7,11,12),mod=12)
#' circular.mean(c(8,11,12),mod=12)
#' 
#' # No identifiable circular mean
#' circular.mean(c(0,pi),mod=2*pi)
#' 
circular.mean <- function(x, mod=2*pi, check=TRUE, na.rm=FALSE){
  nas <- is.na(x)
  if(any(nas)){
    if(na.rm){
      x <- x[!nas]
    }else{
      return(NA)
    }
  }
  
  xscaled <- x/mod*(2*pi)
  xsin <- mean.default(sin(xscaled))
  xcos <- mean.default(cos(xscaled))
  if(check && length(x)>1){
    if(abs(xsin) < .Machine$double.eps && abs(xcos) < .Machine$double.eps){
      return(NA)
    }
  }
  return(atan2(xsin, xcos)/(2*pi)*mod)
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
# TODO: add inequality-balanced accuracy

#' Gini coefficient
#' 
#' Computes the Gini coefficient, an index of inequality.
#'
#' @param x A vector to compute the Gini coefficient of
#' @param unbiased Whether to multiply by n/(n-1) to ensure generalizability to the population-level.
#'
#' @returns The Gini coefficient of the vector. 
#' \code{NA} if the vector contains a negative, zero, or NA value.
#' 
#' @export
#' @seealso \href{https://www.statsdirect.com/help/nonparametric_methods/gini.htm}{Source utilized}
#' 
#' Mills, J.A., Zandvakili, A. (1997). Statistical inference via bootstrapping for measures of inequality. 
#' Journal of Applied Econometrics, 12, 133-150. https://www.jstor.org/stable/2284908
#'
#' @examples
#' gini(state.x77[,"Income"])
#' 
#' # Completely equal gives 0
#' gini(c(100,100,100))
#' 
#' # Outliers increase the value
#' gini(c(100,100,1000))
#' 
gini <- function(x, unbiased=FALSE){
  if(any(x<=0)){
    return(NA)
  }
  x <- sort(x)
  i <- seq_along(x)
  n <- length(x)
  return(sum((2 * i - n - 1) * x) / (n * sum(x)) * ifelse(unbiased,n/(n-1),1))
}

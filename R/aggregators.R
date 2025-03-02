

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

# Add paired/unpaired nonparametric two-var comparison test
# Add correct effect size for paired tests


#' @name faststats
#' @title Statistics functions with concise output
#' @description These were mainly created because the base stats functions
#' produce overly verbose output without effect sizes included.
#' 
#' @param x,y Vectors to compare or correlate
#' @param form A formula where the left size is the values to be compared and
#' the right value is a binary variable defining the groups to compare
#' @param paired whether to perform a paired or unpaired test
#' @param data a \code{data.frame}
#' @param method The correlation type; can be pearson, spearman, or kendall
#' 
#' @md
#' @author Sercan Kahveci
#' 
#' @examples 
#' baseline <- rnorm(50)
#' deviation <- rnorm(50,m=1)
#' 
#' zerodiff(deviation)
#' 
NULL

# Formatted t-test comparing to zero
#' @rdname faststats
#' @export
#' 
zerodiff <- function(x){
  x <- na.omit(x)
  test <- t.test(x)
  cat("t (", test$parameter, ") = ", format_stat(test$statistic, digits=2),
      ", p = ", format_stat(test$p.value, digits=3, type="p"),
      ", g = ", format_stat(hedgesg(x), digits=2),
      ", M = ", format_stat(mean(x), digits=2),
      "\n", sep="")
  return(invisible(test))
}

# Formatted t-test
#' @rdname faststats
#' @export
#' 
twodiff <- function(form, data, paired=FALSE){
  term <- as.character(form)[3]
  outcome <- as.character(form)[2]
  data <- data[,c(term, outcome)]
  data <- na.omit(data)
  
  test <- t.test(form, data, paired=paired)
  x <- data[[outcome]][ data[[term]] == unique(data[[term]])[1] ]
  y <- data[[outcome]][ data[[term]] == unique(data[[term]])[2] ]
  cat("t (", test$parameter, ") = ", format_stat(test$statistic, digits=2),
      ", p = ", format_stat(test$p.value, digits=3, type="p"),
      ", g = ", format_stat(hedgesg(x=x, y=y, paired=paired), digits=2),
      ", Mdiff = ", format_stat(mean(x) - mean(y), digits=2),
      "\n", sep="")
  return(invisible(test))
}

# Formatted t-test for two inputs
#' @rdname faststats
#' @export
#' 
twodiff2 <- function(x, y, paired=FALSE){
  test <- t.test(x=x, y=y, paired=paired)
  cat("t (", test$parameter, ") = ", format_stat(test$statistic, digits=2),
      ", p = ", format_stat(test$p.value, digits=3, type="p"),
      ", g = ", format_stat(hedgesg(x=x, y=y, paired=paired), digits=2),
      ", Mdiff = ", format_stat(mean(x) - mean(y), digits=2),
      "\n", sep="")
  return(invisible(test))
}

# Formatted Wilcoxon signed-rank test
#' @rdname faststats
#' @export
#' 
npr.zerodiff <- function(x){
  x <- na.omit(x)
  test <- coin::wilcoxsign_test(rep(0, length(x))~x, exact=T)
  cat("Z = ", format_stat(test@statistic@teststatistic, digits=2),
      ", p = ", format_stat(test@distribution@pvalue(test@statistic@teststatistic), digits=3, type="p"),
      ", % pos. = ", format_stat(mean(x>0),digits=2),
      ", Mdiff = ", format_stat(mean(x), digits=2),
      "\n", sep="")
  return(invisible(test))
}

# Formatted Wilcoxon-Mann-Whitney test
#' @rdname faststats
#' @export
#' 
npr.twodiff <- function(form, data){
  term <- as.character(form)[3]
  outcome <- as.character(form)[2]
  data <- data[, c(term, outcome)]
  data[[term]] <- factor(data[[term]])
  data <- na.omit(data)
  test <- coin::wilcox_test(form, data, distribution="exact", conf.int=T)
  
  x <- data[[outcome]][ data[[term]] == unique(data[[term]])[1] ]
  y <- data[[outcome]][ data[[term]] == unique(data[[term]])[2] ]
  
  cat("Z = ", format_stat(test@statistic@teststatistic, digits=2),
      ", p = ", format_stat(test@distribution@pvalue(test@statistic@teststatistic), digits=3, type="p"),
      ", Mdiff = ", format_stat(mean(x) - mean(y), digits=2),
      "\n", sep="")
  return(invisible(test))
}

# formatted correlation
#' @rdname faststats
#' @export
#' 
twocor <- function(x, y, method="pearson"){
  test <- cor.test(x, y, method=method)
  cat(method, " r (", test$parameter, ") = ",
      format_stat(test$estimate, digits=2),
      ", p = ",format_stat(test$p.value, digits=3, type="p"),
      "\n", sep="")
  return(invisible(test))
}

#' Multiple correlation
#' Computes the \href{https://en.wikipedia.org/wiki/Multiple_correlation}{multiple correlation coefficient}
#' of variables in \code{ymat} with the variable \code{x}
#' @param x Either a matrix of variables whose multiple correlation with each other is to be estimated; 
#' or a vector of which the multiple correlation with variables in \code{ymat} is to be estimated
#' @param ymat a matrix or data.frame of variables of which 
#' the multiple correlation with \code{x} is to be estimated
#' @param use optional character indicating how to handle missing values (see \link{cor})
#'
#' @return The multiple correlation coefficient
#' @export
#' @seealso https://www.personality-project.org/r/book/chapter5.pdf
#'
#' @examples
#' multiple.cor(mtcars[,1],mtcars[,2:4])
#' 
multiple.cor <- function(x, ymat, use="everything"){
  if(missing(ymat)){
    cv <- cor(x, use=use)
    corvec <- numeric(ncol(x))
    for(i in seq_along(corvec)){
      gfvec <- cv[seq_len(nrow(cv))[-i], i]
      dcm <- cv[seq_len(nrow(cv))[-i], seq_len(ncol(cv))[-i]]
      rsq <- t(gfvec) %*% solve(dcm) %*% gfvec
      corvec[i] <- sqrt(as.vector(rsq))
    }
    names(corvec) <- colnames(cv)
    return(corvec)
  }else{
    cv <- cor(cbind(x,ymat), use=use)
    gfvec <- cv[2:nrow(cv),1]
    dcm <- cv[2:nrow(cv),2:ncol(cv)]
    rsq <- t(gfvec) %*% solve(dcm) %*% gfvec
    return(sqrt(as.vector(rsq)))
  }
}

#' Generate a correlation table
#' 
#' This computes correlations from the imputed data, yielding a correlation matrix
#' as well as matrices for sample size, p-values, and holdout values (computed with
#' [cor.holdout()]). The latter is a robustness check for the obtained correlation,
#' denoting how many observations need to be removed before a certain goal is achieved,
#' such as non-significance, or a flip of the sign of the correlation coefficient
#' (represented by the argument \code{holdout.goal}).
#'
#' @param x A matrix or data.frame of which the columns will be correlated with each other
#' @param method Whether to use pearson or spearman correlations
#' @param holdout.goal When computing the holdout statistic with [cor.holdout()], what should
#' the goal be - non-significance ("nsig") or a flip of the sign ("flip")?
#' @param alpha Alpha level for non-significance testing in case \code{holdout.goal="nsig"}
#' 
#' @details The output of this object can be printed with [print.CorTable()].
#'
#' @return A list with 5 items - 4 matrices 
#' (r - correlation, n = sample size, p = p-value, h = holdout value) 
#' and a list of parameters, such as the correlation type.
#' NA values among the holdout values indicate that computation failed, e.g. 
#' due to insufficient variance or an inability to achieve the holdout goal
#' without depleting the sample size.
#' 
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' test<-CorTable(mtcars,holdout.goal="flip")
#' print(test)
#' 
#' test<-CorTable(mtcars,holdout.goal="nsig")
#' print(test,type="r/h")
#' 
CorTable <- function(x, method=c("pearson","spearman"),
                     holdout.goal = c("nsig","flip"), alpha=.05){
  method <- match.arg(method)
  holdout.goal <- match.arg(holdout.goal)
  
  emptymat <- matrix(NA, nrow=ncol(x), ncol=ncol(x),
                     dimnames=list(colnames(x),colnames(x)))
  coefs<-c("r","n","p","h")
  output <- setNames(lapply(coefs,function(x) { emptymat }),coefs)
  
  testpairs <- allpairs(nval=ncol(x))
  for(i in seq_len(NROW(testpairs))){
    trow <- testpairs[i,1]
    tcol <- testpairs[i,2]
    currvars <- na.omit(x[,c(trow,tcol)])
    output$r[trow,tcol] <- tcor <- cor(currvars,method=method)[1,2]
    output$n[trow,tcol] <- tn <- NROW(currvars)
    output$p[trow,tcol] <- 2*pt(-abs(r2t(tcor,tn-2)), df=tn-2)
    hobj <- cor.holdout(x=currvars[,1,drop=T], y=currvars[,2,drop=T],
                        goal=holdout.goal, method=method, alpha=alpha)
    output$h[trow,tcol] <- ifelse(hobj$success,hobj$h,NA)
  }
  
  # Form output object
  output <- lapply(output, function(x){
    x[lower.tri(x)] <- t(x)[lower.tri(x)]
    return(x)
  })
  output$parameters <- list(method=method,alpha=alpha,holdout.goal=holdout.goal)
  output <- structure(output,class=c("CorTable","list"))
  
  return(output)
}


#' Print a \code{CorTable} object
#'
#' @param x a \code{CorTable} object
#' @param type The type of table to print. 
#' * \code{"full"} Separately print the correlation, sample size, p-value, and holdout matrices
#' * \code{"r/p"} Print a single matrix with p-values in the upper triangle and 
#' correlations in the lower triangle 
#' * \code{"r/h"} Print a single matrix with holdout values in the upper triangle and
#' correlations in the lower triangle 
#' @param alpha Correlations with a p-value below this value will be highlighted 
#' with an asterisk.
#' @param digits Digits to round the values by. p-values will be rounded by this value +1
#' @param ... Ignored
#'
#' @return The printed character matrix is silently returned.
#' @author Sercan Kahveci
#' 
#' @md
#' @export
#' @method print CorTable
#'
#' @examples
#' test<-CorTable(mtcars,holdout.goal="flip")
#' print(test)
#' 
#' test<-CorTable(mtcars,holdout.goal="nsig")
#' print(test,type="r/h")
#' 
print.CorTable <- function(x, type=c("full", "r/p", "r/h"),
                           alpha=0, digits=2, ...){
  type <- match.arg(type)
  
  # Format the values
  printx <- list()
  for(i in 1:4){
    y <- round(x[[i]], digits=digits + ifelse(names(x)[i]=="p", 1, 0))
    y[] <- dropLeadingZero(y)
    diag(y) <- "."
    printx[[ names(x)[i] ]] <- y
  }
  printx$r[which(x$p<alpha)] <- paste0("*",printx$r[which(x$p<alpha)])
  
  hdesc<-paste0("minimum number of cases to be removed to achieve ",
                ifelse(x$parameters$holdout.goal == "nsig",
                       paste0("non-significance with alpha ", x$parameters$alpha),
                       "a sign flip"))
  
  # Print parameters and prepare final printed matrix
  if(type=="full"){
    cat("Correlation type: ", x$parameters$method,
        "\nh coefficient type: ", hdesc,
        "\n", sep="")
  }else{
    newprintx <- printx$r
    cat("Lower triangle: ", x$parameters$method, " correlations\n", sep="")
    if(type=="r/p"){
      cat("Upper triangle: p-values\n")
      newprintx[upper.tri(newprintx)] <- printx$p[upper.tri(newprintx)]
    }else if(type=="r/h"){
      cat("Upper triangle: ", hdesc,
          "\n", sep="")
      newprintx[upper.tri(newprintx)] <- printx$h[upper.tri(newprintx)]
    }
    printx <- newprintx
  }
  
  # Print
  print(printx, quote=F, right=T)
  return(invisible(printx))
}



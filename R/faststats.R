
# Remove rows with OLs from data frame
removeOLs <- function(.tbl,olvars,groups=NULL){
  newtbl <- .tbl %>% group_by(across(all_of(groups))) %>% 
    filter(if_all(.cols=(olvars),.fns=function(x){ abs(vec.scale(x))<3 })) %>% 
    ungroup()
  message("Filtered ",nrow(.tbl)-nrow(newtbl)," rows")
  return(newtbl)
}

# Replace OLs with NA
maskOLs <- function(.tbl,olvars,groups=NULL){
  if(!is.null(groups)){
    groupvar <- interaction(.tbl[groups])
  }else{
    groupvar <- rep(1,nrow(.tbl))
  }
  for(olvar in olvars){
    key <- tapply(.tbl[[olvar]], groupvar,
                  function(x) abs(vec.scale(x)) > 3) %>% 
      unlist() %>% which()
    .tbl[[olvar]][key] <- NA
    message("Masked ", length(key), " outliers from variable ", olvar)
  }
  return(.tbl)
}

# Remove OLs in a vector
vec.removeOLs <- function(x){
  excl<-which(abs(vec.scale(x)) > 3)
  message("Excluding ", length(excl), " observations from vector")
  if(length(excl) > 0){ 
    return(x[-excl])
  }else{ 
    return(x)
  }
}

# Formatted t-test comparing to zero
zerodiff <- function(x){
  x <- na.omit(x)
  tt <- t.test(x)
  cat("t (", tt$parameter, ") = ", tt$statistic,
      ", p = ", tt$p.value,
      ", g = ", format(CohenD(x,correct=T), digits=5),
      ", M = ", format(mean(x), digits=5),
      "\n", sep="")
}

# Formatted t-test
twodiff <- function(form,data,paired=F){
  term <- as.character(form)[3]
  outcome <- as.character(form)[2]
  data <- data[,c(term, outcome)]
  data <- na.omit(data)
  
  tt <- t.test(form, data, paired=paired)
  x <- data[[outcome]][ data[[term]] == unique(data[[term]])[1] ]
  y <- data[[outcome]][ data[[term]] == unique(data[[term]])[2] ]
  cat("t (", tt$parameter, ") = ", tt$statistic, 
      ", p = ", tt$p.value,
      ", g = ", format(CohenD(x=x, y=y, correct=T), digits=5),
      ", Mdiff = ", format(mean(x) - mean(y), digits=5),
      "\n", sep="")
}

# Formatted t-test for two inputs
twodiff2 <- function(x, y, paired=F){
  tt <- t.test(x=x, y=y, paired=paired)
  cat("t (", tt$parameter, ") = ", tt$statistic, 
      ", p = ", tt$p.value,
      ", g = ", format(CohenD(x=x, y=y, correct=T), digits=5),
      ", Mdiff = ", format(mean(x) - mean(y), digits=5),
      "\n", sep="")
}

npr.zerodiff <- function(x){
  x <- na.omit(x)
  test <- coin::wilcoxsign_test(rep(0, length(x))~x, exact=T)
  cat("Z = ", format(test@statistic@teststatistic, digits=5),
      ", p = ", format(test@distribution@pvalue(test@statistic@teststatistic), digits=5),
      ", Mdiff = ", format(mean(x), digits=5),
      "\n", sep="")
}

# Formatted Wilcoxon test
npr.twodiff <- function(form,data){
  term <- as.character(form)[3]
  outcome <- as.character(form)[2]
  data <- data[, c(term,outcome)]
  data[[term]] <- factor(data[[term]])
  data <- na.omit(data)
  test <- coin::wilcox_test(form,data, distribution="exact",conf.int=T)
  
  x <- data[[outcome]][ data[[term]] == unique(data[[term]])[1] ]
  y <- data[[outcome]][ data[[term]] == unique(data[[term]])[2] ]
  
  cat("Z = ", format(test@statistic@teststatistic, digits=5),
      ", p = ", format(test@distribution@pvalue(test@statistic@teststatistic), digits=5),
      ", Mdiff = ", format(mean(x) - mean(y), digits=5),
      "\n", sep="")
}

# formatted correlation
print.cor <- function(x, y, method="pearson"){
  h <- cor.test(x, y, method=method)
  cat(method, " r (", h$parameter, ") = ", dropLeadingZero(format(h$estimate,digits=2)),
      ", p = ",dropLeadingZero(format(h$p.value,digits=3)),
      "\n", sep="")
}





#' Get classification accuracy metrics
#'
#' @param x The true class memberships
#' @param y The predicted class memberships
#'
#' @md
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
#' @examples
#' mymod <- glm(am ~ cyl + disp + hp,
#'              family="binomial", 
#'              data = mtcars)
#' 
#' newpred <- as.numeric(predict(mymod, 
#'                               data=mtcars,
#'                               type="response") > 0.5)
#' 
#' CVMetrics(x=mtcars$am, y=newpred)
#' 
CVMetrics <- function(x, y){
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
      gfvec <- cv[(1:nrow(cv))[-i], i]
      dcm <- cv[(1:nrow(cv))[-i], (1:ncol(cv))[-i]]
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


#' Bonferroni-Holm significance evaluation
#' Evaluate the significance of a vector of p-values according to
#' the Bonferroni-Holm method
#' 
#' @param x A vector of p-values.
#' @param alpha An experiment-wise alpha-level to test against.
#'
#' @return A vector indicating whether values in \code{x} were significant
#' according to the Bonferroni-Holm method.
#' @export
#'
#' @examples
#' bonferroniHolm(c(1,.05,.001,.002,.01))
#' 
bonferroniHolm <- function(x, alpha=.05){
  xord <- rank(x)
  xlen <- length(x)
  sigvec <- logical(xlen)
  for(i in sort(unique(xord)) ){
    idx <- which(xord == i)
    if(x[idx[1]] < alpha/(xlen-sum(sigvec)) ){
      sigvec[idx] <- T
    }else{
      break
    }
  }
  return(sigvec)
}

#' Influence function of the Pearson correlation coefficient
#'
#' @param x,y Numeric vectors
#'
#' @return Influence values of all observations.
#' @export
#'
#' @examples
#' outlier<-numeric(100)
#' outlier[1]<-1000
#' cor.influence(rnorm(100)+outlier,rnorm(100)+outlier)
#' 
cor.influence<-function(x,y){
  x<-x-mean(x)
  y<-y-mean(y)
  x*y-(x^2+y^2)/2*cor(x,y)
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

# CorTable(mtcars) |> print(type="r/h",alpha=0)

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
#' @param ... 
#'
#' @return The printed character matrix is silently returned.
#' 
#' @md
#' @export
#'
#' @examples
#' test<-CorTable(mtcars,holdout.goal="flip")
#' print(test)
#' 
#' test<-CorTable(mtcars,holdout.goal="nsig")
#' print(test,type="r/h")
#' 
print.CorTable<-function(x, type=c("full","r/p","r/h"),
                         alpha=0, digits=2, ...){
  type <- match.arg(type)
  
  # Format the values
  printx <- list()
  for(i in 1:4){
    y <- round(x[[i]],digits=digits+ifelse(names(x)[i]=="p",1,0))
    y[] <- dropLeadingZero(y)
    diag(y)<-"."
    printx[[ names(x)[i] ]]<-y
  }
  printx$r[which(x$p<alpha)]<-paste0("*",printx$r[which(x$p<alpha)])
  
  hdesc<-paste0("minimum number of cases to be removed to achieve ",
                ifelse(x$parameters$holdout.goal=="nsig",
                       paste0("non-significance with alpha ",x$parameters$alpha),
                       "a sign flip"))
  
  # Print parameters and prepare final printed matrix
  if(type=="full"){
    cat("Correlation type: ",x$parameters$method,
        "\nh coefficient type: ",hdesc,
        "\n",sep="")
  }else{
    newprintx<-printx$r
    cat("Lower triangle: ",method," correlations\n",sep="")
    if(type=="r/p"){
      cat("Upper triangle: p-values\n")
      newprintx[upper.tri(newprintx)] <- printx$p[upper.tri(newprintx)]
    }else if(type=="r/h"){
      cat("Upper triangle: ",hdesc,
          "\n",sep="")
      newprintx[upper.tri(newprintx)] <- printx$h[upper.tri(newprintx)]
    }
    printx <- newprintx
  }
  
  # Print
  print(printx,quote=F,right=T)
  return(invisible(printx))
}




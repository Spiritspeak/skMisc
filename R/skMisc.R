#' skMisc: Miscellaneous convenience functions
#'
#' @description Functionality provided
#' 
#' Some of the functions in this package are statistical:
#' * Common outlier exclusion methods and statistical analyses with concise print functions
#' * Functions for manipulating lme4 formulas
#' * Methods for obtaining "holdout statistics"
#' 
#' Some of the functions are geared towards datawrangling:
#' * Functions for merging, subsetting, and transforming data.frames and matrices
#' * Specific date and string manipulations unavailable in common packages
#' * Other small convenience functions that make life easier
#'
#' @name skMisc
#' @useDynLib skMisc
#' @md 
#' @import magrittr dplyr lmerTest ggplot2 coin doParallel parallel foreach iterators MuMIn 
#' Rcpp stringr
#' @importFrom utils install.packages
#' @importFrom methods as Quote
#' @importFrom grDevices rgb col2rgb dev.off png
#' @importFrom graphics abline lines par rect strheight strwidth text
#' @importFrom stats AIC BIC as.formula ave coef cor cor.test deviance
#' formula ks.test lm.influence logLik model.frame na.omit pchisq pt
#' quantile reformulate sd setNames t.test terms update var lm resid
#' optimize pnorm
#' @importFrom car qqp
#' @importFrom lme4 nobars findbars fixef getME lmerControl ranef
#' @importFrom MASS mvrnorm
#' 
#' 
"_PACKAGE"

.onLoad <- function(libname, pkgname){
  #packageStartupMessage("Thank you for loading skMisc v0.02")
  utils::globalVariables(c("currform","x.rowid","y.rowid"))
}

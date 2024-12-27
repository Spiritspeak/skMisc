#' skMisc: Miscellaneous convenience functions
#'
#' @section Functionality provided
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
#' @docType package
#' @name skMisc
#' @md 
#' @import magrittr dplyr lme4 ggplot2
#' 
NULL

.onLoad <- function(libname, pkgname){
  packageStartupMessage("Thank you for loading skMisc v0.02")
}
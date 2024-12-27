#' skMisc: Miscellaneous convenience functions
#'
#' The package provides these categories of functions:
#' * Common statistical analyses with concise print functions
#' * Methods for obtaining "holdout statistics"
#' * Functions for manipulating lme4 formulas
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
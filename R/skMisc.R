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
#' @import qgraph Rcpp deldir 
#' @importFrom utils install.packages
#' @importFrom methods as Quote
#' @importFrom grDevices rgb col2rgb dev.off png
#' @importFrom graphics abline lines par rect strheight strwidth text
#' @importFrom stats AIC BIC as.formula ave coef cor cor.test deviance
#' formula ks.test lm.influence logLik model.frame na.omit pchisq pt
#' quantile reformulate sd setNames t.test terms update var lm resid
#' optimize pnorm
#' @importFrom MASS fitdistr
#' @importFrom car qqp
#' @importFrom lme4 nobars findbars fixef getME lmerControl ranef
#' @importFrom igraph graph_from_adjacency_matrix E V delete_edges
#' layout_with_kk layout_with_fr layout_with_drl layout_with_dh 
#' make_ring subgraph_isomorphisms
#' @importFrom graphlayouts layout_as_dynamic layout_with_stress layout_with_focus
#' @importFrom glmnet glmnet
#' @importFrom arm se.coef
#' @importFrom Matrix Matrix
#' 
#' 
"_PACKAGE"

.onLoad <- function(libname, pkgname){
  #packageStartupMessage("Thank you for loading skMisc v0.02")
  utils::globalVariables(c("currform","x.rowid","y.rowid"))
}

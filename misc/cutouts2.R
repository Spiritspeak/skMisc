
#' Adjusted p-values for false-discovery rates and familywise error rates
#' 
#' @param x A vector of p-values.
#'
#' @return A vector of p-values adjusted in accordance with the chosen method.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' oldp<-c(0.01,0.04,0.02,.1,0.035)
#' bonferroniHolm(oldp)
#' benjaminiHochberg(oldp)
#' benjaminiYekutieli(oldp)
#' 
bonferroniHolm <- function(x){
  xn <- x[!is.na(x)]
  xord <- order(xn)
  xlen <- length(xn)
  newp <- xn[xord] * (xlen-seq_len(xlen)+1)
  newp[newp>1] <- 1
  # Bonferroni Holm monotonicity correction (up-propagated max)
  out <- sapply(seq_len(xlen),function(y)max(newp[1:y]))[order(xord)]
  x[!is.na(x)]<-out
  return(x)
}

#' @rdname bonferroniHolm
#' @export
#' 
benjaminiHochberg <- function(x){
  xn <- x[!is.na(x)]
  xord <- order(xn)
  xlen <- length(xn)
  newp <- xn[xord] * xlen/seq_len(xlen)
  newp[newp>1] <- 1
  # Benjamini Hochberg monotonicity correction (down-propagated min)
  out <- sapply(seq_len(xlen),function(y)min(newp[y:xlen]))[order(xord)]
  x[!is.na(x)]<-out
  return(x)
}

#' @rdname bonferroniHolm
#' @export
#' 
benjaminiYekutieli <- function(x){
  xn <- x[!is.na(x)]
  xord <- rank(xn)
  xlen <- length(xn)
  hf <- sum(1/(1:xlen))
  newp <- xn * hf * xlen/(xord)
  newp[newp>1] <- 1
  # Benjamini Yekutieli monotonicity correction (down-propagated min)
  out <- sapply(seq_len(xlen),function(y)min(newp[y:xlen]))[order(xord)]
  x[!is.na(x)]<-out
  return(x)
}

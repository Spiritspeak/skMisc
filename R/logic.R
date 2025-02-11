
#' Negative \%in\%
#' 
#' Returns which values in the left vector are not in the right vector.
#' 
#' @param x Values whose presence will be checked for in \code{table}.
#' @param table Values that will yield a \code{FALSE} if they exist in \code{x}.
#'
#' @return A logical vector indicating whether each value of \code{x} lacks a match in \code{table}
#' @export
#'
#' @examples
#' (1:5) %nin% (1:3)
#' 
`%nin%` <- function(x,table){
  !(x %in% table)
}

#' @name Case-insensitive-operators
#' @title Case-insensitive inline operators
#' @description 
#' These operators convert both left- and right-hand side to lowercase 
#' before comparing.
#' @param x,y Values to be compared.
#' @return A logical vector.
#' @author Sercan Kahveci
#' @seealso [base::`%in%`], [`%nin%`], [base::tolower()]
#' 
NULL

#' @describeIn Case-insensitive-operators Case-insensitive %in%
#' @export
#' @examples
#' c("apple","pear","coal") %cin% c("Banana","Apple","Pear","Cherry")
#' 
`%cin%` <- function(x, y){
  tolower(x) %in% tolower(y)
}

#' @describeIn Case-insensitive-operators Case-insensitive %nin%
#' @export
#' @examples
#' c("apple","pear","coal") %ncin% c("Banana","Apple","Pear","Cherry")
#' 
`%ncin%` <- function(x, y){
  !(tolower(x) %in% tolower(y))
}

#' @describeIn Case-insensitive-operators Case-insensitive ==
#' @export
#' @examples
#' c("APPLE","COAL") %cis% c("Apple","Pear")
#' 
`%cis%` <- function(x, y){
  tolower(x) == tolower(y)
}

#' @describeIn Case-insensitive-operators Case-insensitive !=
#' @export
#' @examples
#' c("APPLE","COAL") %ncis% c("Apple","Pear")
#' 
`%ncis%` <- function(x, y){
  tolower(x) != tolower(y)
}

#' @name lazylogic
#' @title Logical operators using lazy evaluation
#' @description
#' These functions sequentially evaluate their arguments and return
#' a logical value when sufficient information has been acquired to do so. 
#' Hence, \code{lazy_any()} will not evaluate any arguments 
#' beyond the first \code{TRUE}, since there is already
#' at least one \code{TRUE} value, so \code{TRUE} can be returned.
#' Likewise, \code{lazy_all()} will not evaluate beyond the first \code{FALSE} 
#' it is already clear not all arguments are \code{TRUE}, 
#' so \code{FALSE} is returned. 
#' 
#' This enables logical chains in which you can use logical statements 
#' earlier in the chain to check whether logical statements 
#' further down the chain can be evaluated or would cause errors. 
#' If a checking statement evaluates to the terminating condition 
#' of the used function, 
#' further arguments are not evaluated and hence no error is triggered.
#'
#' @param ... Expressions evaluating to single logical values. 
#'
#' @return A logical value.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' # The final argument is not evaluated because 
#' # the function reaches its termination condition 
#' # (TRUE) before that.
#' lazy_any(FALSE, TRUE, stop())
#' 
#' # Dealing with problematic NULL, NA, and multi-value inputs
#' myvecs <- list(a=NULL,b=NA,c=c(1,5,2),d=10)
#' outcomes <- logical(length(myvecs))
#' for(i in seq_along(myvecs)){
#'   outcomes[i] <- lazy_all(!is.null(myvecs[[i]]),
#'                           !length(myvecs[[i]]) != 1,
#'                           !is.na(myvecs[[i]]),
#'                           myvecs[[i]] == 10)
#' }
#' 
lazy_any <- function(...){
  had_na <- F
  for(i in seq_len(...length())){
    currarg <- ...elt(i)
    if(is.na(currarg)){ 
      had_na <- T
    }else if(currarg){ 
      return(TRUE)
    }
  }
  if(had_na){
    return(NA)
  }else{
    return(F)
  }
}

#' @rdname lazylogic
#' @export
#' 
lazy_all <- function(...){
  had_na <- F
  for(i in seq_len(...length())){
    currarg <- ...elt(i)
    if(is.na(currarg)){ 
      had_na <- T
    }else if(!currarg){ 
      return(FALSE)
    }
  }
  if(had_na){
    return(NA)
  }else{
    return(TRUE)
  }
}


#' @name logical-na-handlers
#' @title Coerce NA values to \code{TRUE} or \code{FALSE}
#' 
#' @param x A logical vector.
#'
#' @return \code{x} with \code{NA} values replaced with 
#' \code{TRUE} or \code{FALSE}.
#' @author Sercan Kahveci
#'
#' @examples
#' fruits <- c("apples","pears",NA,"cherries")
#' NA2TRUE(fruits != "apples")
#' NA2FALSE(fruits == "pears")
#' 
NULL

#' @rdname logical-na-handlers
#' @export
#' 
NA2TRUE <- function(x){
  x[is.na(x)] <- TRUE
  x
}

#' @rdname logical-na-handlers
#' @export
#' 
NA2FALSE <- function(x){
  x[is.na(x)] <- FALSE
  x
}

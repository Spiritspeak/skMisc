
#' Drop leading zeros
#' Remove leading zeroes and return as a character vector.
#'
#' @param x Numeric vector to remove leading zeros from
#'
#' @return A character vector of numbers with leading zeros removed.
#' @export
#' @md
#'
#' @examples
#' dropLeadingZero(c(-1,0,1,0.5,-0.5,1.5,-1.5))
#' 
dropLeadingZero <- function(x){
  gsub("(?<![0-9])0+(?=\\.)", "", x, perl = TRUE)
}


#' Format statistics appropriately
#' This formats numbers according to APA.
#'
#' @param x A numeric vector.
#' @param digits How many digits to round the number to.
#' @param type Type of statistic formatting rule to apply. 
#' "p" prints values rounded to zero as <.001, and 
#' "quotient" prints values <1 as 1/x where x is a value above 1 
#' (e.g. in case of Bayes factors).
#' @param sign.positive If \code{TRUE}, adds a + to every positive number.
#' 
#' @details
#' All leading zeros are also dropped. 
#'
#' @return A character vector of formatted numbers.
#' @export
#'
#' @examples
#' format_stat(0.12345678)
#' 
#' # Proper printing of p-values
#' format_stat(0.0004,type="p")
#' 
#' # Printing of quotients where the range of values between 1 and 0
#' # should be considered equal to that between 1 and infinity
#' format_stat(0.05,type="quotient")
#' 
format_stat <- function(x, digits=2, type=c("default","p","quotient"),
                        sign.positive=FALSE){
  type <- match.arg(type)
  printx <- dropLeadingZero(format(round(x, digits=digits), scientific=F))
  if(type == "quotient"){
    key <- abs(x)<1
    printx[key] <- 
      paste0(ifelse(x[key] < 0,"-",""),"1/",
             dropLeadingZero(format(round(abs(1/x[key]), digits=digits), 
                                    scientific=F)))
  }
  if(type=="p"){
    printx[printx=="0"] <- "<.001"
  }
  if(sign.positive){
    key <- sign(x) == 1
    printx[key] <- paste0("+",printx[key])
  }

  return(printx)
}

#' Convert a vector to an English list
#'
#' @param x A vector of values to convert into a string representing 
#' a grammatically correct English list
#'
#' @return a A string representing a grammatically correct English list
#' @export
#'
#' @examples
#' vec2phrase(c("apples","oranges"))
#' 
#' vec2phrase(c("eggs"))
#' 
#' vec2phrase(c())
#' 
#' vec2phrase(c("cheese","milk","yoghurt","kefir"))
#' 
vec2phrase <- function(x){
  lx <- length(x)
  out <- switch(EXPR=as.character(lx),
                `0`="",
                `1`=as.character(x),
                `2`=paste(x[1], "and", x[2]),
                paste0(paste0(x[-lx], collapse=", "), ", and ", x[lx]))
  return(out)
}



#' Create substrings with a maximal length by splitting at specific characters
#' 
#' This function splits a string into substrings of length \code{width} or shorter.
#' The splitting is done at the characters specified in \code{split}, in order of preference.
#' 
#' This combines the functionality of [base::strwrap()] and [base::strsplit()]; 
#' instead of a string wrapped with newlines, the result is multiple substrings.
#'
#' @param x A character vector of length 1.
#' @param width The maximum character length to break the vector at.
#' @param split A vector of regular expressions to match a character to break the string at.
#' The function will try to break the string at the first value specified in this argument;
#' if that fails, it will move on to the second, then the third, etc.
#'
#' @return A character vector consisting of strings of length \code{width} or shorter, 
#' and split at the characters specified in \code{split}.
#' 
#' @export
#' @md
#'
#' @examples
#' thanks <- paste(readLines(file.path(R.home("doc"), "THANKS")), collapse = "\n")
#' strsplit.wrap(thanks,width=80)
#' 
#' alphabet <- paste0(letters,collapse="")
#' strsplit.wrap(alphabet,width=3)

strsplit.wrap <- function(x, width=2000, split=c("\n"," ",",","")){
  output <- character()
  if(!any(split=="")){ split <- c(split,"") }
  while(nchar(x) > 0){
    cstr <- substr(x,1,width)
    if(nchar(cstr) < width){
      output[length(output)+1] <- x
      x <- ""
    }else{
      for(splitchar in split){
        
        if(nzchar(splitchar)){
          nls <- gregexpr(splitchar,cstr)[[1]]
          end <- nls[length(nls)]
          if(end!=-1){
            output[length(output)+1] <- 
              trimws(substr(cstr,1,end),whitespace=splitchar)
            x <- substr(x, end+1, nchar(x))
            break
          }
        }else{
          output[length(output)+1] <- substr(cstr,1,width)
          x <- substr(x, width+1, nchar(x))
        }
        
      }
    }
  }
  return(output)
}



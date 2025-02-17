#########
# Functions for matching and merging datasets with each other by ID and date

min2<-function(x){ 
  if(length(x) == 0){ 
    NA 
  }else{ 
    min(x) 
  } 
}

#' Match recurring session IDs using their date
#' 
#' Match session IDs from two datasets by taking into account the date of the sessions.
#' This is useful when linking together data from several experiments done by the
#' same participants; often participants accidentally run the same experiment twice, so
#' this should help to distinguish which experiment datasets temporally fit together.
#'
#' @param x.ids,y.ids IDs to match.
#' @param x.dates,y.dates Dates to use in the matching of IDs.
#'
#' @return A \code{data.frame} with two columns, one specifying the indices of x and 
#' one specifying the matched indices of y.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' # In the following example, participant 1 has a single direct match,
#' # one of participant 2's two sessions in x can be matched to a session in y,
#' # both of participant 3's sessions in x have a match in y,
#' # and participant 4 has no match between x and y because the observation in y
#' # occurs before that in x
#' x <- data.frame(id=c(1,2,2,3,3,4),
#'                 date=as.POSIXct(c(1,1,3,1,3,1)))
#' y <- data.frame(id=c(1,2,3,3,4),
#'                 date=as.POSIXct(c(2,2,4,2,0)))
#' match.pps(x.ids=x$id,x.dates=x$date,y.ids=y$id,y.dates=y$date)
#' #>   x.rowid y.rowid
#' #> 1       1       1
#' #> 2       2       2
#' #> 3       5       3
#' #> 4       4       4
#' 
match.pps <- function(x.ids, x.dates, y.ids, y.dates){
  uids <- unique(c(x.ids, y.ids))
  
  matchlist <- list()
  for(uid in uids){
    sm <- expand.grid(x.rowid = which(x.ids == uid), 
                      y.rowid = which(y.ids == uid))
    if(nrow(sm) > 0){
      sm$diff <- as.numeric(y.dates[sm$y.rowid]) - as.numeric(x.dates[sm$x.rowid])
      sm <- sm %>% filter(diff >= 0) %>% 
        group_by(x.rowid) %>% filter(diff == min2(diff)) %>% 
        group_by(y.rowid) %>% filter(diff == min2(diff)) %>% 
        as.data.frame()
      matchlist[[length(matchlist) + 1]] <- sm[, -3]
    }
  }
  matchset <- do.call(rbind,matchlist)
  return(matchset)
}


#' Match and merge two data.frames by ID and date
#' 
#' the rows from \code{y} are matched with rows from \code{x} that 
#' share the same ID and directly precede them in date.
#' 
#' This is useful for matching participant data from different datasets,
#' where one dataset is always collected after another. It is assumed here that
#' \code{x} is always administered before \code{y}.
#'
#' @param x,y The \code{data.frame}s to merge.
#' @param idvar Name of the variable holding participant/session IDs in both x and y.
#' @param datevar Name of the variable holding dates in both x and y.
#' @param keep.unmatched When no match is found for a row, should it be kept or removed?
#'
#' @return A merged \code{data.frame} where rows of \code{x}
#' are merged with rows of \code{y} that match by session id and where
#' the date of each row from \code{x} directly precedes the 
#' date of the matched row from \code{y}.
#' 
#' @details This matches using [match.pps()] and merges using [dplyr::full_join()].
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' x <- data.frame(id=c(1,2,2,3,3,4),
#'                 date=as.POSIXct(c(1,1,3,1,3,1)),
#'                info=letters[1:6])
#' y <- data.frame(id=c(1,2,3,3,4),
#'                date=as.POSIXct(c(2,2,4,2,0)),
#'                 data=LETTERS[1:5])
#' match.merge(x=x,y=y,idvar="id",datevar="date",keep.unmatched="all")
#' match.merge(x=x,y=y,idvar="id",datevar="date",keep.unmatched="none")
#' 
match.merge <- function(x, y, idvar, datevar,
                        keep.unmatched = c("all","x","y","none")){
  keep.unmatched <- match.arg(keep.unmatched)
  matches <- match.pps(x[[idvar]], x[[datevar]], y[[idvar]], y[[datevar]])
  x$.rowmatch <- NA
  y$.rowmatch <- NA
  x$.rowmatch[matches[[1]]] <- seq_along(matches[[1]])
  y$.rowmatch[matches[[2]]] <- seq_along(matches[[2]])
  if(any(keep.unmatched==c("y","none"))){ x<-x[!is.na(x$.rowmatch),] }
  if(any(keep.unmatched==c("x","none"))){ y<-y[!is.na(y$.rowmatch),] }
  
  args <- list(x = x, y = y, 
               by = ".rowmatch",
               na_matches = "never")
  
  z <- do.call(full_join, args)
  z$.rowmatch <- NULL
  return(z)
}


#' Merge Multiple Data Frames
#' 
#' This function makes calls to \code{merge()} to merge every other dataset 
#' with the one next to it, repeating until only one dataset remains. 
#'
#' @param x a list of data frames.
#' @param ... all other arguments for \code{merge} can be provided here.
#'
#' @return A single, merged \code{data.frame}.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' #generate test data
#' testlist<-list()
#' lsize<-50
#' for(i in 1:lsize){
#'   testlist[[i]]<-data.frame(key=sample(1:500,100),
#'                             junk=letters[sample(1:26,100,replace=TRUE)])
#'   colnames(testlist[[i]])[2]<-paste0("info",i)
#' }
#' multimerge(testlist,by="key",all=TRUE)
#' 
multimerge <- function(x, ...){
  while(length(x) > 1){
    out <- list()
    while(length(x) > 0){
      if(length(x) >= 2){
        out[[length(out) + 1]] <- merge(x[[1]], x[[2]], ...)
        x[[1]] <- NULL
        x[[1]] <- NULL
      }else{
        out[[length(out) + 1]] <- x[[1]]
        x[[1]] <- NULL
      }
    }
    x <- out
  }
  return(x)
}


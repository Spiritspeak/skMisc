#########
# Assorted code from the switch project that I haven't yet generalized

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
#' @param x.ids,y.ids IDs to match
#' @param x.dates,y.dates Dates to use in the matching of IDs
#'
#' @return A \code{data.frame} with two columns, one specifying the indices of x and 
#' one specifying the matched indices of y.
#' @export
#'
#' @examples
#' # In the following example, participant 1 has a single direct match,
#' # one of participant 2's two sessions in x can be matched to a session in y,
#' # both of participant 3's sessions in x have a match in y,
#' # and participant 4 has no match between x and y because the observation in y
#' # occurs before that in x
#' x <- data.frame(id=c(1,2,2,3,3,4),
#'                 date=as.POSIXct(c(1,1,3,1,2,1)))
#' y <- data.frame(id=c(1,2,3,3,4),
#'                 date=as.POSIXct(c(2,2,4,3,0)))
#' match.pps(x.ids=x$id,x.dates=x$date,y.ids=y$id,y.dates=y$date)
#' 
match.pps<-function(x.ids, x.dates, y.ids, y.dates){
  uids <- unique(c(x.ids, y.ids))
  
  matchlist <- list()
  for(uid in uids){
    sm <- expand.grid(x.rowid = which(x.ids == uid), 
                      y.rowid = which(y.ids == uid))
    if(nrow(sm) > 0){
      sm$diff <- y.dates[sm$y.rowid] - x.dates[sm$x.rowid]
      sm <- sm %>% group_by(x.rowid) %>% filter(diff >= 0) %>% 
        filter(diff == min2(diff)) %>% as.data.frame()
      matchlist[[length(matchlist) + 1]] <- sm[, -3]
    }
  }
  matchset <- do.call(rbind,matchlist)
  return(matchset)
}
# Fix this code

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
#' @param x.ids,y.ids Vectors of participant/session IDs.
#' @param x.dates,y.dates 
#' @param by Additional variables to merge by
#' @param ... Additional arguments passed to \code{base::merge()}
#'
#' @return A merged \code{data.frame} where rows of \code{x}
#' are merged with rows of \code{y} that match by session id and where
#' the date of each row from \code{x} directly precedes the 
#' date of the matched row from \code{y}.
#' @export
#'
#' @examples
#' 
#' 
match.merge<-function(x, x.ids, x.dates, 
                      y, y.ids, y.dates, 
                      by, ...){
  matches <- match.pps(x.ids, x.dates, y.ids, y.dates)
  x$.rowmatch <- NA
  y$.rowmatch <- NA
  x$.rowmatch[matches[[1]]] <- seq_along(matches[[1]])
  y$.rowmatch[matches[[2]]] <- seq_along(matches[[2]])
  z <- do.call(merge, list(x = x, y = y, by = c(by,".rowmatch"), ...))
  z$.rowmatch <- NULL
  return(z)
}


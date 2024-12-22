#########
# Assorted code from the switch project that I haven't yet generalized

min2<-function(x){ 
  if(length(x) == 0){ 
    NA 
  }else{ 
    min(x) 
  } 
}

match.pps<-function(x.ids, x.dates, y.ids, y.dates){
  uids <- unique(c(x.ids, y.ids))
  
  matchlist <- list()
  for(uid in uids){
    sm <- expand.grid(x.rowid = which(x.ids == uid), 
                      y.rowid = which(y.ids == uid))
    if(nrow(sm) > 0){
      sm$diff <- y.dates[sm$y.rowid] - x.dates[sm$x.rowid]
      sm %<>% group_by(x.rowid) %>% filter(diff >= 0) %>% 
        filter(diff == min2(diff)) %>% as.data.frame()
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


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


#' Match and merge two data.frames by ID and consecutive date
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
#' @param by Any other variables to match by.
#' @param keep.unmatched When no match is found for a row, should it be kept or removed?
#' @param suffixes Character strings to append to the end of variable names that occur in both
#' \code{x} and \code{y}, to make them unique in the resulting \code{data.frame}.
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
#' x <- data.frame(  id=c(1,2,2,3,3,4,5),
#'                 date=c(1,1,3,1,3,1,9),
#'                info=letters[1:7])
#' y <- data.frame(  id=c(1,1,2,3,3,4),
#'                 date=c(2,0,2,4,2,0),
#'                 data=LETTERS[1:6])
#' match.merge(x=x,y=y,idvar="id",datevar="date",keep.unmatched="all")
#' match.merge(x=x,y=y,idvar="id",datevar="date",keep.unmatched="none")
#' 
match.merge <- function(x, y, 
                        idvar, 
                        datevar,
                        by=NULL, 
                        keep.unmatched = c("all","x","y","none"),
                        suffixes=c(".x",".y")){
  keep.unmatched <- match.arg(keep.unmatched)
  
  x.table<-unique(x[c(idvar,datevar)])
  y.table<-unique(y[c(idvar,datevar)])
  tablematches<-match.pps(x.table[[idvar]],x.table[[datevar]],
                          y.table[[idvar]],y.table[[datevar]])
  
  x.table$.rowmatch<-NA
  y.table$.rowmatch<-NA
  x.table$.rowmatch[tablematches$x.rowid]<-seq_len(nrow(tablematches))
  y.table$.rowmatch[tablematches$y.rowid]<-seq_len(nrow(tablematches))
  
  x <- left_join(x,x.table,by=c(idvar,datevar),keep=FALSE,na_matches="never")
  y <- left_join(y,y.table,by=c(idvar,datevar),keep=FALSE,na_matches="never")
  
  args <- list(x = x, y = y, 
               by = c(".rowmatch",by),
               na_matches = "never",
               suffix=suffixes)
  z <- do.call(case_match(keep.unmatched,
                          "all"~"full_join",
                          "x"~"left_join",
                          "y"~"right_join",
                          "none"~"inner_join"), 
               args)
  
  z$.rowmatch <- NULL
  
  z[[paste0(idvar,suffixes[1])]] <- 
    apply(z[paste0(idvar,suffixes)],1,function(x){x[first(which(!is.na(x)))]})
  colnames(z)[colnames(z)==paste0(idvar,suffixes[1])] <- idvar
  z[paste0(idvar,suffixes)] <- NULL
  
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

#' Coalesce mutually exclusive column sets into individual columns
#' 
#' This merges mutually exclusive variables in a \code{data.frame} that 
#' are distinguished from each other by a suffix or some other part of the variable name.
#' 
#' For example, when a survey was administered in multiple languages, and 
#' the responses to the same questions in each language are distinguished by a suffix 
#' like "_EN" or "_FR", this function can be used to 
#' coalesce these values into individual variables.
#' 
#' @param x A \code{data.frame} containing mutually exclusive column sets.
#' @param patterns Multiple regex patterns that only match the part of the variable name 
#' that makes members of the same mutually exclusive variable set unique (i.e. the variant tag).
#'
#' @returns A \code{data.frame} with all mutually exclusive column sets coalesced. 
#' The names of these coalesced columns are based on the original variable names but with 
#' \code{pattern} removed.
#' 
#' @export
#'
#' @examples
#' # Create a dataset where different individuals fill in different columns
#' mydataset <- 
#'    cbind(id=1:5,
#'          help_EN=c(1,NA,0,NA,NA),
#'          help_DE=c(NA,1,NA,1,NA),
#'          help_FR=c(NA,NA,NA,NA,0),
#'          success_EN=c(0,NA,0,NA,NA),
#'          success_DE=c(NA,1,NA,1,NA),
#'          success_FR=c(NA,NA,NA,NA,1)) |> 
#'            as.data.frame()
#'   
#' coalesce.columns(mydataset,c("_EN","_DE","_FR"))
#' 
coalesce.columns <- function(x, patterns){
  # Detect language per row, based on missingness
  nax <- is.na(x)
  langnas <- matrix(0,nrow=nrow(x),ncol=length(patterns)) |> as.data.frame()
  for(i in seq_along(patterns)){
    langnas[,i] <- nax[,str_detect(colnames(nax), patterns[i]), drop=FALSE] |> 
      rowMeans()
  }
  langs <- apply(langnas, 1, which.min)

  
  # Remove unused languages
  newcols <- colnames(x) |> 
    str_subset(paste0(patterns,collapse="|")) |>
    str_remove(paste0(patterns,collapse="|")) |> unique()
  if(any(newcols %in% colnames(x))){
    stop("Column(s) ",paste0(colnames(x)[newcols %in% colnames(x)],collapse=", "),
         " already exist in x. Please rename them.")
  }
  x[newcols] <- NA
  for(i in seq_along(patterns)){
    fromcols <- colnames(x) |> str_subset(patterns[i])
    tocols <- colnames(x) |> str_subset(patterns[i]) |> str_remove(patterns[i])
    x[langs == i, tocols] <- x[langs == i, fromcols]
  }
  oldcols <- colnames(x) |> str_subset(paste0(patterns,collapse="|"))
  x[oldcols] <- NULL
  
  return(x)
}






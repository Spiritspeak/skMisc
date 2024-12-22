#########
# Assorted code from the switch project that I haven't yet generalized


#This function finds out which session b occurred the earliest after each session a
#It gives a unique match out of set b for every item in set a that has a match
#Please note! set "a" must always have been administered before set "b"


min2<-function(x){ 
  if(length(x) == 0){ 
    NA 
  }else{ 
    min(x) 
  } 
}

match.pps<-function(x.ids, x.dates, y.ids, y.dates){
  uids <- unique(c(x.ids,y.ids))
  
  matchlist <- list()
  for(uid in uids){
    sm <- expand.grid(x.rowid=which(x.ids == uid), 
                      y.rowid=which(y.ids == uid))
    if(nrow(sm) > 0){
      sm$diff <- y.dates[sm$y.rowid] - x.dates[sm$x.rowid]
      sm %<>% group_by(x.rowid) %>% filter(diff >= 0) %>% 
        filter(diff == min2(diff)) %>% as.data.frame()
      matchlist[[length(matchlist) + 1]] <- sm[,-3]
    }
  }
  matchset <- do.call(rbind,matchlist)
  return(matchset)
}


match.merge<-function(x, x.ids, x.dates, 
                      y, y.ids, y.dates, 
                      by, ...){
  matches <- match.pps(x.ids, x.dates, y.ids, y.dates)
  x$.rowmatch <- NA
  y$.rowmatch <- NA
  x$.rowmatch[matches[[1]]] <- seq_along(matches[[1]])
  y$.rowmatch[matches[[2]]] <- seq_along(matches[[2]])
  z <- do.call(merge, list(x=x, y=y, by=c(by,".rowmatch"), ...))
  z$.rowmatch <- NULL
  return(z)
}


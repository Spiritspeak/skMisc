#########
# Assorted code from the switch project that I haven't yet generalized


#This function finds out which session b occurred the earliest after each session a
#It gives a unique match out of set b for every item in set a that has a match
#Please note! set "a" must always have been administered before set "b"
match.pps<-function(ia,da,ib,db){
  min2<-function(x){ if(length(x)==0){ NA }else{ min(x) } }
  ma<-rep(NA,length(ia))
  mb<-rep(NA,length(ib))
  
  ua<-unique(ia)
  ub<-unique(ib)
  
  matchset<-data.frame(aid=NA,bid=NA)[F,]
  for(u in ua){
    sm<-expand.grid(aid=which(ia==u),bid=which(ib==u))
    if(nrow(sm)>0){
      sm$diff<-db[sm$bid]-da[sm$aid]
      sm%<>%group_by(aid)%>%filter(diff>=0)%>%filter(diff==min2(diff))%>%as.data.frame()
      matchset<-rbind(matchset,sm[,-3])
    }
  }
  return(matchset)
}

match.merge<-function(a,ia,da,b,ib,db,by,...){
  matches<-match.pps(ia,da,ib,db)
  a$rowmatch<-NA
  b$rowmatch<-NA
  a$rowmatch[matches[[1]]]<-seq_along(matches[[1]])
  b$rowmatch[matches[[2]]]<-seq_along(matches[[2]])
  z<-do.call(merge,list(x=a,y=b,by=c(by,"rowmatch"),...))
  z$rowmatch<-NULL
  return(z)
}



# These are a bunch of functions that shift outliers closer to the mean or downweight them
# Don't really know what to do with this. It's not useful but I don't want to throw away the code!



# "Rein in" outliers by moving them closer to the mean
suppressor<-function(x,soft,hard, strength){
  hard<-hard-soft
  y<-x
  y[x< -soft]<- hard*(y[x< -soft]+soft)/(strength*abs(y[x< -soft]+soft)+hard)-soft
  y[x> soft]<- hard*(y[x> soft]-soft)/(strength*abs(y[x> soft]-soft)+hard)+soft
  return(y)
}

std.suppressor<-function(x,soft=2.5,hard=3,strength=1){
  m<-mean(x)
  s<-sd(x)
  suppressor((x-m)/s,soft,hard,strength)*s+m
}

loop.suppressor<-function(x,soft=2.5,hard=3,strength=1){
  while(any(abs(scale(x))>3)){
    x<-std.suppressor(x,soft,hard,strength)
  }
  x
}




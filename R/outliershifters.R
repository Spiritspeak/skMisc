
# These are a bunch of functions that shift outliers closer to the mean or downweight them
# Don't really know what to do with this. It's not useful but I don't want to throw away the code!

#' Downweight outliers
#' @description Computes weights; trials within certain bounds of the mean receive the maximum weight while trials
#' outside these bounds are downweighted to 0 or an optional minimum.
#'
#' @param x A numeric vector
#' @param mean An optional mean of the vector
#' @param s An optional standard deviation of the vector
#' @param sdist The number of standard deviations beyond which values should be downweighted
#' @param taper A number indicating how strongly values exceeding the standard deviation should taper off
#' @param scale How the weight vector should be scaled: "norm" sets the sum to 1, "max" sets the maximum to 1.
#' @param min A minimum weight. 
#'
#' @return A numeric vector of weights
#' @export
#'
#' @examples
logit.weightfun<-function(x,mean=mean(x),s=sd(x),
                          sdist=3,taper=10,
                          scale=c("max","norm"),min=0){
  scale <- match.arg(scale)
  zx <- (x-m)/s * taper
  out <- inv.logit((zx-sdist*taper)) * inv.logit((-zx-sdist*taper))
  out <- (out/max(out)) * (1-min) + min
  if(scale== "norm"){
    out < -out / sum(out)
  }
  return(out)
}

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







#' Vector transformations
#' 
#' Normalizing transformations for skewed variables
#'
#' @param x Vector to be transformed.
#' @param type Can be "none" (no transform), "log" (\code{log(x)}), "log1p" (\code{log(x+1)}), 
#' "sqrt" (\code{sqrt(x)}), "inv" (\code{1/x}), 
#' "logshift" (\code{log(x+shift)-log(shift)}),
#' "invshift" (\code{(1+shift)/(x+1+shift)}),
#' "bcp" (Box-Cox transform from [car::bcPower()]), 
#' "bcnp" (Shifted Box-Cox transform from [car::bcnPower()]),
#' or "yj" (Yeo-Johnson transformation from [car::yjPower()]). Non-positive values in \code{x} are
#' only allowed in "none", "bcnp" and "yj".
#' @param par A parameter that determines the shape of the transformation. If no parameter is given,
#' the function finds the parameter that gives the most normally distributed values.
#' For \code{"logshift"} and \code{"invshift"}, this is the shift; 
#' for \code{"bcp"} and \code{"yj"}, this is the lambda parameter; 
#' for \code{"bcnp"}, this is two values: the lambda and gamma parameters.
#' For all other transformations, \code{par} is ignored.
#'
#' @returns A transformed vector. 
#' Attribute \code{"transform.method"} indicates the used method, and 
#' \code{"transform.par"} indicates the transformation parameter(s).
#' @export
#'
#' @examples
#' # Normality comparison
#' set.seed(1)
#' u<-rgamma(1000,2)
#' logtrans<-trans(u,"log")
#' shapiro.test(logtrans) # Not normal
#' logshifttrans<-trans(u,"logshift")
#' shapiro.test(logshifttrans) # Not normal but less bad
#' bctrans<-trans(u,"bcp")
#' shapiro.test(bctrans) # Normal
#'
#' u<-rt(100,2) 
#' yjtrans<-trans(u,"yjp")
#' shapiro.test(yjtrans) # Not normal; still leptokurtic
#' 
#' 
trans <- function(x,type=c("none","log","log1p","sqrt","inv","logshift","invshift","bcp","bcnp","yjp"),par=NULL){
  type<-match.arg(type)
  if(type=="none"){
    out<-x
  }else if(type=="log"){
    out<-log(x)
  }else if(type=="log1p"){
    out<-log1p(x)
  }else if(type=="sqrt"){
    out<-sqrt(x)
  }else if(type=="inv"){
    out<-1/x
  }else if(type=="logshift"){
    logshift<-function(x,p){
      log(x+p)-if(p==0){0}else{log(p)}
    }
    if(is.null(par)){
      scorer<-function(p){
        -shapiro.test(logshift(x,exp(p)))$statistic
      }
      par<-exp(optim(par=0,fn=scorer,method="BFGS")$par)
    }
    out<-logshift(x,par[1])
  }else if(type=="invshift"){
    invshift<-function(x,p){
      (p+1)/(x+p+1)
    }
    if(is.null(par)){
      scorer<-function(p){
        -shapiro.test(invshift(x,exp(p)))$statistic
      }
      par<-exp(optim(par=0,fn=scorer,method="BFGS")$par)
    }
    out<-invshift(x,par[1])
  }else if(type=="bcp"){
    if(is.null(par)){
      par<-car::powerTransform(x,family="bcPower")$lambda
    }
    out<-car::bcPower(U=x,lambda=par[1])
  }else if(type=="bcnp"){
    if(is.null(par)){
      res<-car::powerTransform(x,family="bcnPower")
      par<-c(res$lambda,res$gamma)
    }
    out<-car::bcnPower(U=x,lambda=par[1],gamma=par[2])
  }else if(type=="yjp"){
    if(is.null(par)){
      par<-car::powerTransform(x,family="yjPower")$lambda
    }
    out<-car::yjPower(U=x,lambda=par[1])
  }
  attr(out,"transform.method")<-type
  attr(out,"transform.par")<-par
  return(out)
}



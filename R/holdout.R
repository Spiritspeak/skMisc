
# goals right now: 
# add help file to every function
# ensure these functions behave normally with already insignificant models

r2t<-function(r,df){
  sqrt(r^2*df / (1-r^2))
}

cor.influence<-function(x,y){
  x<-x-mean(x)
  y<-y-mean(y)
  x*y-(x^2+y^2)/2*cor(x,y)
}

#' Correlation holdouts
#' 
#' Repeatedly remove the most influential observation using
#' influence functions until a correlation coefficient is 
#' either no longer significant or has changed its sign
#'
#' @param x,y Numeric vectors of variables to correlate
#' @param goal Objective of observation removal: 
#' "nsig" (no longer significant) or "flip" (change of sign).
#' @param method Correlation method, either pearson or spearman
#' @param alpha Alpha level for significance testing
#' @param verbose Logical; if TRUE, will produce more output
#'
#' @return A list with...
#' - h: number of observations needed to reach objective
#' - h.prop: h but as a percentage
#' - final_r: the correlation after removal of observations
#' - final_p: the p value after removal of observations
#' - omit: logical vector denoting which observations were omitted to obtain result
#' @export
#'
#' @examples
#' 
#' xval<-rnorm(100)
#' cor.holdout(x=rnorm(100)+xval,y=rnorm(100)+xval,goal="flip")
#' cor.holdout(x=rnorm(100),y=rnorm(100)+xval,goal="nsig",method="spearman")
#' 
cor.holdout<-function(x,y,
                      goal=c("nsig","flip"),
                      method=c("pearson","spearman"),
                      alpha=.05,verbose=F){
  
  # prep for spearman
  method<-match.arg(method)
  if(method=="spearman"){
    x<-rank(x)
    y<-rank(y)
  }
  
  goal<-match.arg(goal)
  
  
  rval<-origrval<-cor(x,y)
  origsign<-sign(origrval)
  incl<-rep(TRUE,length(x))
  
  # Set goal-reached checker
  if(goal=="nsig"){
    checkGoal<-function(r,df,...){
      t<-r2t(r=r,df=df)
      p<-2*pt(q=abs(t),df=df,lower.tail=F)
      return(p>alpha)
    }
  }else
    if(goal=="flip"){
      checkGoal<-function(r,...){
        return(sign(r)!=origsign)
      }
    }
  
  # Function for finding the next entry to delete
  findNextToDelete<-function(){
    inf<-cor.influence(x[incl],y[incl])
    return(which.max(inf*origsign))
  }
  
  while(T){
    # Evaluate termination conditions
    rval<-cor(x[incl],y[incl])
    currn<-sum(incl)
    if(currn<=2){
      success<-F
      break 
    }else
      if(checkGoal(r=rval,df=sum(incl)-2)){
        success<-T
        break
      }
    
    idx<-findNextToDelete()
    incl[incl][idx]<-F
  }
  
  lastdf<-sum(incl)-2
  
  outlist<-list(final_r=rval,
                final_p=(rval |> r2t(df=lastdf) |> abs() |> pt(df=lastdf,lower.tail=F))*2,
                success=success,
                h=sum(!incl),
                omit=!incl,
                h.prop=sum(!incl)/(length(x)-2))
  out<-structure(.Data=outlist,class="cor.holdout")
  return(out)
}

#' @rdname cor.holdout
#' @param x A cor.holdout object.
#' @param ... Ignored.
#' @export
print.cor.holdout<-function(x,...){
  # cat(sep="",
  #     "h = ",x$h,"; h.prop = ",x$h.prop,"\n",
  #     "final r = ",x$final_r,"; final p = ",x$final_p)
  x[c("h","h.prop","final_r","final_p")] |> unlist() |> print()
}

registerS3method("print","cor.holdout",print.cor.holdout)

##############
# lm.holdout #
##############




#' Linear Regression Holdouts
#' Repeatedly remove the most influential observation using
#' influence functions until a model term is 
#' either no longer significant or has changed its sign
#' 
#' @param model A \code{lm} model
#' @param goal Objective of observation removal: 
#' "nsig" (no longer significant) or "flip" (change of sign).
#' @param terms Names of model terms to compute holdout statistics for. 
#' When NULL, it defaults to all model terms.
#' @param alpha 
#' Maximum p value for significance (only relevant when aiming to make value insignificant)
#' @param verbose Should the function generate verbose output? Defaults to no.
#'
#' @return A list, containing
#' - \code{holdouts}: a \code{data.frame} where each row is a predictor; as for the columns: 
#'  - \code{h} - the number of observations needed to reach the goal
#'  - \code{h.prop} - the percentage of excludable observations that were excluded to reach the goal
#'  - \code{final.beta} - the final beta value of the predictor after exclusions
#'  - \code{final.p} - the final p value of the predictor after exclusions
#' - \code{exclusion.matrix}
#' - \code{model}
#' 
#' @md
#' @export
#'
#' @examples
#' 
#' Sigma <- diag(1,nrow=3)
#' Sigma[1,-1] <- Sigma[-1,1]<- -.5
#' fakedata <- MASS::mvrnorm(n=1000,mu=c(0,0,0),Sigma=Sigma,empirical=TRUE)
#' fakedata <- as.data.frame(fakedata)
#' 
#' lmmod <- lm(V1~V2*V3,data=fakedata)
#' 
#' lmmod.h <- lm.holdout(lmmod,goal="nsig")
#' print(lmmod.h)
#' 
#' lmmod.h2 <- lm.holdout(lmmod,goal="flip")
#' print(lmmod.h2)
#' 
lm.holdout<-function(model,goal=c("nsig","flip"),terms=NULL,alpha=.05,verbose=FALSE){
  
  origcoefs<-coef(model)
  origcoef.signs<-sign(origcoefs)
  
  goal<-match.arg(goal)
  if(goal=="nsig"){
    checkGoal<-function(mod,pred){
      summary(mod)$coefficients[pred,"Pr(>|t|)"] > alpha
    }
  }
  if(goal=="flip"){
    checkGoal<-function(mod,pred){
      sign(coef(mod)[pred]) != origcoef.signs[pred]
    }
  }
  
  if(is.null(terms)){
    terms<-terms(model) |> attr("factors") |> colnames()
  }
  
  mf<-model.frame(model)
  
  lm.holdout.single<-function(pred){
    success<-NA
    currmod<-model
    incl<-rep(T,nrow(mf))
    while(T){
      
      if(sum(incl)<2){
        success<-F
        break
      }
      if(checkGoal(currmod,pred)){
        success<-T
        break
      }
      
      infs<-lm.influence(currmod)$coefficients[,pred]
      idx<-which.max(infs*origcoef.signs[pred])
      incl[incl][idx]<-F
      currmod<-update(object=model,.~.,data=mf[incl,])
      if(verbose){
        message(pred,": ",sum(incl)," remaining")
      }
    }
    
    outlist<-list(final_beta=coef(currmod)[pred],
                  final_p=summary(currmod)$coefficients[pred,"Pr(>|t|)"],
                  success=success,
                  h=sum(!incl),
                  rem=!incl,
                  h.prop=sum(!incl)/(nrow(mf)-2))
  }
  
  
  output<-list()
  for(pred in terms){
    output[[pred]]<-lm.holdout.single(pred)
  }
  out<-list(holdouts=data.frame(h=sapply(output,\(x)x$h),
                                h.prop=sapply(output,\(x)x$h.prop),
                                final_beta=sapply(output,\(x)x$final_beta),
                                final_p=sapply(output,\(x)x$final_p)),
            exclusion.matrix=sapply(output,\(x)x$rem),
            model=model)
  out<-structure(.Data=out,class="lm.holdout")
  return(out)
}

print.lm.holdout<-function(x,...){
  cat("Holdout values for linear regression model\n")
  print(formula(x$model))
  print(x$holdouts)
}

registerS3method("print","lm.holdout",print.lm.holdout)





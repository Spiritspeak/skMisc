
# TODO:
# plot.lm.holdout needs plot titles with the h value



# Remains undocumented, helper function.
r2t<-function(r,df){
  sqrt(r^2*df / (1-r^2))
}

# Seeks most bivariate-outlying pair of values
cor.outlying <- function(x, y){
  scale.vector(x) * scale.vector(y)
}

#' Correlation holdouts
#' 
#' Repeatedly remove the observation whose removal will most strongly 
#' move the correlation coefficient to zero, until it is 
#' either no longer significant or has changed its sign
#'
#' @param x,y Numeric vectors of variables to correlate. 
#' In case of the print function, \code{x} is a \code{cor.holdout} object to print. 
#' @param goal Objective of observation removal: 
#' "nsig" (no longer significant) or "flip" (change of sign).
#' @param method Correlation method, either "pearson" or "spearman"
#' @param alpha Alpha level for significance testing
#' @param verbose Should the function generate verbose output? Defaults to FALSE.
#'
#' @return A list with...
#' - h: number of observations needed to reach objective
#' - h.prop: h but as a percentage of the total degrees of freedom 
#' (total number of observations - 2)
#' - final.r: the correlation after removal of observations
#' - final.p: the p value after removal of observations
#' - omit: integer vector denoting which observations were omitted 
#' at which iteration to obtain the result. A zero value indicates that 
#' the observation was excluded beforehand due to presence of NAs.
#' 
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' xval<-rnorm(100)
#' 
#' # Flip sign
#' cor.holdout(x=rnorm(100)+xval,y=rnorm(100)+xval,goal="flip")
#' 
#' # Make insignificant
#' cor.holdout(x=rnorm(100)+xval,y=rnorm(100)+xval,goal="nsig",method="spearman")
#' 
#' # No true correlation
#' cor.holdout(x=rnorm(100),y=rnorm(100),goal="nsig")
#' 
cor.holdout<-function(x,y,
                      goal=c("nsig","flip"),
                      method=c("pearson","spearman"),
                      alpha=.05,verbose=FALSE){
  
  # prep for spearman
  method<-match.arg(method)
  if(method=="spearman"){
    x<-rank(x)
    y<-rank(y)
  }
  goal<-match.arg(goal)
  
  rval<-origrval<-cor(x,y)
  origsign<-sign(origrval)
  excl<-rep(NA,length(x))
  excl[is.na(x) | is.na(y)]<-0
  incl<-is.na(excl)
  iter<-0
  
  # Set goal-reached checker
  if(goal=="nsig"){
    checkGoal<-function(r,df,...){
      t<-r2t(r=r,df=df)
      p<-2*pt(q=abs(t),df=df,lower.tail=F)
      return(p>alpha)
    }
  }else if(goal=="flip"){
    checkGoal<-function(r,df,...){
      return(sign(r)!=origsign)
    }
  }
  
  # Function for finding the next entry to delete
  findNextToDelete<-function(){
    inf<-cor.outlying(x[incl],y[incl])
    return(which.max(inf*origsign))
  }
  
  while(T){
    iter<-iter+1
    # Evaluate termination conditions
    currn<-sum(incl)
    if(currn<=2 || var(x[incl])==0 || var(y[incl])==0){
      success<-F
      break 
    }else{
      rval<-cor(x[incl],y[incl])
      if(checkGoal(r=rval,df=sum(incl)-2)){
        success<-T
        break
      }
    }
    idx<-findNextToDelete()
    excl[incl][idx]<-iter
    incl[incl][idx]<-F
  }
  
  h<-sum(excl>0,na.rm=T)
  lastdf<-h-2
  
  # Prepare output
  outlist<-list(final.r=rval,
                final.p=(rval |> r2t(df=lastdf) |> abs() |> pt(df=lastdf,lower.tail=F))*2,
                success=success,
                h=h,
                omit=excl,
                h.prop=h/(sum(!is.na(x) & !is.na(y))-2))
  out<-structure(.Data=outlist,class="cor.holdout")
  return(out)
}

#' @rdname cor.holdout
#' @param ... Ignored.
#' @export
#' @method print cor.holdout
#' 
print.cor.holdout<-function(x,...){
  # cat(sep="",
  #     "h = ",x$h,"; h.prop = ",x$h.prop,"\n",
  #     "final r = ",x$final.r,"; final p = ",x$final.p)
  x[c("h","h.prop","final.r","final.p")] |> unlist() |> print()
}

registerS3method("print","cor.holdout",print.cor.holdout)

##############
# lm.holdout #
##############


#' Linear Regression Holdouts
#' 
#' Repeatedly remove the observation whose removal will move a model term 
#' most strongly to zero, until it is either no longer significant 
#' or has changed its sign.
#' 
#' @param model A \code{lm} model.
#' @param goal Objective of observation removal: 
#' "nsig" (no longer significant) or "flip" (change of sign).
#' @param terms Names of model terms to compute holdout statistics for. 
#' When \code{NULL}, it defaults to all model terms.
#' @param alpha 
#' Maximum p value for significance (only relevant when aiming to make value insignificant).
#' @param verbose Should the function generate verbose output? Defaults to \code{FALSE}.
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
#' @details
#' Datapoints to exclude are found with the function \code{[stats::lm.influence()]}.
#' 
#' The plot function displays a scatterplot of each predictor with the dependent variable, 
#' with the excluded datapoints colored differently. 
#' 
#' @md
#' @author Sercan Kahveci
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
#' plot(lmmod.h)
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
    
    outlist<-list(final.beta=coef(currmod)[pred],
                  final.p=summary(currmod)$coefficients[pred,"Pr(>|t|)"],
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
                                final.beta=sapply(output,\(x)x$final.beta),
                                final.p=sapply(output,\(x)x$final.p)),
            exclusion.matrix=sapply(output,\(x)x$rem),
            model=model)
  out<-structure(.Data=out,class="lm.holdout")
  return(out)
}

#' @rdname lm.holdout
#' @param x A \code{lm.holdout} object.
#' @param ... Additional plotting parameters for \code{[base::plot()]}; ignored for \code{print}.
#' @export
#' @method print lm.holdout
#' 
print.lm.holdout<-function(x,...){
  cat("Holdout values for linear regression model\n")
  print(formula(x$model))
  print(x$holdouts)
}

registerS3method("print","lm.holdout",print.lm.holdout)

#' @rdname lm.holdout
#' @export
#' 
plot.lm.holdout<-function(x,...){
  plotterms<-colnames(x$exclusion.matrix)
  modmat<-attr(terms(x$model),"factors")
  dv<-rownames(modmat)[1]
  dat<-x$model$model
  
  oldmfrow<-par("mfrow")
  rowz<-ceiling(sqrt(length(plotterms)))
  colz<-ceiling(length(plotterms)/rowz)
  par(mfrow=c(rowz,colz))
  
  for(cn in plotterms){
    cns<-strsplit(cn,":")[[1]]
    if(length(cns)==1){
      xvar<-dat[[cns]]
    }else{
      xvar<-Reduce(`*`,dat[cns])
    }
    plot(x=xvar,y=dat[[dv]],col=1+x$exclusion.matrix[,cn]*2,
         xlab=cn,ylab=dv,...)
  }
  par(mfrow=oldmfrow)
  return()
}

registerS3method("plot","lm.holdout",plot.lm.holdout)




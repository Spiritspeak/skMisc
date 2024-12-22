

#stolen from stackoverflow
#' Plot highlighted text
#'
#' @param x x position
#' @param y y position
#' @param s text
#' @param bg highlight color
#'
#' @export
#'
#' @examples
#' plot(mtcars$mpg,mtcars$wt,col=mtcars$cyl)
#' hilight(27,2.5,"Light and\nefficient")
#' hilight(17,4.5,"Heavy and\ninefficient")
#' 
hilight<-function(x,y,s, bg="yellow") {
  text.width <- strwidth(s)
  text.height <- strheight(s)
  rect(x,y,x+text.width,y+text.height,col=bg,border=NA)
  text(x,y,s,adj=c(0,0))
}


#FullAutocorPlot
#' Per-subject Autocorrelation Plotting
#'
#' @param ds a dataset
#' @param ppvar name of the variable indicating participant ID
#' @param rtvar name of the variable indicating reaction time
#' @param scope numeric, the maximum lag at which to compute autocorrelation.
#'
#' @export
#'
#' @examples
#' AutocorPlot(ds=ToothGrowth,ppvar="supp",rtvar="len",scope=10)
#' 
AutocorPlot<-function(ds,ppvar,rtvar,scope=64){
  if(missing(ds)){
    ds<-data.frame(ppvar=ppvar,rtvar=rtvar,stringsAsFactors = F)
    ppvar <- "ppvar"
    rtvar <- "rtvar"
  }else{
    ds <- as.data.frame(ds)
  }
  ppvec <- unique(ds[,ppvar]) %>% as.vector()
  npp <- length(ppvec)
  cormat <- matrix(nrow=scope,ncol=npp)
  ct <- 0
  for(pp in ppvec){
    ct <- ct+1
    tds <- ds[ds[[ppvar]] == pp,]
    autocor <- numeric(scope)
    for(i in seq_len(scope)){
      autocor[i] <- cor(tds[,rtvar],dplyr::lag(tds[,rtvar],i),use="complete.obs")
    }
    cormat[,ct] <- autocor
  }
  colnames(cormat) <- ppvec
  plot(rowMeans(cormat),type="l",xlab="lag",ylab="autocorrelation",
       main=paste("Variable",rtvar,"autocorrelation"),
       ylim=c(min(cormat),max(cormat)))
  for(i in seq_len(npp)){ lines(x=cormat[,i],pch=".",col=rgb(0,0,0,0.1),type="l") }
  abline(h=0,col="gray")
  return(invisible(cormat))
}


#TransformPlots
#' Data Transformation Plots
#' 
#' @description Visualize how different transformations of the data 
#' will fit to a normal distribution.
#' @param x A numeric vector.
#'
#' @export
#' 
#' @examples
#' TransformPlots(mtcars$disp)
#' 
TransformPlots<-function(x){
  oldPars<-par("mfrow","mar")
  par(mfrow=c(2,2),mar=c(3,2,3,1))
  inv<-function(x){ 1/x }
  nothing<-function(x){ x }
  titles<-c("Untransformed","Log-transformed","Sqrt-transformed","Inverse-transformed")
  transforms<-c(nothing,log,sqrt,inv)
  for(i in 1:4){
    y<-do.call(transforms[[i]],list(x))
    ks<-try(ks.test(y,"pnorm"),silent=T)
    plottitle<-paste(titles[i],
                     "\nKS-test D = ",round(ks$statistic,digits=3),
                     ", p = ",round(ks$p.value,digits=3))
    car::qqp(y,"norm",main=plottitle)
  }
  par(oldPars)
}




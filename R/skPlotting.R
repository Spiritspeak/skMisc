

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


#' Multi-series Autocorrelation Plotting
#'
#' @param x Variable to compute the autocorrelation of.
#' @param index Grouping ID variable of the same length as \code{X}.
#' @param lag.max The maximum lag at which to compute autocorrelation.
#' @param plot Should a plot be made? Default is \code{TRUE}.
#'
#' @export
#'
#' @examples
#' mycors<-AutocorPlot(x=ToothGrowth$len,index=ToothGrowth$supp,lag.max=10)
#' 
AutocorPlot<-function(x,index,lag.max=64,plot=TRUE){
  mylag<-function(x,i){
    len<-length(x)
    if(i<=len){
      idx<-seq_len(len-i)
    }
    c(rep(NA,i),x[idx])
  }
  indices <- as.vector(unique(index))
  nindices <- length(indices)
  cormat <- matrix(nrow=lag.max,ncol=nindices,dimnames=list(NULL,indices))
  for(currindex in indices){
    currx <- x[index == currindex]
    autocor <- numeric(lag.max)
    for(i in seq_len(lag.max)){
      autocor[i] <- cor(currx,mylag(currx,i),use="complete.obs")
    }
    cormat[,currindex] <- autocor
  }
  if(plot){
    plot(rowMeans(cormat),type="l",xlab="lag",ylab="autocorrelation",
         ylim=c(min(cormat),max(cormat)))
    for(i in seq_len(nindices)){ 
      lines(x=cormat[,i],col=rgb(0,0,0,1/sqrt(nindices)))
    }
    abline(h=0,col="gray")
    lines(x=rowMeans(cormat))
  }
  return(invisible(t(cormat)))
}

# Use for testing.
# library(dplyr)
# newdata <- erotica %>% dplyr::arrange(subject,blocknum,trialnum) %>%
#   mutate(index=paste0(subject,"-",blocknum))
# AutocorPlot(x=newdata$RT,index=newdata$index,lag.max=30)

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




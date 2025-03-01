
colorblind1 <- c("#050D03", "#D5AF00", "#4714D9", "#00EF41", 
                 "#00A2BF", "#D45E00", "#FF29F7")

#' @name ggplot.themes
#' @title extra \code{ggplot2} themes
#' @description A number of ready-to-use ggplot2 themes for scientific manuscripts.
#' @author Sercan Kahveci
#' 
NULL

#' @describeIn ggplot.themes Based on the plot design style of prof. Diane Pecher.
#' @export
#' @examples
#' ggplot(mtcars,aes(x=mpg,y=wt)) + geom_point() + theme_pecher()
#' 
theme_pecher <- function(){
  theme_bw() + 
    theme(text = element_text(size=14, family="serif"),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(-.25,"cm"),
          axis.text = element_text(colour = rgb(0,0,0),size=12),
          axis.text.x=element_text(margin = margin(t=.4,b=.1,r=.1,l=.1,unit="cm")),
          axis.text.y=element_text(margin = margin(t=.1,b=.1,r=.4,l=.1,unit="cm")),
          strip.background=element_blank())
}

#' @describeIn ggplot.themes Theme appropriate for APA manuscripts.
#' @export
#' @examples
#' ggplot(mtcars,aes(x=mpg,y=wt)) + geom_point() + theme_apa()
#' 
theme_apa <- function(){
  theme_bw() + 
    theme(legend.position="bottom",panel.grid=element_blank(),
          panel.border = element_blank(),
          axis.line=element_line(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold",size=unit(14,"pt")),
          axis.title=element_text(face="bold",size=unit(14,"pt"),margin=unit(rep(0,4),"pt")),
          legend.text = element_text(size=unit(14,"pt")),
          axis.ticks.length = unit(-5,"pt"),
          axis.ticks = element_line(color = "black"),
          axis.text=element_text(color="black",size=unit(13,"pt")),
          axis.text.x=element_text(margin=unit(c(8,0,0,0),"pt")),
          axis.text.y=element_text(margin=unit(c(0,8,0,0),"pt")))
}

#stolen from stackoverflow
#' Plot highlighted text
#'
#' @param x x position
#' @param y y position
#' @param s text
#' @param bg highlight color
#'
#' @export
#' @author Sercan Kahveci
#'
#' @examples
#' plot(mtcars$mpg,mtcars$wt,col=mtcars$cyl)
#' hilight(27,2.5,"Light and\nefficient")
#' hilight(17,4.5,"Heavy and\ninefficient")
#' 
hilight <- function(x, y, s, bg="yellow") {
  text.width <- strwidth(s)
  text.height <- strheight(s)
  rect(x ,y, x + text.width, y + text.height, col=bg, border=NA)
  text(x, y, s, adj=c(0, 0))
}


#' Multi-series Autocorrelation Plotting
#'
#' @param x Variable to compute the autocorrelation of.
#' @param index Grouping ID variable of the same length as \code{X}.
#' @param lag.max The maximum lag at which to compute autocorrelation.
#' @param plot Should a plot be made? Default is \code{TRUE}.
#'
#' @export
#' @author Sercan Kahveci
#'
#' @examples
#' mycors<-AutocorPlot(x=ToothGrowth$len,index=ToothGrowth$supp,lag.max=10)
#' 
AutocorPlot <- function(x, index, lag.max=64, plot=TRUE){
  mylag <- function(x,i){
    len <- length(x)
    if(i <= len){
      idx <- seq_len(len-i)
    }
    c(rep(NA,i), x[idx])
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


#' Data Transformation Plots
#' 
#' @description Visualize how different transformations of the data 
#' will fit to a normal distribution.
#' @param x A numeric vector.
#'
#' @export
#' @author Sercan Kahveci
#' 
#' @examples
#' TransformPlots(mtcars$disp)
#' 
TransformPlots <- function(x){
  oldPars <- par("mfrow","mar")
  par(mfrow=c(2,2), mar=c(3,2,3,1))
  inv <- function(x){ 1/x }
  nothing <- function(x){ x }
  titles <- c("Untransformed", "Log-transformed", "Sqrt-transformed", "Inverse-transformed")
  transforms <- c(nothing, log, sqrt, inv)
  for(i in 1:4){
    y <- do.call(transforms[[i]],list(x))
    ks <- try(ks.test(y,"pnorm"), silent=T)
    plottitle <- paste(titles[i],
                       "\nKS-test D = ", round(ks$statistic, digits=3),
                       ", p = ",round(ks$p.value, digits=3))
    car::qqp(y, "norm", main=plottitle)
  }
  par(oldPars)
}




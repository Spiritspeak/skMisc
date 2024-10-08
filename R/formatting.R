#' @rdname ggplot.themes
NULL

#' @describeIn ggplot.themes Based on the plot design style of prof. Diane Pecher.
#' @export
#' @examples
#' ggplot(mtcars,aes(x=mpg,y=wt)) + geom_point() + theme_pecher()
#' 
theme_pecher<-function(){
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

#' @describeIn ggplot.themes Theme appropriate for APA manuscripts
#' @examples
#' ggplot(mtcars,aes(x=mpg,y=wt)) + geom_point() + apatheme()
#' 
apatheme<-function(){
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

#' Drop leading zeros
#' Remove leading zeroes and return as a character vector.
#'
#' @param x Numeric vector to remove leading zeros from
#'
#' @return
#' @author Adapted from [m.evans on StackOverflow](https://stackoverflow.com/a/48431993).
#' @export
#' @md
#'
#' @examples
#' dropLeadingZero(c(-1,0,1,0.5,-0.5,1.5,-1.5))
dropLeadingZero <- function(x){
  gsub("(?<![0-9])0+(?=\\.)", "", x, perl = TRUE)
}

dropLeadingZero_old <- function(x){
  xnew <- c()
  for(i in x){
    if(isTRUE(i==0)){
      xnew <- c(xnew,"0")
    } else if (isTRUE(i>=1) | isTRUE(i<=-1)){
      xnew <- c(xnew, as.character(i))
    } else
      xnew <- c(xnew, gsub("(?<![0-9])0+(?=\\.)", "", i, perl = TRUE))
  }
  return(xnew)
}


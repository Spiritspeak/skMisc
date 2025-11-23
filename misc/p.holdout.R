
# There's a true correlation
xval<-rnorm(100)
thx<-rnorm(100)+xval
thy<-rnorm(100)+xval

# Outlier-driven correlation
thx <- rnorm(100)
thy <- rnorm(100)
thx[1:2] <- thy[1:2] <- 10

origh<-cor.holdout(x=thx,y=thy,goal="flip")

hhh<-permute.cor.holdout(x=thx,y=thy,goal="flip",iters=10000)
hist(hhh)
mean(hhh >= origh$h)

hhh2<-norm.cor.holdout(r=cor(thx,thy),n=length(thx),goal="flip",iters=1000)
hist(hhh2,breaks=length(unique(hhh2)))
mean(hhh2 <= origh$h)

# Consider mid-p adjustment for p-values, as well as adding 1 to both numerator and denominator

norm.cor.holdout<- function(r, n, iters=1000,
                            goal=c("nsig","flip"),
                            method=c("pearson","spearman"),
                            alpha=.05){
  hvec <- 
  replicate(iters,{
    fd <- MASS::mvrnorm(n=n,mu=c(0,0),Sigma=matrix(c(1,r,r,1),nrow=2),empirical=T)
    cor.holdout(x=fd[,1],y=fd[,2],goal=goal,method=method,alpha=alpha)$h
  })
  
}

permute.cor.holdout<-function(x, y,
                              iters=1000,
                              goal=c("nsig","flip"),
                              method=c("pearson","spearman"),
                              alpha=.05){
  unmissing<-!is.na(x) & !is.na(y)
  xt <- x[unmissing]
  yt <- y[unmissing]
  
  hvec <- replicate(iters,
                    cor.holdout(x=sample(xt),y=yt,goal=goal,method=method,alpha=alpha)$h)
  return(hvec)
}

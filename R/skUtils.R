.onLoad<-function(libname, pkgname){
  packageStartupMessage("Thank you for loading skMisc v0.01")
}

#I don't want to import rlang, so it will be done this way instead.
args2strings <- function(...) sapply(substitute({ ... })[-1], deparse)

#' clamp
#'
#' @param val The vector/matrix to clamp
#' @param minval Minimum value; all lower values are clamped to this value
#' @param maxval Maximum value; all higher values are clamped to this value
#'
#' @return Clamped vector.
#' @export 
#'
#' @examples clamp(0:10,2,8)
clamp <- function(val,minval,maxval){
  val[val<minval]<-minval
  val[val>maxval]<-maxval
  val
}

#' coerce a vector to contain only TRUE and FALSE
#'
#' @param x Numeric/logical vector/matrix to coerce into TRUE/FALSE
#' @param default default returned value if NULL or NA is encountered
#'
#' @return logical vector or matrix with only T and F
#' @export
#'
#' @examples 
#' coerce(NULL)
#' # FALSE
#' 
#' coerce(c(T,F,NA,NA,T))
#' # T F F F T
#' 
#' coerce(matrix(c(T,T,F,F,NA,NA),nrow=2))
#' #     [,1]  [,2]  [,3]
#' #[1,] TRUE FALSE FALSE
#' #[2,] TRUE FALSE FALSE
coerce<-function(x,default=FALSE){
  if(length(x)==0){ return(default) }
  x<-x!=0
  x[is.na(x)]<-default
  x
}

#' Crunch Outliers
#'
#' @param x Numeric vector to remove outliers from
#' @param DS A positive numeric value. If value exceeds this many standard deviations, it is counted as an outlier
#' @param hardlimit A numeric vector with two values. If set, values below the first value and above the second will be counted as outliers, and means/standard deviations will be computed from values within these bounds only.
#'
#' @return Vector with outlying values set to NA
#' @export
OLcrunch<-function(x,DS=3,hardlimit=NULL){
  if(!missing(hardlimit)){
    x[x<hardlimit[1] | x>hardlimit[2]]<-NA
  }
  m<-mean(x,na.rm=T)
  s<-sd(x,na.rm=T)
  x[(x>m+s*DS) | (x<m-s*DS)]<-NA
  return(x)
}

#' Smooth a numeric vector using a moving window algorithm
#'
#' @param vect 
#' @param width Over how many values should the vector be averaged?  
#' @param both.sides If TRUE (default), takes the mean of \code{width} values before and after the current index. If FALSE, only takes values ahead of the current index.
#'
#' @return Smoothed numeric vector
#' @export
#'
#' @examples temp<- smoothvect(beaver1$temp)
#' plot(temp,type="l")
smoothvect<-function(vect,width=2,both.sides=T,alg=c("mean","gauss")){
  output<-numeric()
  
  normsum<-function(x){x/sum(x)}
  
  if(alg=="mean"){
    for(i in seq_len(length(vect))){
      output[i]<-mean(vect[ max(i-width*both.sides,1):min(i+width,length(vect))],na.rm=T)
    }
  }
  if(alg=="gauss"){
    winvect<- (-width*both.sides):(width)
    window<-dnorm(winvect/length(vect)*3)
    for(i in seq_len(length(vect))){
      output[i]<-sum(
        vect[max(i-width*both.sides,1):min(i+width,length(vect))] * 
          normsum(window[ (max(i-width,1)-i):(min(i+width,length(vect))-i) +width+1])
        ,na.rm=T)
    }
  }
  output
}


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

#' Initiate an empty data frame
#'
#' @param namelist A character vector of column names.
#'
#' @return A data.frame with 0 rows.
#' @export
df.init<-function(namelist){
  setNames(data.frame(matrix(ncol = length(namelist), nrow = 0)), namelist)
}

#' Change classes of columns in a data.frame
#' @description \code{retype()} changes the class of specific columns; \code{retype_all()} changes the class of all columns of a given class.
#'
#' @param df a data frame
#' @param ... Unquoted column names, paired with the desired class, e.g. 
#' 
#' \code{age = numeric(), language = character()}
#'
#' @export
#'
#' @examples 
#' sapply(ToothGrowth,class)
#' #      len      supp      dose 
#' #"numeric"  "factor" "numeric" 
#' NewToothGrowth <- retype(ToothGrowth, supp = character(), dose = factor())
#' sapply(NewToothGrowth,class)
#' #      len        supp        dose 
#' #"numeric" "character"    "factor" 
retype<-function(df, ...){
  args<-list(...)

  varnames<-names(args)
  vartypes<-sapply(args, class)

  effcols<-names(df)[names(df) %in% varnames]

  for(effcol in effcols){
    df[,effcol]<-as(df[,effcol], vartypes[which(varnames==effcol)])
  }
  return(df)
}
#' @rdname retype
#' @param df A data.frame
#' @param from An empty vector of the class to convert from, or a string. Columns sharing the class of argument \code{from} will be converted to the class of argument \code{to}.
#' @param to An empty vector of the class to convert to, or a string. Columns sharing the class of argument \code{from} will be converted to the class of argument \code{to}.
#'
#' @export
#'
#' @examples 
#'
#' sapply(mtcars,class)
#' #       mpg       cyl      disp        hp      drat        wt
#' # "numeric" "numeric" "numeric" "numeric" "numeric" "numeric"
#' #      qsec        vs        am      gear      carb 
#' # "numeric" "numeric" "numeric" "numeric" "numeric" 
#' 
#' newmtcars <- retype_all(mtcars,"numeric","character")
#' sapply(newmtcars,class)
#' #         mpg         cyl        disp          hp        drat
#' # "character" "character" "character" "character" "character"
#' #          wt        qsec          vs          am        gear        carb 
#' # "character" "character" "character" "character" "character" "character" 
retype_all<-function(df,from,to){
  for(i in which(sapply(df,class)==from)){
    df[[i]]<-as(df[[i]],to)
  }
  df
}

#' Read and merge all .csv files in a folder
#'
#' @param folder path to a folder
#' @param readfunc list of functions that will be used to read the files; if the first function fails, the second function will be used, etc.
#'
#' @return A data.frame containing all merged .csv files 
#' @export
read.csv.folder<-function(folder="./", readfunc=list(read.csv,read.csv2,read.table)){
  flist<-list.files(folder)
  flist<-flist[grepl(".csv",flist)]
  datlist<-list()
  ct<-0
  for(file in flist){
    ct<-ct+1
    datlist[[ct]]<-"fail"
    for(i in seq_len(length(readfunc))){
      #try(
        datlist[[ct]]<-do.call(readfunc[[i]],list(file=paste0(folder,file),stringsAsFactors=F))
      #,silent=T)
      if(any(datlist[[ct]]!="fail")){ break; }
    }
    if(any(datlist[[ct]]=="fail")){
      warning("Failed to read ",file)
    }
  }
  combodat<-datlist[[1]]
  if(length(datlist)>1){
    for(i in 2:length(datlist)){
      combodat<-rbind(combodat,datlist[[i]])
    }
  }
  return(combodat)
}

#' Install packages if neccesary, then load them.
#' @param ... Unquoted names of packages to try loading, and if unable, install and load.
#'
#' @examples trypackages(stats,utils,compiler)
#' @export
trypackages<-function(...){
  packs<-args2strings(...)
  for(pack in packs){
    if(!require(pack,character.only=T)){
      install.packages(pack)
      require(pack,character.only=T)
    }
  }
}

#' Get all possible combinations of strings
#' @description \code{combobulate()} returns all possible combinations of the provided character strings, each combination merged into a single string.
#' @param ... Character vectors to combobulate.
#'
#' @return A character vector.
#' @export
#'
#' @examples combobulate("Hello ",c("Sir","Madam"),", ",c("may I take your order?","what shall it be?"))
#' # [1] "Hello Sir, may I take your order?"    
#' # [2] "Hello Madam, may I take your order?"   
#' # [3] "Hello Sir, what shall it be?"  
#' # [4] "Hello Madam, what shall it be?"
combobulate<-function(...){
  args<-list(...)
  if(length(args)==0){ return("") }
  output<-args[[1]]
  for(i in min(length(args),2):length(args)){
    currarg<-args[[i]]
    newoutput<-character(0)
    for(j in min(length(currarg),1):length(currarg)){
      newoutput<-c(newoutput,paste0(output,currarg[[j]]))
    }
    output<-newoutput
  }
  return(output)
}

#' Test if two correlation coefficients significantly differ
#' @description Uses Fisher's r to z transformation, then performs a z-test on the resulting z-scores
#' @param cor1,cor2 Correlation values being compared
#' @param n1,n2 Sample sizes of the correlation coefficients
#'
#' @return List containing the z-score and p-value 
#' @export
#' 
#' @references http://vassarstats.net/rdiff.html
compcorr<-function(cor1,cor2,n1,n2){
  r2z<-function(r){  z <- .5 * (log(1+r) - log(1-r)) }
  zval<-abs(r2z(cor1)-r2z(cor2)) / sqrt((1/(n1-3)) + (1/(n2-3)))
  pval<-pnorm(zval)
  cat("Z =",zval,"\np =",pval,"\n")
  return(invisible(list(zscore=zval,pvalue=pval)))
}

# CorrCrunch ####

#' Analyse the robustness of a correlation
#' @description \code{CorrCrunch()} computes the minimum number of cases that need to be removed from a dataset to flip the sign of a correlation coefficient.
#' This can be useful in distinguishing genuine correlations from spurious findings that hinge on one or two outliers.
#' Cases are removed iteratively; in each iteration the case that maximally shrinks the correlation coefficient is removed.
#' @param x,y Numeric vectors to correlate.
#' @param verbose if TRUE, prints verbose output.
#'
#' @return A list containing the number of cases that need to be removed to flip the sign of the correlation coefficient; 
#' the proportion removed cases in the data; and a data.frame without these cases.
#' @export
#'
#' @examples CorrCrunch(mtcars$mpg,mtcars$wt)
#' #Holdout needed to flip the sign: 19 (63.33%)
#' #Final r: 0.01181141
CorrCrunch<-function(x,y,verbose=F){
  iter<-0
  delidx<-0
  iterdf<-data.frame(x=x,y=y)
  iterdf<-iterdf[!is.na(rowSums(iterdf)),]
  origrval<-cor(iterdf$x,iterdf$y)
  rval<-origrval
  corsign<-sign(origrval)
  if(verbose){ cat(sep="","Holdout: ",0,", r: ",rval,"\n") }
  while(sign(rval)==corsign & iter<1000){
    iter<-iter+1
    for(i in 1:nrow(iterdf)){
      loopdf<-iterdf[-i,]
      if(corsign*cor(loopdf$x,loopdf$y)<corsign*rval){
        delidx<-i
        rval<-cor(loopdf$x,loopdf$y)
      }
    }
    if(verbose){ cat(sep="","Holdout: ",iter,", r: ",rval,"\n") }
    iterdf<-iterdf[-delidx,]
  }
  structure(.Data=list(h=iter,h.prop=iter/(length(x)-2),lastdf=iterdf),class="CorrCrunch")
}

print.CorrCrunch<-function(x){
  cat("Holdout needed to flip the sign: ",x$h,
      " (",round(x$h.prop*100,digits=2),"%)\n",sep="")
  cat("Final r: ",cor(x$lastdf$x,x$lastdf$y),"\n",sep="")
}
registerS3method("print", "CorrCrunch", print.CorrCrunch)


# Asymmetric CorMat with CorrCrunch ####
#' Create a Correlation Table
#'
#' @param df A data.frame. 
#' @param rowids,columnids character vectors containing column names from \code{df} that need to be correlated.  
#' @param rowdf,columndf data.frames whose columns need to be correlated. 
#' Either \code{df, rowids, & columnids} or \code{rowdf & columndf} are required.
#'
#' @return A formatted markdown table containing correlation coefficients, p-values, and 
#' the number and percentage of cases that need to be removed to flip the sign of each correlation coefficient.
#' @export
#'
#' @examples CorTable(mtcars,rowids=c("mpg","disp","hp"),columnids=c("drat","wt","qsec"))
#' 
#' CorTable(rowdf=mtcars[,c(1,3,4)],columndf=mtcars[,5:7])
CorTable<-function(df,rowids,columnids,rowdf,columndf){
  if(missing(df) | missing(rowids) | missing(columnids)){
    df<-cbind(rowdf,columndf)
    rowids<-colnames(rowdf)
    columnids<-colnames(columndf)
  }

  cormat<-matrix(NA,nrow=length(rowids),ncol=length(columnids))
  colnames(cormat)<-columnids
  rownames(cormat)<-rowids
  dfmat<-pmat<-hmat<-cormat
  
  for(i in rowids){
    for(j in columnids){
      corobj<-cor.test(df[,i],df[,j])
      cormat[i,j]<-corobj$estimate
      pmat[i,j]<-corobj$p.value
      hmat[i,j]<-CorrCrunch(df[,i],df[,j])$h
      dfmat[i,j]<-corobj$parameter
    }
  }
  
  outmat<-matrix("",ncol=length(columnids),nrow=length(rowids)*5-1)
  colnames(outmat)<-abbreviate(columnids)
  outrows<-rep("",length(rowids)*5-1)
  outrows[5*(0:(length(rowids)-1))+1]<-rowids
  rownames(outmat)<-outrows
  
  outmat[5*(0:(length(rowids)-1))+1,]<-gsub("0\\.","\\.",paste0("r=  ",format(cormat,digits=0,nsmall=2)))
  outmat[5*(0:(length(rowids)-1))+2,]<-gsub("0\\.","\\.",paste0("p=  ",format(pmat,digits=0,nsmall=3)))
  outmat[5*(0:(length(rowids)-1))+3,]<-gsub("0\\.","\\.",paste0("h=    ",format(hmat,digits=0)))
  outmat[5*(0:(length(rowids)-1))+4,]<-gsub("0\\.","\\.",paste0("h/df=",format(hmat/dfmat,digits=0,nsmall=2)))
  knitr::kable(outmat,digits=2,align="r")
}


#stolen from stackoverflow
#' @export
hilight<-function(x,y,s, bg="yellow") {
  text.width <- strwidth(s)
  text.height <- strheight(s)
  rect(x,y,x+text.width,y+text.height,col=bg,border=NA)
  text(x,y,s,adj=c(0,0))
}


#FullAutocorPlot
AutocorPlot<-function(ds,ppvar,rtvar,scope=64){
  if(missing(ds)){
    ds<-data.frame(ppvar=ppvar,rtvar=rtvar,stringsAsFactors = F)
    ppvar<-"ppvar"
    rtvar<-"rtvar"
  }
  ds%<>%as.data.frame()
  pplist<-unique(ds[,ppvar])%>%as.vector
  npp<-length(pplist)
  cormat<-matrix(nrow=scope,ncol=npp)
  ct<-0
  for(pp in pplist){
    ct<-ct+1
    tds<-ds%>%filter((!!sym(ppvar)) ==pp)
    autocor<-numeric()
    for(i in seq_len(scope)){
      autocor[i]<-cor(tds[,rtvar],lag(tds[,rtvar],i),use="complete.obs")
    }
    cormat[,ct]<-autocor
  }
  colnames(cormat)<-pplist
  plot(rowMeans(cormat),type="s",xlab="lag",ylab="autocorrelation",main=paste("Variable",rtvar,"autocorrelation"),
       ylim=c(min(cormat),max(cormat)))
  for(i in seq_len(npp)){ lines(x=cormat[,i],pch=".",col=rgb(0,0,0,0.1),type="l") }
  abline(h=0,col="gray")
  return(invisible(cormat))
}


#TransformPlots
#' Title
#' @description Visualize how different transformations of the data will fit to a normal distribution.
#' @param x A numeric vector.
#'
#' @export
#' 
#' @examples
#' TransformPlots(mtcars$disp)
TransformPlots<-function(x){
  par(mfrow=c(2,2),mar=c(3,2,3,1))
  inv<-function(x){ 1/x }
  nothing<-function(x){ x }
  titles<-c("Untransformed","Log-transformed","Sqrt-transformed","Inverse-transformed")
  transforms<-c(nothing,log,sqrt,inv)
  for(i in 1:4){
    y<-do.call(transforms[[i]],list(x))
    ks<-try(ks.test(y,"pnorm"),silent=T)
    car::qqp(y,"norm",main=paste(titles[i],
                                 "\nKS-test D = ",round(ks$statistic,digits=3),
                                 ", p = ",round(ks$p.value,digits=3)))
  }
}

#LEGACY
DistroPlots<-function(x){
  par(mfrow=c(2,3),mar=c(3,2,3,1))
  car::qqp(x,"norm",main="Untransformed")
  car::qqp(log(x),"norm",main="Log-transformed")
  car::qqp(sqrt(x),"norm",main="Sqrt-transformed")
  car::qqp(1/x,"norm",main="Inverse-transformed")
  MASS::fitdistr(x, "gamma") %$%
    car::qqp(x, "gamma", shape = .$estimate[[1]], scale = .$estimate[[2]],main="Gamma-distribution")
  x%>%density%>%plot(main="Density")
}

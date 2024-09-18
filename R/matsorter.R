

# For viewing the matrix as image
showmat<-function(x){
  image(x[,rev(seq_len(dim(x)[2]))])
}

mat<-matrix(rnorm(144),ncol=12,nrow=12)
mat<-mat+t(mat)

showmat(mat)
newmat<-sortmat(mat,"bottomright")
showmat(newmat)

sortmat<-function(mat,sorttype=c("middle","bottomright")){
  dims<-dim(mat)
  k<-dims[1]
  
  if(sorttype=="middle"){
    wt <- (mean(dims)/2-abs(row(mat)-col(mat)))
  }else if(sorttype=="bottomright"){
    wt <- row(mat)+col(mat)
  }

  allswaps<-expand.grid(row=seq_len(k),col=seq_len(k))
  tryswaps<-seq_len(nrow(allswaps)) |> sample()
  oldscore<-sum(mat*wt)
  for(i in 1:nrow(allswaps)){
    key<-seq_len(k)
    key[allswaps[tryswaps[i],1]]<-allswaps[tryswaps[i],2]
    key[allswaps[tryswaps[i],2]]<-allswaps[tryswaps[i],1]
    propmat<-mat[key,key]
    newscore<-sum(propmat*wt)
    if(newscore>oldscore){
      #message("Gain: i=",i,", +",newscore-oldscore)
      oldscore<-newscore
      mat<-propmat
    }
  }

  return(mat)
}


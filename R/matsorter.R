

#' Sort a square matrix
#' This uses an iterative algorithm that swaps the rows and columns of a matrix 
#' to ensure the highest values are either in the middle or the bottom-right.
#'
#' @param mat Square numeric matrix to be sorted
#' @param sorttype Where the highest values should appear - 
#' "diag", "center", or "bottomright".
#'
#' @return A sorted square matrix
#' @export
#'
#' @examples
#' 
#' mat<-matrix(rnorm(144),ncol=12,nrow=12)
#' mat<-mat+t(mat)
#' 
#' newmat1<-sortmat(mat,"bottomright")
#' 
#' newmat2<-sortmat(mat,"diag")
#' 
#' newmat3<-sortmat(mat,"center")
#' 
sortmat<-function(mat,sorttype=c("diag","center","bottomright")){
  dims<-dim(mat)
  k<-dims[1]
  
  if(sorttype=="diag"){
    wt <- (mean(dims)/2-abs(row(mat)-col(mat)))
  }else if(sorttype=="center"){
    wt <- prod(dims/2)-abs((row(mat)-(nrow(mat)+1)/2)*(col(mat)-(ncol(mat)+1)/2))
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

# For viewing the matrix as image
showmat<-function(x){
  image(x[,rev(seq_len(dim(x)[2]))])
}




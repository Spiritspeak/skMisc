

#' Generate unique pairs or N-tuplets
#'
#' @param nval Number of values to arrange into unique tuplets.
#' @param ntuplet N-tuplets to arrange the values uniquely into.
#' @param incl.self Determines whether a value can be paired with itself.
#'
#' @return A matrix where each row is a unique N-tuplet.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' allpairs(nval=5,ntuplet=3)
#' 
allpairs <- function(nval, ntuplet=2, incl.self=FALSE){
  currmat <- matrix(seq_len(nval), ncol=1)
  offset <- ifelse(incl.self, 0, 1)
  for(tuple in seq_len(ntuplet)[-1]){
    newmats <- list()
    for(i in seq_len(nrow(currmat))){
      startval <- currmat[i,tuple-1] + offset
      if(startval <= nval){
        itervec <- startval:nval
        newmats[[i]] <- do.call(cbind, c(as.list(currmat[i,]), list(itervec)))
      }
    }
    currmat <- do.call(rbind, newmats)
  }
  return(currmat)
}


#' Convert between a matrix and a long-format data.frame
#'
#' @param x In case of \code{unwrap.matrix()}, a matrix to unwrap; 
#' in case of \code{rewrap.matrix}, a data.frame with three columns,
#' respectively representing the row name, column name, and value.
#' @param na.value Which value to use in the matrix for elements 
#' not provided in \code{x}
#'
#' @return \code{unwrap.matrix()} returns a \code{data.frame} with three columns: 
#' \code{row} and \code{col} indicating the row and column names, and 
#' \code{value} indicating the respective value in the matrix. 
#' If no row or column names are available, 
#' the row or column number is used instead.
#' 
#' \code{rewrap.matrix()} returns a matrix.
#' 
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' mymatrix <- matrix(1:40,ncol=8,nrow=5)
#' unwrap.matrix(mymatrix)
#' rewrap.matrix(unwrap.matrix(mymatrix))
#' 
#' carmatrix <- as.matrix(mtcars)
#' unwrap.matrix(carmatrix)
#' rewrap.matrix(unwrap.matrix(carmatrix))
#' 
unwrap.matrix <- function(x){
  dn <- dimnames(x)
  unwrap <- expand.grid(row=if(is.null(dn[[1]])){ seq_len(nrow(x)) }else{ dn[[1]] },
                        col=if(is.null(dn[[2]])){ seq_len(ncol(x)) }else{ dn[[2]] })
  unwrap[["value"]] <- as.vector(x)
  return(unwrap)
}

#' @rdname unwrap.matrix
#' @export
#' 
rewrap.matrix <- function(x, na.value=NA){
  rn <- unique(x[,1,drop=T])
  cn <- unique(x[,2,drop=T])
  out <- matrix(NA,
                nrow=length(rn),
                ncol=length(cn),
                dimnames=list(rn, cn))
  out[as.matrix(x[,c(1,2)])] <- x[,3,drop=T]
  out[is.na(out)] <- na.value
  return(out)
}

#' Sort a square matrix
#' This uses an iterative algorithm that swaps the rows and columns of a matrix 
#' to ensure the highest values are either in the middle or the bottom-right.
#'
#' @param mat Square numeric matrix to be sorted.
#' @param sorttype Where the highest values should appear - 
#' "diag", "center", or "bottomright".
#'
#' @return A sorted square matrix.
#' @author Sercan Kahveci
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
sortmat <- function(mat, sorttype=c("diag", "center", "bottomright")){
  dims <- dim(mat)
  k <- dims[1]
  
  if(sorttype == "diag"){
    wt <- (mean(dims) / 2 - abs(row(mat) - col(mat)))
  }else if(sorttype == "center"){
    wt <- prod(dims / 2) - abs((row(mat) - (nrow(mat) + 1) / 2) * (col(mat) - (ncol(mat) + 1) / 2))
  }else if(sorttype == "bottomright"){
    wt <- row(mat) + col(mat)
  }
  
  allswaps <- expand.grid(row=seq_len(k), col=seq_len(k))
  tryswaps <- seq_len(nrow(allswaps)) |> sample()
  oldscore <- sum(mat * wt)
  for(i in seq_len(nrow(allswaps))){
    key <- seq_len(k)
    key[ allswaps[tryswaps[i], 1] ] <- allswaps[tryswaps[i], 2]
    key[ allswaps[tryswaps[i], 2] ] <- allswaps[tryswaps[i], 1]
    propmat <- mat[key,key]
    newscore <- sum(propmat * wt)
    if(newscore > oldscore){
      oldscore <- newscore
      mat <- propmat
    }
  }
  
  return(mat)
}

#' Set column and row names of an object
#' 
#' These are convenience functions that return an object with its column or row names changed.
#' Use it in pipes.
#' 
#' @param x an object.
#' @param names column or row names to be assigned to the object.
#' 
#' @export
#' @author Sercan Kahveci
#' @examples 
#' setColNames(ToothGrowth,c("length","supplement","dosage"))
#' setRowNames(BOD,BOD$Time)
#' 
setColNames <- function(x, names){ colnames(x) <- names; return(x) }
#' @export
#' @rdname setColNames
setRowNames <- function(x, names){ rownames(x) <- names; return(x) }


#' Plot matrix as heatmap
#'
#' @param x A matrix.
#' @param text Whether to print the values of the matrix as text (defaults to FALSE).
#' @param plot Whether to plot the resulting ggplot object (defaults to TRUE).
#' @param ... Ignored.
#'
#' @return Invisibly returns the ggplot object for further modification.
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' # Create a matrix to plot
#' mymat <- matrix(rnorm(100), nrow=10)
#' mymat <- mymat + t(mymat)
#' mymat[cbind(sample(c(1:10)), sample(c(1:10)))] <- NA
#' colnames(mymat) <- rownames(mymat) <- sample(letters[1:10])
#' 
#' plotmat(mymat)
#' 
plotmat <- function(x, text=FALSE, plot=TRUE, ...){
  out <- x |> unwrap.matrix() |> ggplot() + 
    aes(y=.data[["row"]],x=.data[["col"]],fill=.data[["value"]]) + 
    geom_tile() + 
    scale_fill_gradient2(na.value="grey25") + 
    coord_cartesian(expand=0) + 
    theme_bw() + 
    theme(axis.ticks=element_blank(),
          panel.grid=element_blank(),
          panel.border=element_blank(),
          axis.title=element_blank(),
          legend.title=element_blank())
  if(text){
    out <- out + 
      geom_text(aes(label=dropLeadingZero(round(.data[["value"]], digits=2))),
                size=min(2, 2*10/nrow(x)))
  }
  if(plot){
    plot(out)
  }
  return(invisible(out))
}


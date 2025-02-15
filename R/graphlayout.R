
# library(igraph)
# library(graphlayouts)
# mylist<-list(matrix(rnorm(9),nrow=3,
#                     dimnames=list(letters[1:3],letters[1:3])),
#              matrix(rnorm(4),nrow=2,
#                     dimnames=list(letters[2:3],letters[2:3])))
# alignLayouts(mylist,FillInNodes=T)

alignLayouts <- function(matlist, alpha=.5, FillInNodes=FALSE){
  
  # Get a full list of node names
  allnames <- lapply(matlist, colnames) |> unlist() |> unique()
  
  # Fill in missing nodes in each matrix
  filledmats <- list()
  for(i in seq_along(matlist)){
    itermat <- matrix(0, nrow=length(allnames), ncol=length(allnames), 
                      dimnames=list(allnames, allnames))
    itermat[rownames(matlist[[i]]), colnames(matlist[[i]])] <- matlist[[i]]
    filledmats[[i]] <- itermat
  }
  
  # Exponentiate weights, convert to igraph, and apply layout_as_dynamic
  layouts <- filledmats |> 
    lapply(function(x){x[x!=0] <- exp(-x[x!=0]); x}) |>
    lapply(igraph::graph_from_adjacency_matrix) |> 
    layout_as_dynamic(weights=NULL, alpha=alpha) |> 
    lapply(function(x){rownames(x) <- allnames; x})
  
  # Remove filled-in nodes if requested
  if(!FillInNodes){
    for(i in seq_along(filledmats)){
      currnodes <- colnames(matlist[[i]])
      filledmats[[i]] <- filledmats[[i]][currnodes,currnodes]
      layouts[[i]] <- layouts[[i]][currnodes,]
    }
  }
  
  # Return output
  output <- list(mats=filledmats, layouts=layouts)
  return(output)
}


altlayout <- function(x, type=c("stress", "kk", "fr", "drl", "dh", "focus"),
                      repulsion=1, negrepulsion=1, ...){
  
  # Parse arguments
  args <- list(...)
  type <- match.arg(type)
  
  # Convert to igraph
  qgr <- igraph::graph_from_adjacency_matrix(x)
  oldwts <- igraph::E(qgr)$weight
  
  # Compute exponentiated weights for the algorithms that 
  # cannot handle negative values
  if(type!="focus"){
    newwts <- exp(-oldwts*repulsion*ifelse(oldwts<0,negrepulsion,1))
  }
  
  # Run layout functions
  if(type=="kk"){ 
    igraph::E(qgr)$weight <- newwts
    out<-igraph::layout_with_kk(qgr) 
  }else if(type=="stress"){
    out <- graphlayouts::layout_with_stress(qgr,newwts) 
  }else if(type=="fr"){
    igraph::E(qgr)$weight <- 1/newwts
    out <- igraph::layout_with_fr(qgr)
  }else if(type=="drl"){
    igraph::E(qgr)$weight <- 1/newwts
    out <- igraph::layout_with_drl(qgr)
  }else if(type=="dh"){
    igraph::E(qgr)$weight <- 1/newwts
    out <- igraph::layout_with_dh(qgr)
  }else if(type=="focus"){
    qgr %<>% igraph::delete_edges(which(oldwts<0))
    out <- graphlayouts::layout_with_focus(qgr, v=args$focus)$xy
  }
  
  # Return output
  return(out)
}


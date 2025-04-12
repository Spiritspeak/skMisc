

#' Extended Bayesian Information Criterion for \code{glmnet} objects
#'
#' @param x The \code{glmnet} model.
#' @param gamma The gamma parameter. This should be a value between 0 and 1.
#' At a gamma of 0, EBIC is equivalent to BIC, while at higher values, it is
#' more strict on the number of predictors modelled.
#'
#' @returns The EBIC of the model.
#' @export
#'
#' @examples
#' 
#' 
EBIC.glmnet <- function(x, gamma){
  tLL <- -deviance(x)
  k <- x$df
  n <- x$nobs
  p <- nrow(x$beta)
  EBIC <- -tLL + log(n)*k + 2*gamma*k*log(p)
  return(EBIC)
}

#################################
# Predictor generator functions #
#################################

onsets2idmat <- function(pat, preds){
  idmat <- vector(mode="list", length=length(preds))
  names(idmat) <- preds
  for(s in preds){
    idmat[[s]] <- Matrix(as.numeric(pat[["state"]]==s), ncol=1)
  }
  idmat <- do.call(cbind, idmat)
  colnames(idmat) <- preds
  return(idmat)
}

stairmat2predmat <- function(x, stairmat, 
                             type=c("step", "flat", "power", 
                                    "inverse", "accrual")){
  type <- match.arg(type)
  nonzero_stairmat <- stairmat!=0
  
  # Predictor level is x to the power of N new states entered since the state
  if(type == "power"){
    out <- Matrix(data=0,
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)))
    out[nonzero_stairmat] <- x^(stairmat[nonzero_stairmat]-1)
  }
  
  # Predictor level is 1 divided by the number of states entered since that state
  if(type == "inverse"){
    out <- Matrix(data=0,
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)))
    out[nonzero_stairmat] <- (x+1)/(x+stairmat[nonzero_stairmat])
  }
  
  # Predictor level is set to 1 for a set number of rows after a state is entered
  if(type=="flat"){
    out <- Matrix(data=as.numeric(nonzero_stairmat),
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)),
                  sparse=T)
    out[stairmat > x] <- 0
  }
  
  # Predictor level goes down in discrete steps until 0 after a state is entered
  if(type=="step"){
    nonzero_stairmat[stairmat>x] <- F
    out <- Matrix(data=0,
                  nrow=nrow(stairmat),
                  ncol=ncol(stairmat),
                  dimnames=list(c(), colnames(stairmat)))
    out[nonzero_stairmat] <- (x+1-stairmat[nonzero_stairmat])/x
  }
  
  # Once you enter a state, the predictor is set to 1 indefinitely
  if(type=="accrual"){
    out <- nonzero_stairmat+0
  }
  
  return(out)
}

get_recencymat <- function(data, preds, seq.onset=TRUE, verbose=FALSE){
  quickave <- function(X, INDEX, FUN){
    split(X, INDEX) <- lapply(split(X, INDEX), FUN)
    X
  }
  
  data <- dplyr::arrange(data, sequence, date)
  onsets <- !duplicated(paste(data$sequence, data$state, sep="+"))
  pat <- data[onsets, c("sequence", "state", "date")]
  out <- matrix(0, ncol=length(preds), nrow=nrow(pat), dimnames=list(NULL, preds))
  
  for(currstate in preds){
    if(verbose){ cat("\r", currstate, "                        ") }
    indices <- 1:nrow(data)
    indices[data$state != currstate] <- NA
    indices <- quickave(X=indices, INDEX=data$sequence, FUN=carryforward_numeric)
    lastreldate <- data$date[indices]
    out[,currstate] <- data$date[onsets] - lastreldate[onsets]
  }
  rm(indices, lastreldate, onsets)
  out[is.na(out)] <- 0
  
  if(seq.onset){
    if(verbose){ cat("\r(Start)                        ") }
    indices <- 1:nrow(pat)
    sequenceonset <- lag(pat$sequence) != pat$sequence
    sequenceonset[1] <- T
    indices[!sequenceonset] <- NA
    indices <- carryforward_numeric(indices)
    out <- cbind(`(Start)`=pat$date - pat$date[indices] + 1, out)
    rm(indices, sequenceonset)
  }
  
  out <- Matrix(out, dimnames=list(NULL, colnames(out)))
  return(list(pat=pat, stairmat=out))
}

onsets2stairmat <- function(pat,
                            preds,
                            by=c("index", "date"),
                            direction=1,
                            seq.onset=TRUE,
                            verbose=FALSE){
  seqchains <- split(pat[["state"]], pat[["sequence"]])
  by <- match.arg(by)
  if(by == "index"){
    timechains <- split(pat[["idx"]], pat[["sequence"]])
  }else{
    timechains <- split(pat[["date"]], pat[["sequence"]])
  }
  
  allpreds <- preds
  if(end){ allpreds <- c("(end)", allpreds) }
  if(start){ allpreds <- c("(start)", allpreds) }
  
  stairs <- vector(length(allpreds), mode="list")
  names(stairs) <- allpreds
  for(s in preds){
    if(verbose){ cat("\r", s, "                        ") }
    currvec <- direction *
      unsplit(EnumerateFrom(sequences=seqchains,
                            times=timechains,
                            target=s),
              pat[["sequence"]])
    currvec[currvec < 0] <- 0
    stairs[[s]] <- currvec |> Matrix(ncol=1)
  }
  
  if(seq.onset & direction == 1){
    if(verbose){ cat("\r", "(Start)                        ") }
    stairs[["(Start)"]] <- 
      timechains |> 
      lapply(function(x){ x-x[1]+1 }) |> 
      unsplit(pat[["sequence"]])
  }  
  if(seq.onset & direction == -1){
    if(verbose){ cat("\r", "(End)                        ") }
    stairs[["(End)"]] <- 
      timechains |> 
      lapply(function(x){ x[length(x)]-x+1 }) |> 
      unsplit(pat[["sequence"]])
  }
  if(verbose){ cat("\r") }
  stairs <- do.call(cbind, stairs)
  colnames(stairs) <- allpreds
  return(stairs)
}

###################################
# Helper functions for regressors #
###################################

getStateMask <- function(currstate, pat, sampletype, direction){
  pat$idx <- pat$idx * direction
  if(sampletype == "all"){
    include_nonvisitors <- T
  }else if(sampletype == "visitors"){
    include_nonvisitors <- F
  }
  mask <- pat |> 
    group_by(sequence) |> 
    mutate(limited=any(state == currstate)) |>
    transmute(mask= ifelse(limited, 
                           idx <= idx[state == currstate], 
                           include_nonvisitors))
  return(mask$mask)
}

# Assemble betas and ebics into matrices
assembleCoefficients <- function(statefits, extrapreds=NULL){
  # Extract ebics
  ebicmat <- sapply(statefits, \(x)x$ebics)
  
  # Add zeroes for each fit's own state and for start
  statenames <- names(statefits)
  betalist <- statefits |> lapply(\(x)x[["betas"]])
  for(statename in statenames){
    # Add empty row for own state
    rowid <- nrow(betalist[[statename]])+1
    betalist[[statename]] %<>% rbind(0)
    rownames(betalist[[statename]])[rowid] <- statename
    
    # Reorder
    betalist[[statename]] <- 
      betalist[[statename]][c("(Intercept)", extrapreds, statenames),]
  }
  
  # Extract betas by gamma
  gammanames <- rownames(ebicmat)
  gamma_coeflist <-
    vector(length=length(gammanames), mode="list") |> 
    setNames(gammanames)
  for(currgammaname in gammanames){
    currbetas <- betalist |> sapply(\(x)x[, currgammaname])
    colnames(currbetas) <- names(betalist)
    gamma_coeflist[[currgammaname]] <- currbetas
  }
  
  # Return output
  out <- list(gammacoefs=gamma_coeflist, ebicmat=ebicmat)
  return(out)
}

extract_coefs <- function(fit, gammas){
  iterminebics <- iterebics <- numeric(length(gammas))
  for(gammidx in seq_along(gammas)){
    ebicvec <- EBIC.glmnet(fit, gamma=gammas[gammidx])
    iterminebics[gammidx] <- which.min(ebicvec)
    iterebics[gammidx] <- min(ebicvec)
  }
  betamat <- coef(fit)[,iterminebics]
  colnames(betamat) <- names(iterebics) <- paste0("gamma", gammas)
  
  out <- list(betas=betamat, ebics=iterebics)
  return(out)
}

################################
# TRANSITION NETWORK FUNCTIONS #
################################

#' Model a transition network from sequence data
#'
#' @param data A data.frame with 3 columns: sequence (character), 
#' state (character), and date (numeric).
#' @param predictors Which states should be included in the network?
#' @param decaytype Once a state has occurred, how should its predictor look?
#' 
#' * If "flat", the predictor is 1 after the state has occurred until it 
#' no longer satisfies the condition set by the parameter argument.
#' 
#' * If "step", the predictor is 1 immediately after the state has occurred,
#' and decays to 0 until it no longer satisfies the condition set by 
#' the parameter argument.
#' 
#' * If "accrual", the predictor is 1 after the state has occurred 
#' until sequence end.
#' @param parameter This determines the largest temporal distance from 
#' the occurrence of a state that is transformed into a nonzero value 
#' as a predictor. The exact behavior of this argument is determined by
#' the \code{decaytype} argument. See Details for clarification.
#' @param predtype 
#' @param direction In which direction should predictions be made?
#' 1 means current states predict future states, -1 means current states predict
#' past states.
#' @param sampletype When the occurrence of a state is predicted, should every 
#' sequence be included ("all") or only those that contain the state ("visitors")?
#' @param gammas Which EBIC gamma values should be used for model selection?
#' Higher values mean less edges are retained in the network.
#' @param alpha This is the elasticnet mixing parameter.
#' 1 (default) means lasso regularization, 0 means ridge regularization, and 
#' values in-between mix the two. See [glmnet::glmnet()].
#' @param force.positive Should model beta values be restricted to 
#' the positive range?
#' @param ncores Number of models to run in parallel. 
#' Be careful not to overload system memory.
#' @param verbose Produce verbose output in console?
#' 
#' @details
#' # Rationale
#' 
#' This is a method for analyzing sequences of states. The first occurrence 
#' of each state is predicted with the preceding first or most recent occurrences 
#' of other states using lasso-regularized logistic regressions. 
#' You set with the \code{preds} argument which states are predicted and 
#' used as predictors. All other states are kept in the data, 
#' but do not contribute to prediction nor are predicted.
#' 
#' The level of analysis (and the actual data analyzed) is the _first_ 
#' occurrences of all states. That is, each row represents 
#' the first occurrence of a state.
#' For each model, the dependent variable is coded 
#' as 0 if the first occurrence of a state on the current row
#' is not the state being predicted, and 1 if it is that state; 
#' prediction within a single sequence continues until 
#' this state or the end of the sequence is reached. 
#' Hence, if a state is reached before the end of the sequence, 
#' all states after it are ignored in the prediction of that state, 
#' since there can only be one _first_ occurrence of a state within a sequence,
#' and therefore prediction of a first occurrence would not make sense
#' after the first occurrence hadd already occurred.
#' 
#' 
#' @md
#' @returns
#' @author Sercan Kahveci
#' @export
#'
#' @examples
#' 
#' 
transitionNet <- function(data,
                                  predictors,
                                  decaytype=c("step", "flat", "accrual"),
                                  parameter,
                                  predtype=c("onset-to-onset",
                                             "recent-onset-to-onset",
                                             "recency-to-onset"),
                                  direction=1,
                                  sampletype="all",
                                  gammas=seq(0, 1, .1),
                                  alpha=1,
                                  force.positive=FALSE,
                                  ncores=NULL,
                                  verbose=TRUE){
  predtype <- match.arg(predtype)
  mattype <- match.arg(mattype)
  
  if(predtype=="recency-to-onset"){
    if(verbose){ message("Generating predictor matrix precursor") }
    stopifnot(direction==1)
    recencydata <- get_recencymat(data=data,
                                  preds=predictors,
                                  seq.onset=TRUE,
                                  verbose=verbose)
    pat <- recencydata$pat
    stairmat <- recencydata$stairmat
    rm(recencydata)
    gc()
  }
  
  if(predtype %in% c("onset-to-onset","recency-onset-to-onset")){
    stopifnot(!is.null(predictors))
    if(verbose){ message("Generating predictor matrix precursor") }
    pat <- data[!duplicated(paste(data$sequence, data$state, sep="+")), 
                c("sequence", "state", "date")]
    stairmat <- onsets2stairmat(pat=pat,
                                preds=predictors,
                                by=switch(type,
                                          `onset-to-onset`="idx",
                                          `recency-onset-to-onset`="date"),
                                direction=direction,
                                seq.onset=TRUE,
                                verbose=verbose)
  }
  
  # Form predictor matrix
  if(verbose){ message("Generating predictor matrix") }
  predmat <- stairmat2predmat(x=parameter, stairmat=stairmat, type=mattype)
  rm(stairmat, data)
  gc()
  
  # Generate id matrix
  stopifnot(!is.null(predictors))
  if(verbose){ message("Generating DV matrix") }
  idmat <- onsets2idmat(pat=pat, preds=predictors)
  
  # add idx column to pat
  pat$idx <- ave(seq_along(pat$sequence),pat$sequence,FUN=function(x){seq_along(x)})
  
  if(force.positive){
    minweights <- ifelse(colnames(predmat) %in% c("(Start)","(End)","(Intercept)"),
                         -Inf, 0)
  }else{
    minweights <- rep(-Inf, ncol(predmat))
  }
  
  # Prepare cluster
  if(verbose){ message("Preparing cluster") }
  if(is.null(ncores)){ ncores <- detectCores() }
  clust <- makeCluster(ncores)
  registerDoParallel(clust, cores=ncores)
  gc(reset=T, full=T)
  on.exit({
    stopCluster(clust)
    registerDoSEQ()
  })
  
  # Run regressions
  if(verbose){ message("Running regressions on ",ncores," cores") }
  fits <- 
    foreach(currstate=colnames(idmat),
            .packages=c("glmnet","dplyr","Matrix"),
            .export=c("getStateMask","extract_coefs","EBIC.glmnet"),
            .multicombine=T) %dopar% {
              gc(reset=T, full=T)
              mask<-getStateMask(currstate=currstate, pat=pat,
                                 sampletype=sampletype, direction=direction)
              
              iterfit<-glmnet(x=predmat[mask, colnames(predmat) != currstate],
                              y=idmat[mask, currstate, drop=F],
                              family="binomial",
                              alpha=alpha,
                              lower.limits=minweights[colnames(predmat) != currstate])
              extract_coefs(iterfit, gammas)
            }
  if(verbose){ message("Assembling results") }
  names(fits) <- colnames(idmat)
  extrapreds <- ifelse(direction==1, "(Start)", "(End)")
  out <- assembleCoefficients(fits, extrapreds=extrapreds)
  return(out)
}

##########################
# Co-occurrence networks #
##########################

# This is a more efficient clone of IsingFit
cooccurrenceNet <-
  function(coocmat,
           gammas=seq(0, 1, .1),
           ncores=NULL,
           dfmax=Inf){
  coocmat <- as.matrix(coocmat)
  stopifnot(is.numeric(coocmat) & all(unique(as.vector(coocmat)) %in% 0:1))
    
  # Prepare cluster
  if(is.null(ncores)){  ncores <- detectCores() }
  clust <- makeCluster(ncores)
  registerDoParallel(clust, cores=ncores)
  on.exit({
    stopCluster(clust)
    registerDoSEQ()
  })
  
  # Run regressions
  fits <- foreach(currsub=colnames(coocmat),
                  .packages=c("glmnet", "dplyr"),
                  .export=c("extract_coefs", "EBIC.glmnet"),
                  .multicombine=T) %dopar% {
                     
    iterfit<-glmnet(x=coocmat[,colnames(coocmat) != currsub, drop=F],
                    y=coocmat[,currsub, drop=F],
                    family="binomial",
                    alpha=1,
                    dfmax=dfmax)
    extract_coefs(iterfit, gammas)
  }
  names(fits) <- colnames(coocmat)
  out <- assembleCoefficients(fits)
  return(out)
}



EBIC.glmnet <- function(md,gamma){
  tLL <- -deviance(md)
  k <- md$df
  n <- md$nobs
  p <- nrow(md$beta)
  EBIC <- -tLL + log(n)*k + 2*gamma*k*log(p)
  return(EBIC)
}

extract_coefs <- function(fit,gammas){
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

stairmat2predmat <- function(x, stairmat, type=c("step", "flat", "power",
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

get_recencymat <- function(data, maxdate, preds, 
                           type=c("flat", "step"), 
                           seq.onset=TRUE, verbose=FALSE){
  type <- match.arg(type)
  
  quickave <- function(X, INDEX, FUN){
    split(X, INDEX) <- lapply(split(X, INDEX), FUN)
    X
  }
  
  data %<>% arrange(sequence, date)
  onsets <- !duplicated(paste(data$sequence, data$state, sep="+"))
  pat <- data[onsets, c("sequence", "state", "date")]
  out <- matrix(0, ncol=length(preds), nrow=nrow(pat), dimnames=list(NULL, preds))
  
  for(currstate in preds){
    if(verbose){ cat("\r", currstate, "                        ") }
    indices <- 1:nrow(data)
    indices[data$state != currstate] <- NA
    indices <- quickave(X=indices, INDEX=data$sequence, FUN=carryforward)
    lastreldate <- data$date[indices]
    out[,currstate] <- data$date[onsets] - lastreldate[onsets]
  }
  rm(indices, lastreldate, onsets)
  out[is.na(out)] <- 0
  
  if(seq.onset){
    if(verbose){ cat("\r(Start)                        ") }
    indices <- 1:nrow(pat)
    authoronset <- lag(pat$sequence) != pat$sequence
    authoronset[1] <- T
    indices[!authoronset] <- NA
    indices <- carryforward(indices)
    out <- cbind(`(Start)`=pat$date - pat$date[indices] + 1, out)
    rm(indices, authoronset)
  }
  
  if(type=="flat"){
    out <- 0 + (out > 0 & out < maxdate)
  }else if(type == "step"){
    out <- maxdate - out
    out[out < 0 | out == maxdate] <- 0
  }
  out <- Matrix(out, dimnames=list(NULL, colnames(out)))
  return(list(pat=pat, predmat=out))
}


################################
# Predictor generators in Rcpp #
################################

cppFunction(
  "
List EnumerateFrom(List sequences, List times, String target){
  int ll = sequences.length();
  List out (ll);
  for(int i = 0; i < ll; ++i){
    CharacterVector currauth = sequences[i];
    int cl = currauth.length();
    IntegerVector currtimes = times[i];
    
    bool found = false;
    int flippoint = cl+1;
    for(int j = 0; j < cl; ++j){
      if(currauth[j] == target){
        flippoint = currtimes[j];
        found = true;
        break;
      }
    }
    
    NumericVector newtimes (cl);
    if(found){
      newtimes = currtimes - flippoint;
    }else{
      newtimes = rep(0, cl);
    }
    
    out[i] = newtimes;
  }
  return out;
}
")

onsets2stairmat <- function(pat,
                            preds,
                            by=c("index", "date"),
                            direction=1,
                            start=TRUE,
                            end=TRUE,
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
  
  if(start & direction == 1){
    if(verbose){ cat("\r", "(start)                        ") }
    stairs[["(start)"]] <- 
      timechains |> 
      lapply(function(x){ x-x[1]+1 }) |> 
      unsplit(pat[["sequence"]])
  }  
  if(end & direction == -1){
    if(verbose){ cat("\r", "(end)                        ") }
    stairs[["(end)"]] <- 
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
    group_by(author) |> 
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
  
  # Add zeroes for own subreddit and for start
  statenames <- names(statefits)
  betalist <- statefits |> lapply(\(x)x[["betas"]])
  for(statename in statenames){
    # Add empty row for own state
    rowid <- nrow(betalist[[statename]])+1
    betalist[[statename]] %<>% rbind(0)
    rownames(betalist[[statename]])[rowid] <- statename
    
    # Reorder
    betalist[[statename]] <- 
      betalist[[statename]][c("(Intercept)", extrapreds, subnames),]
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

################################
# TRANSITION NETWORK FUNCTIONS #
################################

nodewiseTransitionNet <- function(pat,
                                  mattype,
                                  parameter,
                                  by="index",
                                  direction=1,
                                  sampletype="all",
                                  gammas=seq(0, 1, .1),
                                  alpha=1,
                                  force.positive=FALSE,
                                  idmat=NULL,
                                  stairmat=NULL,
                                  predictors=NULL,
                                  ncores=NULL,
                                  verbose=TRUE){
  
  # Generate id/stair matrices if missing
  if(is.null(idmat)){
    stopifnot(!is.null(predictors))
    if(verbose){ message("Generating DV matrix") }
    idmat <- onsets2idmat(pat=pat, preds=predictors)
  }
  if(is.null(stairmat)){
    stopifnot(!is.null(predictors))
    if(verbose){ message("Generating predictor matrix precursor") }
    stairmat <- onsets2stairmat(pat=pat,
                                preds=predictors,
                                by=by,
                                direction=direction,
                                start=direction==1,
                                end=direction==-1,
                                verbose=verbose)
  }
  
  # Form predictor matrix
  if(verbose){ message("Generating predictor matrix") }
  predmat <- stairmat2predmat(x=parameter, stairmat=stairmat, mattype=mattype)
  rm(stairmat)
  
  if(force.positive){
    minweights <- ifelse(colnames(predmat) %in% c("(start)","(end)","(Intercept)"),
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
                              y=idmat[mask, currstate,drop=F],
                              family="binomial",
                              alpha=alpha,
                              lower.limits=minweights[colnames(predmat) != currstate])
              extract_coefs(iterfit, gammas)
            }
  if(verbose){ message("Assembling results") }
  names(subfits) <- colnames(idmat)
  extrapreds <- ifelse(direction==1, "(start)", "(end)")
  out <- assembleCoefficients(fits, extrapreds=extrapreds)
  return(out)
}

##########################
# Co-occurrence networks #
##########################

# This is a more efficient clone of IsingFit
nodewiseCooccurrenceNet<-
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
  fits<-foreach(currsub=colnames(coocmat),
                .packages=c("glmnet","dplyr"),
                .export=c("extract_coefs","EBIC.glmnet"),
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

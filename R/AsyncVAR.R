

# TODO: add support for alpha-based predictor exclusion
detrendVAR <- function(data,vars,idvar,dayvar=NULL,beepvar=NULL,verbose=FALSE){
  if(is.null(dayvar) && is.null(beepvar)){ 
    stop("At least one of dayvar and beepvar must be provided")
  }
  
  vars <- unlist(vars)
  for(cpred in vars){
    if(verbose){ message("Detrending ",cpred) }
    for(pp in unique(data[[idvar]])){
      key<-data[[idvar]]==pp
      form<-as.formula(paste0(cpred," ~ ",paste0(c(dayvar,beepvar),collapse="+")))
      lmmod<-lm(form, data=data[key,])
      data[key,cpred]<-mean(data[key,cpred])+resid(lmmod)
    }
  }
  return(data)
}

checkAsyncVARComponents <- function(modlist,idvar){
  modlist %>% sapply(function(x) x@optinfo$conv$opt)
  modlist %>% sapply(function(x) x@optinfo$conv$lme4$messages)
  modlist %>% sapply(car::vif)
  modlist %>% sapply(function(x)apply(ranef(x)[[idvar]],2,sd)) %>% 
    {.<0.001} %>% rowMeans()
}

extractAsyncVARComponents <- function(modlist,allpreds,alpha){
  Beta_fixed <- do.call(rbind,lapply(modlist,function(x)fixef(x)[allpreds]))
  Beta_SE <- do.call(rbind,lapply(modlist,function(x)arm::se.fixef(x)[allpreds]))
  Beta_P <- 2*(1-pnorm(abs(Beta_fixed/Beta_SE)))
  rownames(Beta_fixed) <- rownames(Beta_SE) <- rownames(Beta_P) <- allpreds
  
  # Nullify non-sig
  betasig<-Beta_fixed
  betasig[Beta_P>alpha]<-0
  
  mat<-t(betasig)
  return(mat)
}

# TODO: reintroduce support for beepvar
# TODO: remove hardcoded reliance on alpha
# TODO: make contemp and between-subject skippable
# TODO: add t-based significance testing
asyncVAR<-function(data,vars,idvar,dayvar,beepvar,covar=NULL,
                   temporal = c("default", "correlated", 
                                "orthogonal", "fixed"), 
                   contemporaneous = c("default", "correlated",
                                       "orthogonal", "fixed", "skip"), 
                   between = c("default", "skip"),
                   nCores = 1, scale = TRUE, alpha=0.05){
  
  allpreds <- unlist(vars)
  predrank <- sapply(allpreds,function(x){ 
    which(sapply(seq_along(vars),function(y){any(x == vars[[y]])}))
  })
    
  cl<-makeCluster(spec=nCores)
  registerDoParallel(cl=cl,cores=nCores)
  
  
  temporal.models<-#list()
    #for(cpred in allpreds){
    foreach(cpred=allpreds,.packages=c("lmerTest")) %dopar% {

      currdata<-data[,c(idvar,dayvar,beepvar,allpreds,covar)]
      currdata$dv_lagged<-NA
      
      # Lag the DV and all same-row observations that theoretically occur at a preceding moment
      moveback <- c(cpred,allpreds[predrank < predrank[allpreds==cpred]])
      movebacknames <- c("dv_lagged",allpreds[predrank < predrank[allpreds==cpred]])
      for(u in unique(currdata[[idvar]])){
        key<-currdata[[idvar]]==u
        currdata[key,movebacknames]<-
          currdata[key,moveback,drop=F][match(currdata[key,dayvar,drop=T]+1,
                                              currdata[key,dayvar,drop=T]),]
      }
      currdata<-na.omit(currdata)
      
      # Standardize
      if(scale){
        for(currvar in c("dv_lagged",allpreds)){
          currdata[[currvar]]<- (currdata[[currvar]]-mean(currdata[[currvar]]))/
            sd(currdata[[currvar]])
        }
      }
      
      # Run mod
      form<-paste0("dv_lagged ~ ",paste0(c(allpreds,covar),collapse=" + "),
                   " + (",paste0(allpreds,collapse=" + "),
                   ifelse(temporal=="correlated","|","||"),idvar,")") %>% as.formula()
      
      mlvmod<-lmer(formula=form,data=currdata,
                   control=lmerControl(optimizer="bobyqa",calc.derivs=F,
                                       optCtrl=list(maxfun=1e6)))
      
      # Return
      #temporal.models[[cpred]]<-
      mlvmod
    }
  names(temporal.models) <- allpreds
  temporal.check <- checkAsyncVARComponents(modlist=temporal.models,idvar=idvar)
  temporal.mat <- extractAsyncVARComponents(modlist=temporal.models,allpreds=allpreds,
                                            alpha=alpha)
  
  #################
  # Contemp. nets #
  #################
  
  # Extract all residuals
  tempresids<-list()
  for(i in seq_along(allpreds)){
    tempresids[[ allpreds[i] ]] <- model.frame(temporal.models[[i]]) %>% 
      mutate(resid=resid(temporal.models[[i]]))
  }
  residframe <- tempresids[[1]][,c(idvar,covar)]
  for(i in seq_along(allpreds)){
    residframe[[ allpreds[i] ]] <- tempresids[[i]]$resid
  }
  
  if(scale){
    for(currvar in allpreds){
      residframe[[currvar]]<-(residframe[[currvar]]-mean(residframe[[currvar]]))/sd(residframe[[currvar]])
    }
  }
  
  contemp.models <- foreach(cpred=allpreds,.packages=c("lmerTest")) %dopar% {
    
    # Contemp. net is limited to members of own "group"
    currpreds <- allpreds[predrank==predrank[allpreds==cpred]]
    
    # Run mod
    form<-paste0(cpred," ~ ",paste0(c(currpreds,covar),collapse=" + "),
                 " + (",paste0(currpreds,collapse=" + "),
                 ifelse(contemporaneous=="correlated","|","||"),idvar,")") %>% as.formula()
    
    mlvmod<-lmer(formula=form,data=residframe,
                 control=lmerControl(optimizer="bobyqa",calc.derivs=F,
                                     optCtrl=list(maxfun=1e6)))
    
    mlvmod
  }
  names(contemp.models) <- allpreds
  registerDoSEQ()
  stopCluster(cl=cl)
  
  contemp.check <- checkAsyncVarComponent(modlist=contemp.models,idvar=idvar)
  contemp.mat <- extractAsyncVARComponents(modlist=contemp.models,allpreds=allpreds,
                                           alpha=alpha)
  diag(contemp.mat) <- 0
  contemp.mat[is.na(contemp.mat)] <- 0
  
  # Symmetrize
  zeromask <- contemp.mat==0 | t(contemp.mat)==0
  contemp.mat <- (contemp.mat + t(contemp.mat))/2
  contemp.mat[zeromask] <- 0
  
  ############################
  # Between-subjects network #
  ############################
  
  personmeans <- as.data.frame(lapply(allpreds,function(x){
    tapply(X=data[[x]],INDEX=data[[idvar]],FUN=mean,na.rm=TRUE)
  }))
  personmeans %<>% scale()
  personcovars <- as.data.frame(lapply(covar,function(x){
    tapply(X=data[[x]],INDEX=data[[idvar]],FUN=function(y){y[[1]]})
  }))
  personmeans<-cbind(personmeans,personcovars)
  
  between.mat <- matrix(NA,
                        nrow=length(allpreds),
                        ncol=length(allpreds),
                        dimnames=list(allpreds,allpreds))
  
  between.models <- list()
  for(cpred in allpreds){
    form <- as.formula(paste(cpred,"~",paste0(c(covar,allpreds[allpreds!=cpred]),
                                              collapse="+")))
    mod <- lm(data=personmeans,formula=form)
    between.models[[cpred]] <- mod
    coefs <- summary(mod)$coefficients
    sigpreds <- coefs[allpreds[allpreds!=cpred],1]
    sigpreds[coefs[allpreds[allpreds!=cpred],4]>alpha] <- 0
    sigpreds[[cpred]] <- 0
    between.mat[,cpred] <- sigpreds[allpreds]
  }
  
  # Symmetrize
  zeromask <- between.mat==0 | t(between.mat)==0
  between.mat <- (between.mat + t(between.mat))/2
  between.mat[zeromask] <- 0

  ##########
  # Output #
  ##########
  
}

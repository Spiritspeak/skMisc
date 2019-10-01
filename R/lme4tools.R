# lme4tools

#' Extract random terms from a lme4 formula
#'
#' @param form A formula
#'
#' @return A named list containing character vectors with random terms; names are group variables.
#' @export
#'
#' @examples ExtractRandomTerms(grade ~ ChildIQ * TeacherSkill * SchoolType + 
#'                               (ChildIQ * TeacherSkill | School))
#' #$School
#' #[1] "ChildIQ"              "TeacherSkill"         "ChildIQ:TeacherSkill"
ExtractRandomTerms<-function(form){
  bars<-findbars(form)
  terms<-lapply(bars,FUN=function(x){
    x%<>%as.character()
    if(x[2]!="1"){
      x%<>%magrittr::extract(2)%>%reformulate%>%terms%>%attr("term.labels")
    }else{
      x<-"1"
    }
    x
  })
  names(terms)<-lapply(bars,FUN=function(x){
    x%<>%as.character()
    if(length(x)>1){
      x%<>%extract2(3)
    }else{
      x<-gsub(".*\\|","",x)%>%trimws
    }
    x
  })
  terms
}

#' Find all model terms that are not moderated by a higher-order interaction
#'
#' @param form a formula
#'
#' @return A character vector containing all model terms that are not moderated by a higher-order interaction.
#' @export
#'
#' @examples FindTopTerms(speed ~ skill + weight * friction)
#' #[1] "skill"           "weight:friction"
FindTopTerms<-function(form){
  #Do my 1's fit in another column's 1's?
  form%<>%terms()%>%attr("factors")
  dep<-rep(0,ncol(form))
  for(i in seq_len(ncol(form))){
    for(j in seq_len(ncol(form))[-i]){
      #dep[i]<- dep[i] + all(which(as.logical(form[,i])) %in% which(as.logical(form[,j])))
      dep[i] <- dep[i] + all(and(form[,i], form[,j])==form[,i])
    }
  }
  return(colnames(form)[!dep])
}

#' Parse a lme4 formula and return all main effects and interactions as separate terms
#' @param form 
#'
#' @return The same formula, but with all interactions and mai neffects as separate terms
#' @export
#'
#' @examples ExpandFormula(rt ~ pull * target + (pull * target | subjectid))
#' #rt ~ pull + target + pull:target + (pull + target + pull:target | subjectid)
ExpandFormula<-function(form){
  labs<-form%>%terms%>%attr("term.labels")
  norandos<- labs[!grepl("\\|",labs)]
  #norandos<-form%>%nobars%>%terms%>%attr("term.labels")
  randos<- ExtractRandomTerms(form)
  randlist<-character()
  for(i in seq_len(length(randos))){
    randlist[i]<-paste0("(",paste(randos[[i]],collapse=" + "),"|",names(randos)[[i]],")")
  }
  rhs<-paste(paste(norandos,collapse=" + "),paste(randlist,collapse=" + "),sep=" + ")
  formstring<-as.character(form)
  fullform<-paste(formstring[2],formstring[1],rhs) %>%as.formula
  return(fullform)
}


#' Remove all possible models with one unmoderated term removed
#'
#' @param form A formula
#' @param randeff The name of the group from which unmoderated terms should be removed. To remove from fixed effects, use \code{""} (the default).
#'
#' @return A list of formulas which have one unmoderated term removed each. The name of each list item is the term which was removed.
#' @export
#'
#' @examples RemoveTopTerms(a ~ b * c + d + (1|e))
#' #$d
#' #a ~ b + c + b:c + (1 | e)
#' #$`b:c`
#' #a ~ b + c + d + (1 | e)
RemoveTopTerms<-function(form,randeff=""){
  if(randeff==""){
    remform<-form %>% nobars %>% as.character() %>% extract(3) %>% reformulate %>% terms %>% attr("term.labels")
    remcomps<-form %>% nobars %>%FindTopTerms()
    redform<-paste(as.character(form)[2],as.character(form)[1], as.character(form)[3] %>% 
                     paste(remcomps,sep="-")) %>% sapply(FUN=as.formula,USE.NAMES=F)
    redform%<>%lapply(ExpandFormula)
    names(redform)<-remcomps
  }else{
    remform<-ExtractRandomTerms(form)[[randeff]]
    remcomps<-FindTopTerms(reformulate(paste(remform,collapse="+")))
    revcomp<-character()
    for(i in seq_len(length(remcomps))){
      revcomp[i]<- paste0("(",paste(remform[remform != remcomps[i]],collapse=" + ")," | ",randeff,")")
    }
    nonremform<-ExtractRandomTerms(form)
    nonremform<-nonremform[names(nonremform) != randeff]
    miscforms<-character()
    for(i in seq_len(length(nonremform))){
      miscforms[i]<- paste0("(",paste(nonremform[[i]],collapse=" + ")," | ",names(nonremform[i]),")")
    }
    redform<- paste((form %>% nobars %>% as.character() %>% extract(3)),revcomp,paste(miscforms,collapse=" + "),sep=" + ")
    redform<-paste(as.character(form)[2],as.character(form)[1], redform) %>% sapply(FUN=as.formula,USE.NAMES=F)
    names(redform)<-paste0("(",remcomps," | ",randeff,")")
  }
  
  return(redform)
}

ComputeLowerModels<-function(form,data,group="",...){
  args<-list(...)
  cluster<-makeCluster(detectCores()-1)
  registerDoParallel(cluster)
  testforms<-RemoveTopTerms(form,group)
  
  results<-
    foreach(currform=testforms,.packages="lmerTest") %dopar% {
      do.call(lmer,c(list(formula=currform,data=data),args))
    }
  
  stopCluster(cluster)
  results
}

ComputeLowerModels2<-function(model,data,group="",...){
  for(pack in c("magrittr","dplyr","tidyr","lme4","doParallel")){ require(pack,character.only=T) }
  form<-formula(model)
  data<-model.frame(model)
  #data<-model@call$data
  args<-list(...)
  testforms<-RemoveTopTerms(form,group)
  
  cluster<-makeCluster(detectCores()-1)
  registerDoParallel(cluster)
  results<-
    foreach(currform=testforms,.packages="lmerTest") %dopar% {
      do.call(lmer,c(list(formula=currform,data=data,REML=F),args))
    }
  stopCluster(cluster)
  names(results)<-names(testforms)
  
  warns<-sapply(results,function(x){ paste(x@optinfo$warnings,collapse="\n") })
  message(paste0("Warning in model ",names(results)[warns!=""],": ",warns[warns!=""]))
  
  anovatable<-AnovaTable(model,results)
  
  print(anovatable)
  return(invisible(list(anovatable=anovatable,models=results)))
}

#' Compare multilevel models
#'
#' @param ... Model objects to be compared
#' @param fullmodel A model to which all other models are to be compared; only use if \code{...} is not specified.
#' @param models Models to compare to \code{fullmodel}. Only use if \code{...} is not specified.
#' @param serial If TRUE, models are compared serially; if false, all models will be compared to the first.
#' @param suppress Character vector of column names to suppress in printed output.
#'
#' @return A data.frame containing model fit metrics such as AIC, BIC, marginal R-squared (the effect size of fixed effects only),
#' conditional R-squared (the effect size of all model terms), loglikelihood, deviance, and a likelihood ratio test.
#' @export 
#'
#' @examples 
#' 
AnovaTable<-function(...,fullmodel,models,serial=F,suppress=c("AIC","deviance","logLik")){
  if(missing(models) & missing(fullmodel)){    
    models<-list(...)
    fullmodel<-models[[1]]
    models<-models[-1]
  }
  mumin<-MuMIn::r.squaredGLMM(fullmodel)
  anovatable<-data.frame(Df=logLik(fullmodel)%>%attr("df"),
                         AIC=AIC(fullmodel),BIC=BIC(fullmodel),dBIC=0,
                         logLik=logLik(fullmodel),deviance=deviance(fullmodel),
                         R2m=mumin[1],R2c=mumin[2],dR2m=0,dR2c=0,
                         Chisq=NA,ChiDf=NA,P=NA)
  i<-1
  for(mod in models){
    i<-i+1
    mumin<-MuMIn::r.squaredGLMM(mod)
    anovatable%<>%rbind(
      data.frame(Df=logLik(mod)%>%attr("df"),
                 AIC=AIC(mod),
                 BIC=BIC(mod),
                 dBIC=BIC(mod) - ifelse(serial, anovatable[i-1,]$BIC, anovatable[1,]$BIC),
                 logLik=logLik(mod),
                 deviance=deviance(mod),
                 R2m=mumin[1],
                 R2c=mumin[2],
                 dR2m=mumin[1] - ifelse(serial, anovatable[i-1,]$R2m, anovatable[1,]$R2m),
                 dR2c=mumin[2] - ifelse(serial, anovatable[i-1,]$R2c, anovatable[1,]$R2c),
                 Chisq=NA,
                 ChiDf=NA,
                 P=NA))
  }
  for(i in 2:nrow(anovatable)){
    anovatable[i,]$Chisq<-anovatable[ifelse(serial,i-1,1),]$deviance-anovatable[i,]$deviance
    anovatable[i,]$ChiDf<-anovatable[ifelse(serial,i-1,1),]$Df - anovatable[i,]$Df
    anovatable[i,]$P<-pchisq(q = -anovatable[i,]$Chisq,df = anovatable[i,]$ChiDf)
  }
  
  modnames<-c("Full Model",names(models))
  if(length(modnames)<length(models)){ 
    modnames<-args2strings(...)
  }
  rownames(anovatable)<-modnames
  
  formulas<-c(fullmodel,models) %>% sapply(function(x){ x@call$formula })
  header<-paste0(modnames,": ",formulas,"\n",collapse="") %>%paste0("\n")
  
  anovatable<-structure(.Data=anovatable,class=c("anovatable","data.frame"),suppress=suppress,
                        header=header)
  return(anovatable)
}

print.anovatable<-function(x){
  attr(x,"header") %>% cat
  x<-x[,which(!(colnames(x) %in% attr(x,"suppress")))]
  print.data.frame(x,digits=3)
}
registerS3method("print","anovatable",print.anovatable)

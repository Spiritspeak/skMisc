#AATtools
# splithalf engine ####

#multicore splithalf
#' Compute bootstrapped split-half reliability for approach-avoidance task data
#' @description Compute bootstrapped split-half reliability for approach-avoidance task data. 
#' \code{aat_splithalf()} uses multicore computation, which is fast, but provides no clear output when there are errors.
#' \code{aat_splithalf_singlecore()} is much slower, but more easily debugged.
#' @param ds a longformat data.frame
#' @param subjvar Name of the participant identifier column
#' @param pullvar Name of the column indicating pull trials. Pull trials should either be represented by 1, or by the second level of a factor.
#' @param targetvar Name of the column indicating trials featuring the target stimulus. Target stimuli should either be represented by 1, or by the second level of a factor.
#' @param rtvar Name of the reaction time column.
#' @param iters Total number of desired iterations. At least 200 are recommended for reasonable confidence intervals; If you want to see plots of your data, 1 iteration is enough.
#' @param plot Create a scatterplot of the AAT scores computed from each half of the data from the last iteration. This is highly recommended, as it helps to identify outliers that can inflate or diminish the reliability.
#' @param algorithm Function (without brackets or quotes) to be used to compute AAT scores. See \link{aat_doublemeandiff} for a list of usable algorithms.
#' @param trialdropfunc Function (without brackets or quotes) to be used to exclude outlying trials in each half. The way you handle outliers for the reliability computation should mimic the way you do it in your regular analyses. 
#' It is recommended to exclude outlying trials when computing AAT scores using the mean double-dfference scores and multilevel scoring approaches, but not when using d-scores or median double-difference scores.
#' \code{prune_nothing} excludes no trials, while \code{trial_prune_3SD} excludes trials deviating more than 3SD from the mean per participant.
#' @param casedropfunc Function (without brackets or quotes) to be used to exclude outlying participant scores in each half. The way you handle outliers here should mimic the way you do it in your regular analyses.
#' \code{prune_nothing} excludes no participants, while \code{case_prune_3SD} excludes participants deviating more than 3SD from the sample mean.
#' @param ... Other arguments, to be passed on to the algorithm functions (see \code{algorithm} above)
#'
#' @return A list, containing the mean bootstrapped split-half reliability, bootstrapped 95% confidence intervals, a list of data.frames used over each iteration, and a vector containing the split-half reliability of each iteration.
#' @export
#'
#' @examples #Not Run
#' aat_splithalf(ds=ds2,subjvar="subjectid",pullvar="is_pull",targetvar="is_food",
#'               rtvar="rt",iters=1000,trialdropfunc=trial_prune_3SD,
#'               casedropfunc=case_prune_3SD,plot=T,algorithm=aat_dscore)
#' #Mean reliability: 0.521959
#' #Spearman-Brown-corrected r: 0.6859041
#' #95%CI: [0.4167018, 0.6172474]
#' 
#' #Multilevel Splithalf
#' aat_splithalf(ds=ds2,subjvar="subjectid",pullvar="is_pull",targetvar="is_food",
#'               rtvar="rt",iters=100,trialdropfunc=trial_prune_3SD,
#'               casedropfunc=case_prune_3SD,plot=T,algorithm=aat_multilevelscore,
#'               formula = "rt ~ is_pull * is_food + (is_pull * is_food | subjectid)",
#'               aatterm = "is_pull:is_food")
#' #Mean reliability: 0.5313939
#' #Spearman-Brown-corrected r: 0.6940003
#' #95%CI: [0.2687186, 0.6749176]
aat_splithalf<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,plot=T,
                        algorithm=c(aat_doublemeandiff,aat_doublemediandiff,aat_dscore,aat_multilevelscore),
                        trialdropfunc=c(prune_nothing,trial_prune_3SD),
                        casedropfunc=c(prune_nothing,case_prune_3SD),
                        ...){
  for(pack in c("magrittr","dplyr","tidyr","lme4","doParallel")){ require(pack,character.only=T) }
  args<-list(...)
  
  if(deparse(substitute(algorithm))=="aat_multilevelscore_single" & !any(c("formula","aatterm") %in% names(args))){
    args$formula<-paste(rtvar,"~",1,"+(",pullvar,"*",targetvar,"|",subjvar,")")
    args$aatterm<-paste(pullvar,":",targetvar)
    warning("No multilevel formula or AAT-term provided. Defaulting to formula ",args$formula," and AAT-term ",args$aatterm)
  }

  cluster<-makeCluster(detectCores()-1)#,outfile="splithalfmessages.txt")
  registerDoParallel(cluster)
  results<-
    foreach(iter = seq_len(iters), .packages=c("magrittr","dplyr","tidyr","lme4")) %dopar% {
      iterds<-ds%>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
        mutate(key=sample(n())%%2)%>%ungroup()
      iterds<-do.call(trialdropfunc,list(ds=iterds,subjvar=subjvar,rtvar=rtvar))
      #abds<-do.call(algorithm,c(list(ds=iterds,subjvar=subjvar,pullvar=pullvar,
      #                               targetvar=targetvar,rtvar=rtvar),args)) #Use for the deprecated algorithm functions
      
      half0set<-iterds%>%filter(key==0)
      half1set<-iterds%>%filter(key==1)
      abds<-merge(
        do.call(algorithm,c(list(ds=half0set,subjvar=subjvar,pullvar=pullvar,
                                 targetvar=targetvar,rtvar=rtvar),args)),
        do.call(algorithm,c(list(ds=half1set,subjvar=subjvar,pullvar=pullvar,
                                 targetvar=targetvar,rtvar=rtvar),args)),
        by=subjvar,suffixes=c("half0","half1"))
      
      abds<-do.call(casedropfunc,list(ds=abds))
      currcorr<-cor(abds$abhalf0,abds$abhalf1,use="complete.obs") 
      cat("\rCorr for iter ",iter," is ", round(currcorr,digits=2),rep(".",iter/iters*10)," ",sep="")
      list(corr=currcorr,abds=abds)
    }
  stopCluster(cluster)
  cors<-sapply(results,FUN=function(x){x$corr})
  cors%<>%sort
  cat("\nMean reliability: ",mean(cors),
      "\nSpearman-Brown-corrected r: ",SpearmanBrown(mean(cors)),
      "\n95%CI: [", cors[round(iters*0.025)], ", ", cors[round(iters*0.975)],"]\n",
      sep="")
  
  if(plot){
    abds<-results[[length(results)]]$abds
    plot(abds$abhalf0,abds$abhalf1,pch=20,xlab="Half 1 computed bias",ylab="Half 2 computed bias")
    text(abds$abhalf0,abds$abhalf1,abds[[subjvar]],cex= 0.7, pos=3, offset=0.3)
  }
  #cat(scan("splithalfmessages.txt",what=character(),quiet=TRUE))
  output<-list(rsplithalf=mean(cors),
               lowerci=cors[round(iters*0.025)],
               upperci=cors[round(iters*0.975)],
               itercors=sapply(results,function(x){ x$corr }),
               iterdata=lapply(results,function(x){ x$abds }))
  return(invisible(output))
}

#Singlecore splithalf (slower but produces output)
#' @rdname aat_splithalf
aat_splithalf_singlecore<-function(ds,subjvar,pullvar,targetvar,rtvar,iters,plot=F,
                                   algorithm=c(aat_doublemeandiff,aat_doublemediandiff,aat_dscore,aat_multilevelscore),
                                   trialdropfunc=c(prune_nothing,trial_prune_3SD),
                                   casedropfunc=c(prune_nothing,case_prune_3SD),
                                   ...){
  for(pack in c("magrittr","dplyr","tidyr","lme4")){ require(pack,character.only=T) }
  args<-list(...)
  
  cors<-vector()
  for(iter in 1:iters){
    iterds<-ds%>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
      mutate(key=sample(n())%%2)%>%ungroup()
    
    iterds<-do.call(trialdropfunc,list(ds=iterds,subjvar=subjvar,rtvar=rtvar))
    
    abds<-do.call(algorithm,c(list(ds=iterds,subjvar=subjvar,pullvar=pullvar,
                                   targetvar=targetvar,rtvar=rtvar),args))
    
    abds<-do.call(casedropfunc,list(ds=abds))
    
    cors[iter]<-cor(abds$abhalf0,abds$abhalf1,use="complete.obs")
    cat("\rCorr for iter ",iter," is ", round(cors[iter],digits=2),rep(".",iter/iters*10)," ",sep="")
    #cat("\r[",rep("/",iter/iters*20),rep(" ",(iters-iter)/iters*20),"]",sep="")
  }
  cors%<>%sort
  cat("\nMean reliability: ",mean(cors),
      "\nSpearman-Brown-corrected r: ",SpearmanBrown(mean(cors)),
      "\n95%CI: [", cors[round(iters*0.025)], ", ", cors[round(iters*0.975)],"]\n",
      sep="")
  
  if(plot){
    plot(abds$abhalf0,abds$abhalf1,pch=20,xlab="Half 1 computed bias",ylab="Half 2 computed bias",
         main="Scatterplot for last split-half iteration")
    text(abds$abhalf0,abds$abhalf1,abds[[subjvar]],cex= 0.7, pos=3, offset=0.3)
  }
  invisible(list(iterdf=abds,rsplithalf=mean(cors),
                 lowerci=cors[round(iters*0.025)],
                 upperci=cors[round(iters*0.975)]))
}


# Outlier removing algorithms ####

prune_nothing<-function(ds,...){
  ds
}

trial_prune_3SD<-function(ds,subjvar,rtvar){
  ds$ol.rt.var<-ds[[rtvar]]
  ds%>%group_by(!! sym(subjvar),key)%>%
    mutate( prune.ol.mean=mean(ol.rt.var,na.rm=T),prune.ol.sd=sd(ol.rt.var,na.rm=T))%>%
    dplyr::filter((ol.rt.var<prune.ol.mean+3*prune.ol.sd) &
                    (ol.rt.var>prune.ol.mean-3*prune.ol.sd) ) %>%
    dplyr::select(-prune.ol.mean,-prune.ol.sd,-ol.rt.var)
}

#to add: error treatment. remove vs replace (for D score)

case_prune_3SD<-function(ds){
  ds%>%filter((abhalf0 < mean(abhalf0,na.rm=T)+3*sd(abhalf0,na.rm=T) &
                     abhalf0 > mean(abhalf0,na.rm=T)-3*sd(abhalf0,na.rm=T)) &
                    (abhalf1 < mean(abhalf1,na.rm=T)+3*sd(abhalf1,na.rm=T) &
                       abhalf1 > mean(abhalf1,na.rm=T)-3*sd(abhalf1,na.rm=T)))
}



# Score computation algorithms ####

#' AAT score computation algorithms
#' @description
#' \itemize{
#' \item \code{aat_doublemeandiff} computes a mean-based double-difference score: \code{(mean(push_target)-mean(pull_target)) - (mean(push_control) - mean(pull_control))}
#' \item \code{aat_doublemediandiff} computes a median-based double-difference score: \code{(median(push_target)-median(pull_target)) - (median(push_control) - median(pull_control))}
#' \item \code{aat_dscore} computes D-scores (see Greenwald, Nosek, and Banaji, 2003): \code{((mean(push_target)-mean(pull_target)) - (mean(push_control) - mean(pull_control)))/sd(participant_reaction_times)}
#' \item \code{aat_multilevelscore} fits a multilevel model using lme4 and extracts a random effect serving as AAT score. When using this function, additional arguments must be provided: 
#' \itemize {
#' \item \code{formula} - a quoted formula to fit to the data;
#' \item \code{aatterm} the quoted random effect within the subject variable that indicates the approach bias; this is usually the interaction of the pull and target terms.
#' }
#' }
#' @param ds A long-format data.frame
#' @param subjvar Column name of the participant identifier variable
#' @param pullvar Column name of the movement variable (0: avoid; 1: approach)
#' @param targetvar Column name of the stimulus category variable (0: control stimulus; 1: target stimulus)
#' @param rtvar Column name of the reaction time variable
#' @param ... Other arguments passed on by functions (ignored)
#'
#' @return A data.frame containing participant number and computed AAT score.
#' @export
aat_doublemeandiff<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  ds%<>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
    summarise(means =mean(!! sym(rtvar),na.rm=T))%>%group_by()
  abds<-ds%>%mutate(groupcol=paste0("pull",!!sym(pullvar),"target",!!sym(targetvar)))%>%
    dplyr::select(!! subjvar, groupcol,means)%>%tidyr::spread("groupcol",value="means",drop=F)
  abds%<>%
    mutate(ab=(pull0target1-pull1target1)-(pull0target0-pull1target0))%>%
    dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @rdname aat_doublemeandiff
#' @export
aat_doublemediandiff<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  ds%<>%group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar))%>%
    summarise(medians =median(!! sym(rtvar),na.rm=T))%>%group_by()
  abds<-ds%>%mutate(groupcol=paste0("pull",!!sym(pullvar),"target",!!sym(targetvar)))%>%
    dplyr::select(!! subjvar, groupcol,medians)%>%tidyr::spread("groupcol",value="medians",drop=F)
  abds%<>%
    mutate(ab=(pull0target1-pull1target1)-(pull0target0-pull1target0))%>%
    dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @rdname aat_doublemeandiff
#' @export
aat_dscore<-function(ds,subjvar,pullvar,targetvar,rtvar,...){
  setmeans <- ds %>% group_by(!! sym(subjvar), !! sym(pullvar), !! sym(targetvar)) %>%
    summarise(means =mean(!! sym(rtvar),na.rm=T)) %>% group_by() %>%
    mutate(groupcol=paste0("pull",!! sym(pullvar),"target",!! sym(targetvar))) %>%  
    dplyr::select(-!!sym(pullvar),-!!sym(targetvar)) %>% tidyr::spread(key=groupcol,value=means)
  
  setsds <- ds%>%group_by(!! sym(subjvar)) %>%
    summarise(sds =sd(!! sym(rtvar),na.rm=T))
  
  abds <- merge(setmeans,setsds,by=subjvar) %>% 
    mutate(ab=((pull0target1-pull1target1)-(pull0target0-pull1target0))/sds) %>%
    dplyr::select(!!sym(subjvar),ab)
  return(abds)
}

#' @rdname aat_doublemeandiff
#' @param formula A character string containing a formula to fit to the data and derive multilevel scores from
#' @param aatterm The random term, grouped under the subject variable, which represents the approach bias. Usually this is the interaction of the pull and target terms.
#' @export
aat_multilevelscore<-function(ds,subjvar,formula,aatterm,...){
  fit<- lme4::lmer(as.formula(formula),data=ds,control=
                lme4::lmerControl(optimizer="bobyqa",optCtrl = list(maxfun = 2e6)))
  output<-data.frame(ab=-lme4::ranef(fit,whichel=subjvar)[[subjvar]][[aatterm]])
  testset<<-fit
  output[[subjvar]]<-fit@flist[[subjvar]]%>%levels
  return(output)
}

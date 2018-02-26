evaluate <- function(model, data, stat, as.array=FALSE){
  ## output initialisation
  eval <- NULL
  
  if(inherits(model, "BIOMOD.models.out") | inherits(model,"BIOMOD.EnsembleModeling.out")){
    eval <- .evaluate.biomod2.models.out(mod=model,data=data,stat=stat)
    if(as.array) eval <- LIST_to_ARRAY(eval)
  } else if(inherits(model, "biomod2_model") ){
    eval <- .evaluate.biomod2.formal.models(mod=model,data=data,stat=stat)
  } else {
    cat("\n\n! invalid model input => nothing returned !")
  } 
  return(eval)
}

.evaluate.biomod2.formal.models <- function(mod, data, stat='TSS'){
  obs <- data[,mod@resp_name, drop=T]
  fit <- predict(mod, data[,mod@expl_var_names, drop=F], on_0_1000=T)
  if(stat != 'ROC'){
    thresh <- try(mod@model_evaluation[stat,'Cutoff'],silent=T)
  } else { thresh <- 500  } # no need to threshold
  if(inherits(thresh,"try-error")){
    thresh <- mod@model_evaluation[1,'Cutoff']
    cat("\n! no 'true' threshold defined for",stat,"... ", rownames(mod@model_evaluation)[1], "ones' taken.")      
  }
  eval <- Find.Optim.Stat(Stat=stat, Fit=fit, Obs=obs, Fixed.thresh=thresh)
  return(eval)
}

.evaluate.biomod2.models.out <- function(mod, data, stat='TSS'){
  formal.models.names <- BIOMOD_LoadModels(mod)
  eval <- lapply(formal.models.names, function(x, data, stat){
    xx <- get(x)
    eval <- vapply(stat, 
                   function(s){.evaluate.biomod2.formal.models(mod=xx, data=data, stat=s)},
                   FUN.VALUE = c(Evaluating.data=0, Cutoff=500, Sensitivity=0, Specificity=0))
    return(t(eval))
  }, data=data, stat=stat)
  names(eval) <- formal.models.names
  return(eval)
}

.evaluate.biomod2.ensemble.models.out <- function(mod, data, stat='TSS'){
  formal.models.names <- BIOMOD_LoadModels(mod)
  eval <- lapply(formal.models.names, function(x, data, stat){
    xx <- get(x)
    eval <- vapply(stat, 
                   function(s){.evaluate.biomod2.formal.models(mod=xx, data=data, stat=s)},
                   FUN.VALUE = c(Evaluating.data=0, Cutoff=500, Sensitivity=0, Specificity=0))
    return(t(eval))
  }, data=data, stat=stat)
  names(eval) <- formal.models.names
  return(eval)
}

## TEST ##
# em <- evaluate(mod=mod.out,data=data,stat=c('TSS','ROC'))





######## LOW LEVEL FUCTIONS ##############

##' @name Find.Optim.Stat
##' @title Calculate the best score according to a given evaluation method
##' 
##' @description \code{Find.Optim.Stat} is an internal \pkg{biomod2} function 
##' to find the threshold to convert continuous values into binary ones leading 
##' to the best score for a given evaluation metric.
##' 
##' @usage
##'   Find.Optim.Stat(Stat='TSS',
##'                   Fit,
##'                   Obs,
##'                   Nb.thresh.test = 100,
##'                   Fixed.thresh = NULL)
##' 
##' @param Stat either 'ROC', TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 
##'         'POFD', 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR' or 'ORSS'
##' @param Fit vector of fitted values (continuous)
##' @param Obs vector of observed values (binary)
##' @param Nb.thresh.test integer, the numer of thresholds tested over the 
##'        range of fitted value
##' @param Fixed.thresh integer, if not \code{NULL}, the only threshold value tested
##'                 
##' @details
##'   Please refer to \code{\link[biomod2]{BIOMOD_Modeling}} to get more information about this metrics.
##'   If you give a \code{Fixed.thresh}, no optimisation will be done. Only the score for this threshold will be returned.
##'                 
##' @return
##'   A 1 row x 4 column \code{matrix} :
##'   \itemize{
##'     \item{\code{best.iter}:}{ the best score obtained for chosen statistic}
##'     \item{\code{cutoff}:}{ the associated cut-off used for transform fitted vector into binary}
##'     \item{\code{sensibility}:}{ the sensibility with this threshold}
##'     \item{\code{specificity}:}{ the specificity with this threshold}
##'   }
##'                
##' @author Damien Georges
##' 
##' @seealso 
##'   \code{\link[biomod2]{BIOMOD_Modeling}}, 
##'   \code{\link[biomod2]{getStatOptimValue}}, 
##'   \code{\link[biomod2]{calculate.stat}} 
##'                 
##' @examples
##'   a <- sample(c(0,1),100, replace=TRUE)
##'                 
##'   ##' random drawing
##'   b <- runif(100,min=0,max=1000)
##'   Find.Optim.Stat(Stat='TSS',
##'                   Fit=b,
##'                   Obs=a)
##'                 
##'   ##' biased drawing
##'   BiasedDrawing <- function(x, m1=300, sd1=200, m2=700, sd2=200){
##'     return(ifelse(x<0.5, rnorm(1,m1,sd1), rnorm(1,m2,sd2)))
##'   }
##'                 
##'   c <- sapply(a,BiasedDrawing)
##'                 
##'   Find.Optim.Stat(Stat='TSS',
##'                   Fit=c,
##'                   Obs=a,
##'                   Nb.thresh.test = 100)  
##'                 
##'                 
##' @keywords models, options, evaluate, evaluation
Find.Optim.Stat <- function(Stat = 'TSS',
                            Fit, 
                            Obs,
                            Nb.thresh.test = 100,
                            Fixed.thresh = NULL){
  
  ## remove all uninite values
  to_keep <- ( is.finite(Fit) & is.finite(Obs) )
  Fit <- Fit[to_keep]
  Obs <- Obs[to_keep]

  ## guess fit value scale (e.g. 0-1 for a classic fit or 0-1000 for a biomod2 model fit)
  fit.scale <- .guess.scale(Fit)
  
  ## check some data are still here.
  if(!length(Obs) | !length(Fit)){
    cat("Non finite obs or fit available => model evaluation skipped !")
    eval.out <- matrix(NA,1,4, dimnames = list(Stat, c("best.stat","cutoff","sensitivity","specificity")))
    return(eval.out)
  }
  
  if(length(unique(Obs)) == 1 | length(unique(Fit)) == 1){
#     warning("\nObserved or fited data contains only a value.. Evaluation Methods switched off\n",immediate.=T)
#     best.stat <- cutoff <- true.pos <- sensitivity <- true.neg <- specificity <- NA  
      warning("\nObserved or fited data contains a unique value.. Be carefull with this models predictions\n",immediate.=T)
      #best.stat <- cutoff <- true.pos <- sensitivity <- true.neg <- specificity <- NA    
  } #else {
    if(Stat != 'ROC'){
      StatOptimum <- getStatOptimValue(Stat)
      if(is.null(Fixed.thresh)){ # test a range of threshold to get the one giving the best score
        if(length(unique(Fit)) == 1){
          valToTest <- unique(Fit)
          ## add 2 values to test based on maen with 0 and the guessed max of Fit (1 or 1000) 
          valToTest <- round(c(mean(c(fit.scale["min"],valToTest)), 
                               mean(c(fit.scale["max"],valToTest))))
        } else{
#           mini <- max(min(quantile(Fit,0.05, na.rm=T), na.rm=T),0)
#           maxi <- min(max(quantile(Fit,0.95, na.rm=T), na.rm=T),1000)
          mini <- max(min(Fit, na.rm=T), fit.scale["min"], na.rm = T)
          maxi <- min(max(Fit, na.rm=T), fit.scale["max"], na.rm = T) 
          valToTest <- try(unique( round(c(seq(mini, maxi, 
                                               length.out = Nb.thresh.test), 
                                           mini, maxi)) ))
          if(inherits(valToTest, "try-error")){
            valToTest <- seq(fit.scale["min"], fit.scale["max"], 
                             length.out = Nb.thresh.test)
          }
          # deal with unique value to test case
          if(length(valToTest)<3){
            valToTest <- round(c(mean(fit.scale["min"],mini), valToTest, 
                                 mean(fit.scale["max"], maxi)))
          }
        }
#         valToTest <- unique( c(seq(mini,maxi,by=Precision), mini, maxi) )        
      } else{
        valToTest <- Fixed.thresh
      }
      
      calcStat <- sapply(lapply(valToTest, function(x){return(table(Fit>x,Obs))} ), calculate.stat, stat=Stat)
      
      # scal on 0-1 ladder.. 1 is the best
      calcStat <- 1 - abs(StatOptimum - calcStat)
      
      best.stat <- max(calcStat, na.rm=T)
      
      cutoff <- median(valToTest[which(calcStat==best.stat)]) # if several values are selected

      misc <- table(Fit >= cutoff, Obs)
      misc <- .contagency.table.check(misc)
      true.pos <- misc['TRUE','1']
      true.neg <- misc['FALSE','0']
      specificity <- (true.neg * 100)/sum(misc[,'0'])
      sensitivity <- (true.pos * 100)/sum(misc[,'1'])
    } else{
#       require(pROC,quietly=T)
      roc1 <- pROC::roc(Obs, Fit, percent=T, direction="<")
      roc1.out <- pROC::coords(roc1, "best", ret=c("threshold", "sens", "spec"))
      ## if two optimal values are returned keep only the first one
      if(!is.null(ncol(roc1.out))) roc1.out <- roc1.out[, 1]
      best.stat <- as.numeric(pROC::auc(roc1))/100
      cutoff <- as.numeric(roc1.out["threshold"])
      sensitivity <- as.numeric(roc1.out["sensitivity"])
      specificity <- as.numeric(roc1.out["specificity"])
    }
  #}
  eval.out <- cbind(best.stat,cutoff,sensitivity,specificity)
  rownames(eval.out) <- Stat

  return(eval.out)
}

getStatOptimValue <- function(stat){
  if(stat == 'TSS') return(1)
  if(stat == 'KAPPA') return(1)
  if(stat == 'PBC') return(1)
  if(stat == 'ACCURACY') return(1)
  if(stat == 'BIAS') return(1)
  if(stat == 'POD') return(1)
  if(stat == 'FAR') return(0)
  if(stat == 'POFD') return(0)
  if(stat == 'SR') return(1)
  if(stat == 'CSI') return(1)
  if(stat == 'ETS') return(1)
  if(stat == 'HK') return(1)
  if(stat == 'HSS') return(1)
  if(stat == 'OR') return(1000000)
  if(stat == 'ORSS') return(1)
}

calculate.stat <-
function(Misc, stat='TSS')
{
  # Contagency table checking
  Misc <- .contagency.table.check(Misc)
  
  # Defining Classification index
  hits <- Misc['TRUE','1']
  misses <- Misc['FALSE','1']
  false_alarms <- Misc['TRUE','0']
  correct_negatives <- Misc['FALSE','0']
  
  total <- sum(Misc)
  forecast_1 <- sum(Misc['TRUE',])
  forecast_0 <- sum(Misc['FALSE',])
  observed_1 <- sum(Misc[,'1'])
  observed_0 <- sum(Misc[,'0'])
  
  # Calculating choosen evaluating metric
  if(stat=='TSS'){
    return( (hits/(hits+misses)) + (correct_negatives/(false_alarms+correct_negatives)) -1 )
  }
  
  if(stat=='KAPPA'){
    Po <- (1/total) * (hits + correct_negatives)
    Pe <- ((1/total)^2) * ((forecast_1 * observed_1) + (forecast_0 * observed_0))
    return( (Po - Pe) / (1-Pe) )
  }
  
  if(stat=='ACCURACY'){
    return( (hits + correct_negatives) / total)
  }
  
  if(stat=='BIAS'){
    return( (hits + false_alarms) / (hits + misses))
  }
  
  if(stat=='POD'){
    return( hits / (hits + misses))
  }
  
  if(stat=='FAR'){
    return(false_alarms/(hits+false_alarms))
  }
  
  if(stat=='POFD'){
    return(false_alarms / (correct_negatives + false_alarms))
  }
  
  if(stat=='SR'){
    return(hits / (hits + false_alarms))
  }
  
  if(stat=='CSI'){
    return(hits/(hits+misses+false_alarms))
  }
  
  if(stat=='ETS'){
    hits_rand <- ((hits+misses)*(hits+false_alarms)) / total
    return( (hits-hits_rand) / (hits+misses+false_alarms-hits_rand))
  }
  
  if(stat=='HK'){
    return((hits/(hits+misses)) - (false_alarms/(false_alarms + correct_negatives)))
  }
  
  if(stat=='HSS'){
    expected_correct_rand <- (1/total) * ( ((hits+misses)*(hits+false_alarms)) +
      ((correct_negatives + misses)*(correct_negatives+false_alarms)) )
    return((hits+correct_negatives-expected_correct_rand) / (total - expected_correct_rand))
  }
  
  if(stat=='OR'){
    return((hits*correct_negatives)/(misses*false_alarms))
  }
  
  if(stat=='ORSS'){
    return((hits*correct_negatives - misses*false_alarms) / (hits*correct_negatives + misses*false_alarms))
  }
  
  if(stat=="BOYCE"){
    
  }
  
}

.contagency.table.check <- function(Misc){
  # Contagency table checking
  if(dim(Misc)[1]==1){
    if(row.names(Misc)[1]=="FALSE"){
      Misc <- rbind(Misc, c(0,0))
      rownames(Misc) <- c('FALSE','TRUE')
    } else{
      a <- Misc
    	Misc <- c(0,0)
  		Misc <- rbind(Misc, a)
      rownames(Misc) <- c('FALSE','TRUE')
  	}
  }
  
  if(ncol(Misc) != 2 | nrow(Misc) !=2 ){
    Misc = matrix(0, ncol=2, nrow=2, dimnames=list(c('FALSE','TRUE'), c('0','1')))
  }
  
  if((sum(colnames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) | (sum(rownames(Misc) %in% c('FALSE','TRUE','0','1')) < 2) ){
    stop("Unavailable contagency table given")
  }
  
  if('0' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='0')] <- 'FALSE'
  if('1' %in% rownames(Misc)) rownames(Misc)[which(rownames(Misc)=='1')] <- 'TRUE'  
    
  return(Misc)
}


.guess.scale <- function(Fit){
  min <- 0
  max <- ifelse(max(Fit, na.rm = TRUE) <= 1, 1, 1000)
  out <- c(min, max)
  names(out) <- c("min", "max")
  return(out)
}
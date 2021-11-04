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



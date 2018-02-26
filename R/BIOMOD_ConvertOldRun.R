BIOMOD_ConvertOldRun <- function(savedObj, path = NULL){
  .bmCat("BIOMOD results migration")
  
  compress.arg = TRUE #ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  # 1. Check path exists and all objects needed exists too -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(!file.exists(savedObj)){
    stop("Input object doesn't exist")
  }
  
  if(is.null(path)){
    path = sub(tail(unlist(strsplit(savedObj,'/')),1), '', savedObj)
  } else{ # add / at the end of path
    if(tail(unlist(strsplit(path)),1) != '/'){
      path = paste(path,"/",sep="")
    }
  }
  
  # 2. Keep image of current workSpace -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   save.image("WS_tmp.Rdata")
  
  # 3. Load file and extract information -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  load(paste(savedObj))
  
  if(!exists('DataBIOMOD') | !exists('Biomod.material')){
    stop("DataBIOMOD or Biomod.material Object not in your input object, please check it!")
  }
  
  ### tips to remove some compiling warnings 
  if(!exists('Biomod.material')){ Biomod.material <- NULL }
  if(!exists('DataBIOMOD')){ DataBIOMOD <- NULL }
  if(!exists('Biomod.PA.sample')){ Biomod.PA.sample <- NULL }
  if(!exists('Evaluation.results.Roc')){ Evaluation.results.Roc <- NULL }
  if(!exists('Evaluation.results.TSS')){ Evaluation.results.TSS <- NULL }
  if(!exists('Evaluation.results.Kappa')){ Evaluation.results.Kappa <- NULL }
  if(!exists('VarImportance')){ VarImportance <- NULL }
  
  sp.names <- Biomod.material$species.names
  
  NewModelObj <- lapply(sp.names,function(sp.name){
    cat("\n\n", sp.name, 'run convertion...')
    
    dir.create(sp.name,showWarnings=FALSE)
    dir.create(paste(sp.name,"/.BIOMOD_DATA",sep=""),showWarnings=FALSE)
    
    models.out <- new('BIOMOD.models.out',
                      sp.name = sp.name,
                      expl.var.names = Biomod.material$VarNames,
                      rescal.all.models = FALSE)
    
    #   3.1 BIOMOD.formated.data creation
    cat("\n\tBIOMOD.formated.data creation")
    data <- BIOMOD.formated.data(sp = DataBIOMOD[,sp.name],
                                   env = DataBIOMOD[,Biomod.material$VarNames],
                                   xy = data.frame(),
                                   sp.name = sp.name)
    
    if(Biomod.material$NbRepPA > 0){
      PAtmp <- matrix(FALSE,
                      nrow=nrow(DataBIOMOD), 
                      ncol=length(Biomod.PA.sample[[sp.name]]),
                      dimnames=list(NULL, names(Biomod.PA.sample[[sp.name]])))
      for(i in 1: ncol(PAtmp)){
        PAtmp[Biomod.PA.sample[[sp.name]][[i]],i] <- TRUE
      }
      data <- new('BIOMOD.formated.data.PA',
                    sp.name = data@sp.name,
                    coord = data@coord,
                    data.env.var = data@data.env.var,
                    data.species = data@data.species,
                    PA = data.frame(PAtmp))
    }
    # save Input Data
    save(data, file = paste(models.out@sp.name,"/.BIOMOD_DATA/formated.input.data",sep=""), compress=compress.arg)
    models.out@formated.input.data@inMemory <- FALSE
    models.out@formated.input.data@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/formated.input.data",sep="")
    rm(PAtmp,data)
    
    #   3.2 BIOMOD.Model.Options creation
    cat("\n\tBIOMOD.Model.Options creation")
    models.options <- BIOMOD_ModelingOptions()
    # save Model Options
    save(models.options, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.options",sep=""),  compress=compress.arg)
    models.out@models.options@inMemory <- FALSE
    models.out@models.options@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.options",sep="")
    rm(models.options)
    
    #   3.3 Converting Models Computed
    cat("\n\tConverting Models Computed")
    
    if(!file.exists(paste(path,"models/",sep=""))){
      stop("models directory doesn't exist!")
    }
    dir.create(paste(models.out@sp.name,"/models",sep=""), showWarnings=FALSE)
    if(file.exists(paste(path,"models/scaling_models",sep=""))){
      dir.create(paste(models.out@sp.name,"/models/scaling_models",sep=""), showWarnings=FALSE)
    }
    
    old.mod.computed <- list.files(path = paste(path,"models/",sep=""),
                               pattern = models.out@sp.name,
                               recursive = TRUE)
         
    new.mod.computed <- old.mod.computed
    new.mod.computed <- sapply(new.mod.computed, function(x){
      if(length(grep('PA',x)) > 0){ # pseudo absences done case 
        if(length(grep('rep',x)) > 0){ # repetition
          x <- gsub('rep','RUN',x)
        } else{
          x <- paste(x,'_Full',sep='')
        }
        if(length(grep('Rmod_',x)) > 0){ # scaled models
          x <- paste(gsub('Rmod_','',x),'_scaled',sep='')
          
          x.str <- unlist(strsplit(gsub('scaling_models/','',x),'_'))
          x <- paste(x.str[1], x.str[3], x.str[4], x.str[2], x.str[length(x.str)], sep='_')
          x <- paste('scaling_models/',x,sep='')
        } else{
          x.str <- unlist(strsplit(x,'_'))
          x <- paste(x.str[1], x.str[3], x.str[4], x.str[2], sep='_')
        }
      } else{
        x <- gsub('_full','',x)
        if(length(grep('rep',x)) > 0){ # repetition
          x <- gsub('rep','RUN',x)
        } else{
          x <- paste(x,'_Full',sep='')
        }
        if(length(grep('Rmod_',x)) > 0){ # scaled models
          x <- paste(gsub('Rmod_','',x),'_scaled',sep='')
          
          x.str <- unlist(strsplit(gsub('scaling_models/','',x),'_'))
          x <- paste(x.str[1], x.str[3], x.str[2], x.str[length(x.str)], sep='_')
          x <- paste('scaling_models/',x,sep='')
        } else{
          x.str <- unlist(strsplit(x,'_'))
          x <- paste(x.str[1], x.str[3], x.str[2], sep='_')
        }
      }
      return(x)
    })
    
    
    # coping the files with appropriated names
    lapply(1:length(old.mod.computed), function(x){
      file.copy(from = paste(path, "models/", old.mod.computed[x], sep=""),
                to = paste(models.out@sp.name, "/models/", new.mod.computed[x], sep=""),
                overwrite = TRUE,
                recursive = FALSE,
                copy.mode = TRUE )
    })
      
    models.out@models.computed <- unique(as.character(gsub('_scaled','',
                                                    gsub('scaling_models/','',new.mod.computed))))
#     models.out@models.failed <- Biomod.material$calibration.failures
         
    #   3.3 Models evaluation conversion
    algo.choosen.id <- which( Biomod.material$algo.choice == TRUE)     
    algo.choosen.names <- Biomod.material$algo[algo.choosen.id]
    
    if(Biomod.material$NbRunEval>0){
      run.eval.names <- c(paste('RUN',1:Biomod.material$NbRunEval, sep=''),'Full')
    } else{
      run.eval.names <- 'Full'
    }

    if(Biomod.material$NbRepPA>0){
      pa.data.names <- paste('PA',1:(Biomod.material$NbRun[which( Biomod.material$species.names == sp.name)] / (Biomod.material$NbRunEval+1) ),sep='')
    } else{
      pa.data.names <- 'AllData'
    }
                        
    mod.eval.met <- c()
    if(!is.null(Evaluation.results.Roc)) mod.eval.met <- c(mod.eval.met, 'ROC')
    if(!is.null(Evaluation.results.TSS)) mod.eval.met <- c(mod.eval.met, 'TSS')
    if(!is.null(Evaluation.results.Kappa)) mod.eval.met <- c(mod.eval.met, 'KAPPA')
         
    models.evaluation <- array(data=NA,
                               dim=c(length(mod.eval.met), 4, length(algo.choosen.names),
                                     length(run.eval.names), length(pa.data.names)),
                               dimnames=list(mod.eval.met,
                                             c("Testing.data", "Cutoff", "Sensitivity", "Specificity"),
                                             algo.choosen.names,
                                             run.eval.names,
                                             pa.data.names))
    
#     outTmp <- lapply(pa.data.names, function(pdn){
    for( pdn in pa.data.names){
      pdnOld <- ifelse(pdn!='AllData',pdn,'full')
#       lapply(run.eval.names, function(ren){
      for(ren in run.eval.names){
        renOld <- ifelse(length(grep('RUN',ren)>0), gsub('RUN','rep',ren), '')
        runOldName <- ifelse(renOld!='', 
                             paste(models.out@sp.name, pdnOld, renOld, sep='_'), 
                             paste(models.out@sp.name, pdnOld, sep='_'))
        if(!is.null(Evaluation.results.Roc)){
          models.evaluation['ROC',,,ren,pdn] <- matrix(round(as.numeric(as.matrix(Evaluation.results.Roc[[runOldName]][,c('Cross.validation', 'Cutoff', "Sensitivity", "Specificity")])),digits=3), nrow=4, byrow=T)
        }

        if(!is.null(Evaluation.results.TSS)){
          models.evaluation['TSS',,,ren,pdn] <- matrix(round(as.numeric(as.matrix(Evaluation.results.TSS[[runOldName]][,c('Cross.validation', 'Cutoff', "Sensitivity", "Specificity")])),digits=3), nrow=4, byrow=T)
        } 

        if(!is.null(Evaluation.results.Kappa)){
          models.evaluation['KAPPA',,,ren,pdn] <- matrix(round(as.numeric(as.matrix(Evaluation.results.Kappa[[runOldName]][,c('Cross.validation', 'Cutoff', "Sensitivity", "Specificity")]),digits=3)), nrow=4, byrow=T)
        }

      }#)
    }#)
    
    # save model evaluation
    save(models.evaluation, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.evaluation",sep=""),  compress=compress.arg)
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.evaluation",sep="")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)


    #   3.4 Models variable Importances

    # save model variables importances
    if(!is.null(VarImportance)){
      variables.importances <- t(VarImportance[[sp.name]])
      save(variables.importances, file = paste(models.out@sp.name,"/.BIOMOD_DATA/variables.importances",sep=""),  compress=compress.arg)
      models.out@variables.importances@inMemory <- TRUE
      models.out@variables.importances@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/variables.importances",sep="")
      models.out@variables.importances@val <- variables.importances
      rm(variables.importances)
    }


    #   3.5 Models Predictions
    if(!exists(paste('Pred_', sp.name,sep=''))){
      load(paste(path,'pred/Pred_',sp.name,sep=""))
    }
    if(exists(paste('Pred_', sp.name,sep=''))){

      eval(parse(text=paste('models.prediction <- Pred_', sp.name,sep='')))
#       eval(parse(text=paste('rm(Pred_', sp.name,')',sep='')))

      kept.dim <- dim(models.prediction)
      kept.dim[2] <- length(algo.choosen.names)
      kept.dimnames <- dimnames(models.prediction)
#       kept.dimnames[[1]] <- NULL
      kept.dimnames[[2]] <- algo.choosen.names
      
      models.prediction <- models.prediction[,algo.choosen.names,,]
      dim(models.prediction) <- kept.dim
      dimnames(models.prediction) <- kept.dimnames
      
      if(length(grep('PA',dimnames(models.prediction)[[4]] )) == 0){
        dimnames(models.prediction)[[4]] <- 'AllData'
      } 
      if(dim(models.prediction)[3] > 1){
        vecTmp <- models.prediction[,,1,]
        models.prediction[,,1,] <- models.prediction[,,dim(models.prediction)[3],]
        models.prediction[,,dim(models.prediction)[3],] <- vecTmp
        rm(vecTmp)
        dimnames(models.prediction)[[3]] <- c(paste('RUN',1:(dim(models.prediction)[3] - 1),sep=''),'Full')
      } else{
        dimnames(models.prediction)[[3]] <- 'Full'
      }
        
      # save model predictions
      save(models.prediction, file = paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction",sep=""),  compress=compress.arg)
      models.out@models.prediction@inMemory <- FALSE
      models.out@models.prediction@link <- paste(models.out@sp.name,"/.BIOMOD_DATA/models.prediction",sep="")
      rm(models.prediction)
    }

    return(models.out)

  })
         
         

  
  # xx. Reconstruct original workspace -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   rm(list=ls())
#   load("WS_tmp.Rdata")
#   file.remove("WS_tmp.Rdata")

  names(NewModelObj) <- Biomod.material$species.names

  cat("\n\n-=-=-=- Done -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n")
  return(NewModelObj)
}
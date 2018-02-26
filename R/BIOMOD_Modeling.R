####################################################################################################
# BIOMOD_Modeling
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   Compute user's models selected and evaluate those.

# INPUT :
#   data <- a BIOMOD.formated.data object returned by BIOMOD_FormatingData() 
#   models <- vector of models names desired 
#   models.options <- a BIOMOD.models.options object returned by BIOMOD_ModelingOptions() 
#   NbRunEval <- Nb of Evaluation run 
#   DataSplit <- % of data used for models calibrations stuff
#   Yweights <- response points weights 
#   VarImport <- Nb of permutation done for variable Importances evaluation 
#   models.eval.meth <- vector of names of Models evaluation metrix 
#   SavePredictions <- keep or not array of predictions on hard disk (NOTE: Always TRUE for a posteriori EF)
#   KeepPredIndependent <- Not used
#   DoEnsembleForcasting <- Do or not EnsembleModeling.
#   SaveObj <- keep all results on hard disk or not (NOTE: strongly recommanded)

# OUTPUT : 
#   a BIOMOD.models.out (summary of all that were done) that will be given to others BIOMOD functions


# NOTE : 


####################################################################################################

BIOMOD_Modeling <- function( data, 
                               models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka'), 
                               models.options = NULL, 
                               NbRunEval=1, 
                               DataSplit=100, 
                               Yweights=NULL,
                               Prevalence=NULL,
                               VarImport=0, 
                               models.eval.meth = c('KAPPA','TSS','ROC'), 
                               SaveObj = TRUE,
                               rescal.all.models = FALSE,
                               do.full.models = TRUE,
                               modeling.id=as.character(format(Sys.time(), "%s")),
                               ...){

  # 0. loading required libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  .Models.dependencies(silent=TRUE, models.options=models.options )
  
  # 1. args checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- .Models.check.args(data, models, models.options, NbRunEval, DataSplit,
                             Yweights, VarImport, models.eval.meth, Prevalence, 
                             do.full.models, SaveObj,...)
  # updating Models arguments
  models <- args$models
  models.options <- args$models.options
  NbRunEval <- args$NbRunEval
  DataSplit <- args$DataSplit
  Yweights <- args$Yweights
  VarImport <- args$VarImport
  models.eval.meth <- args$models.eval.meth
  Prevalence <- args$Prevalence
  do.full.models <- args$do.full.models
  DataSplitTable <- args$DataSplitTable
  SaveObj <- args$SaveObj
  compress.arg = TRUE#ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
  
  rm(args)
  models.out <- new('BIOMOD.models.out',
                    sp.name = data@sp.name,
                    modeling.id = modeling.id,
                    expl.var.names = colnames(data@data.env.var),
                    has.evaluation.data = data@has.data.eval,
                    rescal.all.models = rescal.all.models)

#   #To keep track of Yweights state at origin (user's input)
#     if(NbRepPA!=0 && is.null(Yweights)) Yweights <- matrix(NA, nc=Biomod.material$NbSpecies, nr=nrow(DataBIOMOD))
   
  
  # 2. creating simulation directories =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #  
  # create the directories in which various objects will be stored (models, predictions and 
  # projections). Projections' directories are created in the Projection() function.
  .Models.prepare.workdir(data@sp.name, models.out@modeling.id)
  
  
  # 3. Saving Data and Model.option objects -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(SaveObj){
    # save Input Data
    save(data, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"formated.input.data"), compress = compress.arg)
    models.out@formated.input.data@inMemory <- FALSE
    models.out@formated.input.data@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"formated.input.data")
    # save Model Options
    save(models.options, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.options"), compress = compress.arg)
    models.out@models.options@inMemory <- FALSE
    models.out@models.options@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.options")

  }
  

  # 3. rearanging data and determining calib and eval data -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #) 
  mod.prep.dat <- .Models.prepare.data(data, NbRunEval, DataSplit, Yweights, Prevalence, do.full.models, DataSplitTable)
  rm(data)
  
  # keeping calibLines
  calib.lines <- mod.prep.dat[[1]]$calibLines
  if(length(mod.prep.dat) > 1){ ## stack calib lines matrix along 3rd dimention of an array
    for(pa in 2:length(mod.prep.dat)){
      calib.lines <- abind(calib.lines, mod.prep.dat[[pa]]$calibLines, along=3)
    }
#     ## update dimnames
#     dimnames(calib.lines) <- list(dimnames(calib.lines)[[1]], dimnames(calib.lines)[[2]], paste("PA", 1:length(mod.prep.dat) ))
  } # else { ## force calib.line object to be a 3D array
#     dim(calib.lines) <- c(dim(calib.lines),1)
#     ## update dimnames
#     dimnames(calib.lines) <- list(dimnames(calib.lines)[[1]], dimnames(calib.lines)[[2]], paste("PA", 1:length(mod.prep.dat) ))
#   }
  # save calib.lines
  save(calib.lines, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"calib.lines"), compress = compress.arg)
  models.out@calib.lines@inMemory <- FALSE
  models.out@calib.lines@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"calib.lines")
  rm(calib.lines)


  # 4. Print modelling summary in console -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  .Models.print.modeling.summary(mod.prep.dat, models)
  
  # 5. Running models -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  
  # loop on PA
  modeling.out <- lapply(mod.prep.dat,.Biomod.Models.loop,
                          modeling.id = models.out@modeling.id,   
                          Model = models,
                          Options = models.options,
                          VarImport = VarImport, 
                          mod.eval.method = models.eval.meth,
                          SavePred = SaveObj,
                          scal.models = rescal.all.models
                          )
  
  # put outputs in good format and save those

  models.out@models.computed <- .transform.outputs(modeling.out, out='models.run')
  models.out@models.failed <- .transform.outputs(modeling.out, out='calib.failure')
  
  if(SaveObj){
    # save model evaluation
    models.evaluation <- .transform.outputs(modeling.out, out='evaluation')
    save(models.evaluation, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.evaluation"), compress = compress.arg)
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.evaluation")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)
    
    # save model variables importances
    if(VarImport > 0 ){
      variables.importances <- .transform.outputs(modeling.out, out='var.import')

      ## trick to put appropriate dimnames
#       vi.dim.names <- dimnames(variables.importances)
#       vi.dim.names[[1]] <- models.out@expl.var.names
#       dimnames(variables.importances) <- vi.dim.names
#       rm('vi.dim.names')
    
      save(variables.importances, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"variables.importance"), compress = compress.arg)
      models.out@variables.importances@inMemory <- TRUE
      models.out@variables.importances@link <-file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"variables.importance")
      models.out@variables.importances@val <- variables.importances
      rm(variables.importances)
    }

    # save model predictions
    models.prediction <- .transform.outputs(modeling.out, out='prediction')
    save(models.prediction, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction"),  compress=compress.arg)
    models.out@models.prediction@inMemory <- FALSE
    models.out@models.prediction@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction)
    
    # save evaluation model predictions
    models.prediction.eval <- .transform.outputs(modeling.out, out='prediction.eval')
    save(models.prediction.eval, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction.eval"), compress = compress.arg)
    models.out@models.prediction.eval@inMemory <- FALSE
    models.out@models.prediction.eval@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction.eval")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction.eval)

  }
  
  # removing MAXENT.Phillips tmp dir
#   if('MAXENT.Phillips' %in% models){
#     .Delete.Maxent.WorkDir(species.name=models.out@sp.name, modeling.id=models.out@modeling.id)
#   }
  
  rm(modeling.out)
  
  # save model object on hard drive
  models.out@link <- file.path(models.out@sp.name, paste(models.out@sp.name, '.', models.out@modeling.id, '.models.out', sep=""))
  assign(x=paste(models.out@sp.name, '.', models.out@modeling.id, '.models.out', sep=""),
         value=models.out)
  save(list=paste(models.out@sp.name, '.', models.out@modeling.id, '.models.out', sep=""), 
       file=models.out@link)
  

  .bmCat("Done")
  return(models.out)
}

# -=-=-=- Several hidden functions -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


.Models.dependencies <- function(silent=TRUE, models.options = NULL){
  # Loading all required libraries
  cat('\n\nLoading required library...')
#   require(nnet, quietly=silent)
#   require(rpart, quietly=silent)
#   require(MASS, quietly=silent)
#   require(gbm, quietly=silent)
#   require(mda, quietly=silent)
#   require(randomForest, quietly=silent)
#   
#   if(!is.null(models.options)){
#     if(grepl('mgcv', models.options@GAM$algo)){
#       if("package:gam" %in% search() ) detach(package:gam)
#       require(mgcv, quietly=silent)
#     } else{
#       if("package:mgcv" %in% search() ) detach(package:mgcv)
#       require(gam, quietly=silent)
#     }
#   } else {
#     if('mgcv' %in% rownames(installed.packages())){
#       if("package:gam" %in% search() ) detach(package:gam)
#       require(mgcv, quietly=silent)
#     } else{
#       if("package:mgcv" %in% search() ) detach(package:mgcv)
#       require(gam, quietly=silent)
#     }    
#   }
# 
#   require(abind, quietly=silent)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.check.args <- function(data, models, models.options, NbRunEval, DataSplit,
                               Yweights, VarImport, models.eval.meth, Prevalence, do.full.models, SaveObj, ...){
  cat('\n\nChecking Models arguments...\n')
  add.args <- list(...)
  
  # data checking
  if( !(class( data ) %in% c("BIOMOD.formated.data",
                             "BIOMOD.formated.data.PA",
                             "BIOMOD.formated.data.indep",
                             "BIOMOD.formated.data.PA.indep") )){
    stop( "data argument mut be a 'BIOMOD.formated.data' (obtain by running Initial.State function) ")
  }
  
  # models checking
  if( !is.character( models ) ){
    stop("models argument must be a 'character' vector")
  }
  
  models <- unique(models)
  
  if(sum(models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka')) != length(models)){
    stop(paste(models[which( (models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka'))
                             == FALSE) ]," is not a availabe model !",sep=""))
  }
  
  categorial_var <- unlist(sapply(colnames(data@data.env.var), function(x){if(is.factor(data@data.env.var[,x])) return(x) else return(NULL)} ))
  
  if(length(categorial_var)){
    models.fact.unsuprort <- c("SRE", "MAXENT.Tsuruoka")
    models.swich.off <- intersect(models, models.fact.unsuprort)
    if(length(models.swich.off)){
      models <- setdiff(models, models.swich.off)
      cat(paste0("\n\t! ", paste(models.swich.off, collapse = ",", sep = " ")," were switch off because of categorical variables !"))
    }
    
  }
  
  # models.options checking ( peut etre permetre l'utilisation de liste de params )
  if( !is.null(models.options) && class(models.options) != "BIOMOD.Model.Options" ){
    stop("models.options argument must be a 'BIOMOD.Model.Options.object' (obtained by running ... ) ")
  }
  
  if( is.null(models.options)){
    warning("Models will run with 'defaults' parameters", immediate.=T)
    # create a default models.options object
    models.options <- BIOMOD_ModelingOptions() # MAXENT.Phillips = list( path_to_maxent.jar = getwd()) 
    
  }
  
  # MAXENT.Phillips specific checking
  if("MAXENT.Phillips" %in% models){
    if(!file.exists(file.path(models.options@MAXENT.Phillips$path_to_maxent.jar ,"maxent.jar")) ){
      models = models[-which(models=='MAXENT.Phillips')]
      warning("The maxent.jar file is missing. You need to download this file (http://www.cs.princeton.edu/~schapire/maxent) and put the maxent.jar file in your working directory -> MAXENT.Phillips was switched off")
    } else if(!.check.java.installed()){
      models = models[-which(models=='MAXENT.Phillips')]
    } else if(nrow(data@coord)==1){
     # no coordinates
      warning("You must give XY coordinates if you want to run MAXENT.Phillips -> MAXENT.Phillips was switched off")
      models = models[-which(models=='MAXENT.Phillips')]
    } 
  }
  
  ## Data split checks
  if(!is.null(add.args$DataSplitTable)){
    cat("\n! User defined data-split table was given -> NbRunEval, DataSplit and do.full.models argument will be ignored")
    if(!(length(dim(add.args$DataSplitTable)) %in% c(2,3) )) stop("DataSplitTable must be a matrix or a 3D array")
    if(dim(add.args$DataSplitTable)[1] != length(data@data.species) ) stop("DataSplitTable must have as many rows (dim1) than your species as data")
    NbRunEval <- dim(add.args$DataSplitTable)[2]
    DataSplit <- 50
    do.full.models <- FALSE
  }
  
  
  # NbRunEval checking (permetre un nb different par espece?)
  if( !is.numeric(NbRunEval) || NbRunEval <= 0 ){
    stop("NbRunEval argument mus be a non null positive 'numeric'")
  }
  
  # DataSplit checking
  if( !is.numeric(DataSplit) || DataSplit < 0 || DataSplit > 100 ){
    stop("DataSplit argument must be a 0-100 'numeric'")
  }
  
  if(DataSplit < 50){
    warning("You choose to allocate more data to evaluation than to calibration of your model 
            (DataSplit<50)\nMake sure you really wanted to do that. \n", immediate.=T)
  }
  
#   # EM weight checking
#   if(!is.null(EM.weight))
#     if(!any(EM.weight==c("Roc","TSS","Kappa")))
#       stop("The 'EM.weight' parameter must be one of the following: NULL, 'Roc', 'TSS' or 'Kappa'.\n")
#       
  # Check that the weight matrix was entered correctly  
  if(!is.null(Yweights)){
     if(!is.numeric(Yweights))
        stop("Yweights must be a numeric vector")
     if(length(Yweights) != length(data@data.species)) 
       stop("The number of 'Weight' does not match with the input calibration data. 
            Simulation cannot proceed.")
  }
  
  # Defining evaluation runs.
  if(NbRunEval <= 0){
      DataSplit <- 100
      if(!(class(data) %in% c("BIOMOD.formated.data.indep", "BIOMOD.formated.data.PA.indep") )){
        warning("The models will be evaluated on the calibration data only (NbRunEval=0 and no
                independent data) \n\t it could lead to over-optimistic predictive performances.\n",
                immediate.=T)
      }
  }
  if(DataSplit==100) NbRunEval <- 0
  
  # Models evaluation method checking
  models.eval.meth <- unique(models.eval.meth)
  
  if(sum(models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS','KAPPA','ACCURACY','BIAS',
                              'POD','PODFD','CSI','ETS','HK','ROC')) != length(models.eval.meth)){
    stop(paste(models.eval.meth[which( (models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS',
                                                                'KAPPA','ACCURACY','BIAS', 'POD',
                                                                'PODFD','CSI', 'ETS','HK','ROC')) 
                                       == FALSE) ]," is not a availabe models evaluation metric !",sep=""))
  }
  
  # Prevalence checking
  if(!is.null(Prevalence)){
    if(!is.numeric(Prevalence) | Prevalence>=1 | Prevalence <=0){
      stop("Prevalence must be a 0-1 numeric")
    } else {
      # update MAXENT.Phillips default prevalence
      if("MAXENT.Phillips" %in% models){
        cat("\n\t MAXENT.Phillips defaultprevalence option was updated to fit with modeling prevalence (i.e",Prevalence,")")
        models.options@MAXENT.Phillips$defaultprevalence = Prevalence
      }
    }
  }
  
  ##### TO BE CHANGE BUT PREVENT FROM BUGS LATTER
  # Force object saving parameter
  if(!SaveObj){
    cat("\n\t SaveObj param was automaticaly set to TRUE to prevent from bugs.")
    SaveObj <- TRUE
  }
  
#   cat('\nChecking done!\n')
  return(list(models = models,
              models.options = models.options,
              NbRunEval = NbRunEval,
              DataSplit = DataSplit,
              Yweights = Yweights,
              VarImport = VarImport,
              models.eval.meth = models.eval.meth,
              Prevalence = Prevalence,
              do.full.models = do.full.models,
              SaveObj = SaveObj,
              DataSplitTable=add.args$DataSplitTable))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.prepare.workdir <- function(sp.name, modeling.id){
  cat("\nCreating suitable Workdir...\n")
  dir.create(sp.name, showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(sp.name,".BIOMOD_DATA",modeling.id), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(sp.name, "models",modeling.id), showWarnings=FALSE, recursive=T)

#   if(sum(models.list %in% c('MARS', 'FDA', 'ANN')) > 0 ){
#     dir.create(paste(getwd(),"/",sp.name, "/models/scaling_models", sep=""), showWarnings=FALSE, recursive=T)
#   }
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.SampleMat <- function(data.sp, dataSplit, nbRun = 1, data.env = NULL){
  # return a matrix with nbRun columns of boolean (T: calib, F= eval)
  # data.sp is a 0,1 vector
  pres <- which(data.sp == 1)
  abs <- (1:length(data.sp))[-pres]
  
  nbPresEval <- round(length(pres) * dataSplit/100)
  nbAbsEval <- round(length(abs) * dataSplit/100)
  
  mat.out <- matrix(FALSE, 
                    nrow = length(data.sp),
                    ncol = nbRun)
  colnames(mat.out) <- paste('_RUN',1:nbRun, sep='')
  
  for (i in 1:ncol(mat.out)){
    ## force to sample at least one level of each factorial variable for calibration
    fact.cell.samp <- NULL
    if(!is.null(data.env)){
      fact.cell.samp <- sample.factor.levels(data.env)
      mat.out[fact.cell.samp, i] <- TRUE
    }
    mat.out[sample(setdiff(pres, fact.cell.samp), 
                   max(nbPresEval - length(fact.cell.samp), 0)), i] <- TRUE
    mat.out[sample(setdiff(abs, fact.cell.samp),
                   max(nbAbsEval - length(fact.cell.samp), 0)), i] <- TRUE
  }
  return(mat.out)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.Models.print.modeling.summary <- function( mod.prep.dat, models){
  cat("\n\n")
  .bmCat(paste(unlist(strsplit(mod.prep.dat[[1]]$name,'_'))[1], "Modeling Summary"))

  cat("\n",ncol(mod.prep.dat[[1]]$dataBM)-1, " environmental variables (", colnames(mod.prep.dat[[1]]$dataBM)[-1], ")")

  cat("\nNumber of evaluation repetitions :" , ncol(mod.prep.dat[[1]]$calibLines))

  cat("\nModels selected :", models, "\n")

  cat("\nTotal number of model runs :",ncol(mod.prep.dat[[1]]$calibLines) * length(models) * length(mod.prep.dat),"\n")
  
  .bmCat()
}
  
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  
.Models.check.EF.args <- function(models, models.eval.meth, models.options){
  # the models selected for doing EF
  if(models.options@EF$models.selected == 'all'){
    EF.algo <- models[which(models != 'SRE')] # remove SRE technique if it was selected (SRE cannot be used for ensemble forecast)
  } else {
    EF.algo <- models[models %in% models.options@EF$models.selected]
  }
  if(length(EF.algo)==0) stop('No models available selected for Ensemble forcastiong stuff')
  # the weight methods
  if(models.options@EF$weight.method == 'all'){
    EF.weight <- models.eval.meth
  } else {
    EF.weight <- models.eval.meth[models.eval.meth %in% models.options@EF$weight.method]
  }
  if(length(EF.weight)==0) stop('No weighting method available selected for Ensemble forcastiong stuff')
  return(list(EF.algo = EF.algo,
            EF.weight = EF.weight))
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.check.java.installed <- function(){
  
  if(.Platform$OS.type == "unix"){
    java.test <- try( expr = eval(system("command -v java", intern=TRUE )) , silent=TRUE)
  } else if(.Platform$OS.type == "windows"){
    java.test <- try( expr = eval(system( "java", intern=TRUE )), silent=TRUE )
  } else java.test <- ""
  
  if(!is.null(attr(java.test,"class"))){
    cat("\n! java software seems not be corectly installed\n  > MAXENT.Phillips modelling was switched off!")
    return(FALSE)
  } else{ return(TRUE) }
  
}


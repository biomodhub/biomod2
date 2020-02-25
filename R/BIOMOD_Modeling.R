##' @name BIOMOD_Modeling
##' @aliases BIOMOD_Modeling
##' @title Run a range of species distribution models
##' @description
##' This function allows to calibrate and evaluate a range of
##' species distribution models techniques run over a given
##' species. Calibrations are made on the whole sample or a
##' random subpart. The predictive power of the different models
##' is estimated using a range of evaluation metrics.
##' 
##' @param data \code{BIOMOD.formated.data} object returned by
##'   \code{\link[biomod2]{BIOMOD_FormatingData}}
##' @param models character, models to be computed names. To be 
##'   chosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE',
##'   'FDA', 'MARS', 'RF', 'MAXENT.Phillips', 'MAXENT.Phillips.2'
##' @param models.options \code{BIOMOD.models.options} object
##'   returned by \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' @param NbRunEval integer, number of Evaluation run.
##' @param DataSplit numeric, \% of data used to calibrate the
##'   models, the remaining part will be used for testing
##' @param Yweights numeric, vector of weights (one per 
##'   observation)
##' @param Prevalence either \code{NULL} (default) or a 0-1
##'   numeric used to build 'weighted response weights'
##' @param VarImport Number of permutation to estimate variable
##'   importance
##' @param models.eval.meth vector of names of evaluation metric
##'   among 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY',
##'   'BIAS', 'POD', 'CSI' and 'ETS'
##' @param SaveObj keep all results and outputs on hard drive or
##'   not (NOTE: strongly recommended)
##' @param rescal.all.models if true, all model prediction will
##'   be scaled with a binomial GLM
##' @param do.full.models if true, models calibrated and
##'   evaluated with the whole dataset are done
##' @param modeling.id character, the ID (=name) of modeling
##'   procedure. A random number by default.
##' @param \ldots further arguments :
##' 
##'  - \code{DataSplitTable} : a \code{matrix}, \code{data.frame}
##'    or a 3D \code{array} filled with \code{TRUE/FALSE} to
##'    specify which part of data must be used for models
##'    calibration (\code{TRUE}) and for models validation
##'    (\code{FALSE}). Each column corresponds to a 'RUN'. If
##'    filled, args \code{NbRunEval}, \code{DataSplit} and
##'    \code{do.full.models} will be ignored.
##'    
##' @details 
##' 
##' 1. \bold{data}
##' .. If you have decide to add pseudo absences to your
##' original dataset (see 
##' \code{\link[biomod2]{BIOMOD_FormatingData}}), 
##' NbPseudoAbsences * \code{NbRunEval + 1} models will be
##' created.
##' 
##' 2. \bold{models}
##' .. The set of models to be calibrated on the data. 10
##' modeling techniques are currently available:
##' 
##' .. - GLM : Generalized Linear Model 
##' (\code{\link[stats]{glm}})
##' 
##' .. - GAM : Generalized Additive Model (\code{\link[gam]{gam}},
##' \code{\link[mgcv]{gam}} or \code{\link[mgcv]{bam}}, see 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions} for details on 
##' algorithm selection})
##' 
##' .. - GBM : Generalized Boosting Model or usually called Boosted
##' Regression Trees (\code{\link[gbm]{gbm}})
##' 
##' .. - CTA: Classification Tree Analysis (\code{\link[rpart]{rpart}})
##' 
##' .. - ANN: Artificial Neural Network (\code{\link[nnet]{nnet}})
##' 
##' .. - SRE: Surface Range Envelop or usually called BIOCLIM
##' 
##' .. - FDA: Flexible Discriminant Analysis (\code{\link[mda]{fda}})
##' 
##' .. - MARS: Multiple Adaptive Regression Splines 
##' (\code{\link[earth]{earth}})
##' 
##' .. - RF: Random Forest (\code{\link[randomForest]{randomForest}})
##' 
##' .. - MAXENT.Phillips: Maximum Entropy (
##' \url{https://biodiversityinformatics.amnh.org/open_source/maxent})
##' 
##' .. - MAXENT.Phillips.2: Maximum Entropy 
##' (\code{\link[maxnet]{maxnet}})
##' 
##' 3. \bold{NbRunEval & DataSplit}
##' .. As already explained in the \code{\link{BIOMOD_FormatingData}}
##' help file, the common trend is to split the original dataset into 
##' two subsets, one to calibrate the models, and another one to evaluate
##' them. Here we provide the possibility to repeat this process
##' (calibration and evaluation) N times (\code{NbRunEval} times). 
##' The proportion of data kept for calibration is determined by the
##' \code{DataSplit} argument (100\% - \code{DataSplit} will be used to
##' evaluate the model). This sort of cross-validation allows to have a
##' quite robust test of the models when independent data are not
##' available. Each technique will also be calibrated on the complete
##' original data. All the models produced by BIOMOD and their related
##' informations are saved on the hard drive.
##' 
##' 4. \bold{Yweights & Prevalence}
##' .. Allows to give more or less weight to some particular 
##' observations. If these arguments is kept to NULL 
##' (\code{Yweights = NULL}, \code{Prevalence = NULL}), each 
##' observation (presence or absence) has the same weight (independent 
##' of the number of presences and absences). If \code{Prevalence = 0.5} 
##' absences will be weighted equally to the presences (i.e. the 
##' weighted sum of presence equals the weighted sum of absences). If
##' prevalence is set below or above 0.5 absences or presences are given
##' more weight, respectively.
##' .. In the particular case that pseudo-absence data have been
##' generated \code{BIOMOD_FormatingData} (\code{PA.nb.rep > 0}), weights
##' are by default (\code{Prevalence = NULL}) calculated such that
##' prevalence is 0.5, meaning that the presences will have the same
##' importance as the absences in the calibration process of the models.
##' Automatically created \code{Yweights} will be composed of integers to
##' prevent different modeling issues.
##' .. Note that the \code{Prevalence} argument will always be ignored if
##' \code{Yweights} are defined.
##' 
##' 5. \bold{models.eval.meth}
##' .. The available evaluations methods are :
##' 
##' .. - \code{ROC} : Relative Operating Characteristic
##' .. - \code{KAPPA} : Cohen's Kappa (Heidke skill score)
##' .. - \code{TSS} : True kill statistic (Hanssen and Kuipers 
##' discriminant, Peirce's skill score)
##' .. - \code{FAR} : False alarm ratio
##' .. - \code{SR} : Success ratio
##' .. - \code{ACCURANCY} : Accuracy (fraction correct)
##' .. - \code{BIAS} : Bias score (frequency bias)
##' .. - \code{POD} : Probability of detection (hit rate)
##' .. - \code{CSI} : Critical success index (threat score)
##' .. - \code{ETS} : Equitable threat score (Gilbert skill score)
##' 
##' Some of them are scaled to have all an optimum at 1. You can choose
##' one of more (vector) evaluation metric. By Default, only 'KAPPA',
##' 'TSS' and 'ROC' evaluation are done. Please refer to the CAWRC
##' website (\url{http://www.cawcr.gov.au/projects/verification/##'Methods_for_dichotomous_forecasts}) 
##' to get detailed description of each metric.
##' 
##' 6. \bold{SaveObj}
##' If this argument is set to False, it may prevent the evaluation of
##' the \sQuote{ensemble modeled} models in further steps. We strongly
##' recommend to always keep this argument \code{TRUE} even it asks for
##' free space onto the hard drive.
##'
##' 7. \bold{rescal.all.models}
##' \bold{This parameter is quite experimental and we advise not to use
##' it. It should lead to reduction in projection scale amplitude}
##' Some categorical models have to be scaled in every case (
##' \sQuote{FDA}, \sQuote{ANN}). But It may be interesting to scale all
##' model computed to ensure that they will produced comparable 
##' predictions (0-1000 ladder). That's particularly useful when you 
##' do some ensemble forecasting to remove the scale prediction effect
##' (the more extended projections are, the more they influence ensemble
##' forecasting results).
##' 
##' 8. \bold{do.full.models}
##' Building models with all information available may be useful in some
##' particular cases (i.e. rare species with few presences points). The
##' main drawback of this method is that, if you don't give separated
##' data for models evaluation, your models will be evaluated with the
##' same data that the ones used for calibration. That will lead to 
##' over-optimistic evaluation scores. Be careful with this '_Full'
##' models interpretation.
##' 
##' @return
##' A BIOMOD.models.out object
##' See \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} 
##' for details.
##' Additional objects are stored out of R in two different directories
##' for memory storage purposes. They are created by the function
##' directly on the root of your working directory set in R ("models"
##' directory). This one contains each calibrated model for each
##' repetition and pseudo-absence run. A hidden folder 
##' \code{.DATA_BIOMOD} contains some files (predictions, original
##' dataset copy, pseudo absences chosen...) used by other functions like
##' \code{\link[biomod2]{BIOMOD_Projection}} or 
##' \code{\link[biomod2]{BIOMOD_EnsembleModeling}}.
##' 
##' The models are currently stored as objects to be read exclusively in
##' R. To load them back (the same stands for all objects stored on the
##' hard disk) use the \code{\link{load}} function (see examples section below).
##' 
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' @seealso \code{\link{BIOMOD_FormatingData}},  
##'   \code{\link{BIOMOD_ModelingOptions}}, 
##'   \code{\link{BIOMOD_Projection}}
##' 
##' @keywords models
##' @keywords regression
##' @keywords nonlinear
##' @keywords multivariate
##' @keywords nonparametric
##' @keywords tree
##'   
##' @examples 
##' ##' species occurrences
##' DataSpecies <- 
##'   read.csv(
##'     system.file(
##'       "external/species/mammals_table.csv",
##'       package="biomod2"
##'     )
##'   )
##' head(DataSpecies)
##' 
##' ##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' ##' the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##' 
##' 
##' ##' Environmental variables extracted from BIOCLIM (bio_3, 
##' ##' bio_4, bio_7, bio_11 & bio_12)
##' myExpl <- 
##'   raster::stack(
##'     system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio12.grd", package = "biomod2")
##'   )
##'
##' ##' 1. Formatting Data
##' myBiomodData <- 
##'   BIOMOD_FormatingData(
##'     resp.var = myResp,
##'     expl.var = myExpl,
##'     resp.xy = myRespXY,
##'     resp.name = myRespName
##'   )
##' 
##' ##' 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' ##' 3. Doing Modelisation
##' myBiomodModelOut <- 
##'   BIOMOD_Modeling(
##'     myBiomodData,
##'     models = c('SRE','RF'),
##'     models.options = myBiomodOption,
##'     NbRunEval = 2,
##'     DataSplit = 80,
##'     VarImport = 0,
##'     models.eval.meth = c('TSS','ROC'),
##'     do.full.models = FALSE,
##'     modeling.id = "test"
##'   )
##' 
##' ##' print a summary of modeling stuff
##' myBiomodModelOut
##' 
BIOMOD_Modeling <- function(
  data,
  models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Phillips.2'),
  models.options = NULL,
  NbRunEval = 1,
  DataSplit = 100,
  Yweights = NULL,
  Prevalence = NULL,
  VarImport = 0,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = TRUE,
  modeling.id = as.character(format(Sys.time(), "%s")),
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
  # browser()
  mod.prep.dat <- 
    .Models.prepare.data(
      data, 
      NbRunEval, 
      DataSplit, 
      Yweights, 
      Prevalence, 
      do.full.models, 
      DataSplitTable
    )
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
  # browser()
  models.out@models.computed <- .transform.outputs.list(modeling.out, out='models.run')
  models.out@models.failed <- .transform.outputs.list(modeling.out, out='calib.failure')

  if(SaveObj){
    # save model evaluation
    models.evaluation <- .transform.outputs.list(modeling.out, out='evaluation')
    save(models.evaluation, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.evaluation"), compress = compress.arg)
    models.out@models.evaluation@inMemory <- TRUE
    models.out@models.evaluation@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.evaluation")
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)

    # save model variables importances
    if(VarImport > 0 ){
      variables.importances <- .transform.outputs.list(modeling.out, out='var.import')

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
    models.prediction <- .transform.outputs.list(modeling.out, out='prediction')
    save(models.prediction, file = file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction"),  compress=compress.arg)
    models.out@models.prediction@inMemory <- FALSE
    models.out@models.prediction@link <- file.path(models.out@sp.name,".BIOMOD_DATA",models.out@modeling.id,"models.prediction")
#     models.out@models.prediction@val <- .transform.outputs(modeling.out, out='prediction')
    rm(models.prediction)

    # save evaluation model predictions
    models.prediction.eval <- .transform.outputs.list(modeling.out, out='prediction.eval')
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
  if(
    !inherits(
      data, 
      c("BIOMOD.formated.data", "BIOMOD.formated.data.PA", "BIOMOD.formated.data.indep", 
        "BIOMOD.formated.data.PA.indep")
    )
  ){
    stop("data argument must be a 'BIOMOD.formated.data' (obtained by running Initial.State function) ")
  }

  # models checking
  if( !is.character( models ) )
  {
    stop("models argument must be a 'character' vector")
  }

  models <- unique(models)
  models.swich.off <- NULL

  avail.models.list <- 
    c(
      'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips', 
      'MAXENT.Phillips.2'
    )
  
  ## check if model is supported
  purrr::map(models, ~ checkmate::assert_choice(.x, avail.models.list))
  # if(sum(models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka')) != length(models)){
  #   stop(paste(models[which( (models %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips', 'MAXENT.Tsuruoka'))
  #                            == FALSE) ]," is not an availabe model !",sep=""))
  # }

  categorial_var <- unlist(sapply(colnames(data@data.env.var), function(x){if(is.factor(data@data.env.var[,x])) return(x) else return(NULL)} ))

  if(length(categorial_var)){
    models.fact.unsuprort <- c("SRE", "MAXENT.Tsuruoka")
    models.swich.off <- c(models.swich.off, intersect(models, models.fact.unsuprort))
    if(length(models.swich.off)){
      models <- setdiff(models, models.swich.off)
      cat(paste0("\n\t! ", paste(models.swich.off, collapse = ",", sep = " ")," were switched off because of categorical variables !"))
    }
  }
  
  ## disable MAXENT.Tsuruoka because of package maintaining issue (request from B Ripley 03-2019)
  if('MAXENT.Tsuruoka' %in% models){
    models.swich.off <- unique(c(models.swich.off, "MAXENT.Tsuruoka"))
    models <- setdiff(models, models.swich.off)
    warning('MAXENT.Tsuruoka has been disabled because of package maintaining issue (request from cran team 03-2019)')
  }

  # models.options checking ( peut etre permetre l'utilisation de liste de params )
  if(!is.null(models.options) & !inherits(models.options, "BIOMOD.Model.Options")){
    stop("models.options argument must be a 'BIOMOD.Model.Options.object' (obtained by running ... ) ")
  }

  if(is.null(models.options)){
    warning("Models will run with 'defaults' parameters", immediate.=T)
    # create a default models.options object
    models.options <- BIOMOD_ModelingOptions() # MAXENT.Phillips = list( path_to_maxent.jar = getwd())

  }

  # MAXENT.Phillips specific checking
  if("MAXENT.Phillips" %in% models){
    if(!file.exists(file.path(models.options@MAXENT.Phillips$path_to_maxent.jar ,"maxent.jar")) ){
      models = models[-which(models=='MAXENT.Phillips')]
      warning("The maxent.jar file is missing. You need to download this file (http://www.cs.princeton.edu/~schapire/maxent) and put the maxent.jar file in your working directory -> MAXENT.Phillips was switched off")
    ## -- 
    ## The java installation check is temporally disable cause it seems to cause 
    ## issues on some Windows users machine.
    ## --
    # } else if(!.check.java.installed()){
    #   models = models[-which(models=='MAXENT.Phillips')]
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
    warning("You chose to allocate more data to evaluation than to calibration of your model
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
      if(!inherits(data, c("BIOMOD.formated.data.indep", "BIOMOD.formated.data.PA.indep"))){
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
    cat("\n\t SaveObj param was automatically set to TRUE to prevent bugs.")
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
  if(length(EF.algo)==0) stop('No models available selected for Ensemble forecasting stuff')
  # the weight methods
  if(models.options@EF$weight.method == 'all'){
    EF.weight <- models.eval.meth
  } else {
    EF.weight <- models.eval.meth[models.eval.meth %in% models.options@EF$weight.method]
  }
  if(length(EF.weight)==0) stop('No weighting method available selected for Ensemble forecasting stuff')
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
    cat("\n! java software seems not be correctly installed\n  > MAXENT.Phillips modelling was switched off!")
    return(FALSE)
  } else{ return(TRUE) }

}

#####################################################################################################
#' Reshape biomod2 objects
#' 
#' This is an internal function (developper only)
#'
#' @param modOut the model object to transform given as a list
#' @param out character, the type of output to be transformed
#'
#' @return extracted statistics of interest from the model object
#'   as `array`.
#' @export
#'
.transform.outputs.array <-
  function(
    modOut, 
    out = 'evaluation'
  ){
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'var.import', 'calib.failure', 'models.run', 'prediction.eval' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'var.import',
                                                            'calib.failure', 'models.run', 'prediction.eval' ))))
    }
    
    # check dim of input list
    if(length(dim(modOut)) != 4 ){
      cat('\n',dim(modOut),'\n')
      print(dimnames(modOut))
      warning("Not computed .transform.outputs because of an incompatible input list dimension", immediate=T)
      return(NULL)
    }
    
    if(dim(modOut)[4] == 1 & length(unlist(strsplit(unlist(dimnames(modOut)[4]),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(length(dimnames(modOut)[[4]]) > 0){
        dataset.names <- unlist(sapply(unlist(dimnames(modOut)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else {
        dataset.names <- paste('PA', 1:dim(modOut)[4])
      }
    }
    
    run.eval.names <- sub('_','',unlist(dimnames(modOut)[3]))
    mod.names <- unlist(dimnames(modOut)[2])
    
    if (out=='evaluation'){
      if( is.null(modOut['evaluation',1,1,1])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut['evaluation',1,1,1]))
      eval.col.names <- colnames(as.data.frame(modOut['evaluation',1,1,1]))
      
      eval.out <- array(data = unlist(modOut['evaluation',,,]),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                        eval.col.names,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(eval.out)
    }
    
    if (out=='prediction'){
      if( is.null(modOut['pred',1,1,1])){ return(NULL) }
      nb.pts.pred <- length(as.numeric(unlist(modOut['pred',1,1,1])))
      pred.out <- array(data = unlist(modOut['pred',,,]),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(pred.out)
    }
    
    if (out=='prediction.eval'){
      if( is.null(modOut['pred.eval',1,1,1])){ return(NULL) }
      nb.pts.pred.eval <- length(as.numeric(unlist(modOut['pred.eval',1,1,1])))
      pred.eval.out <- array(data = unlist(modOut['pred.eval',,,]),
                             dim = c(nb.pts.pred.eval,
                                     length(mod.names),
                                     length(run.eval.names),
                                     length(dataset.names)),
                             dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
      return(pred.eval.out)
    }
    
    if (out=='var.import'){
      if( is.null(unlist(modOut['var.import',1,1,1]))){ return(NULL) }
      nb.var <- length(as.numeric(unlist(modOut['var.import',1,1,1])))
      
      vi.out <- array(data = unlist(modOut['var.import',,,]),
                      dim = c(nb.var,
                              length(mod.names),
                              length(run.eval.names),
                              length(dataset.names)),
                      dimnames = list(paste('Var',1:nb.var,sep=''), # to change
                                      mod.names,
                                      run.eval.names,
                                      dataset.names))
      
      return(vi.out)
    }
    
    if (out == 'calib.failure'){
      cf.out <- unlist(modOut['calib.failure',,,])
      return(cf.out[!is.null(cf.out)])
    }
    
    if (out == 'models.run'){
      mod.run.out <- unlist(modOut['ModelName',,,])
      return(mod.run.out[!is.null(mod.run.out)])
    }
    
  }

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #


#' Reshape biomod2 objects
#' 
#' This is an internal function (developper only)
#'
#' @param modOut the object to transform given as a list
#' @param out character, the type of input object
#' @param dim.names character, if not `NULL` the resshaped object will be stored on the hard drive
#'
#' @return
#' list, the extracted statistics
#' @export
.transform.outputs.list =
  function(
    modOut, 
    out = 'evaluation', 
    dim.names = NULL
  ){
    
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
                    'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'prediction.eval', 'var.import',
                                                            'calib.failure', 'models.run', 'EF.prediction',
                                                            'EF.PCA.median', 'EF.evaluation'))))
    }
    
    if(length(modOut) == 1 & length(unlist(strsplit(unlist(names(modOut)),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(is.null(dim.names)){
        dataset.names <- unlist(sapply(unlist(names(modOut)), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else{
        dataset.names <- unlist(dim.names[1])
      }
    }
    
    if(is.null(dim.names)){
      run.eval.names <- sub('_','',unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
      
      mod.names <- unlist(names(modOut[[1]][[1]]))
    } else{
      run.eval.names <- unlist(dim.names[2])
      mod.names <- unlist(dim.names[3])
    }
    
    if (out=='evaluation'){
      
      eval.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            eval.tab <- modOut[[i]][[j]][[k]][['evaluation']]
            if(!is.null(eval.tab)){ break }
          }
          if(!is.null(eval.tab)){ break }
        }
        if(!is.null(eval.tab)){ break }
      }
      
      if( is.null(eval.tab)){ return(NULL) }
      
      eval.meth.names <- rownames(as.data.frame(eval.tab))
      eval.col.names <- colnames(as.data.frame(eval.tab))
      
      eval.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            if(is.null(modOut[[d1]][[d2]][[d3]][['calib.failure']])){
              return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
            } else { matrix(NA, ncol=length(eval.col.names), nrow=length(eval.meth.names), dimnames=list(eval.meth.names,eval.col.names))}
          })
        })
      })
      
      eval.out <- array(data = unlist(eval.out),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                        eval.col.names,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(eval.out)
    }
    
    if (out=='prediction'){
      
      pred.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            pred.tab <- modOut[[i]][[j]][[k]][['pred']]
            if(!is.null(pred.tab)){ break }
          }
          if(!is.null(pred.tab)){ break }
        }
        if(!is.null(pred.tab)){ break }
      }
      
      if( is.null(pred.tab)){ return(NULL) }
      
      
      nb.pts.pred <- length(as.numeric(pred.tab))
      
      pred.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            if(is.null(modOut[[d1]][[d2]][[d3]][['pred']])){
              return(rep(NA,nb.pts.pred))
            } else{
              return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred']]))
            }
          })
        })
      })
      
      pred.out <- array(data = unlist(pred.out),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                        mod.names,
                                        run.eval.names,
                                        dataset.names))
      
      return(pred.out)
    }
    
    if (out=='prediction.eval'){
      pred.eval.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            pred.eval.tab <- modOut[[i]][[j]][[k]][['pred.eval']]
            if(!is.null(pred.eval.tab)){ break }
          }
          if(!is.null(pred.eval.tab)){ break }
        }
        if(!is.null(pred.eval.tab)){ break }
      }
      
      if( is.null(pred.eval.tab)){ return(NULL) }
      
      
      nb.pts.pred.eval <- length(as.numeric(pred.eval.tab))
      
      pred.eval.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            if(is.null(modOut[[d1]][[d2]][[d3]][['pred.eval']])){
              return(rep(NA,nb.pts.pred.eval))
            } else{
              return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred.eval']]))
            }
          })
        })
      })
      
      pred.eval.out <- array(data = unlist(pred.eval.out),
                             dim = c(nb.pts.pred.eval,
                                     length(mod.names),
                                     length(run.eval.names),
                                     length(dataset.names)),
                             dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
      return(pred.eval.out)
    }
    
    if (out=='var.import'){
      vi.tab <- NULL
      nb_pa <- length(modOut)
      nb_run <- length(modOut[[1]])
      nb_mod <- length(modOut[[1]][[1]])
      
      for(i in 1:nb_pa){
        for(j in 1:nb_run){
          for(k in 1:nb_mod){
            vi.tab <- modOut[[i]][[j]][[k]][['var.import']]
            if(!is.null(vi.tab)){ break }
          }
          if(!is.null(vi.tab)){ break }
        }
        if(!is.null(vi.tab)){ break }
      }
      
      if( is.null(vi.tab)){ return(NULL) }
      
      nb.var <- length(as.numeric(unlist(vi.tab)))
      
      ef.mod <- grep(pattern="EF.",mod.names) # EF models
      if(length(ef.mod)>0){
        kept.mod <- mod.names[-ef.mod]
      } else{
        kept.mod <- mod.names
      }
      
      vi.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(kept.mod, function(d3){ # models without EF ones
            if(is.null(modOut[[d1]][[d2]][[d3]][['var.import']])){
              return(rep(NA,nb.var))
            } else{
              return(as.matrix(modOut[[d1]][[d2]][[d3]][['var.import']]))
            }
          })
        })
      })
      
      vi.out <- array(data = unlist(vi.out),
                      dim = c(nb.var,
                              length(kept.mod),
                              length(run.eval.names),
                              length(dataset.names)),
                      dimnames = list(names(modOut[[1]][[1]][[1]][['var.import']]), # to change
                                      kept.mod,
                                      run.eval.names,
                                      dataset.names))
      
      return(vi.out)
    }
    
    if (out == 'calib.failure'){
      cf.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            return(modOut[[d1]][[d2]][[d3]][['calib.failure']])
          })
        })
      })
      cf.out <- unlist(cf.out)
      if(length(cf.out)) cf.out <- na.omit(cf.out)
      if(length(cf.out)) cf.out <- cf.out[!is.null(cf.out)]
      if(!length(cf.out)) cf.out <- 'none'
      return(cf.out)
    }
    
    if (out == 'models.run'){
      mod.run.out <- lapply(names(modOut),function(d1){ # data set
        lapply(names(modOut[[d1]]), function(d2){ # run eval
          lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
            return(as.character(modOut[[d1]][[d2]][[d3]][['ModelName']]))
          })
        })
      })
      mod.run.out <- unlist(mod.run.out)
      if(length(mod.run.out)) mod.run.out <- na.omit(mod.run.out)
      if(length(mod.run.out)) mod.run.out <- mod.run.out[!is.null(mod.run.out)]
      if(!length(mod.run.out)) mod.run.out <- 'none'
      return(mod.run.out)
    }
    
    
    if (out == 'EF.prediction'){
      if( is.null(modOut[[1]][[1]][[1]][['EM']])){ return(NULL) }
      
      nb.pts.ef.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['EM']])))
      
      ef.pred.out <- lapply(1:length(modOut),function(d1){ # data set
        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
            return(as.numeric(modOut[[d1]][[d2]][[d3]][['EM']]))
          })
        })
      })
      
      ef.pred.out <- array( data = unlist(ef.pred.out),
                            dim = c(nb.pts.ef.pred,
                                    length(modOut[[1]][[1]]),
                                    length(modOut[[1]]),
                                    length(modOut)),
                            dimnames = list(NULL,
                                            mod.names,
                                            run.eval.names,
                                            dataset.names))
      
      return(ef.pred.out)
    }
    
    if (out == 'EF.PCA.median'){
      if( is.null(modOut[[1]][[1]][[1]][['PCA.median']])){ return(NULL) }
      
      ef.pca.out <- lapply(1:length(modOut),function(d1){ # data set
        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
            return(as.character(modOut[[d1]][[d2]][[d3]][['PCA.median']]))
          })
        })
      })
      
      ef.pca.out <- array( data = unlist(ef.pca.out),
                           dim = c(1,
                                   length(modOut[[1]][[1]]),
                                   length(modOut[[1]]),
                                   length(modOut)),
                           dimnames = list(NULL,
                                           mod.names,
                                           run.eval.names,
                                           dataset.names))
      
      return(ef.pca.out)
    }
    
    if (out == 'EF.evaluation'){
      if( is.null(modOut[[1]][[1]][[1]][['EM.eval']])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      
      ef.eval.out <- lapply(1:length(modOut),function(d1){ # data set
        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
            return(data.frame(modOut[[d1]][[d2]][[d3]][['EM.eval']]))
          })
        })
      })
      
      ef.eval.out <- array(data = unlist(ef.eval.out),
                           dim = c(length(eval.meth.names),
                                   length(eval.col.names),
                                   length(modOut[[1]][[1]]),
                                   length(modOut[[1]]),
                                   length(modOut)),
                           dimnames = list(eval.meth.names,
                                           eval.col.names,
                                           mod.names,
                                           run.eval.names,
                                           dataset.names))
      
      return(ef.eval.out)
    }
    
  }

DF_to_ARRAY <- function(df){
  if(!is.data.frame(df) & !is.matrix(df)){
    if(is.list(df)){
      df.names <- names(df)
      df <- as.data.frame(df)
      names(df) <- df.names
    } else{
      stop("You have to give a data.frame")
    }
  }
  
  a <- sapply(strsplit(colnames(df), '_'), tail, n=3)
  b <- lapply(1:3, function(id) return(unique(a[id,])))
  array.dim.names <- c(list(character(0)),rev(b))
  #   array.dim.names <- c(list(c(NULL)),rev(apply(sapply(strsplit(colnames(df), '_'), tail, n=3),1,unique)))
  
  array.dim <- c(nrow(df),sapply(array.dim.names[-1],length))
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
  
  for(x in colnames(df)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    array.out[,dimTmp[1],dimTmp[2],dimTmp[3]] <- df[,x]
  }
  return(array.out)
}

LIST_to_ARRAY <- function(ll){
  test <- sapply(ll, is.array)
  if(!all(test)) stop("list elements should be arrays")
  test <- sapply(ll,dim)
  test <- apply(test,1,function(x){length(unique(x))==1})
  if(!all(test)) stop("list elements differ in dimension")
  
  formal.dim.names <- dimnames(ll[[1]])
  new.dim.names <- rev(apply(sapply(strsplit(names(ll), '_'), tail, n=3),1,unique))
  array.dim.names <- c(formal.dim.names,new.dim.names)
  array.dim <- sapply(array.dim.names,length)
  
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)
  
  for(x in names(ll)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    dimTmp <- paste( paste(rep(",",length(formal.dim.names)),collapse="") , paste("'",dimTmp,"'",sep="",collapse=","),collapse="")
    
    eval(parse(text=paste("array.out[",dimTmp,"] <-  ll[[x]]",sep="")))
  }
  return(array.out)
}



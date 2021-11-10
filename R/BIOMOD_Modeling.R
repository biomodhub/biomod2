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
##' 



BIOMOD_Modeling <- function(data,
                            models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips', 'MAXENT.Phillips.2'),
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
                            ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .Models.check.args(data, models, models.options, NbRunEval, DataSplit, Yweights
                             , VarImport, models.eval.meth, Prevalence, do.full.models, SaveObj,...)
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
  compress.arg = TRUE #ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz'))
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  models.out <- new('BIOMOD.models.out',
                    sp.name = data@sp.name,
                    modeling.id = modeling.id,
                    expl.var.names = colnames(data@data.env.var),
                    has.evaluation.data = data@has.data.eval,
                    rescal.all.models = rescal.all.models)
  
  ## 2. Create simulation directories -------------------------------------------------------------
  ## Various objects will be stored (models, predictions, projections)
  ## Projections directories are created in Projection() function
  .Models.prepare.workdir(data@sp.name, models.out@modeling.id)
  
  ## 3. Prepare internal function to save elements ------------------------------------------------
  name.BIOMOD_DATA = file.path(models.out@sp.name, ".BIOMOD_DATA", models.out@modeling.id)
  .Models.save.object <- function(objName, objValue)
  {
    save(objValue, file = file.path(name.BIOMOD_DATA, objName), compress = compress.arg)
    eval(parse(text = paste0("models.out@", objName, "@inMemory <- FALSE")))
    eval(parse(text = paste0("models.out@", objName, "@link <- file.path(name.BIOMOD_DATA, objName)")))
  }
  
  ## 3.1 Save input data and models options -----------------------------------
  if (SaveObj) {
    .Models.save.object("formated.input.data", data)
    .Models.save.object("models.options", models.options)
  }
  
  ## 3.2 Get and save calibration lines ---------------------------------------
  mod.prep.dat <- .Models.prepare.data(data, 
                                       NbRunEval, 
                                       DataSplit, 
                                       Yweights, 
                                       Prevalence, 
                                       do.full.models, 
                                       DataSplitTable)
  rm(data)
  
  calib.lines <- mod.prep.dat[[1]]$calibLines
  if (length(mod.prep.dat) > 1) { ## stack calib lines matrix along array 3rd-dimension
    for (pa in 2:length(mod.prep.dat)) {
      calib.lines <- abind(calib.lines, mod.prep.dat[[pa]]$calibLines, along = 3)
    }
  } 
  .Models.save.object("calib.lines", calib.lines)
  rm(calib.lines)
  
  ## 4. Print modeling summary in console ---------------------------------------------------------
  .Models.print.modeling.summary(mod.prep.dat, models)
  
  ## 5. Run models with loop over PA --------------------------------------------------------------
  modeling.out <- lapply(mod.prep.dat,
                         .Biomod.Models.loop,
                         modeling.id = models.out@modeling.id,
                         Model = models,
                         Options = models.options,
                         VarImport = VarImport,
                         mod.eval.method = models.eval.meth,
                         SavePred = SaveObj,
                         scal.models = rescal.all.models)
  
  ## 3.3 Rearrange and save outputs -------------------------------------------
  models.out@models.computed <- .transform.outputs.list(modeling.out, out = 'models.run')
  models.out@models.failed <- .transform.outputs.list(modeling.out, out = 'calib.failure')
  
  ## 3.4 Rearrange and save models outputs : ----------------------------------
  ## models evaluation, variables importance, models prediction, predictions evaluation
  if (SaveObj) {
    models.evaluation <- .transform.outputs.list(modeling.out, out = 'evaluation')
    .Models.save.object("models.evaluation", models.evaluation)
    models.out@models.evaluation@val <- models.evaluation
    rm(models.evaluation)
    
    if (VarImport > 0) {
      variables.importance <- .transform.outputs.list(modeling.out, out = 'var.import')
      .Models.save.object("variables.importance", variables.importance)
      models.out@variables.importance@val <- variables.importance
      rm(variables.importance)
    }
    
    models.prediction <- .transform.outputs.list(modeling.out, out = 'prediction')
    .Models.save.object("models.prediction", models.prediction)
    rm(models.prediction)
    
    models.prediction.eval <- .transform.outputs.list(modeling.out, out = 'prediction.eval')
    .Models.save.object("models.prediction.eval", models.prediction.eval)
    rm(models.prediction.eval)
  }
  rm(modeling.out)
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  name.OUT = paste0(models.out@sp.name, '.', models.out@modeling.id, '.models.out')
  models.out@link <- file.path(models.out@sp.name, name.OUT)
  assign(x = name.OUT, value = models.out)
  save(list = name.OUT, file = models.out@link)
  
  .bmCat("Done")
  return(models.out)
}



###################################################################################################

.Models.check.args <- function(data, models, models.options, NbRunEval, DataSplit, Yweights
                               , VarImport, models.eval.meth, Prevalence, do.full.models, SaveObj, ...)
{
  ## 0. Checking data and models arguments ------------------------------------
  cat('\n\nChecking Models arguments...\n')
  add.args <- list(...)
  
  .fun_testIfInherits(TRUE, "data", data, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  if (!is.character(models)) { stop("models must be a 'character' vector") }
  models <- unique(models)
  models.swich.off <- NULL
  
  ## check if model is supported
  avail.models.list <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS'
                         , 'RF', 'MAXENT.Phillips', 'MAXENT.Phillips.2')
  purrr::map(models, ~ checkmate::assert_choice(.x, avail.models.list))
  
  
  ## 1. Remove models not supporting categorical variables --------------------
  categorial_var <- unlist(sapply(colnames(data@data.env.var), function(x) {
    if (is.factor(data@data.env.var[, x])) { return(x) } else { return(NULL) }
  }))
  if (length(categorial_var))
  {
    models.fact.unsuprort <- c("SRE", "MAXENT.Tsuruoka")
    models.swich.off <- c(models.swich.off, intersect(models, models.fact.unsuprort))
    if (length(models.swich.off)) {
      models <- setdiff(models, models.swich.off)
      cat(paste0("\n\t! ", paste(models.swich.off, collapse = ",")," were switched off because of categorical variables !"))
    }
  }
  
  ## 2.1 Disable MAXENT.Tsuruoka ----------------------------------------------
  ## because of package maintaining issue (request from B Ripley 03-2019)
  if ('MAXENT.Tsuruoka' %in% models) {
    models.swich.off <- unique(c(models.swich.off, "MAXENT.Tsuruoka"))
    models <- setdiff(models, models.swich.off)
    warning('MAXENT.Tsuruoka has been disabled because of package maintaining issue (request from cran team 03-2019)')
  }
  
  ## 3. Check models.options arguments ----------------------------------------
  if (!is.null(models.options)) {
    .fun_testIfInherits(TRUE, "models.options", models.options, "BIOMOD.Model.Options")
  } else {
    warning("Models will run with 'defaults' parameters", immediate. = TRUE)
    models.options <- BIOMOD_ModelingOptions()
  }
  
  ## 2.2 Specific check for MAXENT.Phillips -----------------------------------
  if ("MAXENT.Phillips" %in% models)
  {
    if (!file.exists(file.path(models.options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"))) {
      models = models[-which(models == 'MAXENT.Phillips')]
      warning(paste0("MAXENT.Phillips has been disabled because the maxent.jar file is missing. "
                     , "`maxent.jar` file must be downloaded (https://biodiversityinformatics.amnh.org/open_source/maxent/) "
                     , "and put in the working directory."), immediate. = TRUE)
      ## -- 
      ## The java installation check is temporally disabled cause it seems to cause 
      ## issues on some Windows users machine.
      ## --
      # } else if(!.check.java.installed()){
      #   models = models[-which(models=='MAXENT.Phillips')]
    } else if (nrow(data@coord) == 1) {
      warning("MAXENT.Phillips has been disabled because no XY coordinates have been given", immediate. = TRUE)
      models = models[-which(models == 'MAXENT.Phillips')]
    }
  }
  
  ## 4. Check NbRunEval and DataSplitTable arguments --------------------------
  if (!is.null(add.args$DataSplitTable)) {
    cat("\n! User defined data-split table was given -> NbRunEval, DataSplit and do.full.models argument will be ignored")
    if (!(length(dim(add.args$DataSplitTable) %in% c(2, 3))) { stop("DataSplitTable must be a matrix or a 3D array") }
    if (dim(add.args$DataSplitTable)[1] != length(data@data.species)) { stop("DataSplitTable must have as many rows (dim1) than your species as data") }
    NbRunEval <- dim(add.args$DataSplitTable)[2]
    DataSplit <- 50
    do.full.models <- FALSE
  }
  
  .fun_testIfPosInt(TRUE, "NbRunEval", NbRunEval)
  if (DataSplit < 0 || DataSplit > 100) {
    stop("DataSplit argument must be a 0-100 'numeric'")
  } else if (DataSplit < 50) {
    warning("You chose to allocate more data to evaluation than to calibration of your model
            (DataSplit<50)\nMake sure you really wanted to do that. \n", immediate. = TRUE)
  } else if (DataSplit == 100) {
    NbRunEval <- 0
    warning(paste0("The models will be evaluated on the calibration data only "
                   , "(NbRunEval=0 and no independent data) \n\t "
                   , "It could lead to over-optimistic predictive performances.\n")
            , immediate. = TRUE)
  }
  
  ## 5. Check Yweights arguments ----------------------------------------------
  if (!is.null(Yweights)) {
    if (!is.numeric(Yweights)) { stop("Yweights must be a numeric vector") }
    if (length(Yweights) != length(data@data.species)) {
      stop("The number of 'Weight' does not match with the input calibration data. Simulation cannot proceed.")
    }
  }
  
  ## 6. Check models.eval.meth arguments --------------------------------------
  models.eval.meth <- unique(models.eval.meth)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS')
  .fun_testIfIn(TRUE, "models.eval.meth", models.eval.meth, avail.eval.meth.list)
    
  ## 7. Check Prevalence arguments --------------------------------------------
  if (!is.null(Prevalence)) {
    .fun_testIf01(TRUE, "Prevalence", Prevalence)
    if ("MAXENT.Phillips" %in% models) {
      cat("\n\t MAXENT.Phillips default prevalence option was updated to fit with modeling prevalence (i.e", Prevalence, ")")
      models.options@MAXENT.Phillips$defaultprevalence = Prevalence
    }
  }
  
  ##### TO BE CHANGE BUT PREVENT FROM BUGS LATER :  Force object saving parameter
  if (!SaveObj) {
    cat("\n\t SaveObj param was automatically set to TRUE to prevent bugs.")
    SaveObj <- TRUE
  }
  
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

###################################################################################################

.Models.prepare.workdir <- function(sp.name, modeling.id)
{
  cat("\nCreating suitable Workdir...\n")
  dir.create(sp.name, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(sp.name, ".BIOMOD_DATA", modeling.id), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(sp.name, "models", modeling.id), showWarnings = FALSE, recursive = TRUE)
}


###################################################################################################

.Models.print.modeling.summary <- function(mod.prep.dat, models)
{
  cat("\n\n")
  .bmCat(paste(unlist(strsplit(mod.prep.dat[[1]]$name, '_'))[1], "Modeling Summary"))
  cat("\n", ncol(mod.prep.dat[[1]]$dataBM) - 1, " environmental variables (", colnames(mod.prep.dat[[1]]$dataBM)[-1], ")")
  cat("\nNumber of evaluation repetitions :", ncol(mod.prep.dat[[1]]$calibLines))
  cat("\nModels selected :", models, "\n")
  cat("\nTotal number of model runs:", ncol(mod.prep.dat[[1]]$calibLines) * length(models) * length(mod.prep.dat), "\n")
  .bmCat()
}


##' ###############################################################################################
##' 
##' Reshape biomod2 objects
##' 
##' This is an internal function (developper only)
##'
##' @param modOut the object to transform given as a list
##' @param out character, the type of input object
##' @param dim.names character, if not `NULL` the resshaped object will be stored on the hard drive
##'
##' @return list, the extracted statistics
##' 
##' @export
##' 
##' ###############################################################################################


.transform.outputs.list = function(modOut, out = 'evaluation', dim.names = NULL)
{
  out_list = c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
               'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation')
  test = .fun_testIfIn(TRUE, "out", out, out_list)
  
  ## 0.a get dataset names ------------------------------------------------------------------------
  if (length(modOut) == 1 && length(unlist(strsplit(unlist(names(modOut)), '_'))) == 1) {
    dataset.names <- 'AllData'
  } else if (is.null(dim.names)) {
    dataset.names <- unlist(sapply(unlist(names(modOut)), function(name) {
      return(tail(unlist(strsplit(name, '_')), 1))
    }))
  } else {
    dataset.names <- unlist(dim.names[1])
  }
  
  ## 0.b get run.eval and model names -------------------------------------------------------------
  if (is.null(dim.names)) {
    run.eval.names <- sub('_', '', unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
    mod.names <- unlist(names(modOut[[1]][[1]]))
  } else {
    run.eval.names <- unlist(dim.names[2])
    mod.names <- unlist(dim.names[3])
  }
  
  ## 1. CASE evaluation / prediction / prediction.eval / var.import -------------------------------
  
  if (out %in% c("evaluation", "prediction", "prediction.eval", "var.import"))
  {
    nb_pa <- length(modOut)
    nb_run <- length(modOut[[1]])
    nb_mod <- length(modOut[[1]][[1]])
    
    name_slot = out
    if (out == "prediction") { name_slot = "pred" }
    if (out == "prediction.eval") { name_slot = "pred.eval" }
    
    output <- NULL
    for (i in 1:nb_pa) {
      for (j in 1:nb_run) {
        for (k in 1:nb_mod) {
          output <- modOut[[i]][[j]][[k]][[name_slot]]
          if (!is.null(output)) { break }
        }
        if (!is.null(output)) { break }
      }
      if (!is.null(output)) { break }
    }
    if (is.null(output)) { return(NULL) }
    
    if (out == "evaluation") {
      eval.meth.names <- rownames(as.data.frame(output))
      eval.col.names <- colnames(as.data.frame(output))
      dimnames.out = list(eval.meth.names,
                          eval.col.names,
                          mod.names,
                          run.eval.names,
                          dataset.names)
      dim.out = c(length(eval.meth.names),
                  length(eval.col.names),
                  length(mod.names),
                  length(run.eval.names),
                  length(dataset.names))
    } else if (out %in% c("prediction", "prediction.eval", "var.import")) {
      kept.mod = mod.names
      if (out == "var.import") {
        ef.mod <- grep(pattern = "EF.", mod.names) # EF models
        if (length(ef.mod) > 0) {
          kept.mod <- mod.names[-ef.mod]
        }
      }
      
      nb.tmp <- length(as.numeric(output))
      dimnames.out = list(NULL, kept.mod, run.eval.names, dataset.names)
      if (out == "var.import") {
        dimnames.out = list(names(modOut[[1]][[1]][[1]][['var.import']]), # to change
                            , kept.mod, run.eval.names, dataset.names)
      }
      dim.out = c(nb.tmp,
                  length(kept.mod),
                  length(run.eval.names),
                  length(dataset.names))
    }
    
    output <- lapply(names(modOut), function(d1) { # data set
      lapply(names(modOut[[d1]]), function(d2) { # run eval
        lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
          if (out == "evaluation") {
            if (is.null(modOut[[d1]][[d2]][[d3]][['calib.failure']])) {
              return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
            } else { 
              return(matrix(NA, ncol = length(eval.col.names), nrow = length(eval.meth.names)
                            , dimnames = list(eval.meth.names, eval.col.names)))
            }
          } else if (out == "prediction" || out == "prediction.eval" || 
                     (out == "var.import" && d3 %in% kept.mod)) {
            if (is.null(modOut[[d1]][[d2]][[d3]][[name_slot]])) {
              return(rep(NA, nb.pts.pred))
            } else {
              return(as.numeric(modOut[[d1]][[d2]][[d3]][[name_slot]]))
            }
          }
        })
      })
    })
    data.out = unlist(output)
    
    res.out <- array(data = data.out, dim = dim.out, dimnames = dimnames.out)
    return(res.out)
  }
  
  ## 2. CASE calib.failure / models.run -----------------------------------------------------------
  
  if (out %in% c("calib.failure", "models.run"))
  {
    name_slot = out
    if (out == "models.run") { name_slot = "ModelName" }
    
    output <- lapply(names(modOut),function(d1) { # data set
      lapply(names(modOut[[d1]]), function(d2) { # run eval
        lapply(names(modOut[[d1]][[d2]]), function(d3) { # models
          res = modOut[[d1]][[d2]][[d3]][[name_slot]]
          if (out == "calib.failure") {
            return(as.numeric(res))
          } else if (out == "models.run") {
            return(as.character(res))
          }
        })
      })
    })
    
    res.out <- unlist(output)
    if (length(res.out)) { res.out <- na.omit(res.out) }
    if (length(res.out)) { res.out <- res.out[!is.null(res.out)] }
    if (!length(res.out)) { res.out <- 'none' }
    return(res.out)
  }

  ## 3. CASE EF.prediction / EF.PCA.median / EF.evaluation ----------------------------------------
  
  if (out %in% c("EF.prediction", "EF.PCA.median", "EF.evaluation"))
  {
    name_slot = ifelse(out == "EF.prediction", "EM"
                       , ifelse(out == "EF.PCA.median", "PCA.median", "EM.eval"))
    
    if (is.null(modOut[[1]][[1]][[1]][[name_slot]])) { return(NULL) }
    
    nb.tmp = 1
    dimnames.out = list(NULL, mod.names, run.eval.names, dataset.names)
    if (out == "EF.prediction") {
      nb.tmp <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][[name_slot]])))
    } else if (out == "EF.evaluation") {
      eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][[name_slot]]))
      eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][[name_slot]]))
      nb.tmp = c(length(eval.meth.names), length(eval.col.names))
      dimnames.out = list(eval.meth.names, eval.col.names
                          , mod.names, run.eval.names, dataset.names)
    }
    dim.out = c(nb.tmp,
                length(modOut[[1]][[1]]),
                length(modOut[[1]]),
                length(modOut))

    
    output <- lapply(1:length(modOut),function(d1) { # data set
      lapply(1:length(modOut[[d1]]), function(d2) { # run eval
        lapply(1:length(modOut[[d1]][[d2]]), function(d3) { # models
          res = modOut[[d1]][[d2]][[d3]][[name_slot]]
          if (out == "EF.prediction") {
            return(as.numeric(res))
          } else if (out == "EF.PCA.median") {
            return(as.character(res))
          } else if (out == "EF.evaluation") {
            return(data.frame(res))
          }
        })
      })
    })
    data.out = unlist(output)
    
    res.out <- array(data = data.out, dim = dim.out, dimnames = dimnames.out)
    return(res.out)
  }
}

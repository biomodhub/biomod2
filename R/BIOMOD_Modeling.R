##' ###############################################################################################
##' @name BIOMOD_Modeling
##' @aliases BIOMOD_Modeling
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' 
##' @title Run a range of species distribution models
##' 
##' @description This function allows to calibrate and evaluate a range of modeling techniques 
##' for a given species distribution. The dataset can be split up for independent calibration and 
##' validation, and the predictive power of the different models can be estimated using a range 
##' of evaluation metrics (see Details).
##' 
##' @param data a \code{\link{BIOMOD.formated.data} object returned by the 
##' \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}
##' @param models.options a \code{\link{BIOMOD.models.options} object returned by the 
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param NbRunEval an \code{integer} corresponding to the number of repetitions to be done for 
##' calibration/validation splitting (\emph{if specified, \code{DataSplit} and 
##' \code{do.full.models} will be ignored})
##' @param DataSplit a \code{numeric} between \code{0} and \code{1} corresponding to the 
##' percentage of data used to calibrate the models (calibration/validation splitting)
##' @param DataSplitTable (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix} or \code{data.frame} defining for each repetition (in columns) which 
##' observation lines should be used for models calibration (\code{TRUE}) and validation 
##' (\code{FALSE}) (see \code{\link{BIOMOD_CrossValidation}})
##' @param do.full.models (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether models calibrated and evaluated over the whole 
##' dataset must be computed or not
##' @param Yweights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to observation weights (one per 
##' observation, see Details)
##' @param Prevalence (\emph{optional, default} \code{NULL}) \cr 
##' A \code{numeric} between \code{0} and \code{1} corresponding to the species prevalence to 
##' build 'weighted response weights' (see Details)
##' @param VarImport (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param models.eval.meth a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param SaveObj (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether all results and outputs must be saved on hard drive 
##' or not (\emph{! strongly recommended !})
##' @param rescal.all.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions must be scaled with a binomial 
##' GLM or not
##' 
##' 
##' @return
##' 
##' A \code{BIOMOD.models.out} object containing models outputs, or links to saved outputs.
##' Models outputs are stored out of \R (for memory storage reasons) in 2 different folders 
##' created in the current working directory :
##' \enumerate{
##'   \item a \emph{models} folder, named after the \code{resp.name} argument of 
##'   \code{\link{BIOMOD_FormatingData}}, and containing all calibrated models for each 
##'   repetition and pseudo-absence run
##'   \item a \emph{hidden} folder, named \code{.BIOMOD_DATA}, and containing outputs related 
##'   files (original dataset, calibration lines, pseudo-absences selected, predictions, 
##'   variables importance, evaluation values...), that can be retrieved with 
##'   \code{\href{}{get_[...]}} or \code{\link{load}} functions, and used by other 
##'   \pkg{biomod2} functions, like \code{\link{BIOMOD_Projection}} or 
##'   \code{\link{BIOMOD_EnsembleModeling}}
##' }
##' 
##' 
##' @details 
##' 
##' 
##' 1. \bold{data}
##' .. If you have decide to add pseudo absences to your original dataset (see 
##' \code{\link{BIOMOD_FormatingData}}), \code{PA.nb.rep *(NbRunEval + 1)} models will be created.
##' 
##' 
##' 2. \bold{models}
##' .. The set of models to be calibrated on the data. 
##' 10 modeling techniques are currently available :
##' 
##' .. - \code{GLM} : Generalized Linear Model (\code{\link[stats]{glm}})
##' 
##' .. - \code{GAM} : Generalized Additive Model (\code{\link[gam]{gam}}, \code{\link[mgcv]{gam}} 
##' or \code{\link[mgcv]{bam}} (see \code{\link{BIOMOD_ModelingOptions} for details on algorithm 
##' selection})
##' 
##' .. - \code{GBM} : Generalized Boosting Model, or usually called Boosted Regression Trees 
##' (\code{\link[gbm]{gbm}})
##' 
##' .. - \code{CTA} : Classification Tree Analysis (\code{\link[rpart]{rpart}})
##' 
##' .. - \code{ANN} : Artificial Neural Network (\code{\link[nnet]{nnet}})
##' 
##' .. - \code{SRE} : Surface Range Envelop or usually called BIOCLIM
##' 
##' .. - \code{FDA} : Flexible Discriminant Analysis (\code{\link[mda]{fda}})
##' 
##' .. - \code{MARS} : Multiple Adaptive Regression Splines (\code{\link[earth]{earth}})
##' 
##' .. - \code{RF} : Random Forest (\code{\link[randomForest]{randomForest}})
##' 
##' .. - \code{MAXENT.Phillips} : Maximum Entropy 
##' (\url{https://biodiversityinformatics.amnh.org/open_source/maxent})
##' 
##' .. - \code{MAXENT.Phillips.2} : Maximum Entropy (\code{\link[maxnet]{maxnet}})
##' 
##' 
##' 3. \bold{NbRunEval & DataSplit}
##' .. Most simple method in machine learning to calibrate and evaluate a model is to split the 
##' original dataset in two, one to calibrate the model and the other one to evaluate it. The 
##' \code{DataSplit} argument defines the percentage of data that will be randomly selected and 
##' used for the calibration part, the remaining data constituting the evaluation part. This 
##' process is repeated \code{NbRunEval} times, to be sure not to include bias both in the 
##' modeling and evaluation parts.
##' .. Other validation methods are also available to the user :
##' \itemize{
##'   \item evaluation dataset can be directly given to the \code{\link{BIOMOD_FormatingData}} 
##'   function
##'   \item \code{DataSplitTable} argument can be used and obtained from the 
##'   \code{\link{BIOMOD_CrossValidation}} function
##' }
##' 
##' 
##' 4. \bold{Yweights & Prevalence}
##' .. More or less weight can be given to some specific observations.
##' .. If \code{Yweights = Prevalence = NULL}, each observation (presence or absence) will have the 
##' same weight, no matter the total number of presences and absences. If \code{Prevalence = 0.5}, 
##' presences and absences will be weighted equally (i.e. the weighted sum of presences equals the 
##' weighted sum of absences). If \code{Prevalence} is set below (\emph{above}) \code{0.5}, more 
##' weight will be given to absences (\emph{presences}).
##' .. If \code{Yweights} is defined, \code{Prevalence} argument will be ignored, and each 
##' observation will have its own weight.
##' .. If pseudo-absences have been generated (\code{PA.nb.rep > 0} in 
##' \code{\link{BIOMOD_FormatingData}}), weights are by default calculated such that 
##' \code{Prevalence = 0.5}. \emph{Automatically created \code{Yweights} will be \code{integer} 
##' values to prevent some modeling issues.}
##' 
##' 
##' 5. \bold{models.eval.meth}
##' .. The available evaluations methods are :
##' 
##' .. - \code{ROC} : Relative Operating Characteristic
##' .. - \code{KAPPA} : Cohen's Kappa (Heidke skill score)
##' .. - \code{TSS} : True kill statistic (Hanssen and Kuipers discriminant, Peirce's skill score)
##' .. - \code{FAR} : False alarm ratio
##' .. - \code{SR} : Success ratio
##' .. - \code{ACCURANCY} : Accuracy (fraction correct)
##' .. - \code{BIAS} : Bias score (frequency bias)
##' .. - \code{POD} : Probability of detection (hit rate)
##' .. - \code{CSI} : Critical success index (threat score)
##' .. - \code{ETS} : Equitable threat score (Gilbert skill score)
##' 
##' .. Optimal value of each method can be obtained with the \code{\link{get_optim_value}} 
##' function. Several evaluation metrics can be selected. Please refer to the 
##' \href{http://www.cawcr.gov.au/projects/verification/##'Methods_for_dichotomous_forecasts}{CAWRC website} 
##' to get detailed description of each metric.
##' 
##' 
##' 6. \bold{SaveObj}
##' .. \emph{If this argument is set to \code{FALSE}, it may prevent the evaluation of the 
##' ensemble models (see \code{\link{BIOMOD_EnsembleModeling}}) in further steps. Strong 
##' recommandation is to keep \code{SaveObj = TRUE}, even if it requires to have some free 
##' space onto the hard drive.}
##'
##'
##' 7. \bold{rescal.all.models}
##' .. \bold{This parameter is quite experimental and it is recommended not to use it. It may  
##' lead to reduction in projection scale amplitude.} Some categorical models always have to be 
##' scaled (\code{FDA}, \code{ANN}), but it may be interesting to scale all computed models to 
##' ensure comparable predictions (\code{0-1000} range). It might be particularly useful when 
##' doing ensemble forecasting to remove the scale prediction effect (the more extended 
##' projections are, the more they influence ensemble forecasting results).
##' 
##' 
##' 8. \bold{do.full.models}
##' .. Building models with all available information may be useful in some particular cases 
##' (e.g. rare species with few presences points). But calibration and evaluation datasets will 
##' be the same, which might lead to over-optimistic evaluation scores.
##' 
##' 
##' @keywords models, regression, nonlinear, multivariate, nonparametric, tree
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_ModelingOptions}}, 
##' \code{\link{BIOMOD_CrossValidation}}, \code{ \link{bm_VariablesImportance}}, 
##' \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleModeling}}
##' 
##'   
##' @examples
##' 
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv", package="biomod2"), row.names = 1)
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0("external/bioclim/current/bio", c(3, 4, 7, 11, 12), ".grd")
##' myExpl = raster::stack(system.file(myFiles[1], package = "biomod2"),
##'                        system.file(myFiles[2], package = "biomod2"),
##'                        system.file(myFiles[3], package = "biomod2"),
##'                        system.file(myFiles[4], package = "biomod2"),
##'                        system.file(myFiles[5], package = "biomod2"))
##' 
##' # 1. Formating Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' # 3. Doing Modelisation
##' myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
##'                                     models = c('SRE','RF'),
##'                                     models.options = myBiomodOption,
##'                                     NbRunEval = 2,
##'                                     DataSplit = 80,
##'                                     VarImport = 0,
##'                                     models.eval.meth = c('TSS','ROC'),
##'                                     do.full.models = FALSE,
##'                                     modeling.id = "test")
##' 
##' # print a summary of modeling stuff
##' myBiomodModelOut
##' 
##' 
##' ###############################################################################################


BIOMOD_Modeling <- function(data,
                            modeling.id = as.character(format(Sys.time(), "%s")),
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
    eval(parse(text = paste0("models.out@", objName, "@inMemory <<- FALSE")))
    eval(parse(text = paste0("models.out@", objName, "@link <<- file.path(name.BIOMOD_DATA, objName)")))
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
    if (!(length(dim(add.args$DataSplitTable) %in% c(2, 3)))) { stop("DataSplitTable must be a matrix or a 3D array") }
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
  for (i.meth in models.eval.meth) {
    .fun_testIfIn(TRUE, "models.eval.meth", i.meth, avail.eval.meth.list)
  }
  
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


## ###############################################################################################
## 
## Reshape biomod2 objects
## 
## This is an internal function (developper only)
##
## @param modOut the object to transform given as a list
## @param out character, the type of input object
## @param dim.names character, if not `NULL` the resshaped object will be stored on the hard drive
##
## @return list, the extracted statistics
## 
## @export
## 
## ###############################################################################################


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
        dimnames.out = list(names(modOut[[1]][[1]][[1]][['var.import']]) # to change ?
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

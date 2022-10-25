###################################################################################################
##' @name BIOMOD_Modeling
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' 
##' @title Run a range of species distribution models
##' 
##' @description This function allows to calibrate and evaluate a range of modeling techniques 
##' for a given species distribution. The dataset can be split up for independent calibration and 
##' validation, and the predictive power of the different models can be estimated using a range 
##' of evaluation metrics (see Details).
##' 
##' 
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}
##' @param bm.options a \code{\link{BIOMOD.models.options}} object returned by the  
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param nb.rep an \code{integer} corresponding to the number of repetitions to be done for 
##' calibration/validation splitting 
##' @param data.split.perc a \code{numeric} between \code{0} and \code{100} corresponding to the 
##' percentage of data used to calibrate the models (calibration/validation splitting)
##' @param data.split.table (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix} or \code{data.frame} defining for each repetition (in columns) which 
##' observation lines should be used for models calibration (\code{TRUE}) and validation 
##' (\code{FALSE}) (see \code{\link{BIOMOD_CrossValidation}}) \cr (\emph{if specified, 
##' \code{nb.rep}, \code{data.split.perc} and \code{do.full.models} will be ignored})
##' @param do.full.models (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether models calibrated and evaluated over the whole 
##' dataset should be computed or not
##' @param weights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to observation weights (one per 
##' observation, see Details)
##' @param prevalence (\emph{optional, default} \code{NULL}) \cr 
##' A \code{numeric} between \code{0} and \code{1} corresponding to the species prevalence to 
##' build '\emph{weighted response weights}' (see Details)
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param save.output (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether all outputs should be saved on hard drive or not 
##' (\emph{! strongly recommended !})
##' @param scale.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions should be scaled with a 
##' binomial GLM or not
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' @param do.progress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the progress bar is to be rendered or not
##' 
##' 
##' @return
##' 
##' A \code{BIOMOD.models.out} object containing models outputs, or links to saved outputs. \cr
##' Models outputs are stored out of \R (for memory storage reasons) in 2 different folders 
##' created in the current working directory :
##' \enumerate{
##'   \item a \emph{models} folder, named after the \code{resp.name} argument of 
##'   \code{\link{BIOMOD_FormatingData}}, and containing all calibrated models for each 
##'   repetition and pseudo-absence run
##'   \item a \emph{hidden} folder, named \code{.BIOMOD_DATA}, and containing outputs related 
##'   files (original dataset, calibration lines, pseudo-absences selected, predictions, 
##'   variables importance, evaluation values...), that can be retrieved with 
##'   \href{https://biomodhub.github.io/biomod2/reference/getters.out.html}{\code{get_[...]}} or \code{\link{load}} functions, and used by other 
##'   \pkg{biomod2} functions, like \code{\link{BIOMOD_Projection}} or 
##'   \code{\link{BIOMOD_EnsembleModeling}}
##' }
##' 
##' 
##' @details 
##' 
##' \describe{
##'   \item{bm.format}{If you have decided to add pseudo absences to your original dataset (see 
##'   \code{\link{BIOMOD_FormatingData}}), \cr \code{PA.nb.rep *(nb.rep + 1)} models will be 
##'   created.}
##'   
##'   \item{models}{The set of models to be calibrated on the data. 10 modeling techniques 
##'   are currently available :
##'   \itemize{
##'     \item \code{GLM} : Generalized Linear Model (\code{\link[stats]{glm}})
##'     \item \code{GAM} : Generalized Additive Model (\code{\link[gam]{gam}}, \code{\link[mgcv]{gam}} 
##'     or \code{\link[mgcv]{bam}}) \cr 
##'     (see \code{\link{BIOMOD_ModelingOptions} for details on algorithm selection})
##'     \item \code{GBM} : Generalized Boosting Model, or usually called Boosted Regression Trees 
##'     (\code{\link[gbm]{gbm}})
##'     \item \code{CTA} : Classification Tree Analysis (\code{\link[rpart]{rpart}})
##'     \item \code{ANN} : Artificial Neural Network (\code{\link[nnet]{nnet}})
##'     \item \code{SRE} : Surface Range Envelop or usually called BIOCLIM
##'     \item \code{FDA} : Flexible Discriminant Analysis (\code{\link[mda]{fda}})
##'     \item \code{MARS} : Multiple Adaptive Regression Splines (\code{\link[earth]{earth}})
##'     \item \code{RF} : Random Forest (\code{\link[randomForest]{randomForest}})
##'     \item \code{MAXENT.Phillips} : Maximum Entropy 
##'     (\url{https://biodiversityinformatics.amnh.org/open_source/maxent/})
##'     \item \code{MAXENT.Phillips.2} : Maximum Entropy (\code{\link[maxnet]{maxnet}})
##'   }}
##'   
##'   \item{nb.rep & data.split.perc}{
##'   \itemize{
##'     \item Most simple method in machine learning to calibrate and evaluate a model is to 
##'     split the original dataset in two, one to calibrate the model and the other one to 
##'     evaluate it. The \code{data.split.perc} argument defines the percentage of data that will be 
##'     randomly selected and used for the calibration part, the remaining data constituting the 
##'     evaluation part. This process is repeated \code{nb.rep} times, to be sure not to 
##'     include bias both in the modeling and evaluation parts.
##'     \item Other validation methods are also available to the user :
##'     \itemize{
##'       \item evaluation dataset can be directly given to the 
##'       \code{\link{BIOMOD_FormatingData}} function
##'       \item \code{data.split.table} argument can be used and obtained from the 
##'       \code{\link{BIOMOD_CrossValidation}} function
##'     }
##'   }}
##'   
##'   \item{weights & prevalence}{More or less weight can be given to some specific observations.
##'   \itemize{
##'     \item If \code{weights = prevalence = NULL}, each observation (presence or absence) will 
##'     have the same weight, no matter the total number of presences and absences.
##'     \item If \code{prevalence = 0.5}, presences and absences will be weighted equally 
##'     (\emph{i.e. the weighted sum of presences equals the weighted sum of absences}). 
##'     \item If \code{prevalence} is set below (\emph{above}) \code{0.5}, more weight will be 
##'     given to absences (\emph{presences}).
##'     \item If \code{weights} is defined, \code{prevalence} argument will be ignored, and each 
##'     observation will have its own weight.
##'     \item If pseudo-absences have been generated (\code{PA.nb.rep > 0} in 
##'     \code{\link{BIOMOD_FormatingData}}), weights are by default calculated such that 
##'     \code{prevalence = 0.5}. \emph{Automatically created \code{weights} will be \code{integer} 
##'     values to prevent some modeling issues.}
##'   }}
##' 
##'   \item{metric.eval}{
##'   \itemize{
##'     \item \code{ROC} : Relative Operating Characteristic
##'     \item \code{KAPPA} : Cohen's Kappa (Heidke skill score)
##'     \item \code{TSS} : True kill statistic (Hanssen and Kuipers discriminant, Peirce's skill 
##'     score)
##'     \item \code{FAR} : False alarm ratio
##'     \item \code{SR} : Success ratio
##'     \item \code{ACCURANCY} : Accuracy (fraction correct)
##'     \item \code{BIAS} : Bias score (frequency bias)
##'     \item \code{POD} : Probability of detection (hit rate)
##'     \item \code{CSI} : Critical success index (threat score)
##'     \item \code{ETS} : Equitable threat score (Gilbert skill score)
##'   }
##'   Optimal value of each method can be obtained with the \code{\link{get_optim_value}} 
##'   function. Several evaluation metrics can be selected. \emph{Please refer to the 
##'   \href{https://www.cawcr.gov.au/projects/verification/}{CAWRC website (section "Methods for dichotomous forecasts")} 
##'   to get detailed description of each metric.}
##'   }
##'   
##'   \item{save.output}{\emph{If this argument is set to \code{FALSE}, it may prevent the evaluation 
##'   of the ensemble models (see \code{\link{BIOMOD_EnsembleModeling}}) in further steps. Strong 
##'   recommandation is to keep \code{save.output = TRUE}, even if it requires to have some free 
##'   space onto the hard drive.}
##'   }
##'   
##'   \item{scale.models}{\bold{This parameter is quite experimental and it is recommended 
##'   not to use it. It may lead to reduction in projection scale amplitude.} Some categorical 
##'   models always have to be scaled (\code{FDA}, \code{ANN}), but it may be interesting to 
##'   scale all computed models to ensure comparable predictions (\code{0-1000} range). It might 
##'   be particularly useful when doing ensemble forecasting to remove the scale prediction effect 
##'   (\emph{the more extended projections are, the more they influence ensemble forecasting 
##'   results}).
##'   }
##'   
##'   \item{do.full.models}{Building models with all available information may be useful in some 
##'   particular cases (\emph{e.g. rare species with few presences points}). But calibration and 
##'   evaluation datasets will be the same, which might lead to over-optimistic evaluation scores.
##'   }
##' }
##' 
##' 
##' @keywords models regression nonlinear multivariate nonparametric tree
##' 
##' 
##' @seealso \code{\link[stats]{glm}}, \code{\link[gam]{gam}}, \code{\link[mgcv]{gam}}, 
##' \code{\link[mgcv]{bam}}, \code{\link[gbm]{gbm}}, \code{\link[rpart]{rpart}}, 
##' code{\link[nnet]{nnet}}, \code{\link[mda]{fda}}, \code{\link[earth]{earth}}, 
##' \code{\link[randomForest]{randomForest}}, \code{\link[maxnet]{maxnet}},
##' \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_ModelingOptions}}, 
##' \code{\link{BIOMOD_CrossValidation}}, \code{ \link{bm_VariablesImportance}}, 
##' \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleModeling}},
##' \code{\link{bm_PlotEvalMean}}, \code{\link{bm_PlotEvalBoxplot}}, 
##' \code{\link{bm_PlotVarImpBoxplot}}, \code{\link{bm_PlotResponseCurves}}
##' @family Main functions
##' 
##'   
##' @examples
##' 
##' # Load species occurrences (6 species available)
##' myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- 'GuloGulo'
##' 
##' # Get corresponding presence/absence data
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # Get corresponding XY coordinates
##' myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
##' 
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl <- raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExpl <- raster::stack(raster::crop(myExpl, myExtent))
##' }
##' 
##' # ---------------------------------------------------------------
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # Create default modeling options
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##' 
##' # ---------------------------------------------------------------
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     models = c('RF', 'GLM'),
##'                                     bm.options = myBiomodOptions,
##'                                     nb.rep = 2,
##'                                     data.split.perc = 80,
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 2,
##'                                     do.full.models = FALSE,
##'                                     seed.val = 42)
##' myBiomodModelOut
##' 
##' # Get evaluation scores & variables importance
##' get_evaluations(myBiomodModelOut)
##' get_variables_importance(myBiomodModelOut, as.data.frame = TRUE)
##' 
##' # Represent evaluation scores 
##' bm_PlotEvalMean(bm.out = myBiomodModelOut)
##' bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
##' 
##' # Represent variables importance 
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'dataset'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'dataset'))
##' 
##' # Represent response curves 
##' # bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
##' #                       models.chosen = get_built_models(myBiomodModelOut)[c(1:2)],
##' #                       fixed.var = 'median')
##' # bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
##' #                       models.chosen = get_built_models(myBiomodModelOut)[c(1:2)],
##' #                       fixed.var = 'min')
##' # bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
##' #                       models.chosen = get_built_models(myBiomodModelOut)[3],
##' #                       fixed.var = 'median',
##' #                       do.bivariate = TRUE)
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_Modeling <- function(bm.format,
                            modeling.id = as.character(format(Sys.time(), "%s")),
                            models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS'
                                       , 'RF', 'MAXENT.Phillips', 'MAXENT.Phillips.2'),
                            bm.options = NULL,
                            nb.rep = 1,
                            data.split.perc = 100,
                            data.split.table = NULL,
                            do.full.models = TRUE,
                            weights = NULL,
                            prevalence = NULL,
                            metric.eval = c('KAPPA', 'TSS', 'ROC'),
                            var.import = 0,
                            save.output = TRUE,
                            scale.models = FALSE,
                            nb.cpu = 1,
                            seed.val = NULL,
                            do.progress = TRUE)
{
  .bm_cat("Build Single Models")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_Modeling.check.args(bm.format, modeling.id, models, bm.options, nb.rep
                                      , data.split.perc, data.split.table
                                      , do.full.models, weights, prevalence, metric.eval, var.import
                                      , save.output, scale.models, nb.cpu, seed.val, do.progress)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  models.out <- new('BIOMOD.models.out',
                    dir.name = bm.format@dir.name,
                    sp.name = bm.format@sp.name,
                    modeling.id = modeling.id,
                    expl.var.names = colnames(bm.format@data.env.var),
                    has.evaluation.data = bm.format@has.data.eval,
                    scale.models = scale.models)
  
  ## 2. Create simulation directories -------------------------------------------------------------
  ## Various objects will be stored (models, predictions, projections)
  ## Projections directories are created in Projection() function
  .BIOMOD_Modeling.prepare.workdir(bm.format@dir.name, bm.format@sp.name, models.out@modeling.id)
  name.BIOMOD_DATA = file.path(models.out@dir.name, models.out@sp.name, ".BIOMOD_DATA", models.out@modeling.id)
  
  ## 3.1 Save input data and models options -----------------------------------
  if (save.output) {
    models.out = .fill_BIOMOD.models.out("formated.input.data", bm.format, models.out
                                         , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
    models.out = .fill_BIOMOD.models.out("models.options", bm.options, models.out
                                         , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  }
  
  ## 3.2 Get and save calibration lines ---------------------------------------
  mod.prep.dat <- .BIOMOD_Modeling.prepare.data(bm.format, 
                                                nb.rep, 
                                                data.split.perc, 
                                                weights, 
                                                prevalence, 
                                                do.full.models, 
                                                data.split.table,
                                                seed.val)
  rm(bm.format)
  
  calib.lines <- mod.prep.dat[[1]]$calib.lines
  if (length(mod.prep.dat) > 1) { ## stack calib lines matrix along array 3rd-dimension
    for (pa in 2:length(mod.prep.dat)) {
      calib.lines <- abind(calib.lines, mod.prep.dat[[pa]]$calib.lines, along = 3)
    }
  } 
  models.out = .fill_BIOMOD.models.out("calib.lines", calib.lines, models.out
                                       , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  rm(calib.lines)
  
  ## 4. Print modeling summary in console ---------------------------------------------------------
  .BIOMOD_Modeling.summary(mod.prep.dat, models)
  
  ## 5. Run models with loop over PA --------------------------------------------------------------
  mod.out <- lapply(mod.prep.dat,
                    FUN = bm_RunModelsLoop,
                    modeling.id = models.out@modeling.id,
                    model = models,
                    bm.options = bm.options,
                    var.import = var.import,
                    metric.eval = metric.eval,
                    save.output = save.output,
                    scale.models = scale.models,
                    nb.cpu = nb.cpu,
                    seed.val = seed.val,
                    do.progress = do.progress)
  
  ## 3.3 Rearrange and save outputs -------------------------------------------
  models.out@models.computed <- .transform_outputs_list(mod.out, out = 'models.run')
  models.out@models.failed <- .transform_outputs_list(mod.out, out = 'calib.failure')
  
  ## 3.4 Rearrange and save models outputs : ----------------------------------
  ## models evaluation, variables importance, models prediction, predictions evaluation
  if (save.output) {
    models.evaluation <- .transform_outputs_list(mod.out, out = 'evaluation')
    models.out = .fill_BIOMOD.models.out("models.evaluation", models.evaluation, models.out
                                         , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
    rm(models.evaluation)
    
    if (var.import > 0) {
      variables.importance <- .transform_outputs_list(mod.out, out = 'var.import')
      models.out = .fill_BIOMOD.models.out("variables.importance", variables.importance, models.out
                                           , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
      rm(variables.importance)
    }
    
    models.prediction <- .transform_outputs_list(mod.out, out = 'prediction')
    models.out = .fill_BIOMOD.models.out("models.prediction", models.prediction, models.out
                                         , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
    rm(models.prediction)
    
    models.prediction.eval <- .transform_outputs_list(mod.out, out = 'prediction.eval')
    models.out = .fill_BIOMOD.models.out("models.prediction.eval", models.prediction.eval, models.out
                                         , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
    rm(models.prediction.eval)
  }
  rm(mod.out)
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  name.OUT = paste0(models.out@sp.name, '.', models.out@modeling.id, '.models.out')
  models.out@link <- file.path(models.out@dir.name, models.out@sp.name, name.OUT)
  assign(x = name.OUT, value = models.out)
  save(list = name.OUT, file = models.out@link)
  
  .bm_cat("Done")
  return(models.out)
}


###################################################################################################

.BIOMOD_Modeling.prepare.workdir <- function(dir.name, sp.name, modeling.id)
{
  cat("\nCreating suitable Workdir...\n")
  dir.create(file.path(dir.name, sp.name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, ".BIOMOD_DATA", modeling.id), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, "models", modeling.id), showWarnings = FALSE, recursive = TRUE)
}


###################################################################################################

.BIOMOD_Modeling.check.args <- function(bm.format, modeling.id, models, bm.options, nb.rep
                                        , data.split.perc, data.split.table
                                        , do.full.models, weights, prevalence, metric.eval, var.import
                                        , save.output, scale.models, nb.cpu, seed.val, do.progress)
{
  ## 0. Check bm.format and models arguments ----------------------------------
  cat('\n\nChecking Models arguments...\n')
  
  if (!is.character(modeling.id) || length(modeling.id) > 1) { stop("modeling.id must be a 'character' of length 1") }
  
  .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  if (!is.character(models)) { stop("models must be a 'character' vector") }
  models <- unique(models)
  models.switch.off <- NULL
  
  ## check if model is supported
  avail.models.list <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS'
                         , 'RF', 'MAXENT.Phillips', 'MAXENT.Phillips.2')
  .fun_testIfIn(TRUE, "models", models, avail.models.list)
  
  
  ## 1.1 Remove models not supporting categorical variables --------------------
  categorical_var <- unlist(sapply(colnames(bm.format@data.env.var), function(x) {
    if (is.factor(bm.format@data.env.var[, x])) { return(x) } else { return(NULL) }
  }))
  if (length(categorical_var) > 0) {
    models.fact.unsuprort <- c("SRE", "MAXENT.Tsuruoka")
    models.switch.off <- c(models.switch.off, intersect(models, models.fact.unsuprort))
    if (length(models.switch.off) > 0) {
      models <- setdiff(models, models.switch.off)
      cat(paste0("\n\t! ", paste(models.switch.off, collapse = ",")," were switched off because of categorical variables !"))
    }
  }
 
  
  ## 2.1 Disable MAXENT.Tsuruoka ----------------------------------------------
  ## because of package maintaining issue (request from B Ripley 03-2019)
  if ('MAXENT.Tsuruoka' %in% models) {
    models.switch.off <- unique(c(models.switch.off, "MAXENT.Tsuruoka"))
    models <- setdiff(models, models.switch.off)
    warning('MAXENT.Tsuruoka has been disabled because of package maintaining issue (request from cran team 03-2019)')
  }
  
  ## 3. Check bm.options arguments --------------------------------------------
  if (!is.null(bm.options)) {
    .fun_testIfInherits(TRUE, "bm.options", bm.options, "BIOMOD.models.options")
  } else {
    warning("Models will run with 'defaults' parameters", immediate. = TRUE)
    bm.options <- BIOMOD_ModelingOptions()
  }
  
  ## 2.2 Specific check for MAXENT.Phillips -----------------------------------
  if ("MAXENT.Phillips" %in% models)
  {
    if (!file.exists(file.path(bm.options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"))) {
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
    } else if (nrow(bm.format@coord) == 1) {
      warning("MAXENT.Phillips has been disabled because no XY coordinates have been given", immediate. = TRUE)
      models = models[-which(models == 'MAXENT.Phillips')]
    }
  }
  
  ## 4. Check nb.rep and data.split.table arguments ---------------------------
  if (!is.null(data.split.table)) {
    cat("\n! User defined data-split table was given -> nb.rep, data.split.perc and do.full.models argument will be ignored")
    if (!(length(dim(data.split.table) %in% c(2, 3)))) { stop("data.split.table must be a matrix or a 3D array") }
    if (dim(data.split.table)[1] != length(bm.format@data.species)) { stop("data.split.table must have as many rows (dim1) than your species as data") }
    nb.rep <- dim(data.split.table)[2]
    data.split.perc <- 50
    do.full.models <- FALSE
  }
  
  .fun_testIfPosInt(TRUE, "nb.rep", nb.rep)
  if (data.split.perc < 0 || data.split.perc > 100) {
    stop("data.split.perc argument must be a 0-100 'numeric'")
  } else if (data.split.perc < 50) {
    warning("You chose to allocate more data to evaluation than to calibration of your model
            (data.split.perc<50)\nMake sure you really wanted to do that. \n", immediate. = TRUE)
  } else if (data.split.perc == 100) {
    nb.rep <- 0
    warning(paste0("The models will be evaluated on the calibration data only "
                   , "(nb.rep=0 and no independent data) \n\t "
                   , "It could lead to over-optimistic predictive performances.\n")
            , immediate. = TRUE)
  }
  
  ## 5. Check weights arguments -----------------------------------------------
  if (!is.null(weights)) {
    if (!is.numeric(weights)) { stop("weights must be a numeric vector") }
    if (length(weights) != length(bm.format@data.species)) {
      stop("The number of 'Weight' does not match with the input calibration data. Simulation cannot proceed.")
    }
  }
  
  ## 6. Check metric.eval arguments -------------------------------------------
  metric.eval <- unique(metric.eval)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC')
  .fun_testIfIn(TRUE, "metric.eval", metric.eval, avail.eval.meth.list)
  
  ## 7. Check prevalence arguments --------------------------------------------
  if (!is.null(prevalence)) {
    .fun_testIf01(TRUE, "prevalence", prevalence)
    if ("MAXENT.Phillips" %in% models) {
      cat("\n\t MAXENT.Phillips default prevalence option was updated to fit with modeling prevalence (i.e", prevalence, ")")
      bm.options@MAXENT.Phillips$defaultprevalence = prevalence
    }
  } else {
    prevalence = 0.5
  }
  
  ##### TO BE CHANGE BUT PREVENT FROM BUGS LATER :  Force object saving parameter
  if (!save.output) {
    cat("\n\t save.output param was automatically set to TRUE to prevent bugs.")
    save.output <- TRUE
  }
  
  if (!is.null(seed.val)) {
    set.seed(seed.val)
  }
  
  if (is.null(var.import)) {
    var.import = 0
  }
  
  return(list(models = models,
              bm.options = bm.options,
              nb.rep = nb.rep,
              data.split.perc = data.split.perc,
              weights = weights,
              var.import = var.import,
              metric.eval = metric.eval,
              prevalence = prevalence,
              do.full.models = do.full.models,
              save.output = save.output,
              data.split.table = data.split.table,
              seed.val = seed.val,
              do.progress = do.progress))
}


###################################################################################################

.BIOMOD_Modeling.summary <- function(mod.prep.dat, models)
{
  cat("\n\n")
  .bm_cat(paste(unlist(strsplit(mod.prep.dat[[1]]$name, '_'))[1], "Modeling Summary"))
  cat("\n", ncol(mod.prep.dat[[1]]$dataBM) - 1, " environmental variables (", colnames(mod.prep.dat[[1]]$dataBM)[-1], ")")
  cat("\nNumber of evaluation repetitions :", ncol(mod.prep.dat[[1]]$calib.lines))
  cat("\nModels selected :", models, "\n")
  cat("\nTotal number of model runs:", ncol(mod.prep.dat[[1]]$calib.lines) * length(models) * length(mod.prep.dat), "\n")
  .bm_cat()
}


###################################################################################################

## # #For ecospat package
## @export
## 

setGeneric(".BIOMOD_Modeling.prepare.data", def = function(bm.format, ...) { standardGeneric(".BIOMOD_Modeling.prepare.data") })

setMethod('.BIOMOD_Modeling.prepare.data', signature('BIOMOD.formated.data'),
          function(bm.format, nb.rep, data.split.perc, weights = NULL, prevalence = NULL
                   , do.full.models = TRUE, data.split.table = NULL, seed.val = NULL)
          {
            list.out <- list()
            name <- paste0(bm.format@sp.name, '_AllData')
            xy <- bm.format@coord
            dataBM = cbind(bm.format@data.species, bm.format@data.env.var)
            colnames(dataBM)[1] = bm.format@sp.name
            
            ## dealing with evaluation data
            if (bm.format@has.data.eval) {
              eval.data <- data.frame(cbind(bm.format@eval.data.species, bm.format@eval.data.env.var[, , drop = FALSE]))
              colnames(eval.data)[1] <- bm.format@sp.name
              eval.xy <- bm.format@eval.coord
            } else {
              eval.data <- eval.xy <- NULL
            }
            
            ## Calib/Valid lines
            if (!is.null(data.split.table)) {
              calib.lines <- data.split.table
              colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            } else {
              if (nb.rep == 0) { # take all available data
                calib.lines <- matrix(rep(TRUE, length(bm.format@data.species)), ncol = 1)
                colnames(calib.lines) <- '_Full'
              } else {
                calib.lines <- .sample_mat(data.sp = bm.format@data.species,
                                           data.split = data.split.perc,
                                           nb.rep = nb.rep,
                                           data.env = bm.format@data.env.var,
                                           seed.val = seed.val)
                if (do.full.models) {
                  calib.lines <- cbind(calib.lines, rep(TRUE, length(bm.format@data.species)))
                  colnames(calib.lines)[nb.rep + 1] <- '_Full'
                }
              }
            }
            ## force calib.lines object to be 3D array
            if (length(dim(calib.lines)) < 3) {
              dn_tmp <- dimnames(calib.lines) ## keep track of dimnames
              dim(calib.lines) <- c(dim(calib.lines), 1)
              dimnames(calib.lines) <- list(dn_tmp[[1]], dn_tmp[[2]], "_AllData")
            }
            
            if (is.null(weights)) { # 1 for all points
              if (!is.null(prevalence)) {
                cat("\n\t> Automatic weights creation to rise a", prevalence, "prevalence")
                weights <- .automatic_weights_creation(bm.format@data.species , prev = prevalence)
              } else {
                cat("\n\t> No weights : all observations will have the same weight")
                weights <- rep(1, length(bm.format@data.species))
              }
            }
            
            list.out[[name]] <- list(name = name,
                                     dir.name = bm.format@dir.name,
                                     xy = xy,
                                     dataBM = dataBM, 
                                     calib.lines = calib.lines,
                                     weights = weights,
                                     eval.data = eval.data,
                                     eval.xy = eval.xy)
            return(list.out)
          }
)

setMethod('.BIOMOD_Modeling.prepare.data', signature('BIOMOD.formated.data.PA'),
          function(bm.format, nb.rep, data.split.perc, weights = NULL, prevalence = NULL
                   , do.full.models = TRUE, data.split.table = NULL, seed.val = NULL)
          {
            list.out <- list()
            formal_weights <- weights
            for (pa in 1:ncol(bm.format@PA.table))
            {
              weights <- formal_weights
              name <- paste0(bm.format@sp.name, "_", colnames(bm.format@PA.table)[pa])
              xy <- bm.format@coord[bm.format@PA.table[, pa], ]
              resp <- bm.format@data.species[bm.format@PA.table[, pa]] # response variable (with pseudo absences selected)
              resp[is.na(resp)] <- 0
              dataBM <- data.frame(cbind(resp, bm.format@data.env.var[bm.format@PA.table[, pa], , drop = FALSE]))
              colnames(dataBM)[1] <- bm.format@sp.name
              
              ## Calib/Valid lines
              if (!is.null(data.split.table))
              {
                if (length(dim(data.split.table)) == 2) {
                  calib.lines <- data.split.table
                } else {
                  calib.lines <- asub(data.split.table, pa, 3, drop = TRUE)
                }
                colnames(calib.lines) <- paste0('_RUN', 1:ncol(calib.lines))
                calib.lines[which(!bm.format@PA.table[, pa]), ] <- NA
              } else {
                if (nb.rep == 0) { # take all available data
                  calib.lines <- matrix(NA, nrow = length(bm.format@data.species), ncol = 1)
                  calib.lines[bm.format@PA.table[, pa], 1] <- TRUE
                  colnames(calib.lines) <- '_Full'
                } else {
                  calib.lines <- matrix(NA, nrow = length(bm.format@data.species), ncol = nb.rep)
                  sampled.mat <- .sample_mat(
                    data.sp = bm.format@data.species[bm.format@PA.table[, pa]],
                    data.split = data.split.perc,
                    nb.rep = nb.rep,
                    data.env = bm.format@data.env.var[bm.format@PA.table[, pa], , drop = FALSE],
                    seed.val = seed.val
                  )
                  calib.lines[bm.format@PA.table[, pa], ] <- sampled.mat
                  colnames(calib.lines) <- colnames(sampled.mat)
                  if (do.full.models) {
                    calib.lines <- cbind(calib.lines, rep(NA, length(bm.format@data.species)))
                    calib.lines[bm.format@PA.table[, pa], nb.rep + 1] <- TRUE
                    colnames(calib.lines)[nb.rep + 1] <- '_Full'
                  }
                }
              }
              
              ## force calib.lines object to be 3D array
              if (length(dim(calib.lines)) < 3) {
                dn_tmp <- dimnames(calib.lines) ## keep track of dimnames
                dim(calib.lines) <- c(dim(calib.lines), 1)
                dimnames(calib.lines) <- list(dn_tmp[[1]], dn_tmp[[2]], paste0("_PA", pa))
              }
              
              # dealing with evaluation data
              if (bm.format@has.data.eval) {
                eval.data <- data.frame(cbind(bm.format@eval.data.species, bm.format@eval.data.env.var))
                colnames(eval.data)[1] <- bm.format@sp.name
                eval.xy <- bm.format@eval.coord
              } else {
                eval.data <- eval.xy <- NULL
              }
              
              if (is.null(weights)) { # prevalence of 0.5... may be parametrize
                if (is.null(prevalence)) { prevalence <- 0.5 }
                cat("\n\t\t\t! Weights where automatically defined for", name, "to rise a", prevalence, "prevalence !")
                weights <- rep(NA, length(bm.format@data.species))
                weights[bm.format@PA.table[, pa]] <- .automatic_weights_creation(as.numeric(dataBM[, 1]) , prev = prevalence)
              } else { # remove useless weights
                weights[!bm.format@PA.table[, pa]] <- NA
              }
              
              list.out[[name]] <- list(name = name,
                                       dir.name = bm.format@dir.name,
                                       xy = xy, 
                                       dataBM = dataBM,
                                       calib.lines = calib.lines,
                                       weights = weights,
                                       eval.data = eval.data,
                                       eval.xy = eval.xy)
            }
            return(list.out)
          }
)


###################################################################################################

.automatic_weights_creation <- function(resp, prev = 0.5, subset = NULL)
{
  if (is.null(subset)) { subset <- rep(TRUE, length(resp)) }
  
  nbPres <- sum(resp[subset], na.rm = TRUE)
  # The number of true absences + pseudo absences to maintain true value of prevalence
  nbAbsKept <- sum(subset, na.rm = TRUE) - sum(resp[subset], na.rm = TRUE)
  weights <- rep(1, length(resp))
  
  if (nbAbsKept > nbPres) {
    # code absences as 1
    weights[which(resp > 0)] <- (prev * nbAbsKept) / (nbPres * (1 - prev))
  } else {
    # code presences as 1
    weights[which(resp == 0 | is.na(resp))] <- (nbPres * (1 - prev)) / (prev * nbAbsKept)
  }
  weights = round(weights[])
  weights[!subset] <- 0
  
  return(weights)
}

.sample_mat <- function(data.sp, data.split, nb.rep = 1, data.env = NULL, seed.val = NULL)
{
  # data.sp is a 0, 1 vector
  # return a matrix with nb.rep columns of boolean (T: calib, F= eval)
  
  pres <- which(data.sp == 1)
  abs <- (1:length(data.sp))[-pres]
  
  nbPresEval <- round(length(pres) * data.split / 100)
  nbAbsEval <- round(length(abs) * data.split / 100)
  
  mat.out <- matrix(FALSE, nrow = length(data.sp), ncol = nb.rep)
  colnames(mat.out) <- paste0('_RUN', 1:nb.rep)
  
  set.seed(seed.val)
  for (i in 1:ncol(mat.out)) {
    ## force to sample at least one level of each factorial variable for calibration
    fact.cell.samp <- NULL
    if (!is.null(data.env)) {
      fact.cell.samp <- bm_SampleFactorLevels(expl.var = data.env)
      mat.out[fact.cell.samp, i] <- TRUE ## in fact.cell.samp
    }
    mat.out[sample(setdiff(pres, fact.cell.samp), ## in pres, not in fact.cell.samp
                   max(nbPresEval - length(fact.cell.samp), 0)), i] <- TRUE
    mat.out[sample(setdiff(abs, fact.cell.samp), ## in abs, not in fact.cell.samp
                   max(nbAbsEval - length(fact.cell.samp), 0)), i] <- TRUE
  }
  return(mat.out)
}


## ###############################################################################################
## 
## Reshape biomod2 objects
## 
## This is an internal function (developper only)
##
## @param mod.out the object to transform given as a list
## @param out character, the type of input object
## @param dim.names character, if not `NULL` the resshaped object will be stored on the hard drive
##
## @return list, the extracted statistics
## 
## @export
## 
## ###############################################################################################

.transform_outputs_list = function(mod.out, out = 'evaluation', dim.names = NULL)
{
  out_list = c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
               'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation')
  .fun_testIfIn(TRUE, "out", out, out_list)
  
  ## 0.a get dataset names ------------------------------------------------------------------------
  if (length(mod.out) == 1 && length(unlist(strsplit(unlist(names(mod.out)), '_'))) == 1) {
    dataset.names <- 'AllData'
  } else if (is.null(dim.names)) {
    dataset.names <- unlist(sapply(unlist(names(mod.out)), function(name) {
      return(tail(unlist(strsplit(name, '_')), 1))
    }))
  } else {
    dataset.names <- unlist(dim.names[1])
  }
  
  ## 0.b get run.eval and model names -------------------------------------------------------------
  if (is.null(dim.names)) {
    run.eval.names <- sub('_', '', unlist(names(mod.out[[1]]))) # may be good here to test that all names are identics
    mod.names <- unlist(names(mod.out[[1]][[1]]))
  } else {
    run.eval.names <- unlist(dim.names[2])
    mod.names <- unlist(dim.names[3])
  }
  
  ## 1. CASE evaluation / prediction / prediction.eval / var.import -------------------------------
  
  if (out %in% c("evaluation", "prediction", "prediction.eval", "var.import"))
  {
    nb_pa <- length(mod.out)
    nb_run <- length(mod.out[[1]])
    nb_mod <- length(mod.out[[1]][[1]])
    
    name_slot = out
    if (out == "prediction") { name_slot = "pred" }
    if (out == "prediction.eval") { name_slot = "pred.eval" }
    
    output <- NULL
    for (i in 1:nb_pa) {
      for (j in 1:nb_run) {
        for (k in 1:nb_mod) {
          output <- mod.out[[i]][[j]][[k]][[name_slot]]
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
        dimnames.out = list(names(mod.out[[1]][[1]][[1]][['var.import']]) # to change ?
                            , kept.mod, run.eval.names, dataset.names)
      }
      dim.out = c(nb.tmp,
                  length(kept.mod),
                  length(run.eval.names),
                  length(dataset.names))
    }
    
    output <- lapply(names(mod.out), function(d1) { # data set
      lapply(names(mod.out[[d1]]), function(d2) { # run eval
        lapply(names(mod.out[[d1]][[d2]]), function(d3){ # models
          if (out == "evaluation") {
            if (is.null(mod.out[[d1]][[d2]][[d3]][['calib.failure']])) {
              return(data.frame(mod.out[[d1]][[d2]][[d3]][['evaluation']]))
            } else { 
              return(matrix(NA, ncol = length(eval.col.names), nrow = length(eval.meth.names)
                            , dimnames = list(eval.meth.names, eval.col.names)))
            }
          } else if (out == "prediction" || out == "prediction.eval" || 
                     (out == "var.import" && d3 %in% kept.mod)) {
            if (is.null(mod.out[[d1]][[d2]][[d3]][[name_slot]])) {
              return(rep(NA, nb.tmp))
            } else {
              return(as.numeric(mod.out[[d1]][[d2]][[d3]][[name_slot]]))
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
    if (out == "models.run") { name_slot = "model" }
    
    output <- lapply(names(mod.out),function(d1) { # data set
      lapply(names(mod.out[[d1]]), function(d2) { # run eval
        lapply(names(mod.out[[d1]][[d2]]), function(d3) { # models
          res = mod.out[[d1]][[d2]][[d3]][[name_slot]]
          if (out == "calib.failure") {
            return(as.character(res))
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
    
    if (is.null(mod.out[[1]][[1]][[1]][[name_slot]])) { return(NULL) }
    
    nb.tmp = 1
    dimnames.out = list(NULL, mod.names, run.eval.names, dataset.names)
    if (out == "EF.prediction") {
      nb.tmp <- length(as.numeric(unlist(mod.out[[1]][[1]][[1]][[name_slot]])))
    } else if (out == "EF.evaluation") {
      eval.meth.names <- rownames(as.data.frame(mod.out[[1]][[1]][[1]][[name_slot]]))
      eval.col.names <- colnames(as.data.frame(mod.out[[1]][[1]][[1]][[name_slot]]))
      nb.tmp = c(length(eval.meth.names), length(eval.col.names))
      dimnames.out = list(eval.meth.names, eval.col.names
                          , mod.names, run.eval.names, dataset.names)
    }
    dim.out = c(nb.tmp,
                length(mod.out[[1]][[1]]),
                length(mod.out[[1]]),
                length(mod.out))
    
    
    output <- lapply(1:length(mod.out),function(d1) { # data set
      lapply(1:length(mod.out[[d1]]), function(d2) { # run eval
        lapply(1:length(mod.out[[d1]][[d2]]), function(d3) { # models
          res = mod.out[[d1]][[d2]][[d3]][[name_slot]]
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

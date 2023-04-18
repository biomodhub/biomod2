# BIOMOD_Modeling ---------------------------------------------------------
##' @name BIOMOD_Modeling
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' 
##' @title Run a range of species distribution models
##' 
##' @description This function allows to calibrate and evaluate a range of modeling techniques 
##' for a given species distribution. The dataset can be split up in calibration/validation parts,
##' and the predictive power of the different models can be estimated using a range of evaluation 
##' metrics (see Details).
##' 
##' 
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT}, \code{MAXNET}
##' @param models.pa (\emph{optional, default} \code{NULL}) \cr 
##' A \code{list} containing for each model a \code{vector} defining which pseudo-absence datasets 
##' are to be used, must be among \code{colnames(bm.format@PA.table)}
##' @param bm.options a \code{\link{BIOMOD.models.options}} object returned by the  
##' \code{\link{BIOMOD_ModelingOptions}} function
##' 
##' @param CV.strategy a \code{character} corresponding to the cross-validation selection strategy, 
##' must be among \code{random}, \code{kfold}, \code{block}, \code{strat}, \code{env} or 
##' \code{user.defined}
##' @param CV.nb.rep (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'random'} or \code{strategy = 'kfold'}, an \code{integer} corresponding 
##' to the number of sets (repetitions) of cross-validation points that will be drawn
##' @param CV.perc (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'random'}, a \code{numeric} between \code{0} and \code{1} defining the 
##' percentage of data that will be kept for calibration
##' @param CV.k (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'kfold'} or \code{strategy = 'strat'} or \code{strategy = 'env'}, an 
##' \code{integer} corresponding to the number of partitions 
##' @param CV.balance (\emph{optional, default} \code{'presences'}) \cr
##' If \code{strategy = 'strat'} or \code{strategy = 'env'}, a \code{character} corresponding 
##' to how data will be balanced between partitions, must be either \code{presences} or
##' \code{absences} 
##' @param CV.env.var (\emph{optional}) \cr If \code{strategy = 'env'}, a
##'   \code{character} corresponding to the environmental variables used to
##'   build the partition. \code{k} partitions will be built for each
##'   environmental variables. By default the function uses all environmental
##'   variables available.

##' @param CV.strat (\emph{optional, default} \code{'both'}) \cr
##' If \code{strategy = 'env'}, a \code{character} corresponding to how data will partitioned 
##' along gradient, must be among \code{x}, \code{y}, \code{both}
##' @param CV.user.table (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'user.defined'}, a \code{matrix} or \code{data.frame} defining for each 
##' repetition (in columns) which observation lines should be used for models calibration 
##' (\code{TRUE}) and validation (\code{FALSE})
##' @param CV.do.full.models (\emph{optional, default} \code{TRUE}) \cr  
##' A \code{logical} value defining whether models should be also calibrated and validated over 
##' the whole dataset (and pseudo-absence datasets) or not
##' @param nb.rep  \emph{deprecated}, now called \code{CV.nb.rep}
##' @param data.split.perc \emph{deprecated}, now called \code{CV.perc}
##' @param data.split.table \emph{deprecated}, now called \code{CV.user.table}
##' @param do.full.models \emph{deprecated}, now called \code{CV.do.full.models}
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
##'     \item \code{MAXENT} : Maximum Entropy 
##'     (\url{https://biodiversityinformatics.amnh.org/open_source/maxent/})
##'     \item \code{MAXNET} : Maximum Entropy (\code{\link[maxnet]{maxnet}})
##'   }}
##'   
##'   \item{models.pa}{Different models might respond differently to different numbers of 
##'   pseudo-absences. It is possible to create sets of pseudo-absences with different numbers 
##'   of points (see \code{\link{BIOMOD_FormatingData}}) and to assign only some of these 
##'   datasets to each single model.
##'   }
##'   
##'   \item{CV.[...] parameters}{Different methods are available to calibrate/validate the 
##'   single models (see \code{\link{bm_CrossValidation}}.)}
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
##'   \item{scale.models}{\bold{This parameter is quite experimental and it is recommended 
##'   not to use it. It may lead to reduction in projection scale amplitude.} Some categorical 
##'   models always have to be scaled (\code{FDA}, \code{ANN}), but it may be interesting to 
##'   scale all computed models to ensure comparable predictions (\code{0-1000} range). It might 
##'   be particularly useful when doing ensemble forecasting to remove the scale prediction effect 
##'   (\emph{the more extended projections are, the more they influence ensemble forecasting 
##'   results}).
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
##' \code{\link{bm_CrossValidation}}, \code{ \link{bm_VariablesImportance}}, 
##' \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleModeling}},
##' \code{\link{bm_PlotEvalMean}}, \code{\link{bm_PlotEvalBoxplot}}, 
##' \code{\link{bm_PlotVarImpBoxplot}}, \code{\link{bm_PlotResponseCurves}}
##' @family Main functions
##' 
##'   
##' @examples
##' library(terra)
##' 
##' # Load species occurrences (6 species available)
##' data(DataSpecies)
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
##' data(bioclim_current)
##' myExpl <- terra::rast(bioclim_current)
##' 
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExpl <- terra::crop(myExpl, myExtent)
##' }
##' 
##' # ---------------------------------------------------------------------------- #
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
##' # ---------------------------------------------------------------------------- #
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     models = c('RF', 'GLM'),
##'                                     bm.options = myBiomodOptions,
##'                                     CV.strategy = 'random',
##'                                     CV.nb.rep = 2,
##'                                     CV.perc = 0.8,
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 2,
##'                                     seed.val = 42)
##' myBiomodModelOut
##' 
##' # Get evaluation scores & variables importance
##' get_evaluations(myBiomodModelOut)
##' get_variables_importance(myBiomodModelOut)
##' 
##' # Represent evaluation scores 
##' bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'calibration')
##' bm_PlotEvalMean(bm.out = myBiomodModelOut, dataset = 'validation')
##' bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
##' 
##' # # Represent variables importance 
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'dataset'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'dataset'))
##' 
##' # # Represent response curves 
##' # mods <- get_built_models(myBiomodModelOut, run = 'RUN1')
##' # bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
##' #                       models.chosen = mods,
##' #                       fixed.var = 'median')
##' # bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
##' #                       models.chosen = mods,
##' #                       fixed.var = 'min')
##' # mods <- get_built_models(myBiomodModelOut, full.name = 'GuloGulo_allData_RUN2_RF')
##' # bm_PlotResponseCurves(bm.out = myBiomodModelOut, 
##' #                       models.chosen = mods,
##' #                       fixed.var = 'median',
##' #                       do.bivariate = TRUE)
##' 
##' 
##' @export
##' 
##' 

BIOMOD_Modeling <- function(bm.format,
                            modeling.id = as.character(format(Sys.time(), "%s")),
                            models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE'
                                       , 'FDA', 'MARS', 'RF', 'MAXENT', 'MAXNET'),
                            models.pa = NULL,
                            bm.options = NULL,
                            CV.strategy = 'random',
                            CV.nb.rep = 1,
                            CV.perc = NULL,
                            CV.k = NULL,
                            CV.balance = NULL,
                            CV.env.var = NULL,
                            CV.strat = NULL,
                            CV.user.table = NULL,
                            CV.do.full.models = TRUE,
                            nb.rep,
                            data.split.perc,
                            data.split.table,
                            do.full.models,
                            weights = NULL,
                            prevalence = NULL,
                            metric.eval = c('KAPPA', 'TSS', 'ROC'),
                            var.import = 0,
                            scale.models = FALSE,
                            nb.cpu = 1,
                            seed.val = NULL,
                            do.progress = TRUE)
{
  .bm_cat("Build Single Models")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_Modeling.check.args(
    bm.format = bm.format, 
    modeling.id = modeling.id, 
    models = models, 
    models.pa = models.pa, 
    bm.options = bm.options, 
    weights = weights, 
    prevalence = prevalence, 
    metric.eval = metric.eval, 
    var.import = var.import, 
    scale.models = scale.models,
    nb.cpu = nb.cpu,
    seed.val = seed.val,
    do.progress = do.progress
  )
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  # check for obsolete arguments
  args <- .BIOMOD_Modeling.check.args.obsolete(
    CV.strategy = CV.strategy,
    data.split.perc = data.split.perc,
    CV.perc = CV.perc,
    data.split.table = data.split.table,
    CV.user.table = CV.user.table,
    nb.rep = nb.rep,
    CV.nb.rep = CV.nb.rep,
    do.full.models = do.full.models,
    CV.do.full.models = CV.do.full.models
  )
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
  models.out = .fill_BIOMOD.models.out("formated.input.data", bm.format, models.out
                                       , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  models.out = .fill_BIOMOD.models.out("models.options", bm.options, models.out
                                       , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  
  ## 3.2 Get and save calibration lines ---------------------------------------
  calib.lines <- bm_CrossValidation(bm.format = bm.format,
                                    strategy = CV.strategy,
                                    nb.rep = CV.nb.rep,
                                    perc = CV.perc,
                                    k = CV.k,
                                    balance = CV.balance,
                                    env.var = CV.env.var,
                                    strat = CV.strat,
                                    user.table = CV.user.table,
                                    do.full.models = CV.do.full.models)
  models.out = .fill_BIOMOD.models.out("calib.lines", calib.lines, models.out
                                       , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  
  ## 4. Print modeling summary in console ---------------------------------------------------------
  .BIOMOD_Modeling.summary(bm.format, calib.lines, models, models.pa)
  
  ## 5. Run models with loop over PA --------------------------------------------------------------
  mod.out <- bm_RunModelsLoop(bm.format = bm.format,
                              weights = weights,
                              calib.lines = calib.lines,
                              modeling.id = models.out@modeling.id,
                              models = models,
                              models.pa = models.pa,
                              bm.options = bm.options,
                              var.import = var.import,
                              metric.eval = metric.eval,
                              scale.models = scale.models,
                              nb.cpu = nb.cpu,
                              seed.val = seed.val,
                              do.progress = do.progress)
  
  ## 6. Rearrange and save outputs -------------------------------------------
  models.out@models.computed <- .transform_outputs_list("mod", mod.out, out = "model")
  models.out@models.failed <- .transform_outputs_list("mod", mod.out, out = "calib.failure")
  
  if(length(models.out@models.computed) == 1 && models.out@models.computed == "none"){
    cat("\n! All models failed")
    return(models.out)
  }
  
  ## 3.4 Rearrange and save models outputs : ----------------------------------
  ## models evaluation, variables importance, models prediction, predictions evaluation
  models.evaluation <- .transform_outputs_list("mod", mod.out, out = "evaluation")
  models.out = .fill_BIOMOD.models.out("models.evaluation", models.evaluation, models.out
                                       , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
  rm(models.evaluation)
  
  if (var.import > 0) {
    variables.importance <- .transform_outputs_list("mod", mod.out, out = "var.import")
    models.out = .fill_BIOMOD.models.out("variables.importance", variables.importance, models.out
                                         , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
    rm(variables.importance)
  }
  
  models.prediction <- .transform_outputs_list("mod", mod.out, out = "pred")
  models.out = .fill_BIOMOD.models.out("models.prediction", models.prediction, models.out
                                       , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  rm(models.prediction)
  
  models.prediction.eval <- .transform_outputs_list("mod", mod.out, out = "pred.eval")
  models.out = .fill_BIOMOD.models.out("models.prediction.eval", models.prediction.eval, models.out
                                       , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  rm(models.prediction.eval)
  rm(mod.out)
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE ----------------------------
  name.OUT = paste0(models.out@sp.name, '.', models.out@modeling.id, '.models.out')
  models.out@link <- file.path(models.out@dir.name, models.out@sp.name, name.OUT)
  assign(x = name.OUT, value = models.out)
  save(list = name.OUT, file = models.out@link)
  
  .bm_cat("Done")
  return(models.out)
}


# ---------------------------------------------------------------------------- #

.BIOMOD_Modeling.prepare.workdir <- function(dir.name, sp.name, modeling.id)
{
  cat("\nCreating suitable Workdir...\n")
  dir.create(file.path(dir.name, sp.name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, ".BIOMOD_DATA", modeling.id), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, "models", modeling.id), showWarnings = FALSE, recursive = TRUE)
}


# ---------------------------------------------------------------------------- #

.BIOMOD_Modeling.check.args <- function(bm.format, modeling.id, models, models.pa, bm.options
                                        , weights, prevalence, metric.eval, var.import
                                        , scale.models, nb.cpu, seed.val, do.progress)
{
  ## 0. Check bm.format and models arguments ----------------------------------
  cat('\n\nChecking Models arguments...\n')
  
  if (!is.character(modeling.id) || length(modeling.id) > 1) { stop("modeling.id must be a 'character' of length 1") }
  
  .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  if (!is.character(models)) { stop("models must be a 'character' vector") }
  
  # Support for old names in models
  # Deprecated MAXENT.Phillips/MAXENT.Phillips.2
  if (any(models == "MAXENT.Phillips")) {
    models[which(models == "MAXENT.Phillips")] <- "MAXENT"
    cat(paste0("\n\t! 'MAXENT.Phillips' model name is deprecated, please use 'MAXENT' instead."))
  }
  if (any(models == "MAXENT.Phillips.2")) {
    models[which(models == "MAXENT.Phillips.2")] <- "MAXNET"
    cat(paste0("\n\t! 'MAXENT.Phillips.2' model name is deprecated, please use 'MAXNET' instead."))
  }
  ## Deprecated MAXENT.Tsuruoka 
  ## because of package maintaining issue (request from B Ripley 03-2019)
  if ('MAXENT.Tsuruoka' %in% models) {
    models.switch.off <- unique(c(models.switch.off, "MAXENT.Tsuruoka"))
    models <- setdiff(models, models.switch.off)
    warning('MAXENT.Tsuruoka has been disabled because of package maintaining issue (request from cran team 03-2019)')
  }
  models <- unique(models)
  models.switch.off <- NULL
  
  ## check if model is supported
  avail.models.list <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS'
                         , 'RF', 'MAXENT', 'MAXNET')
  .fun_testIfIn(TRUE, "models", models, avail.models.list)
  
  
  ## 1.1 Remove models not supporting categorical variables --------------------
  categorical_var <- .get_categorical_names(bm.format@data.env.var)
  
  if (length(categorical_var) > 0) {
    models.fact.unsupport <- c("SRE", "MAXENT.Tsuruoka")
    models.switch.off <- c(models.switch.off, intersect(models, models.fact.unsupport))
    if (length(models.switch.off) > 0) {
      models <- setdiff(models, models.switch.off)
      cat(paste0("\n\t! ", paste(models.switch.off, collapse = ",")," were switched off because of categorical variables !"))
    }
  }
  
  ## Check models.pa argument
  if (!is.null(models.pa)) {
    if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
      .fun_testIfInherits(TRUE, "models.pa", models.pa, "list")
      .fun_testIfIn(TRUE, "unlist(models.pa)", unlist(models.pa), colnames(bm.format@PA.table))
      .fun_testIfIn(TRUE, "names(models.pa)", names(models.pa), models)
      if (length(models.pa) != length(models)) {
        mod.miss = models[-which(models %in% names(models.pa))]
        list.miss = rep(list(colnames(bm.format@PA.table)), length(mod.miss))
        names(list.miss) = mod.miss
        models.pa = c(models.pa, list.miss)
        warning(paste0(paste0(mod.miss, collapse = ", ")
                       , " have been assigned to all PA datasets as no information was given")
                , immediate. = TRUE)
      }
    } else {
      warning("models.pa has been disabled because no PA datasets have been given", immediate. = TRUE)
      models.pa = NULL
    }
  }
  
  
  ## 3. Check bm.options arguments --------------------------------------------
  if (!is.null(bm.options)) {
    .fun_testIfInherits(TRUE, "bm.options", bm.options, "BIOMOD.models.options")
  } else {
    warning("Models will run with 'defaults' parameters", immediate. = TRUE)
    bm.options <- BIOMOD_ModelingOptions()
  }
  
  ## 2.2 Specific check for MAXENT -----------------------------------
  if ("MAXENT" %in% models)
  {
    if (!file.exists(file.path(bm.options@MAXENT$path_to_maxent.jar, "maxent.jar"))) {
      models = models[-which(models == 'MAXENT')]
      if (!is.null(models.pa)) {
        models.pa = models.pa[-which(names(models.pa) == 'MAXENT')]
      }
      warning(paste0("MAXENT has been disabled because the maxent.jar file is missing. "
                     , "`maxent.jar` file must be downloaded (https://biodiversityinformatics.amnh.org/open_source/maxent/) "
                     , "and put in the working directory."), immediate. = TRUE)
      ## -- 
      ## The java installation check is temporally disabled cause it seems to cause 
      ## issues on some Windows users machine.
      ## --
      # } else if(!.check.java.installed()){
      #   models = models[-which(models=='MAXENT')]
    } else if (nrow(bm.format@coord) == 1) {
      models = models[-which(models == 'MAXENT')]
      if (!is.null(models.pa)) {
        models.pa = models.pa[-which(names(models.pa) == 'MAXENT')]
      }
      warning("MAXENT has been disabled because no XY coordinates have been given", immediate. = TRUE)
    }
  }
  
  ## 5. Check prevalence arguments --------------------------------------------
  if (!is.null(prevalence)) {
    .fun_testIf01(TRUE, "prevalence", prevalence)
    if ("MAXENT" %in% models) {
      cat("\n\t MAXENT default prevalence option was updated to fit with modeling prevalence (i.e", prevalence, ")")
      bm.options@MAXENT$defaultprevalence = prevalence
    }
  } else {
    prevalence = 0.5
  }
  
  ## 6. Check weights arguments -----------------------------------------------
  if (is.null(weights)) {
    if (!is.null(prevalence)) {
      cat("\n\t> Automatic weights creation to rise a", prevalence, "prevalence")
      data.sp <- as.numeric(bm.format@data.species)
      if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
        weights.pa <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
          {
            ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
            data.sp_pa <- data.sp[ind.PA]
            data.sp_pa[which(is.na(data.sp_pa))] <- 0
            weights <- .automatic_weights_creation(data.sp_pa, prev = prevalence)
            
            wei <- rep(NA, length(data.sp))
            wei[ind.PA] <- weights
            return(matrix(wei, ncol = 1))
          }
        weights.pa <- cbind(weights.pa, rep(1, nrow(weights.pa)))
        colnames(weights.pa) <- c(colnames(bm.format@PA.table), "allData")
        weights <- weights.pa
      } else {
        weights <- .automatic_weights_creation(data.sp, prev = prevalence)
        weights <- matrix(weights, nrow = length(weights), ncol = 1)
        colnames(weights) <- "allData"
      }
    }
    # else { ## NEVER OCCURRING NO ??
    #   cat("\n\t> No weights : all observations will have the same weight")
    #   weights <- rep(1, length(bm.format@data.species))
    # }
  } else {
    if (!is.numeric(weights)) { stop("weights must be a numeric vector") }
    if (length(weights) != length(bm.format@data.species)) {
      stop("The number of 'Weight' does not match with the input calibration data. Simulation cannot proceed.")
    }
    if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
      weights.pa <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
        {
          wei <- weights
          wei[which(bm.format@PA.table[, pa] == FALSE | is.na(bm.format@PA.table[, pa]))] <- NA
          return(matrix(wei, ncol = 1))
        }
      weights.pa <- cbind(weights.pa, rep(1, nrow(weights.pa)))
      colnames(weights.pa) <- c(colnames(bm.format@PA.table), "allData")
      weights <- weights.pa
    } else {
      weights <- matrix(weights, nrow = length(weights), ncol = 1)
      colnames(weights) <- "allData"
    }
  }
  
  ## 7. Check metric.eval arguments -------------------------------------------
  metric.eval <- unique(metric.eval)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC')
  .fun_testIfIn(TRUE, "metric.eval", metric.eval, avail.eval.meth.list)
  
  
  if (!is.null(seed.val)) {
    set.seed(seed.val)
  }
  
  if (is.null(var.import)) {
    var.import = 0
  }
  
  return(list(models = models,
              models.pa = models.pa,
              bm.options = bm.options,
              weights = weights,
              var.import = var.import,
              metric.eval = metric.eval,
              prevalence = prevalence,
              seed.val = seed.val,
              do.progress = do.progress))
}


# Obsolete argument Check -------------------------------------------------------


.BIOMOD_Modeling.check.args.obsolete <- function(CV.strategy,
                                                 data.split.perc,
                                                 CV.perc,
                                                 data.split.table,
                                                 CV.user.table,
                                                 nb.rep,
                                                 CV.nb.rep,
                                                 do.full.models,
                                                 CV.do.full.models)
{
  ## do.full.models ------------------------------
  if (!missing(do.full.models)) {
    if (!is.null(CV.do.full.models)) {
      cat("\n! ignored obsolete argument 'do.full.models' as 'CV.do.full.models' was also given")
    } else {
      CV.do.full.models <- do.full.models
      cat("\n!!! argument 'do.full.models' is obsolete, please use 'CV.do.full.models' instead")
    }
  }
  
  ## data.split.perc --------------------------
  if (!missing(data.split.perc)) {
    if (CV.strategy != "random") {
      cat("\n! ignored obsolete argument 'data.split.perc' as 'CV.strategy' was not set to 'random'")
    } else if (!is.null(CV.perc)) {
      cat("\n! ignored obsolete argument 'data.split.perc' as 'CV.perc' was also given")
    } else {
      CV.perc <- data.split.perc/100
      cat("\n!!! argument 'do.full.models' is obsolete, please use 'CV.perc' instead.
          \n /!\ 'CV.perc' is on a scale 0-1 and was set to data.split.perc/100 =",CV.perc)
    }
  }
  
  ## nb.rep --------------------------
  if (!missing(nb.rep)) {
    if (! CV.strategy %in% c("random", "kfold")) {
      cat("\n! ignored obsolete argument 'nb.rep' as 'CV.strategy' was not set to 'random' or 'kfold'")
    } else if (!is.null(CV.nb.rep)) {
      cat("\n! ignored obsolete argument 'nb.rep' as 'CV.nb.rep' was also given")
    } else {
      CV.nb.rep <- nb.rep
      cat("\n!!! argument 'nb.rep' is obsolete, please use 'CV.nb.rep' instead.")
    }
  }
  
  
  ## data.split.table --------------------------
  if (!missing(data.split.table)) {
    if (! CV.strategy %in% c("user.defined")) {
      cat("\n! ignored obsolete argument 'data.split.table' as 'CV.strategy' was not set to 'user.defined'")
    } else if (!is.null(CV.user.table)) {
      cat("\n! ignored obsolete argument 'data.split.table' as 'CV.user.table' was also given")
    } else {
      CV.user.table <- data.split.table
      cat("\n!!! argument 'data.split.table' is obsolete, please use 'CV.user.table' instead.")
    }
  }
  return(list(
    CV.perc = CV.perc,
    CV.user.table = CV.user.table,
    CV.nb.rep = CV.nb.rep,
    CV.do.full.models = CV.do.full.models))
}

## 0. Check bm.format and models arguments ----------------------------------
cat('\n\nChecking Models arguments...\n')

# ---------------------------------------------------------------------------- #

.BIOMOD_Modeling.summary <- function(bm.format, calib.lines, models, models.pa = NULL)
{
  cat("\n\n")
  .bm_cat(paste(bm.format@sp.name, "Modeling Summary"))
  cat("\n", ncol(bm.format@data.env.var), " environmental variables (", colnames(bm.format@data.env.var), ")")
  cat("\nNumber of evaluation repetitions :", ncol(calib.lines))
  cat("\nModels selected :", models, "\n")
  if (is.null(models.pa)) {
    nb.runs = ncol(calib.lines) * length(models)
  } else {
    nb.runs = length(unlist(models.pa))
    nb.runs = ncol(calib.lines) * nb.runs
  }
  cat("\nTotal number of model runs:", nb.runs, "\n")
  .bm_cat()
}

# ---------------------------------------------------------------------------- #

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


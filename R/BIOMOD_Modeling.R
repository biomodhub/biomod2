###################################################################################################
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
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{ANN}, \code{CTA}, \code{DNN}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' @param models.pa (\emph{optional, default} \code{NULL}) \cr 
##' A \code{list} containing for each model a \code{vector} defining which pseudo-absence datasets 
##' are to be used, must be among \code{colnames(bm.format@PA.table)}
##' 
##' @param CV.strategy (\emph{optional, default} \code{NULL}) \cr
##' A \code{character} corresponding to the cross-validation selection strategy, 
##' must be among \code{random}, \code{kfold}, \code{block}, \code{strat}, \code{env} or 
##' \code{user.defined}
##' @param CV.nb.rep (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'random'} or \code{strategy = 'kfold'}, an \code{integer} corresponding 
##' to the number of sets (repetitions) of cross-validation points that will be drawn
##' @param CV.perc (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'random'}, a \code{numeric} between \code{0} and \code{1} defining the 
##' percentage of data that will be kept for calibration
##' @param CV.k (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'kfold'} or \code{strategy = 'strat'} or \code{strategy = 'env'}, an 
##' \code{integer} corresponding to the number of partitions 
##' @param CV.balance (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'strat'} or \code{strategy = 'env'}, a \code{character} corresponding 
##' to how data will be balanced between partitions, must be either \code{presences} or
##' \code{absences} 
##' @param CV.env.var (\emph{optional, default} \code{NULL}) \cr 
##' If \code{strategy = 'env'}, a \code{character} corresponding to the environmental variables 
##' used to build the partition (all available variables by default), and for which \code{CV.k} 
##' partitions will be built
##' @param CV.strat (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'strat'}, a \code{character} corresponding to how data will partitioned 
##' along gradient, must be among \code{x}, \code{y}, \code{both}
##' @param CV.user.table (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'user.defined'}, a \code{matrix} or \code{data.frame} defining for each 
##' repetition (in columns) which observation lines should be used for models calibration 
##' (\code{TRUE}) and validation (\code{FALSE})
##' @param CV.do.full.models (\emph{optional, default} \code{NULL}) \cr  
##' A \code{logical} value defining whether models should be also calibrated and validated over 
##' the whole dataset (and pseudo-absence datasets) or not
##' 
##' @param OPT.strategy (\emph{default} \code{'default'}) \cr
##' A \code{character} corresponding to the method to select models' parameters values, must be 
##' either \code{default}, \code{bigboss}, \code{user.defined}, \code{tuned}
##' @param OPT.user.val (\emph{optional, default} \code{NULL}) \cr
##' A \code{list} containing parameters values for some (all) models
##' @param OPT.user.base (\emph{optional, default} \code{NULL}) \cr A character, 
##' \code{default} or \code{bigboss} used when \code{OPT.strategy = 'user.defined'}. 
##' It sets the bases of parameters to be modified by user defined values.
##' @param OPT.user (\emph{optional, default} \code{NULL}) \cr  
##' A \code{\link{BIOMOD.models.options}} object returned by the \code{\link{bm_ModelingOptions}} 
##' function
##' 
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{AUCroc}, \code{AUCprg}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{OR}, \code{ORSS}, 
##' \code{BOYCE}, \code{MPA} (\emph{binary data}), 
##' \code{RMSE}, \code{MAE}, \code{MSE}, \code{Rsquared}, \code{Rsquared_aj}, \code{Max_error} 
##' (\emph{abundance / count / relative data}), 
##' \code{Accuracy}, \code{Recall}, \code{Precision}, \code{F1} (\emph{multiclass / ordinal data})
##' @param var.import (\emph{optional, default} \code{0}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' 
##' @param weights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to observation weights (one per 
##' observation, see Details)
##' @param prevalence (\emph{optional, default} \code{0.5}) \cr 
##' A \code{numeric} between \code{0} and \code{1} corresponding to the species prevalence to 
##' build '\emph{weighted response weights}' (see Details)
##' @param scale.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions should be scaled with a 
##' binomial GLM or not
##' 
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
##' A \code{\link{BIOMOD.models.out}} object containing models outputs, or links to saved outputs. \cr
##' Models outputs are stored out of \R (for memory storage reasons) in 2 different folders 
##' created in the current working directory :
##' \enumerate{
##'   \item a \emph{models} folder, named after the \code{resp.name} argument of 
##'   \code{\link{BIOMOD_FormatingData}}, and containing all calibrated models for each 
##'   repetition and pseudo-absence run
##'   \item a \emph{hidden} folder, named \code{.BIOMOD_DATA}, and containing outputs related 
##'   files (original dataset, calibration lines, pseudo-absences selected, predictions, 
##'   variables importance, evaluation values...), that can be retrieved with 
##'   \href{https://biomodhub.github.io/biomod2/reference/getters.out.html}{\code{get_[...]}} 
##'   or \code{\link{load}} functions, and used by other \pkg{biomod2} functions, like 
##'   \code{\link{BIOMOD_Projection}} or \code{\link{BIOMOD_EnsembleModeling}}
##' }
##' 
##' 
##' @details 
##' 
##' \describe{
##'   \item{bm.format}{If pseudo-absences have been added to the original dataset (see 
##'   \code{\link{BIOMOD_FormatingData}}), \cr \code{PA.nb.rep *(nb.rep + 1)} models will be 
##'   created.}
##'   
##'   \item{models}{The set of models to be calibrated on the data. 12 modeling techniques 
##'   are currently available :
##'   \itemize{
##'     \item \code{ANN} : Artificial Neural Network (\code{\link[nnet]{nnet}})
##'     \item \code{CTA} : Classification Tree Analysis (\code{\link[rpart]{rpart}})
##'     \item \code{DNN} : Deep Neural Network (\code{\link[cito]{cito}})
##'     \item \code{FDA} : Flexible Discriminant Analysis (\code{\link[mda]{fda}})
##'     \item \code{GAM} : Generalized Additive Model (\code{\link[gam]{gam}}, \code{\link[mgcv]{gam}} 
##'     or \code{\link[mgcv]{bam}}) \cr 
##'     (see \code{\link{bm_ModelingOptions}} for details on algorithm selection)
##'     \item \code{GBM} : Generalized Boosting Model, or usually called Boosted Regression Trees 
##'     (\code{\link[gbm]{gbm}})
##'     \item \code{GLM} : Generalized Linear Model (\code{\link[stats]{glm}})
##'     \item \code{MARS} : Multiple Adaptive Regression Splines (\code{\link[earth]{earth}})
##'     \item \code{MAXENT} : Maximum Entropy 
##'     (\href{https://biodiversityinformatics.amnh.org/open_source/maxent/}{see Maxent website})
##'     \item \code{MAXNET} : Maximum Entropy (\code{\link[maxnet]{maxnet}})
##'     \item \code{RF} : Random Forest (\code{\link[randomForest]{randomForest}})
##'     \item \code{RFd} : Random Forest downsampled (\code{\link[randomForest]{randomForest}})
##'     \item \code{SRE} : Surface Range Envelop or usually called BIOCLIM (\code{\link{bm_SRE}})
##'     \item \code{XGBOOST} : eXtreme Gradient Boosting Training (\code{\link[xgboost]{xgboost}})
##'   }
##'   \tabular{rcccccccccccccc}{
##'     \tab \strong{ANN} \tab \strong{CTA} \tab \strong{DNN} \tab \strong{FDA} \tab \strong{GAM} \tab \strong{GBM} 
##'     \tab \strong{GLM} \tab \strong{MARS} \tab \strong{MAXENT} \tab \strong{MAXNET} 
##'     \tab \strong{RF} \tab \strong{RFd} \tab \strong{SRE} \tab \strong{XGBOOST} \cr
##'    binary \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \tab x \cr
##'    multiclass \tab  \tab x \tab x \tab x \tab  \tab  \tab  \tab x \tab  \tab  \tab x \tab  \tab  \tab x \cr
##'    ordinal \tab  \tab x \tab x \tab x \tab x \tab  \tab x \tab x \tab  \tab  \tab x \tab  \tab  \tab x \cr
##'    abundance / count / relative \tab  \tab x \tab x \tab  \tab x \tab x \tab x \tab x \tab  \tab  \tab x \tab  \tab  \tab x 
##'   }
##'   }
##'   
##'   \item{models.pa}{Different models might respond differently to different numbers of 
##'   pseudo-absences. It is possible to create sets of pseudo-absences with different numbers 
##'   of points (see \code{\link{BIOMOD_FormatingData}}) and to assign only some of these 
##'   datasets to each single model.
##'   }
##'   
##'   \item{CV.[...] parameters}{Different methods are available to calibrate/validate the 
##'   single models (see \code{\link{bm_CrossValidation}}).}
##'   
##'   \item{OPT.[...] parameters}{Different methods are available to parameterize the 
##'   single models (see \code{\link{bm_ModelingOptions}} and 
##'   \code{\link{BIOMOD.options.dataset}}). 
##'   \itemize{
##'     \item \code{default} : only default parameter values of default parameters of the single 
##'     models functions are retrieved. Nothing is changed so it might not give good results.
##'     \item \code{bigboss} : uses parameters pre-defined by \pkg{biomod2} team and that are 
##'     available in the dataset \code{\link{OptionsBigboss}}. \cr 
##'     \emph{to be optimized in near future}
##'     \item \code{user.defined} : updates default or bigboss parameters with some parameters 
##'     values defined by the user (but matching the format of a 
##'     \code{\link{BIOMOD.models.options}} object)
##'     \item \code{tuned} : calling the \code{\link{bm_Tuning}} function to try and optimize 
##'     some default values
##'   }
##'   }
##' 
##'   \item{metric.eval}{
##'   \emph{Please refer to  
##'   \href{https://www.cawcr.gov.au/projects/verification/}{CAWRC website ("Methods for 
##'   dichotomous forecasts")} to get detailed description (simple/complex metrics).} \cr
##'   Several evaluation metrics can be selected. \cr
##'   Optimal value of each method can be obtained with the \code{\link{get_optim_value}} 
##'   function.
##'   \describe{
##'     \item{simple}{
##'     \itemize{
##'       \item \code{POD} : Probability of detection (hit rate)
##'       \item \code{FAR} : False alarm ratio
##'       \item \code{POFD} : Probability of false detection (fall-out)
##'       \item \code{SR} : Success ratio
##'       \item \code{ACCURACY} : Accuracy (fraction correct)
##'       \item \code{BIAS} : Bias score (frequency bias)
##'     }
##'     }
##'     \item{complex}{
##'     \itemize{
##'       \item \code{AUCroc} : Area Under Curve of Relative operating characteristic
##'       \item \code{AUCprg} : Area Under Curve of Precision-Recall-Gain curve
##'       \item \code{TSS} : True skill statistic (Hanssen and Kuipers discriminant, Peirce's 
##'       skill score)
##'       \item \code{KAPPA} : Cohen's Kappa (Heidke skill score)
##'       \item \code{OR} : Odds Ratio
##'       \item \code{ORSS} : Odds ratio skill score (Yule's Q)
##'       \item \code{CSI} : Critical success index (threat score)
##'       \item \code{ETS} : Equitable threat score (Gilbert skill score)
##'     }
##'     }
##'     \item{presence-only}{
##'     \itemize{
##'       \item \code{BOYCE} : Boyce index
##'       \item \code{MPA} : Minimal predicted area (cutoff optimizing MPA to predict 90\% of 
##'       presences)
##'     }
##'     }
##'     \item{abundance / count / relative data}{
##'     \itemize{
##'       \item \code{RMSE} : Root Mean Square Error
##'       \item \code{MSE} : Mean Square Error
##'       \item \code{MAE} : Mean Absolute Error
##'       \item \code{Rsquared} : R squared
##'       \item \code{Rsquared_aj} : R squared adjusted
##'       \item \code{Max_error} : Maximum error
##'     }
##'     }
##'     \item{multiclass/ordinal data}{
##'     \itemize{
##'       \item \code{Accuracy} : Accuracy
##'       \item \code{Recall} : Macro average Recall
##'       \item \code{Precision} : Macro average Precision
##'       \item \code{F1} : Macro F1 score
##'     }
##'     }
##'   }
##'   Results after modeling can be obtained through the \code{\link{get_evaluations}} function. \cr 
##'   Evaluation metric are calculated on the calibrating data (column \code{calibration}), on 
##'   the cross-validation data (column \code{validation}) or on the evaluation data 
##'   (column \code{evaluation}). \cr 
##'   \emph{For cross-validation data, see \code{CV.[...]} parameters in 
##'   \code{\link{BIOMOD_Modeling}} function. \cr For evaluation data, see 
##'   \code{eval.[...]} parameters in \code{\link{BIOMOD_FormatingData}}.}
##'   }
##'   
##'   \item{var.import}{A value characterizing how much each variable has an impact on each model 
##'   predictions can be calculated by randomizing the variable of interest and computing the 
##'   correlation between original and shuffled variables (see \code{\link{bm_VariablesImportance}}).}
##'   
##'   \item{weights & prevalence}{
##'   More or less weight can be given to some specific observations. \cr Automatically created 
##'   \code{weights} will be \code{integer} values to prevent some modeling issues. \cr 
##'   \emph{Note that \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd} and \code{SRE} models 
##'   do not take weights into account.}
##'   \itemize{
##'     \item If \code{prevalence = 0.5} (the default), presences and absences will be weighted equally 
##'     (\emph{i.e. the weighted sum of presences equals the weighted sum of absences}). 
##'     \item If \code{prevalence} is set below (\emph{above}) \code{0.5}, more weight will be 
##'     given to absences (\emph{presences}).
##'     \item If \code{weights} is defined, \code{prevalence} argument will be ignored 
##'     (\emph{EXCEPT for \code{MAXENT}}).
##'   }}
##'   
##'   \item{scale.models}{A binomial GLM is created to scale predictions from 0 to 1. \cr
##'   \code{SRE} is never scaled, and \code{ANN} and \code{FDA} categorical models always are. \cr
##'   \emph{Note that it may lead to reduction in projected scale amplitude.} \cr
##'   \bold{This parameter is quite experimental and it is recommended not to use it.} It was 
##'   developed in the idea to ensure comparable predictions by removing the scale prediction 
##'   effect (\emph{the more extended projections are, the more they influence ensemble 
##'   forecasting results}).
##'   }
##' }
##' 
##' 
##' @keywords models regression nonlinear multivariate nonparametric tree
##' 
##' 
##' @seealso \code{\link[stats]{glm}}, \code{\link[gam]{gam}},
##'   \code{\link[mgcv]{gam}}, \code{\link[mgcv]{bam}}, \code{\link[gbm]{gbm}},
##'   \code{\link[rpart]{rpart}}, \code{\link[nnet]{nnet}}, \code{\link[cito]{cito}},
##'   \code{\link[mda]{fda}}, \code{\link[earth]{earth}},
##'   \code{\link[randomForest]{randomForest}}, \code{\link[maxnet]{maxnet}},
##'   \code{\link[xgboost]{xgboost}}, \code{\link{BIOMOD_FormatingData}},
##'   \code{\link{bm_ModelingOptions}}, \code{\link{bm_Tuning}}, 
##'   \code{\link{bm_CrossValidation}},
##'   \code{ \link{bm_VariablesImportance}}, \code{\link{BIOMOD_Projection}},
##'   \code{\link{BIOMOD_EnsembleModeling}}, \code{\link{bm_PlotEvalMean}},
##'   \code{\link{bm_PlotEvalBoxplot}}, \code{\link{bm_PlotVarImpBoxplot}},
##'   \code{\link{bm_PlotResponseCurves}}
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
##' myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
##'                                      resp.var = myResp,
##'                                      resp.xy = myRespXY,
##'                                      expl.var = myExpl)
##' 
##' 
##' # ---------------------------------------------------------------------------- #
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     models = c('RF', 'GLM'),
##'                                     CV.strategy = 'random',
##'                                     CV.nb.rep = 2,
##'                                     CV.perc = 0.8,
##'                                     OPT.strategy = 'bigboss',
##'                                     metric.eval = c('TSS','AUCroc'),
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
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))
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
##' @importFrom foreach foreach %do%
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_Modeling <- function(bm.format,
                            modeling.id,
                            models = c('ANN', 'CTA', 'DNN', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS'
                                       , 'MAXENT', 'MAXNET', 'RF', 'RFd', 'SRE', 'XGBOOST'),
                            models.pa = NULL,
                            CV.strategy = NULL,
                            CV.nb.rep = NULL,
                            CV.perc = NULL,
                            CV.k = NULL,
                            CV.balance = NULL,
                            CV.env.var = NULL,
                            CV.strat = NULL,
                            CV.user.table = NULL,
                            CV.do.full.models = NULL,
                            OPT.strategy = 'default',
                            OPT.user.val = NULL,
                            OPT.user.base = NULL,
                            OPT.user = NULL,
                            metric.eval = NULL,
                            var.import = 0,
                            weights = NULL,
                            prevalence = 0.5,
                            scale.models = FALSE,
                            nb.cpu = 1,
                            seed.val = NULL,
                            do.progress = TRUE)
{
  .bm_cat("[BIOMOD] Build Single Models")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  cat("\nChecking arguments...")
  args <- .BIOMOD_Modeling.check.args(
    bm.format = bm.format, 
    modeling.id = modeling.id, 
    models = models, 
    models.pa = models.pa, 
    OPT.user = OPT.user,
    CV.strategy = CV.strategy,
    CV.user.table = CV.user.table,
    CV.do.full.models = CV.do.full.models,
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
  cat("\n")
  
  ## 1. Create output object ----------------------------------------------------------------------
  models.out <- new('BIOMOD.models.out',
                    dir.name = bm.format@dir.name,
                    sp.name = bm.format@sp.name,
                    modeling.id = modeling.id,
                    data.type = bm.format@data.type,
                    expl.var.names = colnames(bm.format@data.env.var),
                    has.evaluation.data = bm.format@has.data.eval,
                    scale.models = scale.models)
  
  ## 2. Create simulation directories -------------------------------------------------------------
  ## Various objects will be stored (models, predictions, projections)
  ## Projections directories are created in Projection() function
  .BIOMOD_Modeling.prepare.workdir(bm.format@dir.name, bm.format@sp.name, models.out@modeling.id)
  name.BIOMOD_DATA <- file.path(models.out@dir.name, models.out@sp.name, ".BIOMOD_DATA", models.out@modeling.id)
  
  ## 3.1 Save input data ------------------------------------------------------
  models.out <- .fill_BIOMOD.models.out("formated.input.data", bm.format, models.out
                                        , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  
  ## 3.2 Get and save calibration lines ---------------------------------------
  calib.lines <- bm_CrossValidation(bm.format = bm.format,
                                    strategy = CV.strategy,
                                    nb.rep = CV.nb.rep,
                                    perc = CV.perc,
                                    k = CV.k,
                                    balance = ifelse(!is.null(CV.balance), CV.balance, "presences"),
                                    env.var = CV.env.var,
                                    strat = ifelse(!is.null(CV.strat), CV.strat, "both"),
                                    user.table = CV.user.table,
                                    do.full.models = CV.do.full.models)
  models.out <- .fill_BIOMOD.models.out("calib.lines", calib.lines, models.out
                                        , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  
  ## 3.3 Get and save models options ------------------------------------------
  if (!is.null(OPT.user)) {
    ## Check for model names -----------
    if (sum(!(models %in% sapply(OPT.user@models, function(x) strsplit(x, "[.]")[[1]][1]))) > 0) {
      stop("OPT.user must contain information for ", toString(models), " models. Please check.")
    }
    ## Check data.type coherence 
    data.type.options <- strsplit(OPT.user@models[1],".", fixed = TRUE)[[1]][2]
    if ((bm.format@data.type == "binary" & data.type.options != "binary") ||
        (bm.format@data.type != "binary" & data.type.options == "binary")) {
      stop("data.type of OPT.user should match bm.format@data.type. Please check.")
    }
    ## Check for calib.lines names -----
    for (mod in OPT.user@models) {
      nam <- names(OPT.user@options[[mod]]@args.values)
      vals <- colnames(calib.lines)
      if (any(!(vals %in% nam))) {
        if (length(nam) == 1 && nam == "_allData_allRun") {
          val <- OPT.user@options[[mod]]@args.values[[nam]]
          for (ii in vals) {
            OPT.user@options[[mod]]@args.values[[ii]] <- val
          }
        } else if (all(grepl("_allRun", nam))) { #Check if the user create the options just for PA dataset
          sep.name <- unlist(strsplit(mod, split = '[.]'))
          .message(sep.name[1], " options for '_PAx_allRun' will be given to all PAx runs ('_PAx_RUN1', '_PAx_RUN2', ...).")
          for (run in vals) {
            PA.set <- grep("PA|allData", unlist(strsplit(run, "_")), value = TRUE)
            if (is.null(OPT.user@options[[mod]]@args.values[[paste0("_", PA.set, "_allRun")]])) {
              opt.default <- BIOMOD.options.dataset(mod = sep.name[1], typ = sep.name[2], pkg = sep.name[3], fun = sep.name[4]
                                                    , strategy = "default", bm.format = bm.format, calib.lines = calib.lines)
              OPT.user@options[[mod]]@args.values[[run]] <- opt.default@args.values[["_allData_allRun"]]
            } else {
              OPT.user@options[[mod]]@args.values[[run]] <- OPT.user@options[[mod]]@args.values[[paste0("_", PA.set, "_allRun")]]
            }
            
          }
        } else {
          .fun_testIfIn(paste0("names(OPT.user@options[['", mod, "']]@args.values)"), nam, vals, exact = TRUE)
          stop("names(OPT.user@options[['", mod, "']]@args.values) must be ", toString(vals))
        }
      }
    }
    bm.options <- OPT.user
  } else {
    bm.options <- bm_ModelingOptions(data.type = bm.format@data.type,
                                     models = models,
                                     strategy = OPT.strategy,
                                     user.val = OPT.user.val,
                                     user.base = OPT.user.base,
                                     bm.format = bm.format,
                                     calib.lines = calib.lines)
  }
  if (!is.null(prevalence) && "MAXENT" %in% models) {
    for (nam in names(bm.options@options$MAXENT.binary.MAXENT.MAXENT@args.values)) {
      bm.options@options$MAXENT.binary.MAXENT.MAXENT@args.values[[nam]][['defaultprevalence']] <- prevalence
    }
  }
  models.out <- .fill_BIOMOD.models.out("models.options", bm.options, models.out
                                        , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
  
  ## 4. Print modeling summary in console ---------------------------------------------------------
  .BIOMOD_Modeling.summary(bm.format, calib.lines, models, models.pa)
  
  
  ## 5. Run models with loop over PA --------------------------------------------------------------
  mod.out <- bm_RunModelsLoop(bm.format = bm.format,
                              modeling.id = models.out@modeling.id,
                              models = models,
                              models.pa = models.pa,
                              calib.lines = calib.lines,
                              bm.options = bm.options,
                              metric.eval = metric.eval,
                              var.import = var.import,
                              weights = weights,
                              scale.models = scale.models,
                              nb.cpu = nb.cpu,
                              seed.val = seed.val,
                              do.progress = do.progress)
  
  ## 6. Rearrange and save outputs -------------------------------------------
  models.out@models.computed <- .transform_outputs_list("mod", mod.out, out = "model")
  models.out@models.failed <- .transform_outputs_list("mod", mod.out, out = "calib.failure")
  
  if (length(models.out@models.computed) == 1 && models.out@models.computed == "none") {
    .message("*** Error: all models failed")
    return(models.out)
  }
  
  ## 3.4 Rearrange and save models outputs : ----------------------------------
  ## models evaluation, variables importance, models prediction, predictions evaluation
  if (length(metric.eval) > 0) {
    models.evaluation <- .transform_outputs_list("mod", mod.out, out = "evaluation")
    models.out <- .fill_BIOMOD.models.out("models.evaluation", models.evaluation, models.out
                                          , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
    rm(models.evaluation)
  }
  
  if (var.import > 0) {
    variables.importance <- .transform_outputs_list("mod", mod.out, out = "var.import")
    models.out <- .fill_BIOMOD.models.out("variables.importance", variables.importance, models.out
                                          , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
    rm(variables.importance)
  }
  
  models.prediction <- .transform_outputs_list("mod", mod.out, out = "pred")
  models.out <- .fill_BIOMOD.models.out("models.prediction", models.prediction, models.out
                                        , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  rm(models.prediction)
  
  models.prediction.eval <- .transform_outputs_list("mod", mod.out, out = "pred.eval")
  models.out <- .fill_BIOMOD.models.out("models.prediction.eval", models.prediction.eval, models.out
                                        , inMemory = FALSE, nameFolder = name.BIOMOD_DATA)
  rm(models.prediction.eval)
  rm(mod.out)
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE ----------------------------
  name.OUT <- paste0(models.out@sp.name, '.', models.out@modeling.id, '.models.out')
  models.out@link <- file.path(models.out@dir.name, models.out@sp.name, name.OUT)
  models.out@call <- match.call()
  assign(x = name.OUT, value = models.out)
  save(list = name.OUT, file = models.out@link)
  
  # if (.getOS() == "windows" && "MAXENT" %in% models){
  #   env <- foreach:::.foreachGlobals
  #   rm(list=ls(name=env), pos=env)
  # }
  
  .bm_cat("Done")
  return(models.out)
}


###################################################################################################

.BIOMOD_Modeling.check.args <- function(bm.format, modeling.id, models, models.pa, OPT.user
                                        , CV.strategy, CV.user.table, CV.do.full.models
                                        , weights, prevalence, metric.eval, var.import
                                        , scale.models, nb.cpu, seed.val, do.progress)
{
  ## 0. Check bm.format argument ----------------------------------------------
  .fun_testIfInherits("bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  
  ## 1. Check modeling.id argument --------------------------------------------
  if (missing(modeling.id)) {
    modeling.id <- as.character(format(Sys.time(), "%s"))
    .message("modeling.id set to ", modeling.id)
  }
  modeling.id <- as.character(modeling.id)
  .fun_testIfLength("modeling.id", modeling.id)
  
  ## 2. Check models argument -------------------------------------------------
  models <- unique(as.character(models))
  avail.models.list <- .avail.models.list(bm.format@data.type)
  .fun_testIfIn(paste0("models with ", bm.format@data.type, " data type"), models, avail.models.list)
  
  ## 3.a Specific case of cito
  if ("DNN" %in% models) {
    if (!requireNamespace('torch', quietly = TRUE)) stop("Package 'torch' not found")
  }
  
  ## 3.b Specific case of one variable with GBM / MAXNET
  if ('GBM' %in% models && ncol(bm.format@data.env.var) == 1) {
    .message("GBM might have issues when only one variable is used. "
             , "Please be sure to install the following version : devtools::install_github('rpatin/gbm')")
  }
  if ('MAXNET' %in% models && ncol(bm.format@data.env.var) == 1) {
    .message("MAXNET might have issues when only one variable is used. "
             , "Please be sure to install the following version : devtools::install_github('mrmaxent/maxnet')")
  }
  
  ## 3.c Remove models not supporting categorical variables
  categorical_var <- .get_categorical_names(bm.format@data.env.var)
  models.switch.off <- NULL
  if (length(categorical_var) > 0) {
    models.fact.unsupport <- c("SRE")
    models.switch.off <- c(models.switch.off, intersect(models, models.fact.unsupport))
    if (length(models.switch.off) > 0) {
      models <- setdiff(models, models.switch.off)
      .message(toString(models.switch.off), "switched off (categorical variables)")
    }
  }
  
  ## 4. Check models.pa argument ----------------------------------------------
  if (!is.null(models.pa)) {
    if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
      .fun_testIfInherits("models.pa", models.pa, "list")
      .fun_testIfIn("unlist(models.pa)", unlist(models.pa), colnames(bm.format@PA.table))
      .fun_testIfIn("names(models.pa)", names(models.pa), models)
      if (length(models.pa) != length(models)) {
        mod.miss <- models[-which(models %in% names(models.pa))]
        list.miss <- rep(list(colnames(bm.format@PA.table)), length(mod.miss))
        names(list.miss) <- mod.miss
        models.pa <- c(models.pa, list.miss)
        .message("No models.pa provided for ", toString(mod.miss)
                 , ". They have been assigned to all datasets.")
      }
    } else {
      models.pa <- NULL
      .message("models.pa set to NULL (no PA dataset provided)")
    }
  }
  
  ## 5. Check OPT.user argument -----------------------------------------------
  if (!is.null(OPT.user)) {
    .fun_testIfInherits("OPT.user", OPT.user, "BIOMOD.models.options")
    .message("OPT.user provided, all other OPT.[..] arguments will be ignored.")
  } 
  # else {
  # warning("OPT.user set to BIOMOD_ModelingOptions()", immediate. = TRUE)
  #   bm.options <- BIOMOD_ModelingOptions()
  # }
  
  ## 6. Check CV.user.table / CV.do.full.models arguments ---------------------
  if (is.null(CV.do.full.models)) {
    CV.do.full.models <- FALSE
    .message("CV.do.full.models set to FALSE (no '_allData_allRun' set computed)")
  }
  if (CV.strategy == "user.defined" && !is.null(CV.user.table)) {
    if (!("_allData_allRun" %in% colnames(CV.user.table)) && CV.do.full.models == TRUE) { 
      CV.do.full.models <- FALSE
      .message("CV.do.full.models set to FALSE (no '_allData_allRun' column provided in CV.user.table)")
    }
  } else if (missing(CV.strategy) || is.null(CV.strategy)) {
    CV.do.full.models <- FALSE
  } else if (missing(CV.do.full.models) || is.null(CV.do.full.models)) {
    CV.do.full.models <- FALSE
  }
  
  ## 7. Check prevalence argument ---------------------------------------------
  if (!is.null(prevalence)) {
    .fun_testIf0X("prevalence", prevalence, 1)
  } else {
    prevalence <- 0.5
    .message("prevalence set to 0.5")
  }
  
  ## 7. Check weights argument ------------------------------------------------
  if (is.null(weights)) {
    if (!is.null(prevalence) && !(bm.format@data.type %in% c("ordinal", "multiclass"))) {
      cat("\n > Automatic weights creation to rise a", prevalence, "prevalence")
      data.sp <- as.numeric(bm.format@data.species)
      if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
        weights.pa <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
          {
            ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
            data.sp_pa <- data.sp[ind.PA]
            data.sp_pa[which(is.na(data.sp_pa))] <- 0
            weights <- .automatic_weights_creation(resp = data.sp_pa, prev = prevalence)
            
            wei <- rep(NA, length(data.sp))
            wei[ind.PA] <- weights
            return(matrix(wei, ncol = 1))
          }
        weights.pa <- cbind(weights.pa, rep(1, nrow(weights.pa)))
        colnames(weights.pa) <- c(colnames(bm.format@PA.table), "allData")
        weights <- weights.pa
      } else {
        weights <- .automatic_weights_creation(resp = data.sp, prev = prevalence)
        weights <- matrix(weights, nrow = length(weights), ncol = 1)
        colnames(weights) <- "allData"
      }
    } else { ## NEVER OCCURRING NO ?? --> now happen with the abundance
      .message("All observations will have the same weight (no weights provided).")
    }
  } else {
    .fun_testIfPosNum("weights", weights)
    .fun_testIfSameSize("weights", length(weights), "bm.format@data.species", length(bm.format@data.species))
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
  
  ## 8. Check metric.eval argument --------------------------------------------
  metric.eval <- unique(metric.eval)
  avail.eval.meth.list <- .avail.eval.meth.list(bm.format@data.type)
  .fun_testIfIn(paste0("metric.eval with ", bm.format@data.type, " data type"), metric.eval, avail.eval.meth.list)
  
  ## 9. Check var.import argument ---------------------------------------------
  if (is.null(var.import)) { var.import = 0 }
  
  ## 10. Set the seed (if needed) ----------------------------------------------
  if (!is.null(seed.val)) { set.seed(seed.val) }
  
  return(list(modeling.id = modeling.id,
              models = models,
              models.pa = models.pa,
              CV.do.full.models = CV.do.full.models,
              prevalence = prevalence,
              weights = weights,
              metric.eval = metric.eval,
              var.import = var.import,
              seed.val = seed.val,
              do.progress = do.progress))
}


###################################################################################################

.BIOMOD_Modeling.prepare.workdir <- function(dir.name, sp.name, modeling.id)
{
  cat(" > Creating modeling directory\n")
  dir.create(file.path(dir.name, sp.name), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, ".BIOMOD_DATA", modeling.id), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(dir.name, sp.name, "models", modeling.id), showWarnings = FALSE, recursive = TRUE)
}

# ----------------------------------------------------------------------------------------------- #

.BIOMOD_Modeling.summary <- function(bm.format, calib.lines, models, models.pa = NULL)
{
  .bm_cat(paste(bm.format@sp.name, "Modeling Summary"))
  cat("\nEnvironmental variables :", ncol(bm.format@data.env.var), "(", colnames(bm.format@data.env.var), ")")
  
  if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
    cat("\nNumber of PA datasets :", ncol(bm.format@PA.table))
    nb.PA <- ncol(bm.format@PA.table)
  } else {
    nb.PA <- 0
  }
  
  if ("_allData_allRun" %in% colnames(calib.lines)) {
    nb.full.models <- nb.PA + 1
    nb.eval.rep <- (ncol(calib.lines) - nb.full.models) / ifelse(inherits(bm.format, "BIOMOD.formated.data.PA"), ncol(bm.format@PA.table), 1)
    cat("\nNumber of calibration/validation splits :", nb.eval.rep)
    cat("\n\t > CV.do.full.models activated : +", nb.full.models,  "model", ifelse(nb.full.models > 1, "s", ""))
    nb.eval.rep <- nb.eval.rep +1 #for models.pa
  } else {
    nb.eval.rep <- ncol(calib.lines) / ifelse(inherits(bm.format, "BIOMOD.formated.data.PA"), ncol(bm.format@PA.table), 1)
    cat("\nNumber of calibration/validation splits :", nb.eval.rep)
  }
  
  if (is.null(models.pa)) {
    nb.runs <- ncol(calib.lines) * length(models)
    cat("\n\t >", ncol(calib.lines), ifelse(ncol(calib.lines) > 1, "models", "model"), "for each algorithm")
  } else {
    nb.runs <- length(which(
      sapply(unlist(models.pa), function(x) grepl(colnames(calib.lines), pattern = x))
    ))
    for (algo in names(models.pa)){
      cat("\n\t >", algo, ":", length(models.pa[[algo]]) * nb.eval.rep , "models")
    }
  }
  
  cat("\nAlgorithms selected :", length(models), "(", models, ")")
  cat("\n\nTotal number of model runs :", nb.runs, "\n")
  .bm_cat()
}

# bm_ModelingOptions documentation -----------------------------------------
##' @name bm_ModelingOptions
##' @aliases bm_ModelingOptions
##' @author Damien Georges, Wilfried Thuiller
##' @author Maya Gueguen
##' 
##' @title Configure the modeling options for each selected model
##'
##' @description Parameterize and/or tune \pkg{biomod2}'s single models options.
##'
##' @param data.type a \code{character} corresponding to the data type to 
##' be used, must be either \code{binary}, \code{binary.PA}, \code{abundance}, 
##' \code{compositional}
##' @param models a \code{vector} containing model names to be computed, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT}, \code{MAXNET}, \code{XGBOOST}
##' @param strategy a \code{character} corresponding to the method to 
##' select models' parameters values, must be either \code{default}, 
##' \code{bigboss}, \code{user.defined}, \code{tuned}
##' @param val.list (\emph{optional, default} \code{NULL}) \cr
##' A \code{list} containing parameters values for some (all) models
##' @param bm.format (\emph{optional, default} \code{NULL}) \cr
##' A \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' A \code{data.frame} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{bm_CrossValidation}} functions, to explore the distribution of calibration 
##' and validation datasets
##' 
##'
##'
##' @return 
##' 
##' A \code{\link{BIOMOD.models.options}} of object that can be used to build 
##' species distribution model(s) with the \code{\link{BIOMOD_Modeling}} function.
##' 
##' 
##' @details
##' 
##' This function allows advanced user to change some default parameters of \pkg{biomod2} inner 
##' models. \cr 10 single models are available within the package, and their options can be set 
##' with this function through \code{list} objects.
##' 
##' The \code{\link{bm_DefaultModelingOptions}} function prints all default parameter values for 
##' all available models. \cr This output can be copied and pasted to be used as is (with wanted 
##' changes) as function arguments (see \href{https://biomodhub.github.io/biomod2/reference/bm_ModelingOptions.html#examples}{Examples}).
##' 
##' Below is the detailed list of all modifiable parameters for each available model.
##'
##'
##' @section ANN : (\code{\link[nnet]{nnet}})
##' \itemize{
##'   \item{\code{NbCV = 5}}{ : an \code{integer} corresponding to the number of cross-validation 
##'   repetitions to find best size and decay parameters}
##'   \item{\code{size = NULL}}{ : an \code{integer} corresponding to the number of units in the 
##'   hidden layer. If \code{NULL} then size parameter will be optimized by cross-validation based 
##'   on model AUC (\code{NbCv} cross-validations ; tested size will be the following : 
##'   \code{c(2, 4, 6, 8)}). It is also possible to give a \code{vector} of size values to be tested, 
##'   and the one giving the best model AUC will be kept.}
##'   \item{\code{decay = NULL}}{ : a \code{numeric} corresponding to weight decay. If \code{NULL} 
##'   then decay parameter will be optimized by cross-validation based on model AUC (\code{NbCv} 
##'   cross-validations ; tested size will be the following : \code{c(0.001, 0.01, 0.05, 0.1)}). 
##'   It is also possible to give a \code{vector} of decay values to be tested, and the one giving 
##'   the best model AUC will be kept.}
##'   \item{\code{rang = 0.1}}{ : a \code{numeric} corresponding to the initial random weights on 
##'   \code{[-rang, rang]}}
##'   \item{\code{maxit = 200}}{ : an \code{integer} corresponding to the maximum number of 
##'   iterations}
##' }
##' @section CTA : (\code{\link[rpart]{rpart}})
##' @section FDA : (\code{\link[mda]{fda}})
##' @section GAM : (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}})
##' (see \code{\link[gam]{gam}}, \code{\link[mgcv]{gam}}, \code{\link[mgcv]{bam}})
##' @section GBM : (default \code{\link[gbm]{gbm}})
##' @section GLM : (\code{\link[stats]{glm}})
##' @section MARS : (\code{\link[earth]{earth}})
##' @section MAXENT : (\url{https://biodiversityinformatics.amnh.org/open_source/maxent/})
##' \itemize{
##'   \item{\code{path_to_maxent.jar = getwd()}}{ : a \code{character}
##'   corresponding to \pkg{maxent.jar} file link} 
##'   
##'   \item{\code{memory_allocated = 512}}{ : an \code{integer} corresponding to
##'   the amount of memory (in Mo) reserved for \code{java} to run
##'   \code{MAXENT}, must be \code{64}, \code{128}, \code{256},
##'   \code{512}, \code{1024}... or \code{NULL} to use default \code{java}
##'   memory limitation parameter}
##'   
##'   \item{\code{initial_heap_size = NULL}}{ : a \code{character} initial heap
##'   space (shared memory space) allocated to java. Argument transmitted to
##'   \code{-Xms} when calling java. Used in \code{\link{BIOMOD_Projection}} but
##'   not in \code{\link{BIOMOD_Modeling}}. Values can be \code{1024K},
##'   \code{4096M}, \code{10G} ... or \code{NULL} to use default \code{java}
##'   parameter}
##'   
##'   \item{\code{max_heap_size = NULL}}{ : a \code{character} initial heap
##'   space (shared memory space) allocated to java. Argument transmitted to
##'   \code{-Xmx} when calling java. Used in \code{\link{BIOMOD_Projection}} but
##'   not in \code{\link{BIOMOD_Modeling}}. Must be larger than
##'   \code{initial_heap_size}. Values can be \code{1024K}, \code{4096M},
##'   \code{10G} ... or \code{NULL} to use default \code{java} parameter}
##'   
##'   \item{\code{background_data_dir}}{ : a \code{character} corresponding to
##'   directory path where explanatory variables are stored as \code{ASCII}
##'   files (raster format). If specified, \code{MAXENT} will generate
##'   its own background data from explanatory variables rasters (as usually
##'   done in \code{MAXENT} studies). Otherwise \pkg{biomod2} pseudo-absences
##'   will be used (see \code{\link{BIOMOD_FormatingData}})}
##'   
##'   \item{\code{maximumbackground}}{ : an \code{integer} corresponding to the
##'   maximum number of background data to sample if the
##'   \code{background_data_dir} parameter has been set}
##'   
##'   \item{\code{maximumiterations = 200}}{ : an \code{integer} corresponding
##'   to the maximum number of iterations to do} 
##'   
##'   \item{\code{visible = FALSE}}{ : a \code{logical} to make the
##'   \code{MAXENT} user interface available}
##'   
##'   \item{\code{linear = TRUE}}{ : a \code{logical} to allow linear features
##'   to be used} \item{\code{quadratic = TRUE}}{ : a \code{logical} to allow
##'   quadratic features to be used} 
##'   
##'   \item{\code{product = TRUE}}{ : a \code{logical} to allow product features
##'   to be used}
##'   
##'   \item{\code{threshold = TRUE}}{ : a \code{logical} to allow threshold
##'   features to be used}
##'   
##'   \item{\code{hinge = TRUE}}{ : a \code{logical} to allow hinge features to
##'   be used} \item{\code{lq2lqptthreshold = 80}}{ : an \code{integer}
##'   corresponding to the number of samples at which product and threshold
##'   features start being used} 
##'   
##'   \item{\code{l2lqthreshold = 10}}{ : an
##'   \code{integer} corresponding to the number of samples at which quadratic
##'   features start being used} 
##'   
##'   \item{\code{hingethreshold = 15}}{ : an
##'   \code{integer} corresponding to the number of samples at which hinge
##'   features start being used}
##'   
##'   \item{\code{beta_threshold = -1.0}}{ : a
##'   \code{numeric} corresponding to the regularization parameter to be applied
##'   to all threshold features (\emph{negative value enables automatic
##'   setting})} 
##'   
##'   \item{\code{beta_categorical = -1.0}}{ : a \code{numeric}
##'   corresponding to the regularization parameter to be applied to all
##'   categorical features (\emph{negative value enables automatic setting})}
##'   
##'   \item{\code{beta_lqp = -1.0}}{ : a \code{numeric} corresponding to the
##'   regularization parameter to be applied to all linear, quadratic and
##'   product features (\emph{negative value enables automatic setting})}
##'   
##'   \item{\code{beta_hinge = -1.0}}{ : a \code{numeric} corresponding to the
##'   regularization parameter to be applied to all hinge features
##'   (\emph{negative value enables automatic setting})}
##'   
##'   \item{\code{betamultiplier = 1}}{ : a \code{numeric} to multiply all
##'   automatic regularization parameters \cr (\emph{higher number gives a more
##'   spread-out distribution})} 
##'   
##'   \item{\code{defaultprevalence = 0.5}}{ : a
##'   \code{numeric} corresponding to the default prevalence of the species \cr
##'   (\emph{probability of presence at ordinary occurrence points})}
##' }
##' @section MAXNET : (\code{\link[maxnet]{maxnet}})
##' @section RF : (\code{\link[randomForest]{randomForest}})
##' @section SRE : (\code{\link{bm_SRE}})
##' @section XGBOOST : (default \code{\link[xgboost]{xgboost}})
##' 
##'
##' @keywords models options
##' 
##' 
##' @seealso \code{\link{bm_Tuning}}, \code{\link{BIOMOD_Modeling}}
##' @family Secundary functions
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
##' # ---------------------------------------------------------------#
##' # Print default modeling options
##' bm_DefaultModelingOptions()
##' 
##' # Create default modeling options
##' myBiomodOptions <- bm_ModelingOptions()
##' myBiomodOptions
##' 
##' # # Part (or totality) of the print can be copied and customized
##' # # Below is an example to compute quadratic GLM and select best model with 'BIC' criterium
##' # myBiomodOptions <- bm_ModelingOptions(
##' #   GLM = list(type = 'quadratic',
##' #              interaction.level = 0,
##' #              myFormula = NULL,
##' #              test = 'BIC',
##' #              family = 'binomial',
##' #              control = glm.control(epsilon = 1e-08,
##' #                                    maxit = 1000,
##' #                                    trace = FALSE)))
##' # myBiomodOptions
##' # 
##' # # It is also possible to give a specific GLM formula
##' # myForm <- 'Sp277 ~ bio3 + log(bio10) + poly(bio16, 2) + bio19 + bio3:bio19'
##' # myBiomodOptions <- bm_ModelingOptions(GLM = list(myFormula = formula(myForm)))
##' # myBiomodOptions
##'
##'
##' @importFrom methods as new validObject
##' 
##' 
##' @export
##'
##'

TABLE_MODELS <- data.frame(model = c('ANN', 'CTA', 'FDA', 'GAM', 'GAM', 'GAM', 'GBM', 'GLM'
                                     , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , type = 'binary'
                           , package = c('nnet', 'rpart', 'mda', 'gam', 'mgcv', 'mgcv', 'gbm', 'stats'
                                         , 'earth', 'MAXENT', 'maxnet', 'randomForest', 'biomod2', 'xgboost')
                           , func = c('nnet', 'rpart', 'fda', 'gam', 'bam', 'gam', 'gbm', 'glm'
                                      , 'earth', 'MAXENT', 'maxnet', 'randomForest', 'bm_SRE', 'xgboost')
                           , train = c('avNNet', 'rpart', 'fda', 'gamSpline', 'bam', 'gam', 'gbm', 'glm'
                                       , 'earth', 'ENMevaluate', '', 'rf', 'bm_SRE', 'xgbTree'))



bm_ModelingOptions <- function(data.type
                               , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                          , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                               , strategy, val.list = NULL, bm.format = NULL, calib.lines = NULL)
{
  .bm_cat("Build Modeling Options")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_ModelingOptions.check.args(data.type = data.type, models = models, strategy = strategy
                                         , val.list = val.list, bm.format = bm.format, calib.lines = calib.lines)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Get options -------------------------------------------------------------------------------
  bm.opt <- foreach (model = models, .combine = "c") %do%
    {
      ## Select model / package / function to keep
      if (length(grep("GAM", model)) == 1) {
        tab.model <- TABLE_MODELS[which(TABLE_MODELS$model == "GAM" & TABLE_MODELS$type == data.type &
                                          TABLE_MODELS$package == strsplit(model, "[.]")[[1]][2] & 
                                          TABLE_MODELS$func == strsplit(model, "[.]")[[1]][3]), ]
        model = "GAM"
      } else {
        tab.model <- TABLE_MODELS[which(TABLE_MODELS$model == model & TABLE_MODELS$type == data.type), ]
      }
      
      if (nrow(tab.model) > 0) {
        ## For each kept model : get corresponding options
        BOD.list <- foreach(ii = 1:nrow(tab.model)) %do%
          {
            name_model <- paste0(model, ".", data.type, ".", tab.model$package[ii], ".", tab.model$func[ii])
            val.ii <- NULL
            if (strategy == "user.defined") {
              val.ii <- val.list[[name_model]]
            }
            BOD <- BIOMOD.options.dataset(mod = model
                                          , typ = data.type
                                          , pkg = tab.model$package[ii]
                                          , fun = tab.model$func[ii]
                                          , strategy = strategy
                                          , user.val = val.ii
                                          , tuning.fun = tab.model$train[ii]
                                          , bm.format = bm.format
                                          , calib.lines = calib.lines)
            for (xx in 1:length(BOD@args.values)) { ## SHOULD BE MOVED to place when testing values ??
              if ('...' %in% names(BOD@args.values[[xx]])) {
                BOD@args.values[[xx]][['...']] <- NULL
              }
            }
            return(BOD)
          }
        names(BOD.list) <- paste0(model, ".", data.type, ".", tab.model$package, ".", tab.model$func)
        return(BOD.list)
      }
      # else { warning() } ## BUT SHOULD BE DEALT WITH IN CHECK ?
    }
  bm.options <- new('BIOMOD.models.options', models = names(bm.opt), options = bm.opt)
  cat("\n")
  .bm_cat("Done")
  return(bm.options)
}

# if (!is.null(GBM$perf.method)) { opt@GBM$perf.method <- GBM$perf.method }
# opt@GAM$control <- gam::gam.control()
# opt@GAM$control <- mgcv::gam.control()
# if (!is.null(ANN$NbCV)) { opt@ANN$NbCV <- ANN$NbCV }
# if (!is.null(ANN$size)) { opt@ANN$size <- ANN$size }
# if (!is.null(ANN$decay)) { opt@ANN$decay <- ANN$decay }
# if (!is.null(ANN$rang)) { opt@ANN$rang <- ANN$rang }
# if (!is.null(ANN$maxit)) { opt@ANN$maxit <- ANN$maxit }
# if (!is.null(RF$type)) { opt@RF$type <- RF$type }

# ---------------------------------------------------------------------------- #

.bm_ModelingOptions.check.args <- function(data.type, models, strategy
                                           , val.list = NULL
                                           , bm.format = NULL, calib.lines = NULL)
{
  ## check if type is supported
  avail.types.list <- c('binary', 'binary.PA', 'abundance', 'compositional')
  .fun_testIfIn(TRUE, "data.type", data.type, avail.types.list)

  ## check if model is supported
  avail.models.list <- c('GLM', 'GBM', 'GAM', 'GAM.gam.gam', 'GAM.mgcv.bam', 'GAM.mgcv.gam'
                         , 'CTA', 'ANN', 'SRE', 'FDA', 'MARS'
                         , 'RF', 'MAXENT', 'MAXNET', 'XGBOOST')
  .fun_testIfIn(TRUE, "models", models, avail.models.list)
  if (length(grep('GAM', models)) > 1) {
    stop("Only one GAM model can be activated. Please choose betwen 'GAM', 'GAM.gam.gam', 'GAM.mgcv.bam' or 'GAM.mgcv.gam'")
  } else if ('GAM' %in% models) {
    models[which(models == 'GAM')] = 'GAM.mgcv.gam'
    warning("Only one GAM model can be activated. 'GAM.mgcv.gam' has been set (other available : 'GAM.gam.gam' or 'GAM.mgcv.bam')")
  }

  ## check if strategy is supported
  avail.strategy.list <- c('default', 'bigboss', 'user.defined', 'tuned')
  .fun_testIfIn(TRUE, "strategy", strategy, avail.strategy.list)

  ## USER DEFINED parameterisation --------------
  if (strategy == "user.defined") {
    .fun_testIfInherits(TRUE, "val.list", val.list, c("list"))
    avail.options.list <- paste0(TABLE_MODELS$model, ".", TABLE_MODELS$type, ".", TABLE_MODELS$package, ".", TABLE_MODELS$func)
    .fun_testIfIn(TRUE, "names(val.list)", names(val.list), avail.options.list)
    ## THEN can be directly arguments (all the same for all data) or a list with one set for each dataset (AllData_AllRun, ...)
  }
  
  ## TUNING with bm_Tuning parameterisation -----
  if (strategy == "tuned") {
    .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  }
  
  ## check calib.lines colnames
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun")
    if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
      expected_CVnames <- c(expected_CVnames
                            , sapply(1:ncol(bm.format@PA.table)
                                     , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
                                                           , paste0("_PA", this_PA, "_allRun"))))
    } 
    .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
    
    if (strategy == "user.defined") {
      for (ii in 1:length(val.list)) {
        .fun_testIfIn(TRUE, "names(val.list[[ii]])", names(val.list[[ii]]), expected_CVnames)
      }
    }
  }
  
  return(list(models = models))
}



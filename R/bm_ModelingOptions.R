##----------------------------------------------------------------------------##
##' @name bm_ModelingOptions
##' @aliases bm_ModelingOptions
##' @author Damien Georges, Wilfried Thuiller, Maya Gueguen
##' 
##' @title Configure the modeling options for each selected model
##'
##' @description Parameterize and/or tune \pkg{biomod2}'s single models options.
##'
##' @param data.type a \code{character} corresponding to the data type to be used, must be either 
##' \code{binary}, \code{binary.PA}, \code{abundance}, \code{compositional}. For the moment, only \code{binary} is accept.
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' @param strategy a \code{character} corresponding to the method to select models' parameters 
##' values, must be either \code{default}, \code{bigboss}, \code{user.defined}, \code{tuned}
##' @param user.val (\emph{optional, default} \code{NULL}) \cr
##' A \code{list} containing parameters values for some (all) models
##' @param user.base (\emph{optional, default} \code{bigboss}) \cr A character, 
##' \code{default} or \code{bigboss} used when \code{strategy = 'user.defined'}. 
##' It sets the bases of parameters to be modified by user defined values.
##' @param bm.format (\emph{optional, default} \code{NULL}) \cr
##' A \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} object returned 
##' by the \code{\link{BIOMOD_FormatingData}} function
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' A \code{data.frame} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{bm_CrossValidation}} functions
##' 
##' 
##'
##' @return 
##' 
##' A \code{\link{BIOMOD.models.options}} of object that can be used to build species 
##' distribution model(s) with the \code{\link{BIOMOD_Modeling}} function.
##' 
##' 
##' @details
##' 
##' This function creates a \code{\link{BIOMOD.models.options}} object containing parameter values 
##' for each single model that can be run within \pkg{biomod2} through 
##' \code{\link{BIOMOD_Modeling}} function.
##' 
##' 12 models are currently available, and are listed within the \code{\link{ModelsTable}} dataset.
##' 
##' Different strategies are available to set those parameters, through the \code{strategy} 
##' argument :
##' \describe{
##'   \item{default}{all parameters names and values are directly retrieve from functions to be 
##'   called through \code{\link[methods]{formalArgs}} and \code{\link{formals}} functions respectively}
##'   \item{bigboss}{default parameter values are updated with values predefined by \pkg{biomod2} 
##'   team}
##'   \item{user.defined}{default parameter values are updated with values provided by the user}
##'   \item{tuned}{default parameter values are updated by calling \code{\link{bm_Tuning}} 
##'   function}
##' }
##' 
##' To define the same options for all datasets of a model, you can provide these options as a list in 
##' user.val with the names "for_all_datasets".
##' 
##' @note \code{MAXENT} being the only external model (not called through a \code{R} package), 
##' default parameters, and their values, are the following :
##' 
##' \itemize{
##'   \item \code{path_to_maxent.jar = getwd()} : a \code{character} corresponding to path to 
##'   \code{maxent.jar} file
##'   \item \code{memory_allocated = 512} : an \code{integer} corresponding to the amount of 
##'   memory (in Mo) reserved for \code{java} to run \code{MAXENT}, must be either \code{64}, 
##'   \code{128}, \code{256}, \code{512}, \code{1024}... or \code{NULL} to use default \code{java}
##'   memory limitation parameter
##'   
##'   \item \code{initial_heap_size = NULL} : a \code{character} corresponding to initial heap 
##'   space (shared memory space) allocated to \code{java} (argument \code{-Xms} when calling 
##'   \code{java}), must be either \code{1024K}, \code{4096M}, \code{10G} ... or \code{NULL} to 
##'   use default \code{java} parameter. Used in \code{\link{BIOMOD_Projection}} but not in 
##'   \code{\link{BIOMOD_Modeling}}.
##'   \item \code{max_heap_size = NULL} : a \code{character} corresponding to maximum heap 
##'   space (shared memory space) allocated to \code{java} (argument \code{-Xmx} when calling 
##'   \code{java}), must be either \code{1024K}, \code{4096M}, \code{10G} ... or \code{NULL} to 
##'   use default \code{java} parameter, and must be larger than \code{initial_heap_size}. Used 
##'   in \code{\link{BIOMOD_Projection}} but not in \code{\link{BIOMOD_Modeling}}.
##'   
##'   \item \code{background_data_dir = 'default'} : a \code{character} corresponding to path 
##'   to folder where explanatory variables are stored as \code{ASCII} files (raster format). 
##'   If specified, \code{MAXENT} will generate its own background data from rasters of 
##'   explanatory variables (\code{'default'} value). Otherwise \pkg{biomod2} pseudo-absences
##'   will be used (see \code{\link{BIOMOD_FormatingData}}).
##'   \item \code{visible = FALSE} : a \code{logical} value defining whether \code{MAXENT} 
##'   user interface is to be used or not
##'   
##'   \item \code{linear = TRUE} : a \code{logical} value defining whether linear features are 
##'   to be used or not
##'   \item \code{quadratic = TRUE} : a \code{logical} value defining whether quadratic features are 
##'   to be used or not
##'   \item \code{product = TRUE} : a \code{logical} value defining whether product features are 
##'   to be used or not
##'   \item \code{threshold = TRUE} : a \code{logical} value defining whether threshold features are 
##'   to be used or not
##'   \item \code{hinge = TRUE} : a \code{logical} value defining whether hinge features are 
##'   to be used or not
##'   
##'   \item \code{l2lqthreshold = 10} : an \code{integer} corresponding to the number of 
##'   samples at which quadratic features start being used
##'   \item \code{lq2lqptthreshold = 80} : an \code{integer} corresponding to the number of 
##'   samples at which product and threshold features start being used
##'   \item \code{hingethreshold = 15} : an \code{integer} corresponding to the number of 
##'   samples at which hinge features start being used
##'   
##'   \item \code{beta_lqp = -1.0} : a \code{numeric} corresponding to the regularization 
##'   parameter to be applied to all linear, quadratic and product features (\emph{negative value 
##'   enables automatic setting})
##'   \item \code{beta_threshold = -1.0} : a \code{numeric} corresponding to the regularization 
##'   parameter to be applied to all threshold features (\emph{negative value enables automatic 
##'   setting})
##'   \item \code{beta_hinge = -1.0} : a \code{numeric} corresponding to the regularization 
##'   parameter to be applied to all hinge features (\emph{negative value enables automatic 
##'   setting})
##'   \item \code{beta_categorical = -1.0} : a \code{numeric} corresponding to the regularization 
##'   parameter to be applied to all categorical features (\emph{negative value enables automatic 
##'   setting})
##'   
##'   \item \code{betamultiplier = 1} : a \code{numeric} corresponding to the number by which 
##'   multiply all automatic regularization parameters (\emph{higher number gives a more 
##'   spread-out distribution})
##'   
##'   \item \code{defaultprevalence = 0.5} : a \code{numeric} corresponding to the default 
##'   prevalence of the modelled species (\emph{probability of presence at ordinary occurrence 
##'   points})
##' }
##' 
##'
##' @keywords models options
##' 
##' 
##' @seealso \code{\link{ModelsTable}}, \code{\link{BIOMOD.models.options}}, 
##' \code{\link{bm_Tuning}}, \code{\link{BIOMOD_Modeling}}
##' @family Secondary functions
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
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # k-fold selection
##' cv.k <- bm_CrossValidation(bm.format = myBiomodData,
##'                            strategy = 'kfold',
##'                            nb.rep = 2,
##'                            k = 3)
##' 
##' 
##' # ---------------------------------------------------------------#
##' allModels <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
##'                , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
##' 
##' # default parameters
##' opt.d <- bm_ModelingOptions(data.type = 'binary',
##'                             models = allModels,
##'                             strategy = 'default')
##' 
##' # providing formated data
##' opt.df <- bm_ModelingOptions(data.type = 'binary',
##'                              models = allModels,
##'                              strategy = 'default',
##'                              bm.format = myBiomodData,
##'                              calib.lines = cv.k)
##' 
##' opt.d
##' opt.d@models
##' opt.d@options$ANN.binary.nnet.nnet
##' names(opt.d@options$ANN.binary.nnet.nnet@args.values)
##' 
##' opt.df@options$ANN.binary.nnet.nnet
##' names(opt.df@options$ANN.binary.nnet.nnet@args.values)
##' 
##' 
##' # ---------------------------------------------------------------#
##' # bigboss parameters
##' opt.b <- bm_ModelingOptions(data.type = 'binary',
##'                             models = allModels,
##'                             strategy = 'bigboss')
##' 
##' # user defined parameters
##' user.SRE <- list('_allData_allRun' = list(quant = 0.01))
##' user.XGBOOST <- list('_allData_allRun' = list(nrounds = 10))
##' user.val <- list(SRE.binary.biomod2.bm_SRE = user.SRE
##'                  , XGBOOST.binary.xgboost.xgboost = user.XGBOOST)
##' 
##' opt.u <- bm_ModelingOptions(data.type = 'binary',
##'                             models = c('SRE', 'XGBOOST'),
##'                             strategy = 'user.defined',
##'                             user.val = user.val)
##' 
##' opt.b
##' opt.u
##' 
##' \dontrun{
##' # tuned parameters with formated data
##' opt.t <- bm_ModelingOptions(data.type = 'binary',
##'                             models = c('SRE', 'XGBOOST'),
##'                             strategy = 'tuned',
##'                             bm.format = myBiomodData)
##' opt.t
##' }
##' 

##' 
##'
##' @importFrom foreach foreach %do%
##' @importFrom methods new
##' 
##' @export
##'
##'


bm_ModelingOptions <- function(data.type
                               , models = c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
                                            , 'MARS', 'MAXENT', 'MAXNET', 'RF','RFd', 'SRE', 'XGBOOST')
                               , strategy, user.val = NULL, user.base = "bigboss"
                               , bm.format = NULL, calib.lines = NULL)
{
  .bm_cat("Build Modeling Options")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_ModelingOptions.check.args(data.type = data.type,
                                         models = models,
                                         strategy = strategy, 
                                         user.val = user.val, 
                                         user.base = user.base,
                                         bm.format = bm.format, 
                                         calib.lines = calib.lines)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## Load single models informations
  # data(ModelsTable) # internal data is already readily available
  
  ## 1. Get options -------------------------------------------------------------------------------
  bm.opt <- foreach(model = models, .combine = "c") %do% {
      ## Select model / package / function to keep
      if (length(grep("GAM", model)) == 1) {
        tab.model <- ModelsTable[which(ModelsTable$model == "GAM" & ModelsTable$type == data.type &
                                         ModelsTable$package == strsplit(model, "[.]")[[1]][2] & 
                                         ModelsTable$func == strsplit(model, "[.]")[[1]][3]), ]
        model = "GAM"
      } else {
        tab.model <- ModelsTable[which(ModelsTable$model == model & ModelsTable$type == data.type), ]
      }
      if (nrow(tab.model) > 0) {
        ## For each kept model : get corresponding options
        BOD.list <- foreach(ii = 1:nrow(tab.model)) %do% {
            name_model <- paste0(model, ".", data.type, ".", tab.model$package[ii], ".", tab.model$func[ii])
            val.ii <- NULL
            if (strategy == "user.defined") {
              val.ii <- user.val[[name_model]]
            }
            BOD <- BIOMOD.options.dataset(mod = model
                                          , typ = data.type
                                          , pkg = tab.model$package[ii]
                                          , fun = tab.model$func[ii]
                                          , strategy = strategy
                                          , user.val = val.ii
                                          , user.base = user.base
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



# ---------------------------------------------------------------------------- #

.bm_ModelingOptions.check.args <-
  function(data.type, models, strategy
           , user.val = NULL, user.base = NULL
           , bm.format = NULL, calib.lines = NULL) {
    ## check if type is supported
    avail.types.list <- c('binary', 'binary.PA', 'abundance', 'compositional')
    .fun_testIfIn(TRUE, "data.type", data.type, avail.types.list)
    
    ## check if model is supported
    avail.models.list <- c('ANN', 'CTA', 'FDA', 'GAM', 'GAM.gam.gam', 'GAM.mgcv.bam', 'GAM.mgcv.gam'
                           , 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF','RFd', 'SRE', 'XGBOOST')
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
      avail.user.base <- c('default', 'bigboss')
      .fun_testIfIn(TRUE, "user.base", user.base, avail.user.base)
      
      .fun_testIfInherits(TRUE, "user.val", user.val, c("list"))
      avail.options.list <- paste0(ModelsTable$model, ".", ModelsTable$type, ".", ModelsTable$package, ".", ModelsTable$func)
      .fun_testIfIn(TRUE, "names(user.val)", names(user.val), avail.options.list)
      ## THEN can be directly arguments (all the same for all data) or a list with one set for each dataset (AllData_AllRun, ...)
    }
    
    ## TUNING with bm_Tuning parameterisation -----
    if (strategy == "tuned") {
      .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
    }
    
    ## check calib.lines colnames
    if (!is.null(calib.lines)) {
      if (is.null(bm.format)) {
        stop("`bm.format` must be given along with `calib.lines`")
      }
      .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
      
      expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun", "for_all_datasets")
      if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
        expected_CVnames <- c(expected_CVnames
                              , sapply(1:ncol(bm.format@PA.table)
                                       , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
                                                             , paste0("_PA", this_PA, "_allRun"))))
      }
      .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
      
      if (strategy == "user.defined") {
        for (ii in 1:length(user.val)) {
          .fun_testIfIn(TRUE, "names(user.val[[ii]])", names(user.val[[ii]]), expected_CVnames)
        }
      }
    }
    
    return(list(models = models))
  }



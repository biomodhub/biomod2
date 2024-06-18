
## -------------------------------------------------------------------------- #
## 0. Generic Functions definition ------------------------------------------
## -------------------------------------------------------------------------- #
## Used for different classes 
##    01 = BIOMOD.formated.data, 02 = BIOMOD.formated.data.PA
##    A = BIOMOD.models.out, B = BIOMOD.projection.out, C = BIOMOD.ensemble.models.out

##' @name getters.out
##' @aliases get_species_data
##' @aliases get_eval_data
##' @aliases get_options
##' @aliases get_calib_lines
##' @aliases get_formal_data
##' @aliases get_projected_models
##' @aliases free
##' @aliases get_predictions
##' @aliases get_kept_models
##' @aliases get_built_models
##' @aliases get_evaluations
##' @aliases get_variables_importance
##' @author Damien Georges
##' 
##' @title Functions to extract informations from \code{\link{BIOMOD.models.out}}, 
##' \code{\link{BIOMOD.projection.out}} or \code{\link{BIOMOD.ensemble.models.out}} objects
##' 
##' @description These functions allow the user to easily retrieve informations stored in the 
##' different \pkg{biomod2} objects from the different modeling steps, such as modeling options 
##' and formated data, models used or not, predictions, evaluations, variables importance.
##' 
##' 
##' @param obj a \code{\link{BIOMOD.formated.data}}, \code{\link{BIOMOD.formated.data.PA}}, 
##' \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.projection.out}} or 
##' \code{\link{BIOMOD.ensemble.models.out}} object
##' @param \ldots (\emph{optional, one or several of the following arguments depending on the selected 
##' function}) 
##' @param as.data.frame a \code{logical} defining whether output should be returned as 
##' \code{data.frame} or \code{array} object
##' @param subinfo a \code{character} corresponding to the information to be extracted, must be 
##' among \code{NULL}, \code{expl.var.names}, \code{resp.var}, \code{expl.var}, \code{MinMax}, 
##' \code{eval.resp.var}, \code{eval.expl.var} (see Details)
##' @param evaluation a \code{logical} defining whether evaluation data should be used or not
##' 
##' @param full.name (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing model names to be kept, must be either \code{all} or a 
##' sub-selection of model names that can be obtained with the \code{\link{get_built_models}} 
##' function
##' 
##' @param PA (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing pseudo-absence set to be loaded, must be among \code{PA1}, 
##' \code{PA2}, \code{...}, \code{allData}
##' @param run (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing repetition set to be loaded, must be among \code{RUN1}, 
##' \code{RUN2}, \code{...}, \code{allRun}
##' @param algo (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} containing algorithm to be loaded, must be either 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}
##' 
##' @param merged.by.PA (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing merged pseudo-absence set to be loaded, must be among \code{PA1}, 
##' \code{PA2}, \code{...}, \code{mergedData}
##' @param merged.by.run (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing merged repetition set to be loaded, must be among \code{RUN1}, 
##' \code{RUN2}, \code{...}, \code{mergedRun}
##' @param merged.by.algo (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} containing merged algorithm to be loaded, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}, \code{mergedAlgo}
##' @param filtered.by (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric selected to filter single models to build the 
##' ensemble models, must be among \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, 
##' \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, 
##' \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' 
##' @param metric.eval (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric to be kept, must be among \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, 
##' \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param expl.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing explanatory variables to be kept, that can be obtained with the 
##' \code{\link{get_formal_data}(obj, subinfo = 'expl.var.names')} function
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric selected to transform predictions into binary 
##' values, must be among \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
##' \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, 
##' \code{ETS}, \code{BOYCE}, \code{MPA}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric to filter predictions, must be among \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, \code{ROC}, \code{TSS}, 
##' \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' 
##' @param model.as.col (\emph{optional, default} \code{FALSE}) \cr
##' A \code{boolean} given to \code{\link{get_predictions}}. If \code{TRUE} 
##' prediction are returned as a wide \code{data.frame} with each column containing
##' predictions for a single model and corresponding to the old output given by
##' \pkg{biomod2} in version < 4.2-2. If \code{FALSE} predictions are returned 
##' as a long \code{data.frame} with many additional informations readily 
##' available.
##' 
##' 
##' @return 
##' 
##' \describe{
##'   \item{\code{get_species_data}}{a \code{data.frame} combining \code{data.species}, 
##'   \code{coord}, \code{data.env.var} (and \code{PA.table}) slots of 
##'   \code{\link{BIOMOD.formated.data}} (or \code{\link{BIOMOD.formated.data.PA}}) object}
##'   \item{\code{get_eval_data}}{a \code{data.frame} combining \code{eval.data.species}, 
##'   \code{eval.coord}, \code{eval.data.env.var} slots of 
##'   \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} object}
##' 
##'   \item{\code{get_options}}{a
##'   \code{\link{BIOMOD.stored.options-class}} object from the
##'   \code{models.options} slot of a \code{\link{BIOMOD.models.out-class}}
##'   object} \item{\code{get_calib_lines}}{a
##'   \code{\link{BIOMOD.stored.data.frame-class}} object from the \code{calib.lines}
##'   slot of a \code{\link{BIOMOD.models.out}} object}
##'
##'   \item{\code{get_projected_models}}{a \code{vector} from the
##'   \code{models.projected} slot of a \code{\link{BIOMOD.projection.out}}
##'   object}
##'
##'   \item{\code{get_predictions}}{a \code{\link{BIOMOD.stored.data}} object
##'   from the \code{proj.out} slot of a \code{\link{BIOMOD.models.out}},
##'   \code{\link{BIOMOD.projection.out}} or
##'   \code{\link{BIOMOD.ensemble.models.out}} object}
##'
##'   \item{\code{get_kept_models}}{a \code{vector} containing names of the kept
##'   models of a \code{\link{BIOMOD.ensemble.models.out}} object}
##'
##'   \item{\code{get_formal_data}}{depending on the \code{subinfo} parameter :
##'   \describe{
##'     \item{\code{NULL}}{a \code{\link{BIOMOD.stored.formated.data-class}} (or
##'     \code{\link{BIOMOD.stored.models.out-class}}) object from the
##'     \code{formated.input.data} (or \code{models.out}) slot of a
##'     \code{\link{BIOMOD.models.out}} (or
##'     \code{\link{BIOMOD.ensemble.models.out}}) object}
##'     
##'     \item{\code{expl.var.names}}{a \code{vector} from the
##'     \code{expl.var.names} slot of a \code{\link{BIOMOD.models.out}} or
##'     \code{\link{BIOMOD.ensemble.models.out}} object}

##'     \item{\code{resp.var}}{a \code{vector} from the \code{data.species} slot
##'     of the \code{formated.input.data} slot of a
##'     \code{\link{BIOMOD.models.out}} or
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{expl.var}}{a \code{data.frame} from the \code{data.env.var}
##'     slot of the \code{formated.input.data} slot of a
##'     \code{\link{BIOMOD.models.out}} or
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{MinMax}}{a \code{list} of minimum and maximum values (or
##'     levels if factorial) of variable contained in the \code{data.env.var}
##'     slot of the \code{formated.input.data} slot of a
##'     \code{\link{BIOMOD.models.out}} or
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{eval.resp.var}}{a \code{vector} from the
##'     \code{eval.data.species} slot of the \code{formated.input.data} slot of
##'     a \code{\link{BIOMOD.models.out}} or
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{eval.expl.var}}{a \code{data.frame} from the
##'     \code{eval.data.env.var} slot of the \code{formated.input.data} slot of
##'     a \code{\link{BIOMOD.models.out}} or
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'   }
##'   }
##'   \item{\code{get_built_models}}{a \code{vector} from the
##'   \code{models.computed} slot (or \code{em.computed}) of a
##'   \code{\link{BIOMOD.models.out}} (or
##'   \code{\link{BIOMOD.ensemble.models.out}}) object}
##'   \item{\code{get_evaluations}}{a data.frame from the \code{models.evaluation}
##'    slot (or \code{model_evaluation} of each model in \code{em.computed}) of a
##'   \code{\link{BIOMOD.models.out}} (or \code{\link{BIOMOD.ensemble.models.out}})
##'    object. Contains evaluation metric for different models and dataset. 
##'    Evaluation metric are calculated on the calibrating data (column \code{calibration}),
##'    on the cross-validation data (column \code{validation}) or on the evaluation data 
##'    (column \code{evaluation}). \cr \emph{For cross-validation data, see \code{CV.[...]} 
##'    parameters in \code{\link{BIOMOD_Modeling}} function ; for evaluation data, see 
##'    \code{eval.[...]} parameters in \code{\link{BIOMOD_FormatingData}}.}}
##'   \item{\code{get_variables_importance}}{a
##'   \code{\link{BIOMOD.stored.data.frame-class}} from
##'   the \code{variables.importance} slot (or \code{model_variables_importance}
##'   of each model in \code{em.models}) of a \code{\link{BIOMOD.models.out}}
##'   (or \code{\link{BIOMOD.ensemble.models.out}}) object}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.projection.out}}, 
##' \code{\link{BIOMOD.ensemble.models.out}}
##' @family Toolbox functions
##' 
##' 
##' @importFrom reshape melt.array
##' @importFrom foreach foreach %do%
##' @importFrom abind abind
##' @importFrom terra rast subset
##' 
NULL

setGeneric("get_species_data", function(obj, ...) { standardGeneric("get_species_data") }) ## 012
setGeneric("get_eval_data", function(obj, ...) { standardGeneric("get_eval_data") }) ## 012

setGeneric("get_options", function(obj, ...) { standardGeneric("get_options") }) ## A
setGeneric("get_calib_lines", function(obj, ...) { standardGeneric("get_calib_lines") }) ## A

setGeneric("get_projected_models", function(obj, ...) { standardGeneric("get_projected_models") }) ## B
setGeneric("free", function(obj, ...) { standardGeneric("free") }) ## B

setGeneric("get_predictions", function(obj, ...) { standardGeneric("get_predictions") }) ## ABC

setGeneric("get_kept_models", function(obj, ...) { standardGeneric("get_kept_models") }) ## C

setGeneric("get_formal_data", function(obj, ...) { standardGeneric("get_formal_data") }) ## AC
setGeneric("get_built_models", function(obj, ...) { standardGeneric("get_built_models") }) ## AC
setGeneric("get_evaluations", function(obj, ...) { standardGeneric("get_evaluations") }) ## AC
setGeneric("get_variables_importance", function(obj, ...) { standardGeneric("get_variables_importance") }) ## AC


## get_species_data.BIOMOD.formated.data ------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod('get_species_data', signature('BIOMOD.formated.data'), function(obj) {
  tab.sp <- data.frame(obj@data.species)
  tab.sp <- cbind(tab.sp, obj@coord)
  colnames(tab.sp) <- c(obj@sp.name, "x", "y")
  tab.sp <- cbind(tab.sp, obj@data.env.var)
  return(tab.sp)
})

## get_species_data.BIOMOD.formated.data.PA ---------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod('get_species_data', signature('BIOMOD.formated.data.PA'), function(obj) {
  tab.sp <- data.frame(obj@data.species)
  tab.sp <- cbind(tab.sp, obj@coord)
  colnames(tab.sp) <- c(obj@sp.name, "x", "y")
  tab.sp <- cbind(tab.sp, obj@data.env.var)
  tab.sp <- cbind(tab.sp, obj@PA.table)
  return(tab.sp)
})

## get_eval_data.BIOMOD.formated.data ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod('get_eval_data', signature('BIOMOD.formated.data'), function(obj) {
  if (obj@has.data.eval) {
    tab.sp <- data.frame(obj@eval.data.species)
    tab.sp <- cbind(tab.sp, obj@eval.coord)
    colnames(tab.sp) <- c(obj@sp.name, "x", "y")
    tab.sp <- cbind(tab.sp, obj@eval.data.env.var)
    return(tab.sp)
  } else { return(NULL) }
})


## -------------------------------------------------------------------------- #
## 4. BIOMOD.models.out -----------------------------------------------------
## -------------------------------------------------------------------------- #

##' @name BIOMOD.models.out
##' @aliases BIOMOD.models.out-class
##' @aliases BIOMOD.models.out
## @aliases BIOMOD.stored.models.out
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_Modeling()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_Modeling}}, and used by 
##' \code{\link{BIOMOD_LoadModels}}, \code{\link{BIOMOD_PresenceOnly}}, 
##' \code{\link{BIOMOD_Projection}} and \code{\link{BIOMOD_EnsembleModeling}}
##' 
##' 
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the
##'   simulation set
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory
##'   variables
##' @slot models.computed a \code{vector} containing names of computed models
##' @slot models.failed a \code{vector} containing names of failed models
##' @slot has.evaluation.data a \code{logical} value defining whether evaluation
##'   data is given
##' @slot scale.models a \code{logical} value defining whether models have been
##'   rescaled or not
##' @slot formated.input.data a \code{\link{BIOMOD.stored.formated.data-class}}
##'   object containing informations from \code{\link{BIOMOD_FormatingData}}
##'   object
##' @slot calib.lines a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing calibration lines
##' @slot models.options a \code{\link{BIOMOD.stored.options-class}}
##'   object containing informations from \code{\link{bm_ModelingOptions}}
##'   object
##' @slot models.evaluation a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing models evaluation
##' @slot variables.importance a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing variables importance
##' @slot models.prediction a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing models predictions
##' @slot models.prediction.eval a \code{\link{BIOMOD.stored.data.frame-class}}
##'   object containing models predictions for evaluation data
##' @slot link a \code{character} containing the file name of the saved object
##'
##' @param object a \code{\link{BIOMOD.models.out}} object
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_LoadModels}}, 
##' \code{\link{BIOMOD_PresenceOnly}}, \code{\link{BIOMOD_Projection}}, 
##' \code{\link{BIOMOD_EnsembleModeling}}, \code{\link{bm_VariablesImportance}}, 
##' \code{\link{bm_PlotEvalMean}}, \code{\link{bm_PlotEvalBoxplot}}, 
##' \code{\link{bm_PlotVarImpBoxplot}}, \code{\link{bm_PlotResponseCurves}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.models.out")
##' 
##' ## ----------------------------------------------------------------------- #
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
##' ## ----------------------------------------------------------------------- #
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' ## ----------------------------------------------------------------------- #
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     models = c('RF', 'GLM'),
##'                                     CV.strategy = 'random',
##'                                     CV.nb.rep = 2,
##'                                     CV.perc = 0.8,
##'                                     OPT.strategy = 'bigboss',
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 3,
##'                                     seed.val = 42)
##' myBiomodModelOut
##' 
##' 
NULL

##' @name BIOMOD.models.out-class
##' @rdname BIOMOD.models.out
##' @export
##' 

# 4.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.models.out",
         representation(modeling.id = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.computed = 'character',
                        models.failed = 'character',
                        has.evaluation.data = 'logical',
                        scale.models = 'logical',
                        formated.input.data = 'BIOMOD.stored.formated.data',
                        calib.lines = 'BIOMOD.stored.data.frame',
                        models.options = 'BIOMOD.stored.options',
                        models.evaluation = 'BIOMOD.stored.data.frame',
                        variables.importance = 'BIOMOD.stored.data.frame',
                        models.prediction = 'BIOMOD.stored.data.frame',
                        models.prediction.eval = 'BIOMOD.stored.data.frame',
                        link = 'character'),
         prototype(modeling.id = as.character(format(Sys.time(), "%s")),
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   models.computed = '',
                   models.failed = '',
                   has.evaluation.data = FALSE,
                   scale.models = TRUE,
                   formated.input.data = new('BIOMOD.stored.formated.data'),
                   calib.lines = new('BIOMOD.stored.data.frame'),
                   models.options = new('BIOMOD.stored.options'),
                   models.evaluation = new('BIOMOD.stored.data.frame'),
                   variables.importance = new('BIOMOD.stored.data.frame'),
                   models.prediction = new('BIOMOD.stored.data.frame'),
                   models.prediction.eval = new('BIOMOD.stored.data.frame'),
                   link = ''),
         validity = function(object){ return(TRUE) } )

# BIOMOD.stored.models.out is defined here and not with outher BIOMOD.stored.data
# as its definition require the definition of class BIOMOD.stored.data and files are
# sourced in alphabetical order.
##' @name BIOMOD.stored.models.out-class
##' @rdname BIOMOD.stored.data

setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) } )


# 4.3 Other functions ------------------------------------------------------
## show.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname BIOMOD.models.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.models.out'),
          function(object) {
            .bm_cat("BIOMOD.models.out")
            cat("\nModeling folder :", object@dir.name, fill = .Options$width)
            cat("\nSpecies modeled :", object@sp.name, fill = .Options$width)
            cat("\nModeling id :", object@modeling.id, fill = .Options$width)
            cat("\nConsidered variables :", object@expl.var.names, fill = .Options$width)
            cat("\n\nComputed Models : ", object@models.computed, fill = .Options$width)
            cat("\n\nFailed Models : ", object@models.failed, fill = .Options$width)
            .bm_cat()
          }
)

## get_options.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_options", "BIOMOD.models.out",
          function(obj) {
            model_options <- load_stored_object(obj@models.options)
            return(model_options)
          }
)

## get_calib_lines.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##'

setMethod("get_calib_lines", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, PA = NULL, run = NULL) {
            out <- load_stored_object(obj@calib.lines)
            
            if (!is.null(out) && as.data.frame == TRUE) {
              tmp <- melt(out, varnames = c("points", "PA_run"))
              tmp$PA = strsplit(sub("^_", "", tmp$PA_run), "_")[[1]][1]
              tmp$run = strsplit(sub("^_", "", tmp$PA_run), "_")[[1]][2]
              out <- tmp[, c("PA", "run", "points", "value")]
              colnames(out)[4] = "calib.lines"
              
              keep_lines <- .filter_outputs.df(out, subset.list = list(PA = PA, run = run))
              out <- out[keep_lines, ]
            }
            return(out)
          }
)

## get_formal_data.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_formal_data", "BIOMOD.models.out",
          function(obj, subinfo = NULL) {
            if (is.null(subinfo)) {
              return(load_stored_object(obj@formated.input.data))
            } else if (subinfo == 'MinMax') {
              env = as.data.frame(get_formal_data(obj)@data.env.var)
              MinMax = foreach(i = 1:ncol(env)) %do% {
                x = env[, i]
                if (is.numeric(x)) {
                  return(list(min = min(x, na.rm = TRUE)
                              , max = max(x, na.rm = TRUE)))
                } else if (is.factor(x)) {
                  return(list(levels = levels(x)))
                }
              }
              names(MinMax) = colnames(env)
              return(MinMax)
            } else if (subinfo == 'expl.var') {
              return(as.data.frame(get_formal_data(obj)@data.env.var))
            } else if (subinfo == 'expl.var.names') {
              return(obj@expl.var.names)
            } else if (subinfo == 'resp.var') {
              return(as.numeric(get_formal_data(obj)@data.species))
            } else if (subinfo == 'eval.resp.var') {
              return(as.numeric(get_formal_data(obj)@eval.data.species))
            } else if (subinfo == 'eval.expl.var') {
              return(as.data.frame(get_formal_data(obj)@eval.data.env.var))
            } else { stop("Unknown subinfo tag")}
          }
)



## get_predictions.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "BIOMOD.models.out",
          function(obj, evaluation = FALSE
                   , full.name = NULL, PA = NULL, run = NULL, algo = NULL,
                   model.as.col = FALSE)
          {
            if (evaluation && (!obj@has.evaluation.data)) {
              warning("!   Calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if (evaluation) {
              out <- load_stored_object(obj@models.prediction.eval)
            } else { 
              out <- load_stored_object(obj@models.prediction)
            }
            
            # subselection of models_selected
            keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name, PA = PA
                                                                     , run = run, algo = algo))
            out <- out[keep_lines, ]
            if (model.as.col) {
              out <- .transform_model.as.col(out)
            }
            return(out)
          }
)

## get_built_models.BIOMOD.models.out ---------------------------------------------------
##' @rdname getters.out
##' @export
##' 

setMethod("get_built_models", "BIOMOD.models.out",
          function(obj, full.name = NULL, PA = NULL, run = NULL, algo = NULL) { 
            out <- obj@models.computed
            keep_ind <- .filter_outputs.vec(out, obj.type = "mod", subset.list = list(full.name =  full.name, PA = PA
                                                                                      , run = run, algo = algo))
            out <- out[keep_ind]
            return(out)
          }
)

## get_evaluations.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_evaluations", "BIOMOD.models.out",
          function(obj, full.name = NULL, PA = NULL, run = NULL, algo = NULL, metric.eval = NULL) {
            out <- load_stored_object(obj@models.evaluation)
            if(nrow(out) == 0){
              cat("\n! models have no evaluations\n")
              return(invisible(NULL))
            } else {
              keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name, PA = PA
                                                                       , run = run, algo = algo
                                                                       , metric.eval = metric.eval))
              out <- out[keep_lines, ]
              return(out)
            }
          }
)

## get_variables_importance.BIOMOD.models.out ---------------------------------------------------
##' @rdname getters.out
##' @export
##' 

setMethod("get_variables_importance", "BIOMOD.models.out",
          function(obj, full.name = NULL, PA = NULL, run = NULL, algo = NULL, expl.var = NULL) {
            out <- load_stored_object(obj@variables.importance)
            if(obj@variables.importance@link == ''){
              cat("\n! models have no variables importance\n")
              return(invisible(NULL))
            } else {
              keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name, PA = PA
                                                                       , run = run, algo = algo
                                                                       , expl.var = expl.var))
              out <- out[keep_lines, ]
              return(out)
            }
          }
)


## --------------------------------------------------------------------------- #
## 5. BIOMOD.projection.out --------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.projection.out
##' @aliases BIOMOD.projection.out-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_Projection()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_Projection}}, and used by 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' 
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the simulation set
##' @slot proj.name a \code{character} corresponding to the projection name
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory variables
##' @slot coord a 2-columns \code{matrix} or \code{data.frame} containing the corresponding 
##' \code{X} and \code{Y} coordinates used to project the species distribution model(s)
##' @slot scale.models a \code{logical} value defining whether models have been rescaled or 
##' not
##' @slot models.projected a \code{vector} containing names of projected models
##' @slot models.out a \code{\link{BIOMOD.stored.data}} object
##' @slot type a \code{character} corresponding to the class of the \code{val} slot of the 
##' \code{proj.out} slot
##' @slot proj.out a \code{\link{BIOMOD.stored.data}} object
##' 
##' @param x a \code{\link{BIOMOD.projection.out}} object
##' @param object a \code{\link{BIOMOD.projection.out}} object
##' @param coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' @param plot.output (\emph{optional, default} \code{facet}) a character
##'   determining the type of output: with \code{plot.output = 'list'} the
##'   function will return a list of plots (one plot per model) ; with 'facet' ;
##'   with \code{plot.output = 'facet'} the function will return a single plot
##'   with all asked projections as facet.
##' @param do.plot (\emph{optional, default} \code{TRUE}) a boolean determining
##'   whether the plot should be displayed or just returned.
##' @param std (\emph{optional, default} \code{TRUE}) a boolean controlling the
##'   limits of the color scales. With \code{std = TRUE} color scales are
##'   displayed between 0 and 1 (or 1000). With \code{std = FALSE} color scales
##'   are displayed between 0 and the maximum value observed.
##' @param scales (\emph{optional, default} \code{fixed}) a character
##'   determining whether x and y scales are shared among facet. Argument passed
##'   to \code{\link[ggplot2:facet_wrap]{facet_wrap}}. Possible values: 'fixed', 'free_x',
##'   'free_y', 'free'.
##' @param size (\emph{optional, default} \code{0.75}) a numeric determing the
##'   size of points on the plots and passed to
##'   \code{\link[ggplot2:geom_point]{geom_point}}.
##' @param ... additional parameters to be passed to \code{\link{get_predictions}} 
##' to select the models that will be plotted
##'           
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.projection.out")
##' 
##' ## ----------------------------------------------------------------------- #
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
##' ## ----------------------------------------------------------------------- #
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                        expl.var = myExpl,
##'                                        resp.xy = myRespXY,
##'                                        resp.name = myRespName)
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' }
##' 
##' 
##' ## ----------------------------------------------------------------------- #
##' # Project single models
##' myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                   proj.name = 'Current',
##'                                   new.env = myExpl,
##'                                   models.chosen = 'all',
##'                                   metric.binary = 'all',
##'                                   metric.filter = 'all',
##'                                   build.clamping.mask = TRUE)
##' myBiomodProj
##' plot(myBiomodProj)
##' 
##' 
##' @importFrom grDevices colorRampPalette colors dev.new gray rainbow
##' @importFrom graphics layout legend par points polygon text
##' @importFrom ggplot2 scale_colour_viridis_c scale_fill_viridis_c
##' 
NULL

##' @name BIOMOD.projection.out-class
##' @rdname BIOMOD.projection.out
##' @export
##' 

# 5.1 Class Definition  -----------------------------------

setClass("BIOMOD.projection.out",
         representation(modeling.id = 'character',
                        proj.name = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        coord = 'data.frame',
                        scale.models = 'logical',
                        models.projected = 'character',
                        models.out = 'BIOMOD.stored.data',
                        type = 'character',
                        proj.out = 'BIOMOD.stored.data'),
         prototype(modeling.id = '',
                   proj.name = '',
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   coord = data.frame(),
                   scale.models = TRUE,
                   models.projected = '',
                   type = ''),
         validity = function(object){ return(TRUE) })


# 5.3 Other functions ---------------------------------------------------------
## plot.BIOMOD.projection.out -------------------------------------------------
##' 
##' @rdname BIOMOD.projection.out
##' @export
##' @importFrom terra global
##' @param maxcell maximum number of cells to plot. Argument transmitted to \code{\link[terra]{plot}}.
##' 

setMethod(
  'plot', signature(x = 'BIOMOD.projection.out', y = "missing"),
  function(x,
           coord = NULL,
           plot.output, # list or facet
           do.plot = TRUE, # whether plots are displayed or just returned
           std = TRUE, # limits between 0 and 1000 or between 0 and max
           scales, # transmitted to facet_wrap
           size, # size of points transmitted to geom_point
           maxcell = 5e5, # max number of cells to plot. Transmitted to terra::plot
           ...
  ){
    # extraction of projection happens in argument check
    args <- .plot.BIOMOD.projection.out.check.args(x,
                                                   coord = coord,
                                                   plot.output = plot.output, # list or facet
                                                   do.plot = do.plot,
                                                   std = std,
                                                   scales = scales,
                                                   size = size,
                                                   maxcell = maxcell,
                                                   ...)
    for (argi in names(args)) { 
      assign(x = argi, value = args[[argi]]) 
    }
    rm(args)
    
    
    ### Plot SpatRaster ---------------------------------------------------------
    
    if (inherits(proj,"SpatRaster")) {
      maxi <- ifelse(max(global(proj, "max", na.rm = TRUE)$max) > 1, 1000, 1) 
      if (std) {
        limits <-  c(0,maxi)
      } else {
        limits <- NULL
      }
      
      if (plot.output == "facet") {
        g <- ggplot() +
          tidyterra::geom_spatraster(data = proj,
                                     maxcell = maxcell) +
          scale_fill_viridis_c(NULL, limits = limits) +
          facet_wrap(~lyr)
      } else if (plot.output == "list") {
        g <- lapply(names(proj), function(thislayer){
          ggplot() +
            tidyterra::geom_spatraster(data = subset(proj, thislayer),
                                       maxcell = maxcell) +
            scale_fill_viridis_c(NULL, limits = limits) +
            ggtitle(thislayer)
        })
      }
    } else {
      ### Plot data.frame  -----------------------------------------------------
      maxi <- ifelse(max(proj$pred) > 1, 1000, 1) 
      if (std) {
        limits <-  c(0,maxi)
      } else {
        limits <- NULL
      }
      plot.df <- merge(proj, coord, by = c("points"))
      if(plot.output == "facet"){
        g <- ggplot(plot.df)+
          geom_point(aes(x = x, y = y, color = pred), size = size) +
          scale_colour_viridis_c(NULL, limits = limits) +
          facet_wrap(~full.name)
      } else if (plot.output == "list"){
        g <- lapply(unique(plot.df$full.name), function(thislayer){
          ggplot(subset(plot.df, plot.df$full.name == thislayer)) +
            geom_point(aes(x = x, y = y, color = pred), size = size) +
            scale_colour_viridis_c(NULL, limits = limits) +
            ggtitle(thislayer)
        })
      }
      
    }
    if (do.plot) {
      show(g)
    } 
    return(g)
  }
)

### .plot.BIOMOD.projection.out.check.args ----------------------------------

.plot.BIOMOD.projection.out.check.args <- function(x,
                                                   coord,
                                                   plot.output, # list or facet
                                                   do.plot,
                                                   std,
                                                   scales,
                                                   size,
                                                   ...){
  
  proj <- get_predictions(x, ...)
  
  ## 1 - check for tidyterra ----------------------
  if (inherits(proj, "SpatRaster")) {
    if (!requireNamespace("tidyterra")) {
      stop("Package `tidyterra` is missing. Please install it with `install.packages('tidyterra')`.")
    }
  }
  
  ## 2 - plot.output----------------------
  if (missing(plot.output)) {
    plot.output <- "facet"
  } else {
    .fun_testIfIn(TRUE, "plot.output", plot.output, c("facet","list"))
  }
  
  ## 3 - do.plot ----------------------
  stopifnot(is.logical(do.plot))
  
  ## 4 - std ----------------------
  stopifnot(is.logical(std))
  
  ## 5 - check scales for facet_wrap -------------------------------
  if(missing(scales)){
    scales <- "fixed"
  } else {
    .fun_testIfIn(TRUE, "scales", scales,
                  c("fixed","free","free_x","free_y"))
  }
  
  ## 6 - check coord if x is a data.frame -------------------------------
  if (inherits(proj, 'data.frame')) {
    npred <- length(unique(proj$points))
    
    if (nrow(x@coord) > 0) {
      if(!is.null(coord)){
      cat("! ignoring argument `coord` as coordinates were already given to BIOMOD_Projection")
      }
      coord <- x@coord
    }

    if (nrow(x@coord) == 0 & is.null(coord)) {
        stop("missing coordinates to plot with a data.frame. Either give argument `coord` to plot or argument `new.env.xy` to BIOMOD_Projection")
      
    } else if (!inherits(coord, c("data.frame","matrix"))) {
      stop("`coord` must be a data.frame or a matrix.")
    } else if (ncol(coord) != 2) {
      stop("`coord` must have two columns.")
    } else if (nrow(coord) != npred) {
      stop("`coord` must have as many rows as the number of predictions (", npred, ").")
    } else {
      coord <- as.data.frame(coord)
      colnames(coord) <- c("x","y")
      coord$points <- seq_len(npred)
    }
  }
  
  if(missing(size)){
    size <- 0.75
  } 
  
  ## 7 - check size -------------------------------
  if (inherits(proj, 'data.frame')) {
    .fun_testIfPosNum(TRUE, "size", size)
  }
  
  return(list(proj = proj,
              coord = coord,
              plot.output = plot.output,
              do.plot = do.plot,
              std = std,
              scales = scales, 
              size = size))
}

## show.BIOMOD.projection.out -------------------------------------------------
##' 
##' @rdname BIOMOD.projection.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.projection.out'),
          function(object){
            .bm_cat("BIOMOD.projection.out")
            cat("\nProjection directory :", paste0(object@dir.name, "/", object@sp.name, "/", object@proj.name), fill = .Options$width)
            cat("\n")
            cat("\nsp.name :", object@sp.name, fill = .Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill = .Options$width)
            cat("\n")
            cat("\nmodeling.id :", object@modeling.id , "(", object@models.out@link , ")", fill = .Options$width)
            cat("\nmodels.projected :", toString(object@models.projected), fill = .Options$width)
            df.info <- .extract_projlinkInfo(object)
            if(any(df.info$type == "bin")){
              available.metric <- unique(subset(df.info, df.info$type == "bin")$metric)
              cat("\navailable binary projection :", toString(available.metric), fill = .Options$width)
            }
            if(any(df.info$type == "filt")){
              available.metric <- unique(subset(df.info, df.info$type == "filt")$metric)
              cat("\navailable filtered projection :", toString(available.metric), fill = .Options$width)
            }
            .bm_cat()
          })

## get_projected_models.BIOMOD.projection.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_projected_models", "BIOMOD.projection.out",
          function(obj, full.name = NULL, PA = NULL, run = NULL, algo = NULL
                   , merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL)
          {
            
            out <- obj@models.projected
            if (length(grep("EM|merged", out)) > 0) {
              keep_ind <- .filter_outputs.vec(out, obj.type = "em", subset.list = list(full.name = full.name
                                                                                       , merged.by.PA = merged.by.PA
                                                                                       , merged.by.run = merged.by.run
                                                                                       , merged.by.algo = merged.by.algo
                                                                                       , filtered.by = filtered.by
                                                                                       , algo = algo))
            } else {
              keep_ind <- .filter_outputs.vec(out, obj.type = "mod", subset.list = list(full.name = full.name, PA = PA
                                                                                        , run = run, algo = algo))
            }
            out <- out[keep_ind]
            return(out)
          }
)

## free.BIOMOD.projection.out --------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' @importFrom terra rast

setMethod('free', signature('BIOMOD.projection.out'),
          function(obj) {
            if (inherits(obj@proj.out, "BIOMOD.stored.data.frame")) {
              obj@proj.out@val  <- data.frame()
            } else if (inherits(obj@proj.out, "BIOMOD.stored.SpatRaster")) {
              obj@proj.out@val <- wrap(rast(matrix()))
            } else {
              obj@proj.out@val <- NULL
            }
            obj@proj.out@inMemory <- FALSE
            return(obj)
          })

## get_predictions.BIOMOD.projection.out ---------------------------------------
# (the method is used for EM as well)
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "BIOMOD.projection.out",
          function(obj, metric.binary = NULL, metric.filter = NULL
                   , full.name = NULL, PA = NULL, run = NULL, algo = NULL
                   , merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, 
                   model.as.col = FALSE, ...) {
            
            # extract layers from obj@proj.out@link concerned by metric.filter 
            # or metric.binary
            selected.layers <- .extract_selected.layers(obj, 
                                                        metric.binary = metric.binary,
                                                        metric.filter = metric.filter)
            out <- load_stored_object(obj@proj.out, layer = selected.layers)
            
            # subselection of models_selected
            if (obj@type == "SpatRaster") {
              if (length(grep("EM|merged", names(out))) > 0) {
                keep_layers <- .filter_outputs.vec(names(out), obj.type = "em", 
                                                   subset.list = list(full.name =  full.name
                                                                      , merged.by.PA = merged.by.PA
                                                                      , merged.by.run = merged.by.run
                                                                      , merged.by.algo = merged.by.algo
                                                                      , filtered.by = filtered.by
                                                                      , algo = algo))
              } else {
                keep_layers <- .filter_outputs.vec(names(out), obj.type = "mod",
                                                   subset.list = list(full.name =  full.name, PA = PA
                                                                      , run = run, algo = algo))
              }
              out <- subset(out, keep_layers)
            } else {
              if (length(grep("EM|merged", colnames(out))) > 0) {
                keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name
                                                                         , merged.by.PA = merged.by.PA
                                                                         , merged.by.run = merged.by.run
                                                                         , merged.by.algo = merged.by.algo
                                                                         , filtered.by = filtered.by
                                                                         , algo = algo))
              } else {
                keep_lines <- .filter_outputs.df(out, subset.list = list(full.name =  full.name, PA = PA
                                                                         , run = run, algo = algo))
              }
              out <- out[keep_lines, ]
              if (model.as.col) {
                out <- .transform_model.as.col(out)
              }
            }

            return(out)
          }
)


## --------------------------------------------------------------------------- #
## 6. BIOMOD.ensemble.models.out ---------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.ensemble.models.out
##' @aliases BIOMOD.ensemble.models.out-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_EnsembleModeling()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_EnsembleModeling}}, and used by 
##' \code{\link{BIOMOD_LoadModels}}, \code{\link{BIOMOD_PresenceOnly}} and 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' 
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the
##'   simulation set
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory
##'   variables
##' @slot models.out a \code{\link{BIOMOD.stored.models.out-class}} object
##'   containing informations from \code{\link{BIOMOD_Modeling}} object
##' @slot em.by a \code{character} corresponding to the way kept models have
##'   been combined to build the ensemble models, must be among
##'   \code{PA+run}, \code{PA+algo}, \code{PA},
##'   \code{algo}, \code{all}
##' @slot em.computed a \code{vector} containing names of ensemble models
##' @slot em.failed a \code{vector} containing names of failed ensemble models
# ##' @slot em.models_needed a \code{list} containing single models for each ensemble model
##' @slot em.models_kept a \code{list} containing single models for each ensemble model
##' @slot models.evaluation a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing models evaluation
##' @slot variables.importance a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing variables importance
##' @slot models.prediction a \code{\link{BIOMOD.stored.data.frame-class}} object
##'   containing models predictions
##' @slot models.prediction.eval a \code{\link{BIOMOD.stored.data.frame-class}}
##'   object containing models predictions for evaluation data
##' @slot link a \code{character} containing the file name of the saved object
##'   
##' @param object a \code{\link{BIOMOD.ensemble.models.out}} object
##' 
##' 
##' 
##' @seealso \code{\link{BIOMOD_EnsembleModeling}}, \code{\link{BIOMOD_LoadModels}}, 
##' \code{\link{BIOMOD_PresenceOnly}}, \code{\link{bm_VariablesImportance}}, 
##' \code{\link{bm_PlotEvalMean}}, \code{\link{bm_PlotEvalBoxplot}}, 
##' \code{\link{bm_PlotVarImpBoxplot}}, \code{\link{bm_PlotResponseCurves}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.ensemble.models.out")
##' 
##' ## ----------------------------------------------------------------------- #
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
##' ## ----------------------------------------------------------------------- #
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                        expl.var = myExpl,
##'                                        resp.xy = myRespXY,
##'                                        resp.name = myRespName)
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' }
##' 
##' 
##' ## ----------------------------------------------------------------------- #
##' # Model ensemble models
##' myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
##'                                       models.chosen = 'all',
##'                                       em.by = 'all',
##'                                       em.algo = c('EMmean', 'EMca'),
##'                                       metric.select = c('TSS'),
##'                                       metric.select.thresh = c(0.7),
##'                                       metric.eval = c('TSS', 'ROC'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' myBiomodEM
##' 
##' 
NULL

##' @name BIOMOD.ensemble.models.out-class
##' @rdname BIOMOD.ensemble.models.out
##' @export
##' 

# 6.1 Class Definition ---------------------------------------------------------

setClass("BIOMOD.ensemble.models.out",
         representation(modeling.id = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.out = 'BIOMOD.stored.models.out',
                        em.by = 'character',
                        em.computed = 'character',
                        em.failed = 'character',
                        em.models_kept = 'ANY',
                        models.evaluation = 'BIOMOD.stored.data.frame',
                        variables.importance = 'BIOMOD.stored.data.frame',
                        models.prediction = 'BIOMOD.stored.data.frame',
                        models.prediction.eval = 'BIOMOD.stored.data.frame',
                        link = 'character'),
         prototype(modeling.id = '.',
                   dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   models.out = new('BIOMOD.stored.models.out'),
                   em.by = character(),
                   em.computed = character(),
                   em.failed = character(),
                   em.models_kept = NULL,
                   models.evaluation = new('BIOMOD.stored.data.frame'),
                   variables.importance = new('BIOMOD.stored.data.frame'),
                   models.prediction = new('BIOMOD.stored.data.frame'),
                   models.prediction.eval = new('BIOMOD.stored.data.frame'),
                   link = ''),
         validity = function(object){ return(TRUE) })


# 6.3 Other functions ----------------------------------------------------------
## show.BIOMOD.ensemble.models.out ---------------------------------------------
##' 
##' @rdname BIOMOD.ensemble.models.out
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.ensemble.models.out'),
          function(object){
            .bm_cat("BIOMOD.ensemble.models.out")
            cat("\nsp.name :", object@sp.name, fill = .Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill = .Options$width)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed), fill = .Options$width)
            cat("\nmodels failed:", toString(object@em.failed), fill = .Options$width)
            .bm_cat()
          })
## get_formal_data.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_formal_data", "BIOMOD.ensemble.models.out",
          function(obj, subinfo = NULL) {
            if (is.null(subinfo)) {
              return(load_stored_object(obj@models.out))
            } else {
              bm_form = get_formal_data(obj)
              return(get_formal_data(bm_form, subinfo = subinfo))
            }
          }
)


## get_built_models.BIOMOD.ensemble.models.out ---------------------------------
##'
##' @rdname getters.out
##' @export
##' 

setMethod("get_built_models", "BIOMOD.ensemble.models.out",
          function(obj, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL)
          {
            out <- obj@em.computed
            keep_ind <- .filter_outputs.vec(out, obj.type = "em", subset.list = list(full.name = full.name
                                                                                     , merged.by.PA = merged.by.PA
                                                                                     , merged.by.run = merged.by.run
                                                                                     , merged.by.algo = merged.by.algo
                                                                                     , filtered.by = filtered.by
                                                                                     , algo = algo))
            out <- out[keep_ind]
            return(out)
          }
)


## get_kept_models.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_kept_models", "BIOMOD.ensemble.models.out", function(obj) { return(obj@em.models_kept) })


## get_predictions.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "BIOMOD.ensemble.models.out",
          function(obj, evaluation = FALSE, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL,
                   model.as.col = FALSE)
          {
            # check evaluation data availability
            if (evaluation && (!get_formal_data(obj)@has.evaluation.data)) {
              warning("!   Calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if (evaluation) { 
              out <- load_stored_object(obj@models.prediction.eval)
            } else { 
              out <- load_stored_object(obj@models.prediction)
            }
            
            # subselection of models_selected
            keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name
                                                                     , merged.by.algo = merged.by.algo
                                                                     , merged.by.run = merged.by.run
                                                                     , merged.by.PA = merged.by.PA
                                                                     , filtered.by = filtered.by
                                                                     , algo = algo))
            out <- out[keep_lines, ]
            if (model.as.col) {
              out <- .transform_model.as.col(out)
            }
            return(out)
          }
)


## get_evaluations.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_evaluations", "BIOMOD.ensemble.models.out",
          function(obj, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL, metric.eval = NULL)
          {
            out <- load_stored_object(obj@models.evaluation)
            if(nrow(out) == 0){
              cat("\n! models have no evaluations\n")
              return(invisible(NULL))
            } else {
              keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name
                                                                       , merged.by.algo = merged.by.algo
                                                                       , merged.by.run = merged.by.run
                                                                       , merged.by.PA = merged.by.PA
                                                                       , filtered.by = filtered.by
                                                                       , algo = algo
                                                                       , metric.eval = metric.eval))
              out <- out[keep_lines, ]
              return(out)
            }
          }
)

## get_variables_importance.BIOMOD.ensemble.models.out -------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_variables_importance", "BIOMOD.ensemble.models.out",
          function(obj, full.name = NULL, merged.by.algo = NULL, merged.by.run = NULL
                   , merged.by.PA = NULL, filtered.by = NULL, algo = NULL, expl.var = NULL)
          {
            out <- load_stored_object(obj@variables.importance)
            keep_lines <- .filter_outputs.df(out, subset.list = list(full.name = full.name
                                                                     , merged.by.algo = merged.by.algo
                                                                     , merged.by.run = merged.by.run
                                                                     , merged.by.PA = merged.by.PA
                                                                     , filtered.by = filtered.by
                                                                     , algo = algo
                                                                     , expl.var = expl.var))
            out <- out[keep_lines, ]
            return(out)
          }
)



## -------------------------------------------------------------------------- #
## 0. Generic Functions definition ------------------------------------------
## -------------------------------------------------------------------------- #
## Used for different classes 
##    A = BIOMOD.models.out, B = BIOMOD.projection.out, C = BIOMOD.ensemble.models.out

##' @name getters.out
##' @aliases get_options
##' @aliases get_calib_lines
##' @aliases get_formal_data
##' @aliases get_projected_models
##' @aliases free
##' @aliases get_predictions
##' @aliases get_needed_models
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
##' @param obj a \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.projection.out}} or 
##' \code{\link{BIOMOD.ensemble.models.out}} object
##' @param \ldots (\emph{optional, one or several of the following arguments depending on the selected 
##' function)}) 
##' @param as.data.frame a \code{logical} defining whether output should be returned as 
##' \code{data.frame} or \code{array} object
##' @param subinfo a \code{character} corresponding to the information to be extracted, must be 
##' among \code{NULL}, \code{expl.var.names}, \code{resp.var}, \code{expl.var}, \code{MinMax}, 
##' \code{eval.resp.var}, \code{eval.expl.var} (see \href{getters.out.html#details}{Details})
##' @param evaluation a \code{logical} defining whether evaluation data should be used or not
##' @param full.name a \code{vector} containing model names to be kept, must be either \code{all} 
##' or a sub-selection of model names
##' @param model a \code{character} corresponding to the model name, must be either \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}
##' @param run.eval a \code{vector} containing repetition set to be loaded, must be among 
##' \code{RUN1}, \code{RUN2}, \code{...}, \code{Full}
##' @param data.set a \code{vector} containing pseudo-absence set to be loaded, must be among 
##' \code{PA1}, \code{PA2}, \code{...}
##' @param selected.models a \code{vector} containing names of the needed models of a 
##' \code{\link{BIOMOD.ensemble.models.out}} object
##' 
##' 
##' 
##' @return 
##' 
##' \describe{
##'   \item{\code{get_options}}{a \code{\link{BIOMOD.stored.models.options}} object from the 
##'   \code{models.options} slot of a \code{\link{BIOMOD.models.out}} object}
##'   \item{\code{get_calib_lines}}{a \code{\link{BIOMOD.stored.array}} object from the 
##'   \code{calib.lines} slot of a \code{\link{BIOMOD.models.out}} object}
##'   
##'   \item{\code{get_projected_models}}{a \code{vector} from the \code{models.projected} slot of a 
##'   \code{\link{BIOMOD.projection.out}} object}
##'   
##'   \item{\code{get_predictions}}{a \code{\link{BIOMOD.stored.data}} object from the \code{proj.out} slot 
##'   of a \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.projection.out}} or 
##'   \code{\link{BIOMOD.ensemble.models.out}} object}
##'   
##'   \item{\code{get_needed_models}}{a \code{vector} containing names of the needed models of a 
##'   \code{\link{BIOMOD.ensemble.models.out}} object}
##'   \item{\code{get_kept_models}}{a \code{vector} containing names of the kept models of a 
##'   \code{\link{BIOMOD.ensemble.models.out}} object}
##'   
##'   \item{\code{get_formal_data}}{depending on the \code{subinfo} parameter :
##'   \describe{
##'     \item{\code{NULL}}{a \code{\link{BIOMOD.stored.formated.data}} (or 
##'     \code{\link{BIOMOD.stored.models.out}}) object from the \code{formated.input.data} (or 
##'     \code{models.out}) slot of a \code{\link{BIOMOD.models.out}} (or 
##'     \code{\link{BIOMOD.ensemble.models.out}}) object}
##'     
##'     \item{\code{expl.var.names}}{a \code{vector} from the \code{expl.var.names} slot of a 
##'     \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} object}

##'     \item{\code{resp.var}}{a \code{vector} from the \code{data.species} slot of the 
##'     \code{formated.input.data} slot of a \code{\link{BIOMOD.models.out}} or 
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{expl.var}}{a \code{data.frame} from the \code{data.env.var} slot of the 
##'     \code{formated.input.data} slot of a \code{\link{BIOMOD.models.out}} or 
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{MinMax}}{a \code{list} of minimum and maximum values (or levels if factorial) of 
##'     variable contained in the \code{data.env.var} slot of the 
##'     \code{formated.input.data} slot of a \code{\link{BIOMOD.models.out}} or 
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{eval.resp.var}}{a \code{vector} from the \code{eval.data.species} slot of the 
##'     \code{formated.input.data} slot of a \code{\link{BIOMOD.models.out}} or 
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'     
##'     \item{\code{eval.expl.var}}{a \code{data.frame} from the \code{eval.data.env.var} slot of the 
##'     \code{formated.input.data} slot of a \code{\link{BIOMOD.models.out}} or 
##'     \code{\link{BIOMOD.ensemble.models.out}} object}
##'   }
##'   }
##'   \item{\code{get_built_models}}{a \code{vector} from the \code{models.computed} slot (or 
##'   \code{em.computed}) of a \code{\link{BIOMOD.models.out}} (or 
##'   \code{\link{BIOMOD.ensemble.models.out}}) object}
##'   \item{\code{get_evaluations}}{a \code{\link{BIOMOD.stored.array}} (or \code{matrix}) from the 
##'   \code{models.evaluation} slot (or \code{model_evaluation} of each model in 
##'   \code{em.computed}) of a \code{\link{BIOMOD.models.out}} (or 
##'   \code{\link{BIOMOD.ensemble.models.out}}) object}
##'   \item{\code{get_variables_importance}}{a \code{\link{BIOMOD.stored.array}} from the 
##'   \code{variables.importance} slot (or \code{model_variables_importance} of each model in 
##'   \code{em.models}) of a \code{\link{BIOMOD.models.out}} (or 
##'   \code{\link{BIOMOD.ensemble.models.out}}) object}
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
##' 
NULL

setGeneric("get_options", function(obj, ...) { standardGeneric("get_options") }) ## A
setGeneric("get_calib_lines", function(obj, ...) { standardGeneric("get_calib_lines") }) ## A

setGeneric("get_projected_models", function(obj, ...) { standardGeneric("get_projected_models") }) ## B
setGeneric("free", function(obj, ...) { standardGeneric("free") }) ## B

setGeneric("get_predictions", function(obj, ...) { standardGeneric("get_predictions") }) ## ABC

setGeneric("get_needed_models", function(obj, ...) { standardGeneric("get_needed_models") }) ## C
setGeneric("get_kept_models", function(obj, ...) { standardGeneric("get_kept_models") }) ## C

setGeneric("get_formal_data", function(obj, ...) { standardGeneric("get_formal_data") }) ## AC
setGeneric("get_built_models", function(obj, ...) { standardGeneric("get_built_models") }) ## AC
setGeneric("get_evaluations", function(obj, ...) { standardGeneric("get_evaluations") }) ## AC
setGeneric("get_variables_importance", function(obj, ...) { standardGeneric("get_variables_importance") }) ## AC

## -------------------------------------------------------------------------- #
## 4. BIOMOD.models.out -----------------------------------------------------
## -------------------------------------------------------------------------- #

##' @name BIOMOD.models.out
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
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory variables
##' @slot models.computed a \code{vector} containing names of computed models
##' @slot models.failed a \code{vector} containing names of failed models
##' @slot has.evaluation.data a \code{logical} value defining whether evaluation data is given
##' @slot scale.models a \code{logical} value defining whether models have been rescaled or 
##' not
##' @slot formated.input.data a \code{\link{BIOMOD.stored.formated.data}} object containing 
##' informations from \code{\link{BIOMOD_FormatingData}} object
##' @slot calib.lines a \code{\link{BIOMOD.stored.array}} object containing calibration lines
##' @slot models.options a \code{\link{BIOMOD.stored.models.options}} object containing 
##' informations from \code{\link{BIOMOD_ModelingOptions}} object
##' @slot models.evaluation a \code{\link{BIOMOD.stored.array}} object containing models evaluation
##' @slot variables.importance a \code{\link{BIOMOD.stored.array}} object containing variables 
##' importance
##' @slot models.prediction a \code{\link{BIOMOD.stored.array}} object containing models 
##' predictions
##' @slot models.prediction.eval a \code{\link{BIOMOD.stored.array}} object containing models 
##' predictions for evaluation data
##' @slot link a \code{character} containing the file name of the saved object
##' @slot val a \code{\link{BIOMOD.models.out}} object
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
##' showClass("BIOMOD.stored.models.out")
##' 
##' ## ----------------------------------------------------------------------- #
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
##' 
##' ## ----------------------------------------------------------------------- #
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
##' ## ----------------------------------------------------------------------- #
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     bm.options = myBiomodOptions,
##'                                     nb.rep = 2,
##'                                     data.split.perc = 80,
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 3,
##'                                     do.full.models = FALSE)
##' myBiomodModelOut
##' 
##' 
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
                        calib.lines = 'BIOMOD.stored.array',
                        models.options = 'BIOMOD.stored.models.options',
                        models.evaluation = 'BIOMOD.stored.array',
                        variables.importance = 'BIOMOD.stored.array',
                        models.prediction = 'BIOMOD.stored.array',
                        models.prediction.eval = 'BIOMOD.stored.array',
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
                   calib.lines = new('BIOMOD.stored.array'),
                   models.options = new('BIOMOD.stored.models.options'),
                   models.evaluation = new('BIOMOD.stored.array'),
                   variables.importance = new('BIOMOD.stored.array'),
                   models.prediction = new('BIOMOD.stored.array'),
                   models.prediction.eval = new('BIOMOD.stored.array'),
                   link = ''),
         validity = function(object){ return(TRUE) } )

setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) } )


# 4.3 Other functions ------------------------------------------------------
## show.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname BIOMOD.models.out
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
            if (obj@models.options@inMemory) {
              return(obj@models.options@val)
            } else if (obj@models.options@link != '') {
              return(get(load(obj@models.options@link)))
            } else { return(NA) }
          }
)

## get_calib_lines.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##'

setMethod("get_calib_lines", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, ...) {
            calib_lines <- load_stored_object(obj@calib.lines)
            return(calib_lines)
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
              if (obj@formated.input.data@inMemory) {
                return(obj@formated.input.data@val)
              } else if (obj@formated.input.data@link != '') {
                data <- get(load(obj@formated.input.data@link))
                return(data)
              } else { return(NA) }
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
          function(obj, as.data.frame = FALSE, evaluation = FALSE)
          {
            # check evaluation data availability
            if (evaluation & (!obj@has.evaluation.data)) {
              warning("calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if (evaluation) { 
              pred <- obj@models.prediction.eval 
            } else { 
              pred <- obj@models.prediction 
            }
            
            if (!as.data.frame) {
              if (pred@inMemory) {
                return(pred@val)
              } else if (pred@link != '') {
                return(get(load(pred@link)))
              } else {
                return(NULL)
              }
            } else {
              if (pred@inMemory) {
                mod.pred <- as.data.frame(pred@val)
              } else if (pred@link != '') {
                mod.pred <- as.data.frame(get(load(pred@link)))
              } else { return(NULL) }
              
              names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred), ".", fixed = TRUE), function(x)
              {
                x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo at the end
                data.set.id <- x.rev[1]
                cross.valid.id <- x.rev[2]
                algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
                model.id <- paste(obj@sp.name,
                                  data.set.id,
                                  cross.valid.id,
                                  algo.id, sep = "_")
                return(model.id)
              }))
              return(mod.pred)
            }
          }
)

## get_built_models.BIOMOD.models.out ---------------------------------------------------
##' @rdname getters.out
##' @export
##' 

setMethod("get_built_models", "BIOMOD.models.out", function(obj, ...) { return(obj@models.computed) })

## get_evaluations.BIOMOD.models.out ---------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_evaluations", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, ...) {
            out <- NULL
            if (obj@models.evaluation@inMemory) {
              out <- obj@models.evaluation@val
            } else if (obj@models.evaluation@link != '') {
              out <- get(load(obj@models.evaluation@link))
            }
            
            ## transform into data.frame object if needed
            if (as.data.frame) {
              tmp = melt(out, varnames = c("Eval.metric", "tmp", "Algo", "Run", "Dataset"))
              tmp$Model.name = paste0(tmp$Algo, "_", tmp$Run, "_", tmp$Dataset)
              tmp.split = split(tmp, tmp$tmp)
              tmp.split = lapply(tmp.split, function(x) x[, c("Model.name", "Algo", "Run", "Dataset", "Eval.metric", "value")])
              for (i in 1:length(tmp.split)) { colnames(tmp.split[[i]])[6] = names(tmp.split)[i] }
              out = Reduce(function(x, y) merge(x, y, by = c("Model.name", "Algo", "Run", "Dataset", "Eval.metric")), tmp.split)
            }
            
            return(out)
          }
)

## get_variables_importance.BIOMOD.models.out ---------------------------------------------------
##' @rdname getters.out
##' @export
##' 

setMethod("get_variables_importance", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, ...) {
            out <- NULL
            # if (obj@variables.importance@inMemory) {
            #   out <- obj@variables.importance@val
            # } else if(obj@variables.importance@link != '') {
            #   out <- get(load(obj@variables.importance@link))
            # }
            BIOMOD_LoadModels(bm.out = obj)
            for (mod in get_built_models(obj)) {
              out_tmp <- get(mod)@model_variables_importance
              out <- abind(out, out_tmp, along = 3)
            }
            dimnames(out)[[3]] <- get_built_models(obj)
            
            if (!is.null(out) && as.data.frame == TRUE) {
              tmp <- melt(out, varnames = c("Expl.var", "Rand", "L1"))
              tmp$Model.name <-
                sapply(as.character(tmp$L1),
                       function(x) { 
                         paste(strsplit(x, "_")[[1]][-1], collapse = "_") 
                       })
              tmp$Algo <- 
                sapply(tmp$Model.name, 
                       function(x) { 
                         strsplit(x, "_")[[1]][3] 
                       })
              tmp$Run <-
                sapply(tmp$Model.name, 
                       function(x) {
                         strsplit(x, "_")[[1]][2] 
                       })
              tmp$Dataset <-
                sapply(tmp$Model.name, 
                       function(x) {
                         strsplit(x, "_")[[1]][1]
                       })
              out <- tmp[, c("Model.name", "Algo", "Run",
                             "Dataset", "Expl.var", "Rand", "value")]
              colnames(out)[7] = "Var.imp"
            }
            
            return(out)
          }
)



## --------------------------------------------------------------------------- #
## 5. BIOMOD.projection.out --------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.projection.out
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
##' @param str.grep a \code{character} corresponding to the name (or part of the name(s)) of 
##' models projected
##' @param col a \code{vector} containing colors for plot (default : 
##' \code{colorRampPalette(c("grey90", "yellow4", "green4"))(100)})
##' @param x a \code{\link{BIOMOD.projection.out}} object
##' @param object a \code{\link{BIOMOD.projection.out}} object
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
##'   # Create default modeling options
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       bm.options = myBiomodOptions,
##'                                       nb.rep = 2,
##'                                       data.split.perc = 80,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       do.full.models = FALSE)
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
##' @importFrom raster subset
## @importFrom rasterVis levelplot
##' @importFrom grDevices colorRampPalette colors dev.new gray rainbow
##' @importFrom graphics layout legend par points polygon text
##' 
##' @export
##' 

# 5.1 Class Definition  -----------------------------------

setClass("BIOMOD.projection.out",
         representation(modeling.id = 'character',
                        proj.name = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        coord = 'matrix',
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
                   coord = matrix(),
                   scale.models = TRUE,
                   models.projected = '',
                   type = ''),
         validity = function(object){ return(TRUE) })


# 5.3 Other functions ---------------------------------------------------------
## plot.BIOMOD.projection.out -------------------------------------------------
##' 
##' @rdname BIOMOD.projection.out
##' @export
##' 

setMethod(
  'plot', signature(x = 'BIOMOD.projection.out', y = "missing"),
  function(x, col = NULL, str.grep = NULL){
    models_selected <- x@models.projected
    if (length(str.grep)) { 
      models_selected <- grep(paste(str.grep, collapse = "|"),
                              models_selected, value = TRUE)
    }
    if (!length(models_selected)) { 
      stop("invalid str.grep arg")
    }
    
    if (inherits(x@proj.out, "BIOMOD.stored.raster.stack")) {
      requireNamespace("rasterVis")
      maxi <- 
        try(
          cellStats(
            get_predictions(x, full.name = models_selected), 
            max
          )
        )
      maxi <- max(maxi)
      maxi <- ifelse(maxi <= 1, 1, ifelse(maxi < 1000, 1000, maxi))
      my.at <- seq(0, maxi, by = 100 * maxi / 1000) ## breaks of color key
      my.labs.at <- seq(0, maxi, by = 250 * maxi / 1000) ## labels placed vertically centered
      my.lab <- seq(0, maxi, by = 250 * maxi / 1000) ## labels
      my.col <- colorRampPalette(c("grey90", "yellow4", "green4"))(100) ## colors
      
      ## try to use levelplot function
      try_plot <- try(
        rasterVis::levelplot(
          get_predictions(x, full.name = models_selected),
          at = my.at,
          margin = TRUE,
          col.regions = my.col,
          main = paste(x@sp.name, x@proj.name, "projections"),
          colorkey = list(labels = list(labels = my.lab, 
                                        at = my.labs.at))
        )
      )
      if (!inherits(try_plot, "try-error")) { ## produce plot
        print(try_plot)
      } else { ## try classical plot
        cat("\nrasterVis' levelplot() function failed. Try to call standard raster plotting function.",
            "It can lead to unoptimal representations.",
            "You should try to do it by yourself extracting predicions (see : get_predictions() function).",
            fill = options()$width)
        try_plot <- try(plot(get_predictions(x, full.name = models_selected)))
        if (inherits(try_plot,"try-error")) {
          cat("\n Plotting function failed.. You should try to do it by yourself!")
        }
      }
    } else if (inherits(x@proj.out, "BIOMOD.stored.array")) {
      if (ncol(x@coord) != 2) {
        cat("\n ! Impossible to plot projections because xy coordinates are not available !")
      } else {
        .multiple.plot(Data = get_predictions(x, full.name = models_selected, as.data.frame = TRUE)
                       , coor = x@coord)
      }
    } else {
      cat("\n !  Biomod Projection plotting issue !", fill = .Options$width)
    }
  }
)

## show.BIOMOD.projection.out -------------------------------------------------
##' 
##' @rdname BIOMOD.projection.out
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
            .bm_cat()
          })

## get_projected_models.BIOMOD.projection.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_projected_models", "BIOMOD.projection.out", function(obj){ return(obj@models.projected) })

## free.BIOMOD.projection.out --------------------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod('free', signature('BIOMOD.projection.out'),
          function(obj) {
            if (inherits(obj@proj.out, "BIOMOD.stored.array")) {
              obj@proj.out@val  <- array()
            } else if (inherits(obj@proj.out, "BIOMOD.stored.raster.stack")) {
              obj@proj.out@val <- stack()
            } else {
              obj@proj.out@val <- NULL
            }
            obj@proj.out@inMemory <- FALSE
            return(obj)
          })

## get_predictions.BIOMOD.projection.out ---------------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "BIOMOD.projection.out",
          function(obj, as.data.frame = FALSE, full.name = NULL, 
                   model = NULL, run.eval = NULL, data.set = NULL) {
            models_selected <- get_projected_models(obj)
            if (length(full.name)) {
              models_selected <- intersect(full.name, models_selected)
            } else if (length(model) | length(run.eval) | length(data.set)) {
              # models subselection according to model, run.eval and data.set parameters
              grep_model = grep_run.eval = grep_data.set = "*"
              if (length(model)) { 
                grep_model <- paste0("(", paste(model, collapse = "|"), ")")
              }
              if (length(run.eval)) { 
                grep_run.eval <- paste0("(", paste(run.eval, collapse = "|"), ")") 
              }
              if (length(data.set)) { 
                grep_data.set <- paste0("(", paste(data.set, collapse = "|"), ")") 
              }
              grep_full <- paste0(grep_data.set, "_", grep_run.eval, "_", grep_model, "$")
              models_selected <- grep(pattern = grep_full, models_selected, value = TRUE)
            }
            
            if (length(models_selected)) {
              proj <- load_stored_object(obj@proj.out)
              names(proj) <- obj@models.projected
              if (inherits(proj, 'Raster')) {
                proj <- subset(proj, models_selected, drop = FALSE)
              } else if (length(dim(proj)) == 4) { ## 4D arrays
                proj <- proj[, .extract_modelNamesInfo(model.names = models_selected, info = 'models'),
                             .extract_modelNamesInfo(model.names = models_selected, info = 'run.eval'),
                             .extract_modelNamesInfo(model.names = models_selected, info = 'data.set'), drop = FALSE]
              } else { ## matrix (e.g. from ensemble models projections)
                proj <- proj[, models_selected, drop = FALSE]
              }
              
              if (as.data.frame) {
                proj <- as.data.frame(proj)
                ## set correct names
                if (obj@type == 'array' &
                    sum(!(names(proj) %in% models_selected)) > 0) { ## from array & not valid names
                  names(proj) <- unlist(
                    lapply(strsplit(names(proj), ".", fixed = TRUE), 
                           function(x){
                             x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                             data.set.id <- x.rev[1]
                             cross.valid.id <- x.rev[2]
                             algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
                             model.id <- paste(obj@sp.name,
                                               data.set.id,
                                               cross.valid.id,
                                               algo.id, sep = "_")
                             return(model.id)
                           }))
                }
                proj <- proj[, models_selected, drop = FALSE] # reorder the data.frame
              }
            } else { 
              proj <- NULL 
            }
            return(proj)
          }
)


## --------------------------------------------------------------------------- #
## 6. BIOMOD.ensemble.models.out ---------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.ensemble.models.out
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_EnsembleModeling()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_EnsembleModeling}}, and used by 
##' \code{\link{BIOMOD_LoadModels}}, \code{\link{BIOMOD_PresenceOnly}} and 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' 
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot expl.var.names a \code{vector} containing names of explanatory variables
##' @slot models.out a \code{\link{BIOMOD.stored.models.out}} object containing 
##' informations from \code{\link{BIOMOD_Modeling}} object
##' @slot em.computed a \code{vector} containing names of ensemble models
##' @slot em.by a \code{character} corresponding to the way kept models have been combined to 
##' build the ensemble models, must be among \code{PA_dataset+repet}, \code{PA_dataset+algo}, 
##' \code{PA_dataset}, \code{algo}, \code{all}
##' @slot em.models a \code{list} containing ensemble models
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the simulation set
##' @slot link a \code{character} containing the file name of the saved object
##' 
##' @param object a \code{\link{BIOMOD.ensemble.models.out}} object
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
##'   # Create default modeling options
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       bm.options = myBiomodOptions,
##'                                       nb.rep = 2,
##'                                       data.split.perc = 80,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       do.full.models = FALSE)
##' }
##' 
##' 
##' ## ----------------------------------------------------------------------- #
##' # Model ensemble models
##' myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
##'                                       models.chosen = 'all',
##'                                       em.by = 'all',
##'                                       metric.select = c('TSS'),
##'                                       metric.select.thresh = c(0.7),
##'                                       metric.eval = c('TSS', 'ROC'),
##'                                       var.import = 3,
##'                                       prob.mean = TRUE,
##'                                       prob.median = TRUE,
##'                                       prob.cv = TRUE,
##'                                       prob.ci = TRUE,
##'                                       prob.ci.alpha = 0.05,
##'                                       committee.averaging = TRUE,
##'                                       prob.mean.weight = TRUE,
##'                                       prob.mean.weight.decay = 'proportional')
##' myBiomodEM
##' 
##' 
##' @export
##' 

# 6.1 Class Definition ---------------------------------------------------------

setClass("BIOMOD.ensemble.models.out",
         representation(dir.name = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.out = 'BIOMOD.stored.models.out',
                        em.computed = 'character',
                        em.by = 'character',
                        em.models = 'ANY',
                        modeling.id = 'character',
                        link = 'character'),
         prototype(dir.name = '.',
                   sp.name = '',
                   expl.var.names = '',
                   models.out = new('BIOMOD.stored.models.out'),
                   em.models = NULL,
                   em.computed = character(),
                   modeling.id = '.',
                   link = ''),
         validity = function(object){ return(TRUE) })


# 6.3 Other functions ----------------------------------------------------------
## show.BIOMOD.ensemble.models.out ---------------------------------------------
##' 
##' @rdname BIOMOD.ensemble.models.out
##' @export
##' 

setMethod('show', signature('BIOMOD.ensemble.models.out'),
          function(object){
            .bm_cat("BIOMOD.ensemble.models.out")
            cat("\nsp.name :", object@sp.name, fill = .Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill = .Options$width)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed), fill = .Options$width)
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
              if (obj@models.out@inMemory) {
                return(obj@models.out@val)
              } else if (obj@models.out@link != '') {
                data <- get(load(obj@models.out@link))
                return(data)
              } else { return(NA) }
            } else {
              bm_form = get(load(get_formal_data(obj)@formated.input.data@link))
              if (subinfo == 'MinMax') {
                env = as.data.frame(bm_form@data.env.var)
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
                return(as.data.frame(bm_form@data.env.var))
              } else if (subinfo == 'expl.var.names') {
                return(obj@expl.var.names)
              } else if (subinfo == 'resp.var') {
                return(as.numeric(bm_form@data.species))
              } else if (subinfo == 'eval.resp.var') {
                return(as.numeric(bm_form@eval.data.species))
              } else if (subinfo == 'eval.expl.var') {
                return(as.data.frame(bm_form@eval.data.env.var))
              } else { stop("Unknown subinfo tag")}
            }
          }
)

## get_needed_models.BIOMOD.ensemble.models.out --------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_needed_models", "BIOMOD.ensemble.models.out",
          function(obj, selected.models = 'all', ...) {
            add.args <- list(...)
            if (selected.models[[1]] == "all") {
              selected.index <- c(1:length(obj@em.models))
            } else {
              selected.index <- which(names(obj@em.models) %in% selected.models)
            }
            needed_models <- lapply(obj@em.models[selected.index], function(x) { return(x@model) })
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          })


## get_kept_models.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_kept_models", "BIOMOD.ensemble.models.out",
          function(obj, model, ...) {
            if (is.character(model) | is.numeric(model)) {
              return(obj@em.models[[model]]@model)
            } else {
              kept_mod <- lapply(obj@em.models, function(x) { return(x@model) })
              names(kept_mod) <- names(obj@em.models)
              return(kept_mod)
            }
          })

## get_predictions.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_predictions", "BIOMOD.ensemble.models.out",
          function(obj, ...) {
            ## note: ensemble models predictions are stored within the directory
            ##  <dir.name>/<sp.name>/.BIOMOD_DATA/<modeling.id>/ensemble.models/ensemble.models.projections/
            ##  This function is just a friendly way to load this data
            
            ## get the path to projections files we want to load
            files.to.load <- file.path(obj@dir.name, obj@sp.name, ".BIOMOD_DATA", obj@modeling.id, "ensemble.models",
                                       "ensemble.models.predictions", paste0(obj@em.computed, ".predictions"))
            ## load and merge projection files within a data.frame
            bm.pred <- do.call(cbind, lapply(files.to.load, function(ftl) get(load(ftl))))
            colnames(bm.pred) <- obj@em.computed
            return(bm.pred)
          })

## get_built_models.BIOMOD.ensemble.models.out ---------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_built_models", "BIOMOD.ensemble.models.out", 
          function(obj, ...){ 
            return(obj@em.computed) 
          })

## get_evaluations.BIOMOD.ensemble.models.out ----------------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_evaluations", "BIOMOD.ensemble.models.out",
          function(obj, as.data.frame = FALSE, ...) {
            ## extract evaluation scores as a list
            out <- list()
            models <- obj@em.computed ## list of computed models
            for (mod in models) {
              out[[mod]] <- obj@em.models[[mod]]@model_evaluation[, , drop = FALSE]
            }
            
            ## transform into data.frame object if needed
            if(as.data.frame) {
              tmp = melt(out, varnames = c("Eval.metric", "tmp"))
              tmp$Model.name = sapply(tmp$L1, function(x) { paste(strsplit(x, "_")[[1]][-1], collapse = "_") })
              tmp$Model = sapply(tmp$Model.name, function(x) { strsplit(x, "_")[[1]][1] })
              tmp$Algo = sapply(tmp$Model.name, function(x) { strsplit(x, "_")[[1]][2] })
              tmp$Run = sapply(tmp$Model.name, function(x) { strsplit(x, "_")[[1]][3] })
              tmp$Dataset = sapply(tmp$Model.name, function(x) { strsplit(x, "_")[[1]][4] })
              tmp.split = split(tmp, tmp$tmp)
              tmp.split = lapply(tmp.split, function(x) x[, c("Model.name", "Model", "Algo", "Run", "Dataset", "Eval.metric", "value")])
              for (i in 1:length(tmp.split)) { colnames(tmp.split[[i]])[7] = names(tmp.split)[i] }
              out = Reduce(function(x, y) merge(x, y, by = c("Model.name", "Model", "Algo", "Run", "Dataset", "Eval.metric")), tmp.split)
            }
            
            return(out)
          }
)

## get_variables_importance.BIOMOD.ensemble.models.out -------------------------
##' 
##' @rdname getters.out
##' @export
##' 

setMethod("get_variables_importance", "BIOMOD.ensemble.models.out",
          function(obj, as.data.frame = FALSE, ...) {
            out <- NULL
            for (mod in get_built_models(obj)) {
              out_tmp <- obj@em.models[[mod]]@model_variables_importance
              out <- abind(out, out_tmp, along = 3)
            }
            dimnames(out)[[3]] <- get_built_models(obj)
            
            if(!is.null(out) && as.data.frame == TRUE) {
              tmp <- melt(out, varnames = c("Expl.var", "Rand", "L1"))
              tmp$Model.name <- 
                sapply(as.character(tmp$L1), 
                       function(x) {
                         paste(strsplit(x, "_")[[1]][-1], collapse = "_") 
                       })
              tmp$Model <-  
                sapply(tmp$Model.name, 
                       function(x) { 
                         strsplit(x, "_")[[1]][1]
                       })
              tmp$Algo <-  
                sapply(tmp$Model.name,
                       function(x) { 
                         strsplit(x, "_")[[1]][2] 
                       })
              tmp$Run <- 
                sapply(tmp$Model.name,
                       function(x) { 
                         strsplit(x, "_")[[1]][3] 
                       })
              tmp$Dataset <- 
                sapply(tmp$Model.name, 
                       function(x) { 
                         strsplit(x, "_")[[1]][4] 
                       })
              out <- tmp[, c("Model.name", "Model", "Algo",
                             "Run", "Dataset", "Expl.var",
                             "Rand", "value")]
              colnames(out)[8] = "Var.imp"
            }
            
            return(out)
          }
)


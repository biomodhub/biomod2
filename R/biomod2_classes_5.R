
##' @name predict.em
## @aliases predict
##' @aliases .predict.EMmean_biomod2_model.RasterStack
##' @aliases .predict.EMmean_biomod2_model.data.frame
##' @aliases .predict.EMmedian_biomod2_model.RasterStack
##' @aliases .predict.EMmedian_biomod2_model.data.frame
##' @aliases .predict.EMcv_biomod2_model.RasterStack
##' @aliases .predict.EMcv_biomod2_model.data.frame
##' @aliases .predict.EMci_biomod2_model.RasterStack
##' @aliases .predict.EMci_biomod2_model.data.frame
##' @aliases .predict.EMca_biomod2_model.RasterStack
##' @aliases .predict.EMca_biomod2_model.data.frame
##' @aliases .predict.EMwmean_biomod2_model.RasterStack
##' @aliases .predict.EMwmean_biomod2_model.data.frame
##' @author Damien Georges
##' 
##' @title Functions to get predictions from \code{\link{biomod2_ensemble_model}} objects
##' 
##' @description This function allows the user to predict single models from 
##' \code{\link{biomod2_ensemble_model}} on (new) explanatory variables.
##' 
##' 
##' @param object a \code{\link{biomod2_ensemble_model}} object
##' @param newdata a \code{data.frame} or \code{\link[raster:stack]{RasterStack}} object 
##' containing data for new predictions
##' @param formal_predictions a \code{matrix} containing formal predictions
##' @param \ldots (\emph{optional)}) 
##' 
##' 
##' @seealso \code{\link{biomod2_ensemble_model}}
##' @family Toolbox functions
##' 
##' 
##' @importFrom raster calc reclassify cv
##' 
NULL

#setGeneric("predict", def = function(object, ...) { standardGeneric("predict") })

# ---------------------------------------------------------------------------- #
# 9. biomod2_ensemble_model --------------------------------------------------
# ---------------------------------------------------------------------------- #

##' @name biomod2_ensemble_model
##' @aliases biomod2_ensemble_model-class
##' @aliases biomod2_ensemble_model
##' @aliases EMmean_biomod2_model
##' @aliases EMmedian_biomod2_model
##' @aliases EMcv_biomod2_model
##' @aliases EMci_biomod2_model
##' @aliases EMca_biomod2_model
##' @aliases EMwmean_biomod2_model
##' @author Damien Georges
##' 
##' @title Ensemble model output object class (when running \code{BIOMOD_EnsembleModeling()})
##' 
##' @description Class created by \code{\link{BIOMOD_EnsembleModeling}}
##' 
##' 
##' @slot modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' @slot model_name a \code{character} corresponding to the model name
##' @slot model_class a \code{character} corresponding to the model class
##' @slot model_options a \code{list} containing the model options
##' @slot model the corresponding model object
##' @slot scaling_model the corresponding scaled model object
##' @slot dir_name a \code{character} corresponding to the modeling folder
##' @slot resp_name a \code{character} corresponding to the species name
##' @slot expl_var_names a \code{vector} containing names of explanatory variables
##' @slot expl_var_type a \code{vector} containing classes of explanatory variables
##' @slot expl_var_range a \code{list} containing ranges of explanatory variables
##' @slot model_evaluation a \code{matrix} containing the model evaluations
##' @slot model_variables_importance a \code{matrix} containing the model variables importance
##' 
##' @details 
##' 
##' \code{biomod2_model} is the basic object for \pkg{biomod2} ensemble species distribution models. 
##' \cr All listed classes below are derived from \code{biomod2_model}, and have a 
##' \code{model_class} slot specific value :
##' 
##' \itemize{
##'   \item{\code{biomod2_ensemble_model} : }{\code{model_class} is \code{EM}}
##'   \item{\code{EMmean_biomod2_model} : }{\code{model_class} is \code{EMmean}}
##'   \item{\code{EMmedian_biomod2_model} : }{\code{model_class} is \code{EMmedian}}
##'   \item{\code{EMcv_biomod2_model} : }{\code{model_class} is \code{EMcv}}
##'   \item{\code{EMci_biomod2_model} : }{\code{model_class} is \code{EMci}}
##'   \item{\code{EMca_biomod2_model} : }{\code{model_class} is \code{EMca}}
##'   \item{\code{EMwmean_biomod2_model} : }{\code{model_class} is \code{EMwmean}}
##' }
##' 
##' 
##' @seealso \code{\link{biomod2_model}}, \code{\link{BIOMOD_EnsembleModeling}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("biomod2_ensemble_model")
##' showClass("EMmean_biomod2_model")
##' showClass("EMmedian_biomod2_model")
##' showClass("EMcv_biomod2_model")
##' showClass("EMci_biomod2_model")
##' showClass("EMca_biomod2_model")
##' showClass("EMwmean_biomod2_model")
##' 
##' 
NULL

##' @name biomod2_ensemble_model-class
##' @rdname biomod2_ensemble_model
##' @export
##' 

## 9.1 Class Definition ---------------------------------------------------------
setClass('biomod2_ensemble_model',
         representation(modeling.id = 'character'), ##maybe some additional args should be added here
         contains = 'biomod2_model',
         prototype = list(model_class = 'EM'),
         validity = function(object) { return(TRUE) })

### biomod2_ensemble_model predict2 method
### 
### biomod2_ensemble_model + Raster  -------------------------------------------------
setMethod('predict2', signature(object = 'biomod2_ensemble_model', newdata = "RasterStack"),
          function(object, newdata, predfun, seedval = NULL,
                   data_as_formal_predictions = FALSE, ...) {
            args <- list(...)
            filename <- args$filename
            on_0_1000 <- args$on_0_1000
            # additional arg retrieved for EMci
            sd_prediction <- args$sd_prediction
            mean_prediction <- args$mean_prediction
            side <- args$side
            # additional arg retrived for EMca
            thresh <- args$thresh
            # additional arg retrived for EMwmean
            penalization_scores <- args$penalization_scores
            
            if (is.null(filename)) { 
              filename <- "" 
            }
            if (is.null(on_0_1000)) { 
              on_0_1000 <- FALSE
            }
            
            if (!data_as_formal_predictions) {
              newdata <- .template_predictEM.formal_predictions(object, newdata, on_0_1000 = on_0_1000, seedval = seedval)
            }
            
            out <- predfun(newdata, 
                           on_0_1000 = on_0_1000,
                           mean_prediction = mean_prediction,
                           sd_prediction = sd_prediction, 
                           side = side,
                           thresh = thresh, 
                           penalization_scores = penalization_scores)
            if (!is.null(out)) {
              writeRaster(out, filename = filename, overwrite = TRUE)
            }
            return(out)
            
          })

### biomod2_ensemble_model + data.frame  -------------------------------------------------
setMethod('predict2', signature(object = 'biomod2_ensemble_model', newdata = "data.frame"),
          function(object, newdata, predfun, seedval = NULL, 
                   data_as_formal_predictions = FALSE, ...) {
            
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            if (is.null(on_0_1000)) {
              on_0_1000 <- FALSE 
            }
            # additional arg retrieved for EMci
            sd_prediction <- args$sd_prediction
            mean_prediction <- args$mean_prediction
            side <- args$side
            # additional arg retrived for EMca
            thresh <- args$thresh
            # additional arg retrived for EMwmean
            penalization_scores <- args$penalization_scores
            
            if (!data_as_formal_predictions) {
              newdata = .template_predictEM.formal_predictions(object, newdata, on_0_1000 = on_0_1000, seedval = seedval)
            }
            out <- predfun(newdata,
                           on_0_1000 = on_0_1000,
                           mean_prediction = mean_prediction,
                           sd_prediction = sd_prediction,
                           side = side,
                           thresh = thresh, 
                           penalization_scores = penalization_scores)
            return(out)
            
          })


### -------------------------------------------------------------------------- #
### 10.1 EMmean_biomod2_model ------------------------------------------------
### -------------------------------------------------------------------------- #


setClass('EMmean_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmean'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict2', signature(object = 'EMmean_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, ...){
              calc(newdata,function(x){
                m <- mean(x)
                if (on_0_1000) { 
                  m <- round(m)
                }
                return(m)
              })
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

setMethod('predict2', signature(object = 'EMmean_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, ...){
              out <- rowMeans(newdata, na.rm = TRUE)
              if (on_0_1000) { 
                out <- round(out) 
              }
              out
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


### -------------------------------------------------------------------------- #
### 10.2 EMmedian_biomod2_model ----------------------------------------------
### -------------------------------------------------------------------------- #

setClass('EMmedian_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 


setMethod('predict2', signature(object = 'EMmedian_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, ...){
              calc(newdata,function(x){
                m <- median(x)
                if (on_0_1000) { 
                  m <- round(m)
                }
                return(m)
              })
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

setMethod('predict2', signature(object = 'EMmedian_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, ...){
              out <- apply(newdata, 1, median, na.rm = TRUE)
              if (on_0_1000) { 
                out <- round(out) 
              }
              out
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

### -------------------------------------------------------------------------- #
### 10.3 EMcv_biomod2_model --------------------------------------------------
### -------------------------------------------------------------------------- #

setClass('EMcv_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMcv'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict2', signature(object = 'EMcv_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, ...){
              if (nlayers(newdata) > 1) {
                out <- calc(newdata, cv, na.rm = TRUE, aszero = TRUE)
                return(out)
              } else {
                warning(paste0("\n Model EMcv was not computed because only one single model was kept in ensemble modeling ("
                               , names(newdata), ")"))
                return(NULL)
              }
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

setMethod('predict2', signature(object = 'EMcv_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, ...){
              if (ncol(newdata) > 1) {
                out <- apply(newdata, 1, cv, na.rm = TRUE, aszero = TRUE)
                return(out)
              } else {
                warning(paste0("\n Model EMcv was not computed because only one single model was kept in ensemble modeling ("
                               , colnames(newdata), ")"))
                return(NULL)
              }
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

### -------------------------------------------------------------------------- #
### 10.4 EMci_biomod2_model --------------------------------------------------
### -------------------------------------------------------------------------- #

setClass('EMci_biomod2_model',
         representation(alpha = 'numeric', side = 'character'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMci', alpha = 0.05, side = 'superior'),
         validity = function(object) {
           if (!(object@side %in% c('inferior','superior'))) {
             stop("side arg should be 'inferior' or 'superior")
           } else { return(TRUE) }
         })

##' 
##' @rdname predict.em
##' @export
##' 


setMethod('predict2', signature(object = 'EMci_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, ...){
              args <- list(...)
              mean_prediction <- args$mean_prediction
              sd_prediction <- args$sd_prediction
              side <- args$side
              
              if (is.null(mean_prediction)) { 
                mean_prediction <- calc(newdata, mean) 
              }
              if (is.null(sd_prediction)) { 
                sd_prediction <- calc(newdata, sd) 
              }
              
              ci_prediction <-  switch(
                side,
                inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) )
              )
              
              if (on_0_1000) {
                ci_prediction <- reclassify(round(ci_prediction), c(-Inf, 0, 0, 1000, Inf, 1000))
              } else {
                ci_prediction <- reclassify(ci_prediction, c(-Inf, 0, 0, 1, Inf, 1))
              }
              ci_prediction
            } 
            
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, side = object@side, ...)
          }
)

setMethod('predict2', signature(object = 'EMci_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, ...){
              args <- list(...)
              mean_prediction <- args$mean_prediction
              sd_prediction <- args$sd_prediction
              side <- args$side
              on_0_1000 <- args$on_0_1000
              
              if (is.null(mean_prediction)) { 
                mean_prediction <- round(rowMeans(newdata, na.rm = TRUE)) 
              }
              if (is.null(sd_prediction)) { 
                sd_prediction <- apply(newdata, 1, sd, na.rm = TRUE)
              }
              
              # browser()
              ci_prediction <- switch(side,
                                      inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                                      superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))
              
              # reclassify prediction to prevent from out of bounds prediction
              if (on_0_1000) {
                ci_prediction <- round(ci_prediction * 1000)
                ci_prediction[ci_prediction > 1000] <- 1000
                ci_prediction[ci_prediction < 0] <- 0
              } else {
                ci_prediction[ci_prediction > 1] <- 1
                ci_prediction[ci_prediction < 0] <- 0
              }
              ci_prediction
            }
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, side = object@side, ...)
          }
)


### -------------------------------------------------------------------------- #
### 10.5 EMca_biomod2_model --------------------------------------------------
### -------------------------------------------------------------------------- #

setClass('EMca_biomod2_model',
         representation(thresholds = 'numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMca'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict2', signature(object = 'EMca_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            if (is.null(on_0_1000)) { 
              on_0_1000 <- FALSE
            }
            
            predfun <- function(newdata, on_0_1000, thresh, ...){
              
              out <- calc(bm_BinaryTransformation(newdata, thresh), function(x)
              {
                m <- mean(x)
                if (on_0_1000) { 
                  m <- round(m * 1000)
                }
                return(m)
              })
            }
            if (on_0_1000) {
              thresh <- object@thresholds 
            } else {
              thresh <- object@thresholds / 1000 
            }
            
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, thresh = thresh, ...)
          }
)

setMethod('predict2', signature(object = 'EMca_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            if (is.null(on_0_1000)) { 
              on_0_1000 <- FALSE 
            }
            predfun <- function(newdata, ...){
              out <- rowMeans(as.data.frame(bm_BinaryTransformation(newdata, thresh)), na.rm = TRUE)
              if (on_0_1000) {
                out <- round(out * 1000)
              }
              out
            }
            if (on_0_1000) { 
              thresh <- object@thresholds
            } else {
              thresh <- object@thresholds / 1000
            }
            
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


### -------------------------------------------------------------------------- #
### 10.6 EMwmean_biomod2_model -----------------------------------------------
### -------------------------------------------------------------------------- #

setClass('EMwmean_biomod2_model',
         representation(penalization_scores='numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMwmean'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict2', signature(object = 'EMwmean_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, penalization_scores, ...){
              
              out <- calc(newdata, function(x)
              {
                wm <- sum(x * penalization_scores)
                if (on_0_1000) { 
                  wm <- round(wm) 
                }
                return(wm)
              })
            }
            
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun,
                           penalization_scores = object@penalization_scores, ...)
          }
)

setMethod('predict2', signature(object = 'EMwmean_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, penalization_scores, ...){
              out <- as.vector(
                as.matrix(newdata) %*% penalization_scores
              )
              if (on_0_1000) { 
                out <- round(out) 
              }
              out
            }
            
            # redirect to predict2.biomod2_ensemble_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, 
                           penalization_scores = object@penalization_scores,  ...)
          }
)


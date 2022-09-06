
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

###################################################################################################
## 9. biomod2_ensemble_model
###################################################################################################

##' @name biomod2_ensemble_model
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
##' @export
##' 

# 9.1 Class Definition ----------------------------------------------------------------------------
setClass('biomod2_ensemble_model',
         representation(modeling.id = 'character'), ##maybe some additional args should be added here
         contains = 'biomod2_model',
         prototype = list(model_class = 'EM'),
         validity = function(object) { return(TRUE) })


###################################################################################################
## 10.1 EMmean_biomod2_model
###################################################################################################

setClass('EMmean_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmean'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict', signature(object = 'EMmean_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMmean", object, newdata, formal_predictions, ...))
          })

.predict.EMmean_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  out <- calc(formal_predictions, function(x)
  {
    m <- mean(x)
    if (on_0_1000) { m <- round(m) }
    return(m)
  })
  writeRaster(out, filename = filename, overwrite = TRUE)
  
  return(out)
}

.predict.EMmean_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ... )
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  seedval <- args$seedval
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                   function(mod.name, dir_name, resp_name, modeling.id){
                                     ## check if model is loaded on memory
                                     if (is.character(mod.name)) {
                                       mod <- get(load(file.path(dir_name, resp_name, "models", modeling.id, mod.name)))
                                     }
                                     temp_workdir = NULL
                                     if (length(grep("MAXENT.Phillips$", mod.name)) == 1) {
                                       temp_workdir = mod@model_output_dir
                                     }
                                     return(predict(mod, newdata = newdata, on_0_1000 = on_0_1000
                                                    , temp_workdir = temp_workdir, seedval = seedval))
                                   } , dir_name = object@dir_name, resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  out <- rowMeans(formal_predictions, na.rm = TRUE)

  if (on_0_1000){
    out <- round(out)
  }

  return(out)
}


###################################################################################################
## 10.2 EMmedian_biomod2_model
###################################################################################################

setClass('EMmedian_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict', signature(object = 'EMmedian_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMmedian", object, newdata, formal_predictions, ...))
          })

.predict.EMmedian_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  out <- calc(formal_predictions, function(x)
  {
    m <- median(x)
    if (on_0_1000) { m <- round(m) }
    return(m)
  })
  writeRaster(out, filename = filename, overwrite = TRUE)
  
  return(out)
}

.predict.EMmedian_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  out <- apply(formal_predictions, 1, median, na.rm = TRUE)
  if (on_0_1000) { out <- round(out) }
  return(out)
}


###################################################################################################
## 10.3 EMcv_biomod2_model
###################################################################################################

setClass('EMcv_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMcv'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict', signature(object = 'EMcv_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMcv", object, newdata, formal_predictions, ...))
          })

.predict.EMcv_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  filename <- args$filename
  if (is.null(filename)) { filename <- "" }
  
  if (nlayers(formal_predictions) > 1) {
    out <- calc(formal_predictions, cv, na.rm = TRUE, aszero = TRUE)
    writeRaster(out, filename = filename, overwrite = TRUE)
    return(out)
  } else {
    warning(paste0("\n Model EMcv was not computed because only one single model was kept in ensemble modeling ("
                   , names(formal_predictions), ")"))
    return(NULL)
  }
}

.predict.EMcv_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  if (ncol(formal_predictions) > 1) {
    out <- apply(formal_predictions, 1, cv, na.rm = TRUE, aszero = TRUE)
    return(out)
  } else {
    warning(paste0("\n Model EMcv was not computed because only one single model was kept in ensemble modeling ("
                   , colnames(formal_predictions), ")"))
    return(NULL)
  }
}


###################################################################################################
## 10.4 EMci_biomod2_model
###################################################################################################

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

setMethod('predict', signature(object = 'EMci_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMci", object, newdata, formal_predictions, ...))
          })

.predict.EMci_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving
  if (is.null(mean_prediction)) { mean_prediction <- calc(formal_predictions, mean) }
  if (is.null(sd_prediction)) { sd_prediction <- calc(formal_predictions, sd) }
  
  ci_prediction <- switch(object@side,
                          inferior = mean_prediction -  (sd_prediction * (qt((1 - object@alpha / 2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                          superior = mean_prediction +  (sd_prediction * (qt((1 - object@alpha / 2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))
  
  # reclassify prediction to prevent from out of bounds prediction
  if (on_0_1000) {
    ci_prediction <- reclassify(round(ci_prediction), c(-Inf, 0, 0, 1000, Inf, 1000))
  } else {
    ci_prediction <- reclassify(ci_prediction, c(-Inf, 0, 0, 1, Inf, 1))
  }
  
  return(ci_prediction)
}

.predict.EMci_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving
  if (is.null(mean_prediction)) { mean_prediction <- round(rowMeans(formal_predictions, na.rm = TRUE)) }
  if (is.null(sd_prediction)) { sd_prediction <- apply(formal_predictions, 1, sd, na.rm = TRUE) }
  
  ci_prediction <- switch(object@side,
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
  
  return(ci_prediction)
}


###################################################################################################
## 10.5 EMca_biomod2_model
###################################################################################################

setClass('EMca_biomod2_model',
         representation(thresholds = 'numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMca'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict', signature(object = 'EMca_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMca", object, newdata, formal_predictions, ...))
          })

.predict.EMca_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  if (on_0_1000) { thresh <- object@thresholds } else { thresh <- object@thresholds / 1000 }
  
  out <- calc(bm_BinaryTransformation(formal_predictions, thresh), function(x)
  {
    m <- mean(x)
    if (on_0_1000) { m <- round(m * 1000) }
    return(m)
  })
  writeRaster(out, filename = filename, overwrite = TRUE)
  
  return(out)
}

.predict.EMca_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  if (on_0_1000) { thresh <- object@thresholds } else { thresh <- object@thresholds / 1000 }
  
  out <- rowMeans(as.data.frame(bm_BinaryTransformation(formal_predictions, thresh)), na.rm = TRUE)
  if (on_0_1000) { out <- round(out * 1000) }
  return(out)
}


###################################################################################################
## 10.6 EMwmean_biomod2_model
###################################################################################################

setClass('EMwmean_biomod2_model',
         representation(penalization_scores='numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMwmean'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.em
##' @export
##' 

setMethod('predict', signature(object = 'EMwmean_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMwmean", object, newdata, formal_predictions, ...))
          })

.predict.EMwmean_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  out <- calc(formal_predictions, function(x)
  {
    wm <- sum(x * object@penalization_scores)
    if (on_0_1000) { wm <- round(wm) }
    return(wm)
  })
  writeRaster(out, filename = filename, overwrite = TRUE)

  return(out)
}

.predict.EMwmean_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, on_0_1000 = on_0_1000, ...)
  }
  
  out <- as.vector(as.matrix(formal_predictions) %*% object@penalization_scores)
  if (on_0_1000) { out <- round(out) }
  return(out)
}

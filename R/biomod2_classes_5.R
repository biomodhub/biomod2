# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD models definition(to make it easier to access plot, predict, ...)
# Damien Georges, Maya Gueguen
# 20/11/2012, update 22/10/2021
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



###################################################################################################
## 9. biomod2_ensemble_model
###################################################################################################

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

setMethod('predict', signature(object = 'EMmean_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMmean", object, newdata, formal_predictions, ...))
          })

.predict.EMmean_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  out <- calc(formal_predictions, function(x)
  {
    m <- mean(x)
    if (on_0_1000) { m <- round(m) }
    return(m)
  }, filename = filename, overwrite = TRUE)
  
  return(out)
}

.predict.EMmean_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                   function(mod, resp_name, modeling.id){
                                     ## check if model is loaded on memory
                                     if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                     return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                   } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  out <- rowMeans(formal_predictions, na.rm=T)

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

setMethod('predict', signature(object = 'EMmedian_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMmedian", object, newdata, formal_predictions, ...))
          })

.predict.EMmedian_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  out <- calc(formal_predictions, function(x)
  {
    m <- median(x)
    if (on_0_1000) { m <- round(m) }
    return(m)
  }, filename = filename, overwrite = TRUE)
  
  return(out)
}

.predict.EMmedian_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  out <- apply(formal_predictions, 1, median, na.rm = T)
  if (on_0_1000) { out <- round(out) }
  return(out)
}


###################################################################################################
## 10.3 EMcv_biomod2_model
###################################################################################################

setClass('EMcv_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object) { return(TRUE) })

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
    out <- calc(formal_predictions, cv, filename = filename, overwrite = TRUE, na.rm = TRUE, aszero = TRUE)
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
    out <- apply(formal_predictions, 1, cv, na.rm = T, aszero = T)
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

setMethod('predict', signature(object = 'EMci_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMci", object, newdata, formal_predictions, ...))
          })

.predict.EMci_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
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
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving
  if (is.null(mean_prediction)) { mean_prediction <- round(rowMeans(formal_predictions, na.rm = T)) }
  if (is.null(sd_prediction)) { sd_prediction <- apply(formal_predictions, 1, sd, na.rm =  T) }
  
  ci_prediction <- switch(object@side,
                          inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                          superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))
  
  # reclassify prediction to prevent from out of bounds prediction
  if (on_0_1000) {
    ci_prediction <- round(ci_prediction)
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
         representation(tresholds = 'numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMca'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'EMca_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMca", object, newdata, formal_predictions, ...))
          })

.predict.EMca_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (on_0_1000) { thresh <- object@tresholds } else { thresh <- object@tresholds / 1000 }
  
  out <- calc(bm_BinaryTransformation(formal_predictions, thresh), function(x)
  {
    m <- mean(x)
    if (on_0_1000) { m <- round(m * 1000) }
    return(m)
  }, filename = filename, overwrite = TRUE)
  
  return(out)
}


.predict.EMca_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (on_0_1000) { thresh <- object@tresholds } else { thresh <- object@tresholds / 1000 }
  
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

setMethod('predict', signature(object = 'EMwmean_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            return(.template_predictEM(mod = "EMwmean", object, newdata, formal_predictions, ...))
          })

.predict.EMwmean_biomod2_model.RasterStack <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  out <- calc(formal_predictions, function(x)
  {
    wm <- sum(x * object@penalization_scores)
    if (on_0_1000) { wm <- round(wm) }
    return(wm)
  }, filename = filename, overwrite = TRUE)
  
  return(out)
}

.predict.EMwmean_biomod2_model.data.frame <- function(object, newdata = NULL, formal_predictions = NULL, ...)
{
  if (is.null(formal_predictions)) {
    formal_predictions = .template_predictEM.formal_predictions(object, newdata, formal_predictions, ...)
  }
  
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  out <- as.vector(as.matrix(formal_predictions) %*% object@penalization_scores)
  if (on_0_1000) { out <- round(out) }
  return(out)
}

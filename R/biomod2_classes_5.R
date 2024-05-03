
##' @name predict.em
##' @author Damien Georges
##' 
##' @title Functions to get predictions from \code{\link{biomod2_ensemble_model}} objects
##' 
##' @description This function allows the user to predict single models from 
##' \code{\link{biomod2_ensemble_model}} on (new) explanatory variables.
##' 
##' 
##' @param object a \code{\link{biomod2_ensemble_model}} object
##' @param newdata a \code{data.frame} or \code{\link[terra:rast]{SpatRaster}} object 
##' containing data for new predictions
##' @param \ldots (\emph{optional}) 
##' 
##' 
##' @seealso \code{\link{biomod2_ensemble_model}}
##' @family Toolbox functions
##' 
##' 
##' @importFrom terra app classify nlyr
##' @importFrom stats qt sd
##' 
NULL

#setGeneric("predict", def = function(object, ...) { standardGeneric("predict") })

# ---------------------------------------------------------------------------- #
# 9. biomod2_ensemble_model --------------------------------------------------
# ---------------------------------------------------------------------------- #

##' @name biomod2_ensemble_model
##' @aliases biomod2_ensemble_model-class
##' @aliases EMmean_biomod2_model-class
##' @aliases EMmedian_biomod2_model-class
##' @aliases EMcv_biomod2_model-class
##' @aliases EMci_biomod2_model-class
##' @aliases EMca_biomod2_model-class
##' @aliases EMwmean_biomod2_model-class
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
##' @slot model_evaluation a \code{data.frame} containing the model evaluations
##' @slot model_variables_importance a \code{data.frame} containing the model variables importance
##' 
##' @param object a \code{\link{biomod2_ensemble_model}} object
##' 
##' @details 
##' 
##' \code{biomod2_model} is the basic object for \pkg{biomod2} ensemble species distribution models. 
##' \cr All listed classes below are derived from \code{biomod2_model}, and have a 
##' \code{model_class} slot specific value :
##' 
##' \itemize{
##'   \item \code{biomod2_ensemble_model} : \code{model_class} is \code{EM}
##'   \item \code{EMmean_biomod2_model} : \code{model_class} is \code{EMmean}
##'   \item \code{EMmedian_biomod2_model} : \code{model_class} is \code{EMmedian}
##'   \item \code{EMcv_biomod2_model} : \code{model_class} is \code{EMcv}
##'   \item \code{EMci_biomod2_model} : \code{model_class} is \code{EMci}
##'   \item \code{EMca_biomod2_model} : \code{model_class} is \code{EMca}
##'   \item \code{EMwmean_biomod2_model} : \code{model_class} is \code{EMwmean}
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
         representation(modeling.id = 'character'), ## maybe some additional args should be added here
         contains = 'biomod2_model',
         prototype = list(model_class = 'EM'),
         validity = function(object) { return(TRUE) })


## 9.2 Show method -------------------------------------------------------------
##' @rdname biomod2_ensemble_model
##' @importMethodsFrom methods show
##' @importFrom methods callNextMethod
##' @export
##' 

setMethod('show', signature('biomod2_ensemble_model'),
          function(object) {
            callNextMethod(object)
          })


## 9.3 predict2 method -------------------------------------------------------------

### biomod2_ensemble_model predict2.em doc -------------------------------------

##' @name predict2.em
##' @aliases predict2.biomod2_ensemble_model.SpatRaster
##' @aliases predict2.biomod2_ensemble_model.data.frame
##' @aliases predict2.EMmean_biomod2_model.SpatRaster
##' @aliases predict2.EMmean_biomod2_model.data.frame
##' @aliases predict2.EMmedian_biomod2_model.SpatRaster
##' @aliases predict2.EMmedian_biomod2_model.data.frame
##' @aliases predict2.EMcv_biomod2_model.SpatRaster
##' @aliases predict2.EMcv_biomod2_model.data.frame
##' @aliases predict2.EMci_biomod2_model.SpatRaster
##' @aliases predict2.EMci_biomod2_model.data.frame
##' @aliases predict2.EMca_biomod2_model.SpatRaster
##' @aliases predict2.EMca_biomod2_model.data.frame
##' @aliases predict2.EMwmean_biomod2_model.SpatRaster
##' @aliases predict2.EMwmean_biomod2_model.data.frame
##' @author Remi Patin
##' 
##' @title Functions to get predictions from \code{\link{biomod2_ensemble_model}} objects
##' 
##' @description This function allows the user to predict single models from 
##' \code{\link{biomod2_ensemble_model}} on (new) explanatory variables.
##' 
##' 
##' @param object a \code{\link{biomod2_ensemble_model}} object
##' @param newdata a \code{data.frame} or \code{\link[terra:rast]{SpatRaster}} object 
##' containing data for new predictions
##' @param data_as_formal_predictions (\emph{optional, default} \code{FALSE}). A
##' \code{boolean} describing whether \code{newdata} is given as raw environmental 
##' data (\code{FALSE}) or as formal predictions of the individual models 
##' used to build the ensemble model (\code{TRUE}).
##' 
##' @param \ldots (\emph{optional}) 
##' @inheritParams predict2.bm
##' 
##' @seealso \code{\link{biomod2_ensemble_model}}
##' @family Toolbox functions
##' 
##' 
##' @importFrom terra rast app classify writeRaster
##' @keywords internal

NULL

### biomod2_ensemble_model + SpatRaster  -------------------------------------------------
##' @rdname predict2.em
setMethod('predict2', signature(object = 'biomod2_ensemble_model', newdata = "SpatRaster"),
          function(object, newdata, predfun, seedval = NULL, ...) {
            args <- list(...)
            data_as_formal_predictions <- args$data_as_formal_predictions
            if (is.null(data_as_formal_predictions)) {
              data_as_formal_predictions <- FALSE 
            }
            filename <- args$filename
            overwrite <- args$overwrite
            on_0_1000 <- args$on_0_1000
            na.rm <- args$na.rm
            
            if (is.null(na.rm)) { 
              na.rm <- TRUE 
            }
            
            if (is.null(overwrite)) { 
              overwrite <- TRUE 
            }
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
            mod.name <- args$mod.name
            
            
            if (data_as_formal_predictions) { 
              newdata <- subset(newdata, object@model)
            } else {
              newdata <- .get_formal_predictions(object, newdata, on_0_1000 = on_0_1000, seedval = seedval)
            }
            
            out <- predfun(newdata, 
                           on_0_1000 = on_0_1000,
                           mean_prediction = mean_prediction,
                           sd_prediction = sd_prediction, 
                           side = side,
                           thresh = thresh, 
                           penalization_scores = penalization_scores,
                           mod.name = mod.name,
                           na.rm = na.rm)
            
            if (!is.null(out) & !is.null(filename)) {
              cat("\n\t\tWriting projection on hard drive...")
              if (on_0_1000 & !inherits(object, "EMcv_biomod2_model")) { ## projections are stored as positive integer
                writeRaster(out, filename = filename, overwrite = overwrite, 
                            datatype = "INT2S", NAflag = -9999)
              } else { ## keep default data format for saved raster
                writeRaster(out, filename = filename, overwrite = overwrite)
              }
              out <- rast(filename)
            }
            return(out)
            
          })

### biomod2_ensemble_model + data.frame  -------------------------------------
##' @rdname predict2.em
setMethod('predict2', signature(object = 'biomod2_ensemble_model', newdata = "data.frame"),
          function(object, newdata, predfun, seedval = NULL,  ...) {
            args <- list(...)
            data_as_formal_predictions <- args$data_as_formal_predictions
            if (is.null(data_as_formal_predictions)) {
              data_as_formal_predictions <- FALSE 
            }
            
            on_0_1000 <- args$on_0_1000
            if (is.null(on_0_1000)) {
              on_0_1000 <- FALSE 
            }
            
            na.rm <- args$na.rm
            if (is.null(na.rm)) { 
              na.rm <- TRUE 
            }
            
            # additional arg retrieved for EMci
            sd_prediction <- args$sd_prediction
            mean_prediction <- args$mean_prediction
            side <- args$side
            # additional arg retrived for EMca
            thresh <- args$thresh
            # additional arg retrived for EMwmean
            penalization_scores <- args$penalization_scores
            
            if (data_as_formal_predictions) { 
              newdata <- newdata[ , object@model, drop = FALSE]
            } else  {
              newdata <- .get_formal_predictions(object, newdata, on_0_1000 = on_0_1000, seedval = seedval)
            }
            out <- predfun(newdata,
                           on_0_1000 = on_0_1000,
                           mean_prediction = mean_prediction,
                           sd_prediction = sd_prediction,
                           side = side,
                           thresh = thresh, 
                           penalization_scores = penalization_scores,
                           na.rm = na.rm)
            return(out)
            
          })

# --------------------------------------------------------------------------- #
# 10.1 biomod2_ensemble_model subclass ---------------------------------------
# ---------------------------------------------------------------------------- #

### -------------------------------------------------------------------------- #
### 10.1 EMmean_biomod2_model ------------------------------------------------
### -------------------------------------------------------------------------- #
##' @name EMmean_biomod2_model-class
##' @rdname biomod2_ensemble_model
##' @export


setClass('EMmean_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmean'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict2.em
##' 

setMethod('predict2', signature(object = 'EMmean_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, mod.name, na.rm, ...){
              if (nlyr(newdata) == 1) {
                return(newdata)
              } else {
                return(
                  app(newdata,function(x){
                    m <- mean(x, na.rm = na.rm)
                    if (on_0_1000) { 
                      m <- round(m)
                    }
                    return(m)
                  }, wopt = list(names = mod.name))
                )
              }
            }
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.em
setMethod('predict2', signature(object = 'EMmean_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, na.rm, ...){
              out <- rowMeans(newdata, na.rm = na.rm)
              if (on_0_1000) { 
                out <- round(out) 
              }
              out
            }
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


### -------------------------------------------------------------------------- #
### 10.2 EMmedian_biomod2_model ----------------------------------------------
### -------------------------------------------------------------------------- #
##' @name EMmedian_biomod2_model-class
##' @rdname biomod2_ensemble_model
##' @export

setClass('EMmedian_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict2.em
##' 


setMethod('predict2', signature(object = 'EMmedian_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, mod.name, na.rm, ...){
              if (nlyr(newdata) == 1) {
                return(newdata)
              } else {
                return(
                  app(newdata,function(x){
                    m <- median(x, na.rm = na.rm)
                    if (on_0_1000) { 
                      m <- round(m)
                    }
                    return(m)
                  }, wopt = list(names = mod.name))
                )
              }
            }
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.em
setMethod('predict2', signature(object = 'EMmedian_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, na.rm, ...){
              out <- apply(newdata, 1, median, na.rm = na.rm)
              if (on_0_1000) { 
                out <- round(out) 
              }
              out
            }
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

### -------------------------------------------------------------------------- #
### 10.3 EMcv_biomod2_model --------------------------------------------------
### -------------------------------------------------------------------------- #
##' @name EMcv_biomod2_model-class
##' @rdname biomod2_ensemble_model
##' @export

setClass('EMcv_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMcv'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict2.em
##' 

setMethod('predict2', signature(object = 'EMcv_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata,  ...) {
            predfun <- function(newdata, on_0_1000, mod.name, na.rm, ...){
              if (nlyr(newdata) <= 1) {
                stop(paste0("\n Model EMcv was not computed because only one single model was kept in ensemble modeling (", names(newdata), ")"))
              }
              out <- app(newdata, function(x){
                sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm) * 100
              }, wopt = list(names = mod.name))
              return(out)
            }
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.em
setMethod('predict2', signature(object = 'EMcv_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, na.rm, ...){
              if (ncol(newdata) <= 1) {
                stop(paste0("\n Model EMcv was not computed because only one single model was kept in ensemble modeling ("
                            , colnames(newdata), ")"))
              }
              out <- apply(newdata, 1,
                           function(x) {
                             ifelse(length(x) == 1, 0, 
                                    sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm)*100)
                           })
              return(out)
            } 
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

### -------------------------------------------------------------------------- #
### 10.4 EMci_biomod2_model --------------------------------------------------
### -------------------------------------------------------------------------- #
##' @name EMci_biomod2_model-class
##' @rdname biomod2_ensemble_model
##' @export

setClass('EMci_biomod2_model',
         representation(alpha = 'numeric', side = 'character'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMci', alpha = 0.05, side = 'superior'),
         validity = function(object) {
           if (!(object@side %in% c('inferior','superior'))) {
             stop("side arg should be 'inferior' or 'superior")
           } else { 
             return(TRUE) 
           }
         })

##' 
##' @rdname predict2.em
##' 


setMethod('predict2', signature(object = 'EMci_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(newdata, on_0_1000, mod.name, na.rm, ...){
              args <- list(...)
              mean_prediction <- args$mean_prediction
              sd_prediction <- args$sd_prediction
              side <- args$side

              if (is.null(mean_prediction)) { 
                mean_prediction <- app(newdata, mean, wopt = list(names = mod.name),
                                       na.rm = na.rm)
              }
              if (is.null(sd_prediction)) { 
                sd_prediction <- app(newdata, sd, wopt = list(names = mod.name),
                                     na.rm = na.rm) 
              }
              
              ci_prediction <-  switch(
                side,
                inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) )
              )
              
              if (on_0_1000) {
                ci_prediction <- classify(round(ci_prediction),
                                          matrix(c(-Inf, 0, 0,
                                                   1000, Inf, 1000),
                                                 nrow = 2, byrow = TRUE))
              } else {
                ci_prediction <- classify(round(ci_prediction),
                                          matrix(c(-Inf, 0, 0,
                                                   1, Inf, 1),
                                                 nrow = 2, byrow = TRUE))     
              }
              ci_prediction
            } 
            
            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, side = object@side, ...)
          }
)

##' @rdname predict2.em
setMethod('predict2', signature(object = 'EMci_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(newdata, na.rm, ...){
              args <- list(...)
              mean_prediction <- args$mean_prediction
              sd_prediction <- args$sd_prediction
              side <- args$side
              on_0_1000 <- args$on_0_1000
              
              if (is.null(mean_prediction)) { 
                mean_prediction <- rowMeans(newdata, na.rm = na.rm)
              }
              if (is.null(sd_prediction)) { 
                sd_prediction <- apply(newdata, 1, sd, na.rm = na.rm)
              }
              
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
            # redirect to predict2.biomod2_ensemble_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, side = object@side, ...)
          }
)


### -------------------------------------------------------------------------- #
### 10.5 EMca_biomod2_model --------------------------------------------------
### -------------------------------------------------------------------------- #
##' @name EMca_biomod2_model-class
##' @rdname biomod2_ensemble_model
##' @export

setClass('EMca_biomod2_model',
         representation(thresholds = 'numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMca'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict2.em
##' 

setMethod('predict2', signature(object = 'EMca_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, data_as_formal_predictions = FALSE, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            if (is.null(on_0_1000)) { 
              on_0_1000 <- FALSE
            }
            
            predfun <- function(newdata, on_0_1000, thresh, mod.name, na.rm, ...){
              if (nlyr(newdata) == 1) {
                return(bm_BinaryTransformation(newdata, thresh))
              } else {
                return(
                  app(bm_BinaryTransformation(newdata, thresh),
                      function(x){
                        m <- mean(x, na.rm = na.rm)
                        if (on_0_1000) { 
                          m <- round(m * 1000)
                        }
                        return(m)
                      }, wopt = list(names = mod.name))
                )
              }
            }
            if (on_0_1000) {
              thresh <- object@thresholds 
            } else {
              thresh <- object@thresholds / 1000 
            }
            

            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, thresh = thresh,
                           data_as_formal_predictions = data_as_formal_predictions,
                           ...)
          }
)

##' @rdname predict2.em
setMethod('predict2', signature(object = 'EMca_biomod2_model', newdata = "data.frame"),
          function(object, newdata, data_as_formal_predictions = FALSE, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            if (is.null(on_0_1000)) { 
              on_0_1000 <- FALSE 
            }
            predfun <- function(newdata, na.rm, ...){
              out <- rowMeans(bm_BinaryTransformation(newdata, thresh), 
                              na.rm = na.rm)
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
            

            # redirect to predict2.biomod2_ensemble_model.data.frame
            callNextMethod(object, newdata, predfun = predfun,
                           data_as_formal_predictions = data_as_formal_predictions,
                           ...)
          }
)


### -------------------------------------------------------------------------- #
### 10.6 EMwmean_biomod2_model -----------------------------------------------
### -------------------------------------------------------------------------- #
##' @name EMwmean_biomod2_model-class
##' @rdname biomod2_ensemble_model
##' @export

setClass('EMwmean_biomod2_model',
         representation(penalization_scores='numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMwmean'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict2.em
##' 

setMethod('predict2', signature(object = 'EMwmean_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, data_as_formal_predictions = FALSE, ...) {
            if (ncol(newdata) < 1) {
              stop("Model EMwmean was not computed because no single model was kept in ensemble modeling")
            }
            predfun <- function(newdata, on_0_1000, penalization_scores,
                                mod.name, ...){
              if (nlyr(newdata) == 1) {
                return(newdata)
              } else {
                return(
                  app(newdata, function(x) {
                    wm <- sum(x * penalization_scores)
                    if (on_0_1000) { 
                      wm <- round(wm) 
                    }
                    return(wm)
                  }, wopt = list(names = mod.name))
                )
              }
            }

            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata,
                           predfun = predfun,
                           data_as_formal_predictions = data_as_formal_predictions,
                           penalization_scores = object@penalization_scores, ...)
          }
)

##' @rdname predict2.em
setMethod('predict2', signature(object = 'EMwmean_biomod2_model', newdata = "data.frame"),
          function(object, newdata, data_as_formal_predictions = FALSE, ...) {
            
            if (ncol(newdata) < 1) {
              stop("Model EMwmean was not computed because no single model was kept in ensemble modeling")
            }
            predfun <- function(newdata, on_0_1000, penalization_scores, ...){
              out <- as.vector(
                as.matrix(newdata) %*% penalization_scores
              )
              if (on_0_1000) { 
                out <- round(out) 
              }
              out
            }
            

            # redirect to predict2.biomod2_ensemble_model.SpatRaster
            callNextMethod(object, newdata,
                           predfun = predfun, 
                           penalization_scores = object@penalization_scores,
                           data_as_formal_predictions = data_as_formal_predictions,
                           ...)
          }
)


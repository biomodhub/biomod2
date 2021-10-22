
# EM parent Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('biomod2_ensemble_model',
         representation(modeling.id = 'character'), ##maybe some additional args should be added here
         contains = 'biomod2_model',
         prototype = list(model_class = 'EM'),
         validity = function(object) { return(TRUE) })

# EMmean Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMmean_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmean'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'EMmean_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...){
            args <- list(...)
            do_check <- args$do_check
            if (is.null(do_check)) { do_check <- TRUE }
            if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMmean_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMmean_biomod2_model.data.frame(object, newdata, formal_predictions,  ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMmean_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
          function(mod, resp_name, modeling.id){
            ## check if model is loaded on memory
            if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
            return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
          }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  out <- calc(formal_predictions,
             function(x){
               m <- mean(x)
               if(on_0_1000) m <- round(m)
               return(m)
             },
             filename = filename, overwrite = TRUE)

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



# EMmedian Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMmedian_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'EMmedian_biomod2_model'),
          function(object, newdata = NULL, formal_predictions = NULL, ...)
          {
            args <- list(...)
            do_check <- args$do_check
            if (is.null(do_check)) { do_check <- TRUE }
            if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMmedian_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMmedian_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMmedian_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  out <- calc(formal_predictions,
              function(x){
                m <- median(x)
                if(on_0_1000) m <- round(m)
                return(m)
              },
              filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMmedian_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
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

  out <- apply(formal_predictions, 1, median, na.rm=T)

  if (on_0_1000){
    out <- round(out)
  }

  return(out)

}



# EMcv Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMcv_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmedian'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'EMcv_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)
            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMcv_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMcv_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMcv_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions
  mean_prediction <- args$mean_prediction

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  out <- calc(formal_predictions, cv,
              filename = filename, overwrite = TRUE, na.rm=TRUE, aszero=TRUE)

  return(out)

}

.predict.EMcv_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions
#   mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

#   if(is.null(mean_prediction)){
#     # calculate mean of predictions
#     mean_prediction <- round(rowMeans(formal_predictions, na.rm=T))
#   }
#   # transforming 0 into Inf to produce null cv where mean is null
#   mean_prediction[mean_prediction==0] <- Inf
#
#   # calculate cv of formal models predictions
#   sd_prediction <- apply(formal_predictions,1,sd, na.rm=T)
#
#   return(round(sd_prediction / mean_prediction, 2))

  out <- apply(formal_predictions, 1, cv, na.rm=T, aszero=T)

#   if (on_0_1000){
#     out <- round(out)
#   }

  return(out)

}



# EMci Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
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
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMci_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMci_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMci_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  if(is.null(mean_prediction)){
    mean_prediction <- calc(formal_predictions, mean)
  }

  if(is.null(sd_prediction)){
    sd_prediction <- calc(formal_predictions, sd)
  }

  ci_prediction <- switch(object@side,
                          inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                          superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))

  # reclassify prediction to prevent from out of bounds prediction
  if (on_0_1000){
    ci_prediction <- reclassify(round(ci_prediction), c(-Inf,0,0, 1000,Inf,1000))
  } else {
    ci_prediction <- reclassify(ci_prediction, c(-Inf,0,0, 1,Inf,1))
  }

  return(ci_prediction)
}

.predict.EMci_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  #formal_predictions <- args$formal_predictions
  mean_prediction <- args$mean_prediction # mean of predictions should be given for time saving
  sd_prediction <- args$sd_prediction # mean of predictions should be given for time saving

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- sapply(object@model,
                                 function(mod, resp_name, modeling.id){
                                   ## check if model is loaded on memory
                                   if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                   return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                 } , resp_name = object@resp_name, modeling.id = object@modeling.id)
  }

  if(is.null(mean_prediction)){
    # calculate mean of predictions
    mean_prediction <- round(rowMeans(formal_predictions, na.rm=T))
  }

  if(is.null(sd_prediction)){
    # calculate cv of formal models predictions
    sd_prediction <- apply(formal_predictions,1,sd, na.rm=T)
  }

  ci_prediction <- switch(object@side,
                          inferior = mean_prediction -  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ),
                          superior = mean_prediction +  (sd_prediction * (qt((1-object@alpha/2), df = length(object@model) + 1 ) / sqrt(length(object@model))) ))

  # reclassify prediction to prevent from out of bounds prediction
  if (on_0_1000){
  ci_prediction <- round(ci_prediction)
  ci_prediction[ci_prediction > 1000] <- 1000
  ci_prediction[ci_prediction < 0] <- 0
  } else {
    ci_prediction[ci_prediction > 1] <- 1
    ci_prediction[ci_prediction < 0] <- 0
  }

  return(ci_prediction)
}


# EMca Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMca_biomod2_model',
         representation(tresholds = 'numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMca'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'EMca_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMca_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMca_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMca_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  filename <- args$filename
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

  if (on_0_1000){
    thresh <- object@tresholds
  } else { thresh <- object@tresholds / 1000 }

#   out <- raster::mean(BinaryTransformation(formal_predictions, thresh), na.rm=T)
  out <- calc(BinaryTransformation(formal_predictions, thresh),
                function(x){
                  m <- mean(x)
                  if(on_0_1000) m <- round(m * 1000)
                  return(m)
                  },
              filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMca_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
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

  if (on_0_1000){
    thresh <- object@tresholds
  } else { thresh <- object@tresholds / 1000 }

  out <- apply(as.data.frame(BinaryTransformation(formal_predictions, thresh)), 1, mean, na.rm=T)

  if (on_0_1000){
    out <- round(out * 1000)
  }

  return(out)

}

# EMwmean Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMwmean_biomod2_model',
         representation(penalization_scores='numeric'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMwmean'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'EMwmean_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#             ## check if models are formal loaded
#             if(is.character(object@model)){
#               model_tmp <- lapply(object@model, function(x){
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               })
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

            if(inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')){
              return(.predict.EMwmean_biomod2_model.RasterStack(object, newdata, formal_predictions, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')){
              return(.predict.EMwmean_biomod2_model.data.frame(object, newdata, formal_predictions, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.EMwmean_biomod2_model.RasterStack <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
  args <- list(...)

  on_0_1000 <- args$on_0_1000
  filename <- args$filename
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(filename)) filename <- ''

  #formal_predictions <- args$formal_predictions

  if(is.null(formal_predictions)){
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id){
                                                 ## check if model is loaded on memory
                                                 if(is.character(mod)) mod <- get(load(file.path(resp_name, "models", modeling.id, mod)))
                                                 return(predict(mod, newdata=newdata, on_0_1000=on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }

#   out <- sum(formal_predictions * object@penalization_scores)
  out <- calc(formal_predictions,
              function(x){
                wm <- sum(x * object@penalization_scores)
                if(on_0_1000) wm <- round(wm)
                return(wm)
              },
              filename = filename, overwrite = TRUE)

  return(out)

}

.predict.EMwmean_biomod2_model.data.frame <- function(object, newdata=NULL, formal_predictions=NULL, ... ){
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

  out <- as.vector(as.matrix(formal_predictions) %*% object@penalization_scores)

  if (on_0_1000){
    out <- round(out)
  }

  return(out)

}


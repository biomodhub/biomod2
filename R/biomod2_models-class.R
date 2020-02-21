# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Definition of models build with biomod2 to make it easy to access
# plot, .predict...
# Damien G. - 20/11/2012
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# NOTE -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# This script will induce strong changes in biomod2. The version
# of biomod2 using objects defined there will be biomod2_2.x.x
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Formal Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('biomod2_model',
         representation(model_name = 'character',
                        model_class = 'character',
                        model_options = 'list',
                        model = 'ANY',
                        scaling_model = 'ANY',
                        resp_name = 'character',
                        expl_var_names = 'character',
                        expl_var_type = 'character',
                        expl_var_range = 'list',
                        model_evaluation = 'matrix',
                        model_variables_importance = 'matrix'),
         prototype = list(model_name = 'mySpecies_DataSet_RunName_myModelClass',
                   model_class = 'myModelClass',
                   model_options = list(),
                   model = list(),
                   scaling_model = list(),
                   resp_name = 'mySpecies',
                   expl_var_names = 'myRespVar',
                   expl_var_type = 'unknown',
                   expl_var_range = list(),
                   model_evaluation = matrix(),
                   model_variables_importance = matrix()),
         validity = function(object){

           # check that scaler is a glm if it is defined
           if(length(object@scaling_model))
             if(!inherits(object@scaling_model, c("glm", "lm"))) 
               return(FALSE)

           return(TRUE)
           } )

setMethod('show', signature('biomod2_model'),
          function(object){
            .bmCat("'biomod2_model'")
            cat("\n\t model name :", object@model_name, fill=.Options$width)
            cat("\n\t model class :", object@model_class, fill=.Options$width)
            cat("\n\t This model", ifelse(length(object@scaling_model), "has", "doesn't have"),"its own scaler", fill=.Options$width)

            cat("\n")
            cat("\n\t response modelled :", object@resp_name, fill=.Options$width)
            cat("\n\n\t explanatory variables used:", fill=.Options$width)
            cat("\n\t", "name", "\t", "type", "\t", "range", fill=.Options$width)
            for(i in 1: length(object@expl_var_names)){
              cat("\n\t", object@expl_var_names[i],"\t", object@expl_var_type[i], "\t", object@expl_var_range[[i]], fill=.Options$width)
            }

            cat("\n")
            cat("\n\t NOTE : ")
            cat("\n\t\t You can access 'formal' model with get_formal_model function")
            cat(ifelse(length(object@scaling_model), "\n\t\t You can access scaling model with get_scaling_model function\n", "\n"))

            .bmCat()
          })

# Models getters -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
setGeneric( "get_formal_model",
            def = function(object){
              standardGeneric( "get_formal_model" )
            } )

setMethod('get_formal_model', signature('biomod2_model'),
          function(object){
            return(object@model)
          })

setGeneric( "get_scaling_model",
            def = function(object){
              standardGeneric( "get_scaling_model" )
            } )

setMethod('get_scaling_model', signature('biomod2_model'),
          function(object){
            return(object@scaling_model)
          })

# Fuction to get variables ranges
get_var_type <- function(data){
    return(sapply(data,class))
}

get_var_range <- function(data){
  get_range <- function(x){
    if(is.numeric(x)){
      return(c(min=min(x,na.rm=T), max=max(x,na.rm=T)))
    }
    if(is.factor(x)){
      return(levels(x))
    }
  }
  xx <- lapply(data,get_range)
  names(xx) <- names(data)
  return(xx)
}

# Function to check new data range compatibility with calibrating data #
check_data_range <- function(model, new_data){
  ## TODO : remettre en marche cette fonction

#   # get calibration data caracteristics
#   expl_var_names <- model@expl_var_names
#   expl_var_type <- model@expl_var_type
#   expl_var_range <- model@expl_var_range
#
#   if(inherits(new_data, "Raster")){ ## raster data case =-=-=-=-=-=-=- #
#     # check var names compatibility
#     nd_expl_var_names <- names(new_data)
#     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
#       stop("calibration and projections variables names mismatch")
#     }
#     # reorder the stack
#     new_data <- raster::subset(new_data,expl_var_names)
#     # check var types compatibility (factors)
#     expl_var_fact <- (expl_var_type=='factor')
#     nd_expl_var_fact <- is.factor(new_data)
#     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
#       stop("calibration and projections variables class mismatch")
#     }
#     # check var range compatibility
#     ### remove all new factors
#     if(sum(expl_var_fact)>0){ ## there are factorial variables
#       for(fact_var_id in which(expl_var_fact)){
#         ## check if new factors occurs
#         nd_levels <- levels(raster::subset(new_data,fact_var_id))[[1]]
#         nd_levels <- as.character(nd_levels[,ncol(nd_levels)])
#         names(nd_levels) <- levels(raster::subset(new_data,fact_var_id))[[1]]$ID
#         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
#
#         ## detect new levels
#         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
#
#         if(length(new_levels)){
#           for(n_l in new_levels){
#             # remove points where out of range factors have been detected
#             new_data[subset(new_data,fact_var_id)[]==as.numeric(names(nd_levels)[which(nd_levels==n_l)])] <- NA
#           }
#           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
#         }
#       }
#     }
#     ## convert data to be sure to get RasterStack output
# #     new_data <- stack(new_data)
#   } else{ ## table data case -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#     # check var names compatibility
#     nd_expl_var_names <- colnames(new_data)
#     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
#       stop("calibration and projections variables names mismatch")
#     }
#     # reorder the stack
#     new_data <- new_data[,expl_var_names, drop=F]
#     # check var types compatibility (factors)
#     expl_var_fact <- (expl_var_type=='factor')
#     nd_expl_var_fact <- sapply(new_data,is.factor)
#
#     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
#       stop("calibration and projections variables class mismatch")
#     }
#     # check var range compatibility
#     ### remove all new factors
#     if(sum(expl_var_fact)>0){ ## there are factorial variables
#       for(fact_var_id in which(expl_var_fact)){
#         ## check if new factors occurs
#         nd_levels <- levels(new_data[,fact_var_id])
#         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
#
#         ## detect new levels
#         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
#
#         if(length(new_levels)){
#           # remove points where out of range factors have been detected
# #           new_data <- new_data[- which(new_data[,fact_var_id] %in% new_levels),]
#           new_data[which(new_data[,fact_var_id] %in% new_levels),] <- NA
#           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
#         }
#       }
#     }
#   }

  return(new_data)
}


# ANN Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setClass('ANN_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'ANN'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "nnet")) return(FALSE)
           return(TRUE)
           })

setMethod('predict', signature(object = 'ANN_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.ANN_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
               return(.predict.ANN_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.ANN_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

#   set.seed(555)
  proj <- predict(newdata, get_formal_model(object), type="raw")

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }

  return(proj)
}

.predict.ANN_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE



  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }


  set.seed(555)
  proj <- as.numeric( predict(get_formal_model(object), newdata[not_na_rows,,drop=F], type="raw") )

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# CTA Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('CTA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'CTA'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "rpart")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'CTA_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.CTA_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.CTA_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.CTA_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  set.seed(123) # to be able to refind our trees MAY BE BAD
  proj <- predict(newdata, model=get_formal_model(object), type='prob', index=2)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.CTA_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  set.seed(123)
  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]),type="prob")[,2])

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# FDA Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('FDA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'FDA'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "fda")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'FDA_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.FDA_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.FDA_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.FDA_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- predict(newdata, model=get_formal_model(object), type='posterior', index=2)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.FDA_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000

  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]),type = "posterior")[,2])

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# GAM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GAM_biomod2_model',
         representation(model_subclass='character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GAM',
                   model_subclass = 'GAM_mgcv'),
         validity = function(object){
           # check model class
           if(! (object@model_subclass %in% c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv') ) )
             return(FALSE)

           if(object@model_subclass %in% c('GAM_mgcv','GAM_gam'))
             if(!inherits(object@model, c("gam")))
               return(FALSE)

           if(object@model_subclass == 'BAM_mgcv')
             if(!inherits(object@model, c("bam")))
               return(FALSE)

           return(TRUE)
         })

setMethod('predict', signature(object = 'GAM_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE


            if(object@model_subclass %in% c("GAM_mgcv","BAM_mgcv")){
              # cat("\n*** unloading gam package / loading mgcv package")
              if(isNamespaceLoaded("gam")){unloadNamespace("gam")}
              if(!isNamespaceLoaded("mgcv")){requireNamespace("mgcv", quietly = TRUE)}
            }

            if(object@model_subclass == "GAM_gam"){
              # cat("\n*** unloading mgcv package / loading gam package")
              if(isNamespaceLoaded("mgcv")){
                if(isNamespaceLoaded("caret")){unloadNamespace("caret")} ## need to unload caret before car
                if(isNamespaceLoaded("car")){unloadNamespace("car")} ## need to unload car before mgcv
                unloadNamespace("mgcv")
              }
              if(!isNamespaceLoaded("gam")){requireNamespace("gam", quietly = TRUE)}

            }

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.GAM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GAM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.GAM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.GAM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows,,drop=FALSE])))
#   proj <- as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows,,drop=FALSE])))

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# GBM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GBM_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GBM',
                   n.trees_optim = 1000),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "gbm")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'GBM_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.GBM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.GBM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.GBM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

#   proj <- predict(newdata, model=get_formal_model(object), n.trees = object@n.trees_optim, type = "response")
  proj <- predict(newdata, model=get_formal_model(object), fun=gbm::predict.gbm, n.trees = object@n.trees_optim, type = "response")

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.GBM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), n.trees = object@n.trees_optim, type = "response"))

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# GLM Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('GLM_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GLM'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "glm")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'GLM_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.GLM_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              ## transform matrix into dataframe to be able to use predict.glm fct
              if (inherits(newdata, 'matrix')) newdata <- as.data.frame(newdata)
              return(.predict.GLM_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.GLM_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.GLM_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata[not_na_rows,,drop=FALSE]))

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# MARS Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('MARS_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MARS'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "MARS")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MARS_biomod2_model'),
          function(object, newdata, ...){

            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.MARS_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MARS_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.MARS_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  ##' @note we have to handle separatly rasterstack depending on the presence or not
  ##' of factorial variable.
  fact.var <- is.factor(newdata)
  if(any(fact.var)){
    ## get factor levels
    fact.var.levels <- subset(levels(newdata), fact.var)
    proj <- calc(newdata, function(x) {
      xx <- data.frame(x)
      ## ensure that the data.frame has the right set of levels
      for(i in which(fact.var)){
        xx[[i]] <- factor(xx[[i]], levels = unlist(fact.var.levels[[i]]))
      }
      ## do the projection
      proj.out <- as.numeric(predict(get_formal_model(object), xx, type = 'response'))
      return(proj.out)
    })
  } else {
    proj <- predict(newdata, model = get_formal_model(object), type = 'response')
  }



  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.MARS_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), type = 'response'))


  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# MAXENT.Phillips Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('MAXENT.Phillips_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXENT.Phillips'),
         validity = function(object){
           # check model class
#            if(sum(! ( c("randomForest.formula", "randomForest") %in% class(object@model) ) ) > 0) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'MAXENT.Phillips_biomod2_model'),
          function(object, newdata, silent=TRUE, ...){
            args <- list(...)
            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.MAXENT.Phillips_biomod2_model.RasterStack(object, newdata, silent=TRUE, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.MAXENT.Phillips_biomod2_model.data.frame(object, newdata, silent=TRUE, ... ))
            } else{ stop("invalid newdata input") }

          })


.predict.MAXENT.Phillips_biomod2_model.RasterStack <- function(object, newdata, silent=TRUE,  ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  rm_tmp_files <- args$rm_tmp_files
  temp_workdir <- args$temp_workdir
  split.proj <- args$split.proj

#   if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if(is.null(split.proj)) split.proj <- 1

  MWD <- .Prepare.Maxent.Proj.WorkDir(Data = newdata, species.name = object@resp_name,
                                      silent = TRUE, split.proj = split.proj )

  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if(!file.exists(path_to_maxent.jar)){
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }

  if(!silent) cat("\n\t\tRunning Maxent...")

  for(spl in 1:split.proj){
    maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated),"",paste0("-mx", object@model_options$memory_allocated, "m")),
                            " -cp ", "\"", path_to_maxent.jar, "\"",
                            " density.Project ",
                            "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                            "\"", MWD$m_workdir[[spl]], "\" ",
                            "\"", file.path(MWD$m_workdir[[spl]], "projMaxent.asc"), "\" ",
                            " doclamp=false visible=false autorun nowarnings notooltips")
    system(command = maxent.command, wait = TRUE, intern=TRUE)
  }

  if(!silent) cat("\n\t\tReading Maxent outputs...")

  ## get the list of projections part by part
  proj.list <- lapply(file.path(unlist(MWD$m_workdir),"projMaxent.asc"), raster, RAT = FALSE, crs=projection(newdata))
  ## merge all parts in a single raster
  if(length(proj.list) > 1){
    proj <- do.call(raster::merge, proj.list)
  } else {
    proj <- proj.list[[1]]
  }

#   # keep the coordinates ref system of new data
#   # TO DO => do it in .Prepare.Maxent.Proj.WorkDir()
#   proj.ref <- projection(newdata)

  #

  #   if(length(get_scaling_model(object))){
  #     names(proj) <- "pred"
  #     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  #   }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  } else if(!inMemory(proj)){
    proj <- readAll(proj) # to prevent from tmp files removing
  }




  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
#       unlink(x=file.path(object@resp_name, temp_workdir), recursive=TRUE, force=TRUE )
      .Delete.Maxent.WorkDir(MWD, silent=silent)
    }
  }

  return(proj)
}

.predict.MAXENT.Phillips_biomod2_model.data.frame <- function(object, newdata, silent=TRUE, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  temp_workdir <- args$temp_workdir
  rm_tmp_files <- args$rm_tmp_files
  xy <- args$xy

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
#   if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) rm_tmp_files <- TRUE

#   if( is.null(xy) ){
#     if( sum(c('x','y') %in% colnames(newdata) ) == 2 ){
#       coor_col <- c( which(colnames(newdata) == 'x'), which(colnames(newdata) == 'y') )
#       xy <- newdata[,coor_col]
#       newdata <- newdata[,- coor_col]
#     } else {
#       xy <- data.frame(x=rep(0,nrow(newdata)), y=rep(0,nrow(newdata)))
#     }
#   }

  ## no xy needed for models projections
  xy <- NULL

  ## check if na occurs in newdata cause they are not well supported
  not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})

  MWD <- .Prepare.Maxent.Proj.WorkDir(Data = as.data.frame(newdata[not_na_rows,,drop=FALSE]), xy = xy , species.name = object@resp_name, silent=T)
#   .Prepare.Maxent.Proj.WorkDir(Data = newdata, species.name = object@resp_name, proj.name = temp_workdir )

  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if(!file.exists(path_to_maxent.jar)){
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }

  if(!silent) cat("\n\t\tRunning Maxent...")

  maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated),"",paste0("-mx", object@model_options$memory_allocated, "m")),
                           " -cp ", "\"", path_to_maxent.jar, "\"",
                           " density.Project ",
                           "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                           "\"", file.path(MWD$m_workdir, "Pred_swd.csv"), "\" ",
                           "\"", file.path(MWD$m_workdir, "projMaxent.asc") , "\" ",
                           "doclamp=false visible=false autorun nowarnings notooltips")

  system(command = maxent.command, wait = TRUE, intern=TRUE)


  if(!silent) cat("\n\t\tReading Maxent outputs...")
  proj <- as.numeric(read.asciigrid(file.path(MWD$m_workdir, "projMaxent.asc"))@data[,1])

  ## add original NAs in table
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(!is.null(rm_tmp_files)){
    if(rm_tmp_files){
      .Delete.Maxent.WorkDir(MWD, silent=silent)
    }
  }

#   if(length(get_scaling_model(object))){
#     proj <- data.frame(pred = proj)
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)

}

## MAXENT old Class ------------------------------------------------------------

##' @note Here is define a conversion of the MAXENT models until the biomod 2.xx
##'   This model class was equivqlent to the current 'MAXENT.Phillips_biomod2_model'

setClass('MAXENT_biomod2_model',
         contains = 'MAXENT.Phillips_biomod2_model',
         prototype = list(model_class = 'MAXENT'),
         validity = function(object){
           return(TRUE)
         })

## End MAXENT old Class --------------------------------------------------------

## Maxent Phillips 2 end ===========================================================================
setClass(
  'MAXENT.Phillips.2_biomod2_model',
  representation(),
  contains = 'biomod2_model',
  prototype = list(model_class = 'MAXENT.Phillips.2'),
  validity = function(object) checkmate::test_class(object@model, 'maxnet')
)

setMethod(
  'predict', 
  signature(object = 'MAXENT.Phillips.2_biomod2_model'),
  function(object, newdata, ...){
    args <- list(...)
    do_check <- args$do_check
    if(is.null(do_check)) do_check <- TRUE
    
    ## data checking
    if(do_check){
      newdata <- check_data_range(model=object, new_data=newdata)
    }
    
    if(inherits(newdata, 'Raster')){
      return(.predict.MAXENT.Phillips.2_biomod2_model.RasterStack(object, newdata, ... ))
    } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
      return(.predict.MAXENT.Phillips.2_biomod2_model.data.frame(object, newdata, ... ))
    } else { stop("invalid newdata input") }
    
  })

.predict.MAXENT.Phillips.2_biomod2_model.RasterStack <- function(object, newdata, ...){
  newdata.df <- 
    newdata %>%
    as.matrix()
  
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  
  proj <- 
    predict(
      object = get_formal_model(object), 
      newdata = newdata.df,
      clamp = FALSE,
      type = 'logistic'
    )[, 1]
  
  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  ## convert back to raster file
  proj.ras <- raster(newdata)
  proj.ras[apply(newdata.df, 1, function(.x) all(!is.na(.x)))] <- proj
  proj <- proj.ras
  
  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }
  
  return(proj)
}

.predict.MAXENT.Phillips.2_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na
  
  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE
  
  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }
  
  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), type='logistic')[, 1])
  
  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }
  
  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }
  
  if(on_0_1000) proj <- round(proj*1000)
  
  return(proj)
}


## Maxent Phillips 2 end ===========================================================================

## MAXENT.Tsuruoka Class -------------------------------------------------------

# setClass('MAXENT.Tsuruoka_biomod2_model',
#          representation(),
#          contains = 'biomod2_model',
#          prototype = list(model_class = 'MAXENT.Tsuruoka'),
#          validity = function(object){
#            # check model class
#            if(sum(! ( c("maxent") %in% class(object@model) ) ) > 0) return(FALSE)
#            return(TRUE)
#          })

# setMethod('predict', signature(object = 'MAXENT.Tsuruoka_biomod2_model'),
#           function(object, newdata, ...){
#             args <- list(...)
# 
#             do_check <- args$do_check
#             if(is.null(do_check)) do_check <- TRUE
# 
#             ## data checking
#             if(do_check){
#               newdata <- check_data_range(model=object, new_data=newdata)
#             }
# 
#             if(inherits(newdata, 'Raster')){
#               return(.predict.MAXENT.Tsuruoka_biomod2_model.RasterStack(object, newdata, ... ))
#             } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
#               return(.predict.MAXENT.Tsuruoka_biomod2_model.data.frame(object, newdata, ... ))
#             } else{ stop("invalid newdata input") }
# 
#           })

# .predict.MAXENT.Tsuruoka_biomod2_model.RasterStack <- function(object, newdata, ...){
#   args <- list(...)
#   filename <- args$filename
#   overwrite <- args$overwrite
#   on_0_1000 <- args$on_0_1000
# 
#   if (is.null(overwrite)) overwrite <- TRUE
#   if (is.null(on_0_1000)) on_0_1000 <- FALSE
# 
#   proj <- calc(newdata, function(x) {
#     proj.out <- rep(NA, nrow(x))
#     x.no.na <- na.omit(x)
#     if(nrow(x.no.na)){
#       proj.not.na <- as.numeric(predict.maxent(get_formal_model(object), x.no.na)[, '1'])
#       proj.out[-attr(x.no.na, "na.action")] <- proj.not.na
#     }
#     return(proj.out)
#     })
# 
#   if(length(get_scaling_model(object))){
#     names(proj) <- "pred"
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }
# 
#   if(on_0_1000) proj <- round(proj*1000)
# 
#   # save raster on hard drive ?
#   if(!is.null(filename)){
#     cat("\n\t\tWriting projection on hard drive...")
#     if(on_0_1000){ ## projections are stored as positive integer
#       writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
#     } else { ## keep default data format for saved raster
#       writeRaster(proj, filename=filename, overwrite=overwrite)
#     }
#     proj <- raster(filename,RAT=FALSE)
#   }
# 
#   return(proj)
# }

# .predict.MAXENT.Tsuruoka_biomod2_model.data.frame <- function(object, newdata, ...){
#   args <- list(...)
#   on_0_1000 <- args$on_0_1000
#   omit.na <- args$omit.na
# 
#   if (is.null(on_0_1000)) on_0_1000 <- FALSE
#   if (is.null(omit.na)) omit.na <- FALSE
# 
#   ## check if na occurs in newdata cause they are not well supported
#   if(omit.na){
#     not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
#   } else {
#     not_na_rows <- rep(T, nrow(newdata))
#   }
# 
#   proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]))[,'1'])
# 
#   ## add original NAs in table if it needed
#   if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
#     tmp <- rep(NA,length(not_na_rows))
#     tmp[not_na_rows] <- proj
#     proj <- tmp
#     rm('tmp')
#   }
# 
#   if(length(get_scaling_model(object))){
#     proj <- data.frame(pred = proj)
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }
# 
#   if(on_0_1000) proj <- round(proj*1000)
# 
#   return(proj)
# }

## End MAXENT.Tsuruoka Class ---------------------------------------------------

# RF Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('RF_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'RF'),
         validity = function(object){
           # check model class
           if(!inherits(object@model, "randomForest")) return(FALSE)
           return(TRUE)
         })

setMethod('predict', signature(object = 'RF_biomod2_model'),
          function(object, newdata, ...){
            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.RF_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.RF_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.RF_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- predict(newdata, model=get_formal_model(object), type='prob', index=2)

  if(length(get_scaling_model(object))){
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }

  return(proj)
}

.predict.RF_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na

  if (is.null(on_0_1000)) on_0_1000 <- FALSE
  if (is.null(omit.na)) omit.na <- FALSE

  ## check if na occurs in newdata cause they are not well supported
  if(omit.na){
    not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }

  proj <- as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows,,drop=FALSE]), type='prob')[,'1'])

  ## add original NAs in table if it needed
  if(sum(!not_na_rows) > 0 ){ # some NAs in formal dataset
    tmp <- rep(NA,length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if(length(get_scaling_model(object))){
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}





# SRE Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('SRE_biomod2_model',
         representation(extremal_conditions='data.frame'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'SRE'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'SRE_biomod2_model'),
          function(object, newdata, ...){
            args <- list(...)

            do_check <- args$do_check
            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

            if(inherits(newdata, 'Raster')){
              return(.predict.SRE_biomod2_model.RasterStack(object, newdata, ... ))
            } else if(inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')){
              return(.predict.SRE_biomod2_model.data.frame(object, newdata, ... ))
            } else{ stop("invalid newdata input") }

          })

.predict.SRE_biomod2_model.RasterStack <- function(object, newdata, ...){
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000

  if (is.null(overwrite)) overwrite <- TRUE
  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .sre.projection(NewData=newdata, ExtremCond=object@extremal_conditions)

  if(on_0_1000) proj <- round(proj*1000)

  # save raster on hard drive ?
  if(!is.null(filename)){
    cat("\n\t\tWriting projection on hard drive...")
    if(on_0_1000){ ## projections are stored as positive integer
      writeRaster(proj, filename=filename, overwrite=overwrite, datatype="INT2S", NAflag=-9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename=filename, overwrite=overwrite)
    }
    proj <- raster(filename,RAT=FALSE)
  }


  return(proj)
}

.predict.SRE_biomod2_model.data.frame <- function(object, newdata, ...){
  args <- list(...)
  on_0_1000 <- args$on_0_1000

  if (is.null(on_0_1000)) on_0_1000 <- FALSE

  proj <- .sre.projection(NewData=newdata, ExtremCond=object@extremal_conditions)

  if(on_0_1000) proj <- round(proj*1000)

  return(proj)
}

























# EM parent Class =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('biomod2_ensemble_model',
         representation(modeling.id='character'), ##maybe some additional args should be added here
         contains = 'biomod2_model',
         prototype = list(model_class = 'EM'),
         validity = function(object){
           return(TRUE)
         })

# EMmean Class -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
setClass('EMmean_biomod2_model',
         representation(),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMmean'),
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMmean_biomod2_model'),
          function(object, newdata=NULL, formal_predictions=NULL, ...){
            args <- list(...)

            do_check <- args$do_check

            if(is.null(do_check)) do_check <- TRUE

            ## data checking
            if(do_check){
              newdata <- check_data_range(model=object, new_data=newdata)
            }

#            ## check if models are formal loaded
#             if(is.character(object@model)){
#               cat("\n*** models that will be loaded", object@model)
#               mem.change <- mem_change(model_tmp <- lapply(object@model, function(x){
#                 cat("\n***\t\t loading ", file.path(object@resp_name, "models", object@modeling.id, x))
#                 return(get(load(file.path(object@resp_name, "models", object@modeling.id, x))))
#               }))
#               cat("\n*** mem change:", mem.change)
#               names(model_tmp) <- object@model
#               object@model <- model_tmp
#             }

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
         validity = function(object){
           return(TRUE)
         })

setMethod('predict', signature(object = 'EMmedian_biomod2_model'),
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
         validity = function(object){
           return(TRUE)
         })

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
         representation(alpha = 'numeric',
                        side = 'character'),
         contains = 'biomod2_ensemble_model',
         prototype = list(model_class = 'EMci',
                   alpha = 0.05,
                   side = 'superior'),
         validity = function(object){
           if(!(object@side %in% c('inferior','superior'))) stop("side arg should be 'inferior' or 'superior")
           return(TRUE)
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
         validity = function(object){
           return(TRUE)
         })

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
         validity = function(object){
           return(TRUE)
         })

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


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# BIOMOD models definition(to make it easier to access plot, predict, ...)
# Damien Georges, Maya Guéguen
# 20/11/2012, update 22/10/2021
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



###################################################################################################
## 7. biomod2_model
###################################################################################################

# 7.1 Class Definition ----------------------------------------------------------------------------
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
         validity = function(object) { # check that scaler is a glm if it is defined
           if(length(object@scaling_model) > 0 && !inherits(object@scaling_model, c("glm", "lm"))) { return(FALSE) } else { return(TRUE) }
         })

# 7.2 Getters -------------------------------------------------------------------------------------
setGeneric("get_formal_model", def = function(object) { standardGeneric("get_formal_model") })

setMethod('get_formal_model', signature('biomod2_model'), function(object) { return(object@model) })

setGeneric("get_scaling_model", def = function(object) { standardGeneric("get_scaling_model") })

setMethod('get_scaling_model', signature('biomod2_model'), function(object) { return(object@scaling_model) })


# 7.3 Other Functions -----------------------------------------------------------------------------
setMethod('show', signature('biomod2_model'),
          function(object)
          {
            .bmCat("'biomod2_model'")
            cat("\n\t model name :", object@model_name, fill = .Options$width)
            cat("\n\t model class :", object@model_class, fill = .Options$width)
            cat("\n\t This model", ifelse(length(object@scaling_model), "has", "doesn't have"), "its own scaler", fill = .Options$width)

            cat("\n")
            cat("\n\t response modelled:", object@resp_name, fill = .Options$width)
            cat("\n\n\t explanatory variables used:", fill = .Options$width)
            cat("\n\t", "name", "\t", "type", "\t", "range", fill = .Options$width)
            for (i in 1:length(object@expl_var_names)) {
              cat("\n\t", object@expl_var_names[i],"\t", object@expl_var_type[i], "\t", object@expl_var_range[[i]], fill = .Options$width)
            }

            cat("\n")
            cat("\n\t NOTE : ")
            cat("\n\t\t You can access 'formal' model with get_formal_model function")
            cat(ifelse(length(object@scaling_model), "\n\t\t You can access scaling model with get_scaling_model function\n", "\n"))

            .bmCat()
          })



###################################################################################################
## 8.1 ANN_biomod2_model
###################################################################################################

setClass('ANN_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'ANN'),
         validity = function(object) { if(!inherits(object@model, "nnet")) { return(FALSE) } else { return(TRUE) }})

setMethod('predict', signature(object = 'ANN_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "ANN", object, newdata, ...))
          })

.predict.ANN_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = "predict(newdata, get_formal_model(object), type = 'raw')", object, newdata, ...))
}

.predict.ANN_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = 555, predcommand = "as.numeric(predict(get_formal_model(object), newdata[not_na_rows, , drop = F], type = 'raw'))", object, newdata, ...))
}


###################################################################################################
## 8.2 CTA_biomod2_model
###################################################################################################

setClass('CTA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'CTA'),
         validity = function(object) { # check model class
           if (!inherits(object@model, "rpart")) { return(FALSE) } else { return(TRUE) }
         })

setMethod('predict', signature(object = 'CTA_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "CTA", object, newdata, ...))
          })

.predict.CTA_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = 123, predcommand = "predict(newdata, model = get_formal_model(object), type = 'prob', index = 2)", object, newdata, ...))
}

.predict.CTA_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = 123, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'prob')[, 2])", object, newdata, ...))
}  


###################################################################################################
## 8.3 FDA_biomod2_model
###################################################################################################

setClass('FDA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'FDA'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "fda")) { return(FALSE) } else { return(TRUE) }
         })

setMethod('predict', signature(object = 'FDA_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "FDA", object, newdata, ...))
          })

.predict.FDA_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = "predict(newdata, model = get_formal_model(object), type = 'posterior', index = 2)", object, newdata, ...))
}

.predict.FDA_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'posterior')[, 2])", object, newdata, omit.na = TRUE, ...))
}  


###################################################################################################
## 8.4 GAM_biomod2_model
###################################################################################################

setClass('GAM_biomod2_model',
         representation(model_subclass = 'character'), 
         contains = 'biomod2_model',
         prototype = list(model_class = 'GAM', model_subclass = 'GAM_mgcv'),
         validity = function(object) { ## check model class
           if ((!(object@model_subclass %in% c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv'))) ||
               (object@model_subclass %in% c('GAM_mgcv', 'GAM_gam') && !inherits(object@model, c("gam"))) ||
               (object@model_subclass == 'BAM_mgcv' && !inherits(object@model, c("bam")))) {
             return(FALSE)
           } else {
             return(TRUE)
           }
         })

setMethod('predict', signature(object = 'GAM_biomod2_model'),
          function(object, newdata, ...)
          {
            args <- list(...)
            do_check <- args$do_check
            if (is.null(do_check)) { do_check <- TRUE }

            if (object@model_subclass %in% c("GAM_mgcv", "BAM_mgcv")) {
              # cat("\n*** unloading gam package / loading mgcv package")
              if (isNamespaceLoaded("gam")) { unloadNamespace("gam") }
              if (!isNamespaceLoaded("mgcv")) { requireNamespace("mgcv", quietly = TRUE) }
            }

            if (object@model_subclass == "GAM_gam") {
              # cat("\n*** unloading mgcv package / loading gam package")
              if (isNamespaceLoaded("mgcv")) {
                if (isNamespaceLoaded("caret")) { unloadNamespace("caret")} ## need to unload caret before car
                if (isNamespaceLoaded("car")) { unloadNamespace("car") } ## need to unload car before mgcv
                unloadNamespace("mgcv")
              }
              if (!isNamespaceLoaded("gam")) { requireNamespace("gam", quietly = TRUE) }
            }
            
            ## data checking
            if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }
            
            if (inherits(newdata, 'Raster')) {
              return(.predict.GAM_biomod2_model.RasterStack(object, newdata, ...))
            } else if (inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')) {
              return(.predict.GAM_biomod2_model.data.frame(object, newdata, ...))
            } else {
              stop("invalid newdata input")
            }
          })

.predict.GAM_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = ".testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata)", object, newdata, ...))
}

.predict.GAM_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = 555, predcommand = "as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows, , drop = FALSE])))", object, newdata, ...))
} 


###################################################################################################
## 8.5 GBM_biomod2_model
###################################################################################################

setClass('GBM_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GBM', n.trees_optim = 1000),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "gbm")) { return(FALSE) } else { return(TRUE) }
         })

setMethod('predict', signature(object = 'GBM_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "GBM", object, newdata, ...))
          })

.predict.GBM_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = "predict(newdata, model = get_formal_model(object), fun = gbm::predict.gbm, n.trees = object@n.trees_optim, type = 'response')", object, newdata, ...))
}

.predict.GBM_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), n.trees = object@n.trees_optim, type = 'response'))", object, newdata, ...))
}  


###################################################################################################
## 8.6 GLM_biomod2_model
###################################################################################################

setClass('GLM_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GLM'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "glm")) { return(FALSE) } else { return(TRUE) }
         })

setMethod('predict', signature(object = 'GLM_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "GLM", object, newdata, ...))
          })

.predict.GLM_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = ".testnull(object = get_formal_model(object), Prev = 0.5 , dat = newdata)", object, newdata, ...))
}

.predict.GLM_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(.testnull(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows, , drop = FALSE])))", object, newdata, ...))
}


###################################################################################################
## 8.7 MARS_biomod2_model
###################################################################################################

setClass('MARS_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MARS'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, c('earth', 'MARS', 'mars'))) { return(FALSE) } else { return(TRUE) }
         })

setMethod('predict', signature(object = 'MARS_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "MARS", object, newdata, ...))
          })

.predict.MARS_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) { overwrite <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  ##' @note we have to handle separatly rasterstack depending on the presence or not
  ##' of factorial variable.
  fact.var <- is.factor(newdata)
  if (any(fact.var))
  {
    ## get factor levels
    fact.var.levels <- subset(levels(newdata), fact.var)
    proj <- calc(newdata, function(x)
    {
      xx <- data.frame(x)
      ## ensure that the data.frame has the right set of levels
      for (i in which(fact.var)) {
        xx[[i]] <- factor(xx[[i]], levels = unlist(fact.var.levels[[i]]))
      }
      ## do the projection
      proj.out <- as.numeric(predict(get_formal_model(object), xx, type = 'response'))
      return(proj.out)
    })
  } else {
    proj <- predict(newdata, model = get_formal_model(object), type = 'response')
  }
  
  if (length(get_scaling_model(object))) {
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }
  
  if (on_0_1000) { proj <- round(proj * 1000)}
  
  # save raster on hard drive ?
  if (!is.null(filename)) {
    cat("\n\t\tWriting projection on hard drive...")
    if (on_0_1000) { ## projections are stored as positive integer
      writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename = filename, overwrite = overwrite)
    }
    proj <- raster(filename, RAT = FALSE)
  }
  return(proj)
}

.predict.MARS_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'response'))", object, newdata, ...))
}


###################################################################################################
## 8.8 MAXENT.Phillips_biomod2_model
###################################################################################################

setClass('MAXENT.Phillips_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXENT.Phillips'),
         validity = function(object) { return(TRUE) })

setMethod('predict', signature(object = 'MAXENT.Phillips_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "MAXENT.Phillips", object, newdata, ...))
          })

.predict.MAXENT.Phillips_biomod2_model.RasterStack <- function(object, newdata, silent = TRUE,  ...)
{
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  rm_tmp_files <- args$rm_tmp_files
  temp_workdir <- args$temp_workdir
  split.proj <- args$split.proj

  # if (is.null(temp_workdir)) temp_workdir <- paste("maxentWDtmp", format(Sys.time(), "%s"), sep="")
  if (is.null(rm_tmp_files)) { rm_tmp_files <- TRUE }
  if (is.null(overwrite)) { overwrite <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(split.proj)) { split.proj <- 1 }

  MWD <- .Prepare.Maxent.Proj.WorkDir(Data = newdata, species.name = object@resp_name,
                                      silent = TRUE, split.proj = split.proj )

  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if (!file.exists(path_to_maxent.jar)) {
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }
  
  if (!silent) { cat("\n\t\tRunning Maxent...") }
  for (spl in 1:split.proj) {
    maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated), "", paste0("-mx", object@model_options$memory_allocated, "m")),
                            " -cp ", "\"", path_to_maxent.jar, "\"",
                            " density.Project ",
                            "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                            "\"", MWD$m_workdir[[spl]], "\" ",
                            "\"", file.path(MWD$m_workdir[[spl]], "projMaxent.asc"), "\" ",
                            " doclamp=false visible=false autorun nowarnings notooltips")
    system(command = maxent.command, wait = TRUE, intern = TRUE)
  }
  
  if(!silent) { cat("\n\t\tReading Maxent outputs...") }

  ## get the list of projections part by part
  # check crs is not NA
  if (!is.na(projection(newdata))) {
    proj.list <- lapply(file.path(unlist(MWD$m_workdir), "projMaxent.asc"), raster, RAT = FALSE, crs = projection(newdata))
  } else {
    proj.list <- lapply(file.path(unlist(MWD$m_workdir), "projMaxent.asc"), raster, RAT = FALSE)
  }
  ## merge all parts in a single raster
  if (length(proj.list) > 1) {
    proj <- do.call(raster::merge, proj.list)
  } else {
    proj <- proj.list[[1]]
  }

  if (on_0_1000) { proj <- round(proj * 1000) }

  # save raster on hard drive ?
  if (!is.null(filename)) {
    cat("\n\t\tWriting projection on hard drive...")
    if (on_0_1000) { ## projections are stored as positive integer
      writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename = filename, overwrite = overwrite)
    }
    proj <- raster(filename, RAT = FALSE)
  } else if (!inMemory(proj)) {
    proj <- readAll(proj) # to prevent from tmp files removing
  }

  if (!is.null(rm_tmp_files) && rm_tmp_files) {
    .Delete.Maxent.WorkDir(MWD, silent = silent)
  }
  return(proj)
}

.predict.MAXENT.Phillips_biomod2_model.data.frame <- function(object, newdata, silent = TRUE, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  temp_workdir <- args$temp_workdir
  rm_tmp_files <- args$rm_tmp_files
  xy <- args$xy

  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(rm_tmp_files)) { rm_tmp_files <- TRUE }

  ## no xy needed for models projections
  xy <- NULL

  ## check if na occurs in newdata cause they are not well supported
  not_na_rows <- apply(newdata, 1, function(x){sum(is.na(x))==0})

  MWD <- .Prepare.Maxent.Proj.WorkDir(Data = as.data.frame(newdata[not_na_rows, , drop = FALSE])
                                      , xy = xy , species.name = object@resp_name, silent = T)

  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if (!file.exists(path_to_maxent.jar)) {
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }

  if (!silent) cat("\n\t\tRunning Maxent...")
  maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated), "", paste0("-mx", object@model_options$memory_allocated, "m")),
                           " -cp ", "\"", path_to_maxent.jar, "\"",
                           " density.Project ",
                           "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                           "\"", file.path(MWD$m_workdir, "Pred_swd.csv"), "\" ",
                           "\"", file.path(MWD$m_workdir, "projMaxent.asc") , "\" ",
                           "doclamp=false visible=false autorun nowarnings notooltips")
  system(command = maxent.command, wait = TRUE, intern = TRUE)

  if (!silent) { cat("\n\t\tReading Maxent outputs...") }
  proj <- as.numeric(read.asciigrid(file.path(MWD$m_workdir, "projMaxent.asc"))@data[, 1])

  ## add original NA from formal dataset
  if (sum(!not_na_rows) > 0) {
    tmp <- rep(NA, length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }

  if (!is.null(rm_tmp_files) && rm_tmp_files) {
    .Delete.Maxent.WorkDir(MWD, silent = silent)
  }
  if (on_0_1000) { proj <- round(proj * 1000) }
  return(proj)
}

##' @note Here is define a conversion of the MAXENT models until the biomod 2.xx
##'   This model class was equivqlent to the current 'MAXENT.Phillips_biomod2_model'

setClass('MAXENT_biomod2_model',
         contains = 'MAXENT.Phillips_biomod2_model',
         prototype = list(model_class = 'MAXENT'),
         validity = function(object) { return(TRUE) })


###################################################################################################
## 8.9 MAXENT.Phillips.2_biomod2_model
###################################################################################################

setClass( 'MAXENT.Phillips.2_biomod2_model',
          representation(),
          contains = 'biomod2_model',
          prototype = list(model_class = 'MAXENT.Phillips.2'),
          validity = function(object) { checkmate::test_class(object@model, 'maxnet') })

setMethod('predict', signature(object = 'MAXENT.Phillips.2_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "MAXENT.Phillips.2", object, newdata, ...))
          })

.predict.MAXENT.Phillips.2_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  newdata.df <- newdata %>% raster::as.matrix()
  
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) { overwrite <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  proj <- predict(object = get_formal_model(object), newdata = newdata.df, clamp = FALSE, type = 'logistic')[, 1]
  
  if (length(get_scaling_model(object))) {
    proj.to.scale <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj.to.scale)
  }
  
  if (on_0_1000) { proj <- round(proj * 1000) }
  
  ## convert back to raster file
  proj.ras <- raster(newdata)
  proj.ras[apply(newdata.df, 1, function(.x) { all(!is.na(.x)) })] <- proj
  proj <- proj.ras
  
  # save raster on hard drive ?
  if (!is.null(filename)) {
    cat("\n\t\tWriting projection on hard drive...")
    if (on_0_1000) { ## projections are stored as positive integer
      writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename = filename, overwrite = overwrite)
    }
    proj <- raster(filename, RAT = FALSE)
  }
  
  return(proj)
}

.predict.MAXENT.Phillips.2_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'logistic')[, 1])", object, newdata, ...))
}


###################################################################################################
## 8.9 MAXENT.Tsuruoka_biomod2_model
###################################################################################################

# setClass('MAXENT.Tsuruoka_biomod2_model',
#          representation(),
#          contains = 'biomod2_model',
#          prototype = list(model_class = 'MAXENT.Tsuruoka'),
#          validity = function(object) { ## check model class
#            if( sum(!(c("maxent") %in% class(object@model))) > 0) { return(FALSE) } else { return(TRUE) }
#          })
# 
# setMethod('predict', signature(object = 'MAXENT.Tsuruoka_biomod2_model'),
#           function(object, newdata, ...)
#           {
#             return(.template_predict(mod = "MAXENT.Tsuruoka", object, newdata, ...))
#           })
# 
# .predict.MAXENT.Tsuruoka_biomod2_model.RasterStack <- function(object, newdata, ...)*
# {
#   args <- list(...)
#   filename <- args$filename
#   overwrite <- args$overwrite
#   on_0_1000 <- args$on_0_1000
# 
#   if (is.null(overwrite)) { overwrite <- TRUE }
#   if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
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
#   if (length(get_scaling_model(object))) {
#     names(proj) <- "pred"
#     proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }
# 
#   if (on_0_1000) { proj <- round(proj * 1000) }
# 
#   # save raster on hard drive ?
#   if (!is.null(filename)) {
#     cat("\n\t\tWriting projection on hard drive...")
#     if (on_0_1000) { ## projections are stored as positive integer
#       writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
#     } else { ## keep default data format for saved raster
#       writeRaster(proj, filename = filename, overwrite = overwrite)
#     }
#     proj <- raster(filename, RAT = FALSE)
#   }
#   return(proj)
# }
# 
# .predict.MAXENT.Tsuruoka_biomod2_model.data.frame <- function(object, newdata, ...)
# {
#   return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]))[,'1'])", object, newdata, ...))
# }


###################################################################################################
## 8.10 RF_biomod2_model
###################################################################################################

setClass('RF_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'RF'),
         validity = function(object) { # check model class
           if (!inherits(object@model, "randomForest")) { return(FALSE)} else { return(TRUE) }
         })

setMethod('predict', signature(object = 'RF_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "RF", object, newdata, ...))
          })

.predict.RF_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = "predict(newdata, model = get_formal_model(object), type = 'prob', index = 2)", object, newdata, ...))
}

.predict.RF_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'prob')[, '1'])", object, newdata, ...))
}


###################################################################################################
## 8.11 SRE_biomod2_model
###################################################################################################

setClass('SRE_biomod2_model',
         representation(extremal_conditions='data.frame'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'SRE'),
         validity = function(object){ return(TRUE) })

setMethod('predict', signature(object = 'SRE_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "SRE", object, newdata, ...))
          })

.predict.SRE_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  return(.template_predict.RasterStack(seedval = NULL, predcommand = ".sre.projection(NewData = newdata, ExtremCond = object@extremal_conditions)", object, newdata, ...))
}

.predict.SRE_biomod2_model.data.frame <- function(object, newdata, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  proj <- .sre.projection(NewData = newdata, ExtremCond = object@extremal_conditions)
  if (on_0_1000) { proj <- round(proj * 1000) }
  return(proj)
}



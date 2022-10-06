###################################################################################################
##' @name BIOMOD_EnsembleForecasting
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' 
##' @title Project ensemble species distribution models onto new environment
##' 
##' @description This function allows to project ensemble models built with the 
##' \code{\link{BIOMOD_EnsembleModeling}} function onto new environmental data 
##' (\emph{which can represent new areas, resolution or time scales for example}).
##' 
##' 
##' @param bm.em a \code{\link{BIOMOD.ensemble.models.out}} object returned by the 
##' \code{\link{BIOMOD_EnsembleModeling}} function
##' @param bm.proj a \code{\link{BIOMOD.projection.out}} object returned by the 
##' \code{\link{BIOMOD_Projection}} function
##' @param proj.name (\emph{optional, default} \code{NULL}) \cr 
##' If \code{bm.proj = NULL}, a \code{character} corresponding to the name (ID) of the 
##' projection set (\emph{a new folder will be created within the simulation folder with this 
##' name})
##' @param new.env (\emph{optional, default} \code{NULL}) \cr 
##' If \code{bm.proj = NULL}, a \code{matrix}, \code{data.frame} or 
##' \code{\link[raster:stack]{RasterStack}} object containing the new explanatory variables (in 
##' columns or layers, with names matching the variables names given to the 
##' \code{\link{BIOMOD_FormatingData}} function to build \code{bm.em}) that will be used to 
##' project the ensemble species distribution model(s)
##' @param new.env.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{new.env} is a \code{matrix} or a \code{data.frame}, a 2-columns \code{matrix} or 
##' \code{data.frame} containing the corresponding \code{X} and \code{Y} coordinates that will 
##' be used to project the ensemble species distribution model(s)
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into binary values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{modeling.output}) or \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, 
##' \code{BIAS}, \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, 
##' \code{HK}, \code{HSS}, \code{OR}, \code{ORSS}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into filtered values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{modeling.output}) or \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, 
##' \code{BIAS}, \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, 
##' \code{HK}, \code{HSS}, \code{OR}, \code{ORSS}
##' 
##' @param compress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} or a \code{character} value defining whether and how objects should be 
##' compressed when saved on hard drive. Must be either \code{TRUE}, \code{FALSE}, \code{xz} or 
##' \code{gzip} (see Details)
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' 
##' @param \ldots (\emph{optional, see Details})
##' 
##' 
##' @return
##' 
##' A \code{BIOMOD.projection.out} object containing models projections, or links to saved 
##' outputs. \cr Models projections are stored out of \R (for memory storage reasons) in 
##' \code{proj.name} folder created in the current working directory :
##' \enumerate{
##'   \item the output is a 4-dimensional array if \code{new.env} is a \code{matrix} or a 
##'   \code{data.frame}
##'   \item it is a \code{rasterStack} if \code{new.env} is a \code{rasterStack} (or several 
##'   \code{rasterLayer} objects, if \code{new.env} is too large)
##'   \item raw projections, as well as binary and filtered projections (if asked), are saved in 
##'   the \code{proj.name} folder
##' }
##' 
##' 
##' @details 
##' 
##' If \code{models.chosen = 'all'}, projections are done for all evaluation and pseudo absences 
##' runs if applicable. \cr These projections may be used later by the 
##' \code{\link{BIOMOD_EnsembleForecasting}} function. \cr \cr
##' 
##' If \code{build.clamping.mask = TRUE}, a raster file will be saved within the projection 
##' folder. This mask values will correspond to the number of variables in each pixel that are out 
##' of their calibration / training range, identifying locations where predictions are uncertain. 
##' \cr \cr
##' 
##' \code{...} can take the following values :
##' \itemize{
##'   \item{\code{on_0_1000} : }{a \code{logical} value defining whether \code{0 - 1} 
##'   probabilities are to be converted to \code{0 - 1000} scale to save memory on backup}
##'   \item{\code{do.stack} : }{a \code{logical} value defining whether all projections are to be 
##'   saved as one \code{RasterStack} object or several \code{RasterLayer} files (\emph{the 
##'   default if projections are too heavy to be all loaded at once in memory})}
##'   \item{\code{keep.in.memory} : }{a \code{logical} value defining whether all projections are 
##'   to be kept loaded at once in memory, or only links pointing to hard drive are to be returned}
##'   \item{\code{output.format} : }{a \code{character} value corresponding to the projections 
##'   saving format on hard drive, must be either \code{.grd}, \code{.img} or \code{.RData} (the 
##'   default if \code{new.env} is given as \code{matrix} or \code{data.frame})}
##' }
##' 
##' 
##' @keywords models projection
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_ModelingOptions}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_RangeSize}}
##' @family Main functions
##' 
##'   
##' @examples
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
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExpl <- raster::stack(raster::crop(myExpl, myExtent))
##' }
##' 
##'  
##' # ---------------------------------------------------------------
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
##'                                       models = c('RF', 'GLM'),
##'                                       bm.options = myBiomodOptions,
##'                                       nb.rep = 2,
##'                                       data.split.perc = 80,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       do.full.models = FALSE,
##'                                       seed.val = 42)
##' }
##' 
##' 
##' file.proj <- paste0(myRespName, "/proj_Current/", myRespName, ".Current.projection.out")
##' if (file.exists(file.proj)) {
##'   myBiomodProj <- get(load(file.proj))
##' } else {
##' 
##'   # Project single models
##'   myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                     proj.name = 'Current',
##'                                     new.env = myExpl,
##'                                     models.chosen = 'all',
##'                                     metric.binary = 'all',
##'                                     metric.filter = 'all',
##'                                     build.clamping.mask = TRUE)
##' }
##' 
##' 
##' file.EM <- paste0(myRespName, "/", myRespName, ".AllModels.ensemble.models.out")
##' if (file.exists(file.EM)) {
##'   myBiomodEM <- get(load(file.EM))
##' } else {
##' 
##'   # Model ensemble models
##'   myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
##'                                         models.chosen = 'all',
##'                                         em.by = 'all',
##'                                         metric.select = c('TSS'),
##'                                         metric.select.thresh = c(0.7),
##'                                         metric.eval = c('TSS', 'ROC'),
##'                                         var.import = 3,
##'                                         prob.mean = TRUE,
##'                                         prob.median = FALSE,
##'                                         prob.cv = FALSE,
##'                                         prob.ci = FALSE,
##'                                         prob.ci.alpha = 0.05,
##'                                         committee.averaging = TRUE,
##'                                         prob.mean.weight = FALSE,
##'                                         prob.mean.weight.decay = 'proportional',
##'                                         seed.val = 42)
##' }
##' 
##' 
##' # ---------------------------------------------------------------
##' # Project ensemble models (from single projections)
##' myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, 
##'                                              bm.proj = myBiomodProj,
##'                                              models.chosen = 'all',
##'                                              metric.binary = 'all',
##'                                              metric.filter = 'all')
##' 
##' # Project ensemble models (building single projections)
##' myBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
##'                                              proj.name = 'CurrentEM',
##'                                              new.env = myExpl,
##'                                              models.chosen = 'all',
##'                                              metric.binary = 'all',
##'                                              metric.filter = 'all')
##' myBiomodEMProj
##' plot(myBiomodEMProj)
##' 
##' 
##' @importFrom raster stack subset predict writeRaster
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_EnsembleForecasting <- function(bm.em,
                                       bm.proj = NULL,
                                       proj.name = NULL,
                                       new.env = NULL,
                                       new.env.xy = NULL,
                                       models.chosen = 'all',
                                       metric.binary = NULL,
                                       metric.filter = NULL,
                                       compress = TRUE,
                                       nb.cpu = 1,
                                       ...)
{
  .bm_cat("Do Ensemble Models Projection")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_EnsembleForecasting.check.args(bm.em, bm.proj, proj.name, new.env,
                                                 models.chosen, metric.binary, metric.filter, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Create output object ----------------------------------------------------------------------
  proj_out <- new('BIOMOD.projection.out',
                  proj.name = proj.name,
                  dir.name =  bm.em@dir.name,
                  sp.name =  bm.em@sp.name,
                  expl.var.names = bm.em@expl.var.names,
                  models.projected = models.chosen,
                  coord = new.env.xy,
                  modeling.id = bm.em@modeling.id)
  proj_out@models.out@link = bm.em@link
  
  proj_is_raster <- FALSE
  if (inherits(new.env, 'Raster') || (length(bm.proj) && inherits(bm.proj@proj.out, 'BIOMOD.stored.raster.stack'))) {
    proj_is_raster <- TRUE
  }
  if (proj_is_raster) {
    proj_out@proj.out <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj.out <- new('BIOMOD.stored.array')
    do.stack = TRUE
  }
  
  ## 2. Create simulation directories -------------------------------------------------------------
  nameProj <- paste0("proj_", proj.name)
  nameProjSp <- paste0(nameProj, "_", bm.em@sp.name, "_ensemble")
  namePath <- file.path(bm.em@dir.name, bm.em@sp.name, nameProj)
  indiv_proj_dir <- .BIOMOD_EnsembleForecasting.prepare.workdir(dir.name = bm.em@dir.name
                                                                , sp.name = bm.em@sp.name
                                                                , proj.folder = nameProj)
  
  
  ## 3. Get needed projections --------------------------------------------------------------------
  needed_predictions <- get_needed_models(bm.em, models.chosen = models.chosen)
  if (length(bm.proj)) {
    formal_pred <- get_predictions(bm.proj,
                                   full.name = needed_predictions,
                                   as.data.frame = ifelse(bm.proj@type == 'array', TRUE, FALSE))
  } else {
    # make prediction according to given environment
    tmp_dir <- paste0('Tmp', format(Sys.time(), "%s"))
    formal_pred <- BIOMOD_Projection(bm.mod = load_stored_object(bm.em@models.out),
                                     new.env = new.env,
                                     proj.name = tmp_dir,
                                     new.env.xy = NULL,
                                     models.chosen = needed_predictions,
                                     compress = TRUE,
                                     build.clamping.mask = FALSE,
                                     do.stack = TRUE,
                                     on_0_1000 = on_0_1000,
                                     nb.cpu = nb.cpu)
    # getting the results
    formal_pred <- get_predictions(formal_pred,
                                   full.name = needed_predictions,
                                   as.data.frame = ifelse(inherits(new.env, 'Raster'), FALSE, TRUE))
    # remove tmp directory
    unlink(file.path(bm.em@dir.name, bm.em@sp.name, paste0("proj_", tmp_dir))
           , recursive = TRUE, force = TRUE)
  }
  
  ## 4. MAKING PROJECTIONS ------------------------------------------------------------------------
  ef.out <- NULL
  saved.files <- proj_names <- vector()
  for (em.comp in bm.em@em.computed[which(bm.em@em.computed %in% models.chosen)]) {
    cat("\n\t> Projecting", em.comp, "...")
    file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp,output.format))
    
    model.tmp <- NULL
    BIOMOD_LoadModels(bm.out = bm.em, full.name = em.comp, as = 'model.tmp')
    if (inherits(formal_pred, 'Raster')) {
      ef.tmp <- predict(model.tmp,
                        newdata = subset(formal_pred, subset = model.tmp@model, drop = FALSE),
                        data_as_formal_predictions = TRUE,
                        on_0_1000 = on_0_1000,
                        filename = ifelse(output.format == '.RData', '', file_name_tmp))
    } else {
      ef.tmp <- predict(model.tmp,
                        newdata = formal_pred[, model.tmp@model, drop = FALSE],
                        data_as_formal_predictions = TRUE,
                        on_0_1000 = on_0_1000)
    }
    
    if (!is.null(ef.tmp)) {
      proj_names <- c(proj_names, em.comp)
      if (inherits(ef.tmp, 'Raster')) {
        if (do.stack) {
          if (length(ef.out)) {
            ef.out <- stack(ef.out, ef.tmp)
          } else {
            ef.out <- stack(ef.tmp)
          }
        } else {
          file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp, output.format))
          if (output.format == '.RData') {
            save(ef.tmp, file = file_name_tmp, compress = compress)
          }
          saved.files <- c(saved.files, file_name_tmp)
        }
      } else if (!is.null(ef.tmp)) {
        # if (!is.null(ef.out) && length(ef.tmp) != nrow(ef.out)) {
        #   proj_names = proj_names[-length(proj_names)]
        #   warning(paste0("Predictions for ", em.comp, " do not match initial number of points (", length(ef.tmp), " instead of ", nrow(ef.out), ")"))
        # } else {
          ef.out <- cbind(ef.out, matrix(ef.tmp, ncol = 1))
        # }
      }
    }
  }
  
  # proj_out@models.projected <- bm.em@em.computed[which(bm.em@em.computed %in% models.chosen)]
  proj_out@models.projected <- proj_names
  
  if(do.stack) {
    if (inherits(ef.out, "Raster")) {
      names(ef.out) <- proj_out@models.projected
    } else {
      colnames(ef.out) <- proj_out@models.projected
    }
    # save object
    file_name_tmp <- file.path(namePath, paste0(nameProjSp, output.format))
    if (output.format == '.RData') {
      save(ef.out, file = file_name_tmp, compress = compress)
    } else if (inherits(ef.out, "Raster")) {
      ## TODO : define the raster dataformat (depends if em.cv has been computed)
      writeRaster(ef.out, filename = file_name_tmp, overwrite = TRUE)
    }
    saved.files <- c(saved.files, file_name_tmp)
  } 
  proj_out@proj.out@link <- saved.files #bm.em@em.computed
  
  if (!is.null(ef.out)) {
    proj_out@proj.out@val <- ef.out
    proj_out@proj.out@inMemory <- TRUE
  }
  
  
  ## 5. Compute binary and/or filtered transformation ---------------------------------------------
  if (length(metric.binary) | length(metric.filter))
  {
    cat("\n")
    eval.meth <- unique(c(metric.binary, metric.filter))
    
    ## Get all evaluation thresholds
    thresholds <- sapply(models.chosen, function(x) {
      get_evaluations(bm.em)[[x]][eval.meth, "Cutoff"]
    })
    if (!on_0_1000) { thresholds <- thresholds / 1000 }
    
    ## TODO : define the raster dataformat (depends if em.cv has been computed)
    ## Do binary transformation
    for (eval.meth in metric.binary) {
      cat("\n\t> Building", eval.meth, "binaries")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj.out@link)) {
          file.tmp <- proj_out@proj.out@link[i]
          writeRaster(x = bm_BinaryTransformation(raster(file.tmp, RAT = FALSE), thresholds[i]),
                      filename = sub(output.format,
                                     paste0("_", eval.meth, "bin", output.format),
                                     file.tmp),
                      overwrite = TRUE)
        }
      } else {
        nameBin = paste0(nameProjSp, "_", eval.meth, "bin")
        assign(x = nameBin, value = bm_BinaryTransformation(ef.out, thresholds))
        
        if (output.format == '.RData') {
          save(list = nameBin,
               file = file.path(namePath, paste0(nameBin, output.format)),
               compress = compress)
        } else {
          writeRaster(x = get(nameBin),
                      filename = file.path(namePath, paste0(nameBin, output.format)),
                      overwrite = TRUE,)
        }
      }
    }
    
    ## Do filtered transformation
    for (eval.meth in metric.filter) {
      cat("\n\t> Building", eval.meth, "filtered")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj.out@link)) {
          file.tmp <- proj_out@proj.out@link[i]
          writeRaster(x = bm_BinaryTransformation(raster(file.tmp, RAT = FALSE), thresholds[i], do.filtering = TRUE),
                      filename = sub(output.format,
                                     paste0("_", eval.meth, "filt", output.format),
                                     file.tmp),
                      overwrite = TRUE)
        }
      } else {
        nameFilt = paste0(nameProjSp, "_", eval.meth, "filt")
        assign(x = nameFilt, value = bm_BinaryTransformation(ef.out, thresholds, do.filtering = TRUE))
        if (output.format == '.RData') {
          save(list = nameFilt,
               file = file.path(namePath, paste0(nameFilt, output.format)),
               compress = compress)
        } else {
          writeRaster(x = get(nameFilt),
                      filename = file.path(namePath, paste0(nameFilt, output.format)),
                      overwrite = TRUE)
        }
      }
    }
    cat("\n")
  }
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  ## save a copy of output object without value to be lighter
  nameOut <- paste0(bm.em@sp.name, ".", proj.name, ".ensemble.projection.out")
  if (!keep.in.memory) { proj_out <- free(proj_out) }
  assign(nameOut, proj_out)
  save(list = nameOut, file = file.path(namePath, nameOut))
  
  .bm_cat("Done")
  return(proj_out)
}



###################################################################################################

.BIOMOD_EnsembleForecasting.prepare.workdir <- function(dir.name, sp.name, proj.folder)
{
  cat("\nCreating suitable Workdir...\n")
  dir.create(file.path(dir.name, sp.name, proj.folder), showWarnings = FALSE, recursive = TRUE, mode = "777")
  indiv_proj_dir <- file.path(dir.name, sp.name, proj.folder, "individual_projections")
  dir.create(indiv_proj_dir, showWarnings = FALSE, recursive = TRUE, mode = "777")
  return(indiv_proj_dir)
}

###################################################################################################

.BIOMOD_EnsembleForecasting.check.args <- function(bm.em, bm.proj, proj.name, new.env,
                                                   models.chosen, metric.binary, metric.filter, ...)
{
  args <- list(...)
  
  ## 1. Check bm.em -----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.em", bm.em, "BIOMOD.ensemble.models.out")
  
  ## 2. Check needed data and predictions -------------------------------------
  if ((is.null(bm.proj) && is.null(new.env)) ||
      (!is.null(bm.proj) & !is.null(new.env))) {
    stop("You have to refer at one of 'bm.proj' or 'new.env' argument")
  }
  
  if (!is.null(bm.proj)) {
    .fun_testIfInherits(TRUE, "bm.proj", bm.proj, "BIOMOD.projection.out")
    
    ## check all needed predictions are available
    needed_pred <- get_needed_models(bm.em, models.chosen = models.chosen)  
    missing_pred <- needed_pred[!(needed_pred %in% bm.proj@models.projected)]
    if (length(missing_pred)) {
      stop("Some models predictions missing :", toString(missing_pred))
    }
  }
  
  ## 3. Check models.chosen ---------------------------------------------------
  if (models.chosen[1] == 'all') {
    models.chosen <- get_built_models(bm.em)
  } else {
    models.chosen <- intersect(models.chosen, get_built_models(bm.em))
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  ## 4. Check proj.name -------------------------------------------------------
  if (!length(proj.name) && !length(bm.proj)) {
    stop("You have to give a valid 'proj.name' if you don't work with bm.proj")
  } else if (!length(proj.name)) {
    proj.name <- bm.proj@proj.name
  }
  
  ## 5. Check metric.binary & metric.filter -----------------------------------
  if (!is.null(metric.binary) | !is.null(metric.filter)) {
    models.evaluation <- get_evaluations(bm.em)
    if (is.null(models.evaluation)) {
      warning("Binary and/or Filtered transformations of projection not ran because of models evaluation information missing")
    } else {
      available.evaluation <- unique(dimnames(models.evaluation[[1]])[[1]])
      if (!is.null(metric.binary) && metric.binary[1] == 'all') {
        metric.binary <- available.evaluation
      } else if (!is.null(metric.binary) && sum(!(metric.binary %in% available.evaluation)) > 0) {
        warning(paste0(toString(metric.binary[!(metric.binary %in% available.evaluation)]),
                       " Binary Transformation were switched off because no corresponding evaluation method found"))
        metric.binary <- metric.binary[metric.binary %in% available.evaluation]
      }
      
      if (!is.null(metric.filter) && metric.filter[1] == 'all') {
        metric.filter <- available.evaluation
      } else if (!is.null(metric.filter) && sum(!(metric.filter %in% available.evaluation)) > 0) {
        warning(paste0(toString(metric.filter[!(metric.filter %in% available.evaluation)]),
                       " Filtered Transformation were switched off because no corresponding evaluation method found"))
        metric.filter <- metric.filter[metric.filter %in% available.evaluation]
      }
    }
  }
  
  ## 6. Check output.format ---------------------------------------------------
  output.format <- args$output.format
  if (is.null(output.format)) {
    if (length(bm.proj) > 0) {
      output.format = ifelse(bm.proj@type != 'RasterStack', ".RData", ".grd")
    } else {
      output.format = ifelse(!inherits(new.env, 'Raster'), ".RData", ".grd")
    }
  }
  
  ## 7. Check do.stack --------------------------------------------------------
  do.stack <- args$do.stack
  if (is.null(do.stack)) {
    do.stack <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than bm.proj ones
    if (!is.null(bm.proj) &&
        all(grepl("individual_projections", bm.proj@proj.out@link))) {
      do.stack <- FALSE
    }
  }
  
  ## 8. Check keep.in.memory --------------------------------------------------
  keep.in.memory <- args$keep.in.memory
  if (is.null(keep.in.memory)) {
    keep.in.memory <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than bm.proj ones
    if (!is.null(bm.proj)) {
      keep.in.memory <- bm.proj@proj.out@inMemory
    }
  }
  
  ## 9. Check new.env.xy ------------------------------------------------------
  new.env.xy <- args$new.env.xy
  if (is.null(new.env.xy)) {
    if (!is.null(bm.proj)) {
      new.env.xy <- bm.proj@coord
    } else {
      new.env.xy <- matrix()
    }
  }
  
  return(list(bm.em = bm.em,
              bm.proj = bm.proj,
              models.chosen = models.chosen,
              proj.name = proj.name,
              metric.binary = metric.binary,
              metric.filter = metric.filter,
              output.format = output.format,
              compress = ifelse(is.null(args$compress), FALSE, args$compress),
              on_0_1000 = ifelse(is.null(args$on_0_1000), TRUE, args$on_0_1000),
              do.stack = do.stack,
              keep.in.memory = keep.in.memory,
              new.env.xy = new.env.xy))
}

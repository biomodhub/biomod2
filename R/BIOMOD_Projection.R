###################################################################################################
##' @name BIOMOD_Projection
##' @author Wilfried Thuiller, Damien Georges
##' 
##' @title Project a range of calibrated species distribution models onto new environment
##' 
##' @description This function allows to project a range of models built with the 
##' \code{\link{BIOMOD_Modeling}} function onto new environmental data (\emph{which can 
##' represent new areas, resolution or time scales for example}).
##' 
##' 
##' @param bm.mod a \code{\link{BIOMOD.models.out}} object returned by the 
##' \code{\link{BIOMOD_Modeling}} function
##' @param proj.name a \code{character} corresponding to the name (ID) of the projection set 
##' (\emph{a new folder will be created within the simulation folder with this name})
##' @param new.env A \code{matrix}, \code{data.frame} or
##'  \code{\link[terra:rast]{SpatRaster}} object containing the new 
##' explanatory variables (in columns or layers, with names matching the
##'  variables names given to the \code{\link{BIOMOD_FormatingData}} function to build 
##' \code{bm.mod}) that will be used to project the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} are still supported such as 
##' \code{RasterStack} objects. }
##' 
##' @param new.env.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{new.env} is a \code{matrix} or a \code{data.frame}, a 2-columns \code{matrix} or 
##' \code{data.frame} containing the corresponding \code{X} and \code{Y} coordinates that will be 
##' used to project the species distribution model(s)
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into binary values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{bm.mod}) or \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, 
##' \code{BIAS}, \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, 
##' \code{HK}, \code{HSS}, \code{OR}, \code{ORSS}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into filtered values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_Modeling}} function. Must be among \code{all} (same evaluation metrics than 
##' those of \code{bm.mod}) or \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, 
##' \code{BIAS}, \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, 
##' \code{HK}, \code{HSS}, \code{OR}, \code{ORSS}
##' 
##' @param compress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} or a \code{character} value defining whether and how objects should be 
##' compressed when saved on hard drive. Must be either \code{TRUE}, \code{FALSE}, \code{xz} or 
##' \code{gzip} (see Details)
##' @param build.clamping.mask (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether a clamping mask should be built and saved on hard 
##' drive or not (see Details)
##' 
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' 
##' @param \ldots (\emph{optional, see Details)}) 
##' 
##' 
##' @return
##' 
##' A \code{BIOMOD.projection.out} object containing models projections, or links to saved 
##' outputs. \cr Models projections are stored out of \R (for memory storage reasons) in 
##' \code{proj.name} folder created in the current working directory :
##' \enumerate{
##'   \item the output is a \code{data.frame} if \code{new.env} is a \code{matrix} or a 
##'   \code{data.frame}
##'   \item it is a \code{\link[terra:rast]{SpatRaster}} if \code{new.env} is a 
##'   \code{\link[terra:rast]{SpatRaster}} (or several \code{\link[terra:rast]{SpatRaster}} 
##'   objects, if \code{new.env} is too large)
##'   \item raw projections, as well as binary and filtered projections (if asked), are saved in 
##'   the \code{proj.name} folder
##' }
##' 
##' 
##' @details 
##' 
##' If \code{models.chosen = 'all'}, projections are done for all calibration and pseudo absences 
##' runs if applicable. \cr These projections may be used later by the 
##' \code{\link{BIOMOD_EnsembleForecasting}} function. \cr \cr 
##' 
##' If \code{build.clamping.mask = TRUE}, a raster file will be saved within the projection folder. 
##' This mask values will correspond to the number of variables in each pixel that are out of their 
##' calibration / validation range, identifying locations where predictions are uncertain. \cr \cr
##' 
##' \code{...} can take the following values :
##' \itemize{
##   \item{\code{clamping.level} : }{a \code{logical} value defining whether \code{clamping.mask} 
##   cells with at least one variable out of its calibration range are to be removed from the 
##   projections or not
##'   \item{\code{omit.na} : }{a \code{logical} value defining whether all not fully referenced 
##'   environmental points will get \code{NA} as predictions or not}
##'   \item{\code{on_0_1000} : }{a \code{logical} value defining whether \code{0 - 1} probabilities 
##'   are to be converted to \code{0 - 1000} scale to save memory on backup}
##'   \item{\code{do.stack} : }{a \code{logical} value defining whether all projections are to be 
##'   saved as one \code{\link[terra:rast]{SpatRaster}} object or several 
##'   \code{\link[terra:rast]{SpatRaster}} files (\emph{the default if projections are too heavy to 
##'   be all loaded at once in memory})}
##'   \item{\code{keep.in.memory} : }{a \code{logical} value defining whether all projections are 
##'   to be kept loaded at once in memory, or only links pointing to hard drive are to be returned}
##'   \item{\code{output.format} : }{a \code{character} value corresponding to the projections 
##'   saving format on hard drive, must be either \code{.grd}, \code{.img}, \code{.tif} or \code{.RData} (the 
##'   default if \code{new.env} is given as \code{matrix} or \code{data.frame})}
##' }
##' 
##' 
##' @keywords models projection
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_RangeSize}}
##' @family Main functions
##' 
##'   
##' @examples
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
##' # ---------------------------------------------------------------#
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
##' # ---------------------------------------------------------------#
##' # Project single models
##' file.proj <- paste0(myRespName, "/proj_Current/", myRespName, ".Current.projection.out")
##' if (file.exists(file.proj)) {
##'   myBiomodProj <- get(load(file.proj))
##' } else {
##' myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                   proj.name = 'Current',
##'                                   new.env = myExpl,
##'                                   models.chosen = 'all')
##' }
##' myBiomodProj
##' plot(myBiomodProj)
##' 
##' 
##' @importFrom foreach foreach %dopar% 
## @importFrom doParallel registerDoParallel
##' @importFrom terra rast subset nlyr writeRaster terraOptions wrap mem_info app
##' @importFrom utils capture.output
##' @importFrom abind asub
##' 
##' @export
##' 
##' 
###################################################################################################

BIOMOD_Projection <- function(bm.mod,
                              proj.name,
                              new.env,
                              new.env.xy = NULL,
                              models.chosen = 'all',
                              metric.binary = NULL,
                              metric.filter = NULL,
                              compress = TRUE,
                              build.clamping.mask = TRUE,
                              nb.cpu = 1,
                              seed.val = NULL,
                              ...) {
  .bm_cat("Do Single Models Projection")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_Projection.check.args(bm.mod, proj.name, new.env, new.env.xy
                                        , models.chosen, metric.binary, metric.filter, compress, seed.val, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  proj_out <- new('BIOMOD.projection.out',
                  proj.name = proj.name,
                  dir.name = bm.mod@dir.name,
                  sp.name =  bm.mod@sp.name,
                  expl.var.names = bm.mod@expl.var.names,
                  models.projected = models.chosen,
                  scale.models = bm.mod@scale.models,
                  coord = new.env.xy,
                  modeling.id = bm.mod@modeling.id)
  proj_out@models.out@link = bm.mod@link
  
  proj_is_raster <- FALSE
  if (inherits(new.env, 'SpatRaster')) {
    proj_is_raster <- TRUE
  }
  if (proj_is_raster) {
    proj_out@proj.out <- new('BIOMOD.stored.SpatRaster')
  } else {
    proj_out@proj.out <- new('BIOMOD.stored.data.frame')
  }
  
  ## 2. Create simulation directories -------------------------------------------------------------
  nameProj <- paste0("proj_", proj.name)
  nameProjSp <- paste0(nameProj, "_", bm.mod@sp.name)
  namePath <- file.path(bm.mod@dir.name, bm.mod@sp.name, nameProj)
  dir.create(namePath, showWarnings = FALSE, recursive = TRUE, mode = "777")
  if (!do.stack) {
    dir.create(file.path(namePath, "individual_projections"),
               showWarnings = FALSE, recursive = TRUE, mode = "777")
  }
  
  ## 3. Define the clamping mask ------------------------------------------------------------------
  if (build.clamping.mask) {
    cat("\n\t> Building clamping mask\n")
    nameMask <- paste0(nameProjSp, "_ClampingMask")
    MinMax <- get_formal_data(bm.mod, 'MinMax')
    assign(x = nameMask, value = .build_clamping_mask(new.env, MinMax))
    
    if (output.format == '.RData') {
      if(proj_is_raster){
        save(list = wrap(nameMask),
             file = file.path(namePath, 
                              paste0(nameProj, "_ClampingMask", output.format)),
             compress = compress)
      } else {
        save(list = nameMask,
             file = file.path(namePath, 
                              paste0(nameProj, "_ClampingMask", output.format)),
             compress = compress)
      }
    } else {
      writeRaster(x = get(nameMask),
                  filename = file.path(namePath, paste0(nameProj, "_ClampingMask", output.format)),
                  datatype = "INT2S",
                  NAflag = -9999,
                  overwrite = TRUE)
    }
  }
  
  ## 4. MAKING PROJECTIONS ------------------------------------------------------------------------
  
  
  if (nb.cpu > 1) {
    if (.getOS() != "windows") {
      if (!isNamespaceLoaded("doParallel")) { 
        if(!requireNamespace('doParallel', quietly = TRUE)) stop("Package 'doParallel' not found")
      }
      doParallel::registerDoParallel(cores = nb.cpu)
    } else {
      warning("Parallelisation with `foreach` is not available for Windows. Sorry.")
    }
  }
  
  proj <- foreach(mod.name = models.chosen) %dopar% {
      cat("\n\t> Projecting", mod.name, "...")
      if (do.stack) {
        filename <- NULL
      } else {
        filename <- file.path(namePath, "individual_projections",
                              paste0(nameProj, "_", mod.name, 
                                     ifelse(output.format == ".RData"
                                            , ".tif", output.format)))
      }
      
      mod <- get(BIOMOD_LoadModels(bm.out = bm.mod, full.name = mod.name))
      temp_workdir = NULL
      if (length(grep("MAXENT$", mod.name)) == 1) {
        temp_workdir = mod@model_output_dir
      }
      pred.tmp <- predict(mod, new.env, on_0_1000 = on_0_1000, 
                          filename = filename, omit.na = omit.na, 
                          temp_workdir = temp_workdir, seedval = seed.val, 
                          overwrite = TRUE)
      if (do.stack) {
        if (proj_is_raster) {
          return(wrap(pred.tmp)) 
        } else {
          return(pred.tmp)
        }
      } else {
        return(filename)
      }
    }
  
  ## Putting predictions into the right format
  if (do.stack){
    if (proj_is_raster) {
      proj <- rast(lapply(proj, rast)) # SpatRaster needs to be wrapped before saving
      names(proj) <- models.chosen
      proj <- wrap(proj)
      proj.trans <- proj
    } else {
      proj <- as.data.frame(proj)
      names(proj) <- models.chosen
      proj.trans <- proj
      proj <- .format_proj.df(proj, obj.type = "mod")
    }
    
    if (keep.in.memory) {
      proj_out@proj.out@val <- proj
      proj_out@proj.out@inMemory <- TRUE
    }
  }
  
  ## save projections
  proj_out@type <- class(new.env)
  if (!do.stack){
    saved.files = unlist(proj)
  } else {
    assign(x = nameProjSp, value = proj)
    saved.files <- file.path(namePath, paste0(nameProjSp, output.format))
    if (output.format == '.RData') {
      save(list = nameProjSp, file = saved.files, compress = compress)
    } else  {
      writeRaster(x = rast(get(nameProjSp)), filename = saved.files,
                  overwrite = TRUE, datatype = ifelse(on_0_1000, "INT2S", "FLT4S"), NAflag = -9999)
    }
  }
  proj_out@proj.out@link <- saved.files
  
  # now that proj have been saved, it can be unwrapped if it is a SpatRaster
  if (proj_is_raster && do.stack) {
    proj.trans <- rast(proj.trans) 
  }
  
  ## 5. Compute binary and/or filtered transformation ---------------------------------------------
  if (!is.null(metric.binary) | !is.null(metric.filter)) {
    cat("\n")
    saved.files.binary <- NULL
    saved.files.filtered <- NULL
    thresholds <- get_evaluations(bm.mod, full.name = models.chosen)
    if (!on_0_1000) { thresholds[, "cutoff"]  <- thresholds[, "cutoff"] / 1000 }
    
    ## Do binary/filtering transformation
    for (eval.meth in unique(c(metric.binary, metric.filter))) {
      thres.tmp <- thresholds[which(thresholds$metric.eval == eval.meth), ]
      rownames(thres.tmp) <- thres.tmp$full.name
      thres.tmp <- thres.tmp[models.chosen, "cutoff"]
      
      cat("\n\t> Building", eval.meth, "binaries / filtered")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj.out@link)) {
          file.tmp <- proj_out@proj.out@link[i]
          
          if (eval.meth %in% metric.binary) {
            file.tmp.binary <- sub(output.format,
                                   paste0("_", eval.meth, "bin", output.format),
                                   file.tmp)
            saved.files.binary <- c(saved.files.binary, file.tmp.binary)
            writeRaster(x = bm_BinaryTransformation(rast(file.tmp), thres.tmp[i]),
                        filename = file.tmp.binary,
                        overwrite = TRUE,
                        datatype = "INT2S",
                        NAflag = -9999)
          }
          
          if (eval.meth %in% metric.filter) {
            file.tmp.filtered <- sub(output.format,
                                     paste0("_", eval.meth, "filt", output.format),
                                     file.tmp)
            saved.files.filtered <- c(saved.files.filtered, file.tmp.filtered)
            writeRaster(x = bm_BinaryTransformation(rast(file.tmp), thres.tmp[i],
                                                    do.filtering = TRUE),
                        filename = file.tmp.filtered,
                        overwrite = TRUE,
                        datatype = ifelse(on_0_1000, "INT2S", "FLT4S"),
                        NAflag = -9999)
          }
        }
      } else {
        if (eval.meth %in% metric.binary) {
          nameBin <- paste0(nameProjSp, "_", eval.meth, "bin")
          assign(x = nameBin, value = bm_BinaryTransformation(proj.trans, thres.tmp))
          
          file.tmp.binary <- file.path(namePath, paste0(nameBin, output.format))
          saved.files.binary <- c(saved.files.binary, file.tmp.binary)
          
          if (output.format == '.RData') {
            if (proj_is_raster) {
              assign(x = nameBin, value = wrap(get(nameBin)))
            }  else {
              assign(x = nameBin, value = .format_proj.df((get(nameBin)), obj.type = "mod"))
            }
            save(list = nameBin,
                 file = file.tmp.binary,
                 compress = compress)
          } else {
            writeRaster(x = get(nameBin),
                        filename = file.tmp.binary,
                        overwrite = TRUE,
                        datatype = "INT2S",
                        NAflag = -9999)
          }
        }
        
        
        if (eval.meth %in% metric.filter) {
          nameFilt <- paste0(nameProjSp, "_", eval.meth, "filt")
          assign(x = nameFilt,
                 value = bm_BinaryTransformation(proj.trans, thres.tmp, 
                                                 do.filtering = TRUE))
          file.tmp.filtered <- file.path(namePath, paste0(nameFilt, output.format))
          saved.files.filtered <- c(saved.files.filtered, file.tmp.filtered)
          
          if (output.format == '.RData') {
            if (proj_is_raster) {
              assign(x = nameFilt, value = wrap(get(nameFilt)))
            } else {
              assign(x = nameFilt, value = .format_proj.df((get(nameFilt)), obj.type = "mod"))
            }
            save(list = nameFilt,
                 file = file.path(namePath, paste0(nameFilt, output.format)),
                 compress = compress)
          } else {
            writeRaster(x = get(nameFilt),
                        filename = file.tmp.filtered,
                        overwrite = TRUE ,
                        datatype = ifelse(on_0_1000, "INT2S", "FLT4S"),
                        NAflag = -9999)
          }
        }
      }
    }
    
    
    ### save binary/filtered file link into proj_out ----------------------------
    if (!is.null(metric.binary)) {
      proj_out@proj.out@link <- c(proj_out@proj.out@link, saved.files.binary)
    }
    if (!is.null(metric.filter)) {
      proj_out@proj.out@link <- c(proj_out@proj.out@link, saved.files.filtered)
    }
    cat("\n")
  }
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  ## save a copy of output object without value to be lighter
  nameOut <- paste0(bm.mod@sp.name, ".", proj.name, ".projection.out")
  if (!keep.in.memory) { proj_out <- free(proj_out) }
  assign(nameOut, proj_out)
  save(list = nameOut, file = file.path(namePath, nameOut))
  
  .bm_cat("Done")
  return(proj_out)
}

# .BIOMOD_Projection.check.args---------------------------------------------

.BIOMOD_Projection.check.args <- function(bm.mod, proj.name, new.env, new.env.xy,
                                          models.chosen, metric.binary, metric.filter, compress, seed.val, ...)
{
  args <- list(...)
  
  ## 1. Check bm.mod ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
  
  ## 2. Check proj.name -------------------------------------------------------
  if (is.null(proj.name)) {
    stop("\nYou must define a name for Projection Outputs")
  } else {
    dir.create(paste0(bm.mod@sp.name, '/proj_', proj.name, '/'), showWarnings = FALSE)
  }
  
  ## 3. Check new.env ---------------------------------------------------------
  .fun_testIfInherits(TRUE, "new.env", new.env, c('matrix', 'data.frame', 'SpatRaster','Raster'))
  
  if(inherits(new.env, 'matrix')){
    if (any(sapply(get_formal_data(bm.mod, "expl.var"), is.factor))) {
      stop("new.env cannot be given as matrix when model involves categorical variables")
    }
    new.env <- data.frame(new.env)
  }
  if (inherits(new.env, 'Raster')) {
    # conversion into SpatRaster
    if(any(raster::is.factor(new.env))){
      new.env <- categorical_stack_to_terra(raster::stack(new.env))
    } else {
      new.env <- rast(new.env)
    }
  }
  
  if (inherits(new.env, 'SpatRaster')) {
    .fun_testIfIn(TRUE, "names(new.env)", names(new.env), bm.mod@expl.var.names)
  } else {
    .fun_testIfIn(TRUE, "colnames(new.env)", colnames(new.env), bm.mod@expl.var.names)
  }
  
  ## 4. Check new.env.xy ------------------------------------------------------
  if (!is.null(new.env.xy)  & !inherits(new.env, 'SpatRaster')) {
    new.env.xy = as.data.frame(new.env.xy)
    if (ncol(new.env.xy) != 2 || nrow(new.env.xy) != nrow(new.env)) {
      stop("invalid xy coordinates argument given -- dimensions mismatch !")
    }
  } else {
    new.env.xy = data.frame()
  }
  ## 5. Check models.chosen ---------------------------------------------------
  if (models.chosen[1] == 'all') {
    models.chosen <- bm.mod@models.computed
  } else {
    models.chosen <- intersect(models.chosen, bm.mod@models.computed)
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  ## check that given models exist
  files.check <- paste0(bm.mod@dir.name, "/", bm.mod@sp.name, "/models/",
                        bm.mod@modeling.id, "/", models.chosen)
  
  not.checked.files <- grep('MAXENT|SRE', files.check)
  if (length(not.checked.files) > 0) {
    files.check <- files.check[-not.checked.files]
  }
  missing.files <- files.check[!file.exists(files.check)]
  if (length(missing.files) > 0) {
    stop(paste0("Projection files missing : ", toString(missing.files)))
    if (length(missing.files) == length(files.check)) {
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }
  
  ## 6. Check metric.binary & metric.filter -----------------------------------
  if (!is.null(metric.binary) | !is.null(metric.filter)) {
    models.evaluation <- get_evaluations(bm.mod)
    if (is.null(models.evaluation)) {
      warning("Binary and/or Filtered transformations of projection not ran because of models evaluation information missing")
    } else {
      available.evaluation <- unique(models.evaluation$metric.eval)
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
  
  ## 7. Check compress --------------------------------------------------------
  if (compress == 'xz') {
    compress <- ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }
  
  ## 9. Check output.format ---------------------------------------------------
  output.format <- args$output.format # raster output format
  if (!is.null(output.format)) {
    if (!output.format %in% c(".img", ".grd", ".tif", ".RData")) {
      stop(paste0("output.format argument should be one of '.img','.grd', '.tif' or '.RData'\n"
                  , "Note : '.img','.grd', '.tif' are only available if you give environmental condition as a SpatRaster object"))
    }
    if (output.format %in% c(".img", ".grd", ".tif") && !inherits(new.env, "SpatRaster")) {
      warning("output.format was automatically set to '.RData' because environmental conditions are not given as a raster object")
    }
  } else {
    output.format <- ifelse(!inherits(new.env, "SpatRaster"), ".RData", ".tif")
  }
  
  ## 9. Check do.stack --------------------------------------------------------
  do.stack <- ifelse(is.null(args$do.stack), TRUE, args$do.stack)
  if (!inherits(new.env, 'SpatRaster')) {
    if (!do.stack) { 
      cat("\n\t\t! 'do.stack' arg is always set as TRUE for data.frame/matrix dataset") 
    }
    do.stack <- TRUE
  } else if (do.stack) { 
    # test if there is enough memory to work with multilayer SpatRaster
    capture.output({
      test <-
        mem_info(
          subset(new.env, 1),
          n = 2 * length(models.chosen) + nlyr(new.env)
        )
    })
    if (test["needed"] >= test["available"]) { 
      terraOptions(todisk = TRUE) 
    }
    # option without capture.output but with :::
    # ncopies = 2 * length(models.chosen) + nlyr(new.env)
    # test <- new.env@ptr$mem_needs(terra:::spatOptions(ncopies = ncopies))
    # if (test[1] > test[2]) {
    #   terraOptions(todisk = TRUE) 
    # }
  } else if (output.format == "RData"){
    cat("\n\t\t! 'do.stack' arg is always set as TRUE for .RData output format") 
    do.stack <- TRUE
  }
  
  
  
  return(list(proj.name = proj.name,
              new.env = new.env,
              new.env.xy = new.env.xy,
              models.chosen = models.chosen,
              metric.binary = metric.binary,
              metric.filter = metric.filter,
              compress = compress,
              do.stack = do.stack,
              output.format = output.format,
              omit.na = ifelse(is.null(args$omit.na), TRUE, args$omit.na),
              do.stack = do.stack,
              keep.in.memory = ifelse(is.null(args$keep.in.memory), TRUE, args$keep.in.memory),
              on_0_1000 = ifelse(is.null(args$on_0_1000), TRUE, args$on_0_1000),
              seed.val = seed.val))
}


### .build_clamping_mask -------------------------------------------------------

.build_clamping_mask <- function(env, MinMax)
{
  if (inherits(env, 'SpatRaster')) { 
    ## raster case ------------------------------------------------------------
    env <- subset(env, names(MinMax))
    ## create an empty mask
    clamp.mask <- app(subset(env, 1), function(x) ifelse(is.na(x), NA, 0))
    ref.mask <- app(subset(env, 1), function(x) ifelse(is.na(x), NA, 1))
    
    for (e.v in names(env)) {
      if (!is.null(MinMax[[e.v]]$min)) { # numeric variable
        clamp.mask <- clamp.mask + 
          bm_BinaryTransformation(subset(env, e.v), MinMax[[e.v]]$max) + ## pix with values outside of calib range
          (ref.mask - bm_BinaryTransformation(subset(env, e.v), MinMax[[e.v]]$min)) ## pix with no values (within env[[1]] area)
      } else if (!is.null(MinMax[[e.v]]$levels)) { # factorial variable
        clamp.mask <- 
          clamp.mask + 
          app(subset(env, e.v),
              function(x) {
                ifelse(as.character(x) %in% MinMax[[e.v]]$levels,
                       1,
                       0)
              }) ## pix with values outside of calib range
      }
    }
    
  } else if (is.data.frame(env) | is.matrix(env) | is.numeric(env)) { ## matrix and data.frame case ---------------------------------------------
    env <- as.data.frame(env)
    
    # create an empty mask
    clamp.mask <- rep(0, nrow(env))
    
    for (e.v in names(MinMax)) {
      if (!is.null(MinMax[[e.v]]$min)) { # numeric variable
        clamp.mask <- clamp.mask + 
          bm_BinaryTransformation(env[, e.v], MinMax[[e.v]]$max) + ## pix with values outside of calib range
          (1 - bm_BinaryTransformation(env[, e.v], MinMax[[e.v]]$min)) ## pix with no values (within env[[1]] area)
      } else if (!is.null(MinMax[[e.v]]$levels)) { # factorial variable
        clamp.mask <- clamp.mask + (env[, e.v] %in% MinMax[[e.v]]$levels)
      }
    }
  } else {
    stop("Unsupported env arg") 
  }
  return(clamp.mask)
}


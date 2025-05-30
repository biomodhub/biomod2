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
##'  \code{\link[terra:rast]{SpatRaster}} object containing the new 
##' explanatory variables (in columns or layers, with names matching the
##'  variables names given to the \code{\link{BIOMOD_FormatingData}} function to build 
##' \code{bm.mod}) that will be used to project the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} are still supported such as 
##' \code{RasterStack} objects. }
##' 
##' @param new.env.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{new.env} is a \code{matrix} or a \code{data.frame}, a 2-columns \code{matrix} or 
##' \code{data.frame} containing the corresponding \code{X} and \code{Y} coordinates that will 
##' be used to project the ensemble species distribution model(s)
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function applied to \code{bm.mod}
##' 
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into binary values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_EnsembleModeling}} function. Must be among \code{all} (same evaluation 
##' metrics than those of \code{bm.mod}) or \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, 
##' \code{ACCURACY}, \code{BIAS}, \code{AUCroc}, \code{AUCprg}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, 
##' \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' \cr \emph{Note that this is for binary data only.}
##' @param metric.filter (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric names to be used to transform prediction values 
##' into filtered values based on models evaluation scores obtained with the 
##' \code{\link{BIOMOD_EnsembleModeling}} function. Must be among \code{all} (same evaluation 
##' metrics than those of \code{bm.mod}) or \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, 
##' \code{ACCURACY}, \code{BIAS}, \code{AUCroc}, \code{AUCprg}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, 
##' \code{CSI}, \code{ETS}, \code{BOYCE}, \code{MPA}
##' \cr \emph{Note that this is for binary data only.}
##' 
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr
##' A \code{logical} value defining whether ensemble model projection should ignore missing values 
##' in single model projections or not (\emph{ignored by \code{EMwmean} algorithm})
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computatio
##' 
##' @param \ldots (\emph{optional, see Details})
##' 
##' 
##' @return
##' 
##' A \code{\link{BIOMOD.projection.out}} object containing models projections, or links to saved 
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
##' \code{...} can take the following values :
##' \describe{
##'   \item{digits}{(\emph{optional, default} \code{0}) : \cr
##'   an \code{integer} value corresponding to the number of digits of the predictions
##'   }
##'   \item{on_0_1000}{(\emph{optional, default} \code{TRUE}) : \cr
##'   a \code{logical} value defining whether \code{0 - 1} probabilities are to be converted to 
##'   \code{0 - 1000} scale to save memory on backup
##'   }
##'   \item{keep.in.memory}{(\emph{optional, default} \code{TRUE}) : \cr
##'   a \code{logical} value defining whether all projections are to be kept loaded at once in 
##'   memory, or only links pointing to hard drive are to be returned
##'   }
##'   \item{do.stack}{(\emph{optional, default} \code{TRUE}) : \cr 
##'   a \code{logical} value defining whether all projections are to be saved as one 
##'   \code{\link[terra:rast]{SpatRaster}} object or several \code{\link[terra:rast]{SpatRaster}} 
##'   files (\emph{the default if projections are too heavy to be all loaded at once in memory})
##'   }
##'   \item{output.format}{(\emph{optional, default} \code{.RData} or \code{.tif}) : \cr
##'   a \code{character} value corresponding to the projections saving format on hard drive, must 
##'   be either \code{.grd}, \code{.img}, \code{.tif} or \code{.RData} (the default if 
##'   \code{new.env} is given as \code{matrix} or \code{data.frame})
##'   }
##'   \item{compress}{(\emph{optional, default} \code{TRUE}) : \cr
##'   a \code{logical} or a \code{character} value defining whether and how objects should be 
##'   compressed when saved on hard drive. Must be either \code{TRUE}, \code{FALSE}, \code{gzip} 
##'   (for Windows OS) or \code{xz} (for other OS)
##'   }
##' }
##' 
##' 
##' @keywords models projection
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{bm_ModelingOptions}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_RangeSize}}
##' @family Main functions
##' 
##'   
##' @examples
##' 
##' library(terra)
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
##'  
##' # --------------------------------------------------------------- #
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
##'                                        resp.var = myResp,
##'                                        resp.xy = myRespXY,
##'                                        expl.var = myExpl)
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS', 'AUCroc'),
##'                                       var.import = 3,
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
##'                                         em.algo = c('EMmean', 'EMca'),
##'                                         metric.select = c('TSS'),
##'                                         metric.select.thresh = c(0.7),
##'                                         metric.eval = c('TSS', 'AUCroc'),
##'                                         var.import = 3,
##'                                         seed.val = 42)
##' }
##' 
##' 
##' # --------------------------------------------------------------- #
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
##' @importFrom foreach foreach %dopar%
##' @importFrom terra rast `add<-` wrap writeRaster mask
##' 
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
                                       na.rm = TRUE,
                                       nb.cpu = 1,
                                       ...)
{
  .bm_cat("Do Ensemble Models Projection")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_EnsembleForecasting.check.args(bm.em = bm.em,
                                                 bm.proj = bm.proj, 
                                                 proj.name = proj.name, 
                                                 new.env = new.env, 
                                                 new.env.xy = new.env.xy,
                                                 models.chosen = models.chosen, 
                                                 metric.binary = metric.binary, 
                                                 metric.filter = metric.filter, 
                                                 na.rm = na.rm, 
                                                 ...)
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
                  modeling.id = bm.em@modeling.id,
                  data.type = bm.em@data.type)
  proj_out@models.out@link = bm.em@link
  proj_out@call <- match.call()
  
  proj_is_raster <- FALSE
  if (inherits(new.env, 'SpatRaster') || (length(bm.proj) && inherits(bm.proj@proj.out, 'BIOMOD.stored.SpatRaster'))) {
    proj_is_raster <- TRUE
  }
  if (proj_is_raster) {
    proj_out@proj.out <- new('BIOMOD.stored.SpatRaster')
  } else {
    proj_out@proj.out <- new('BIOMOD.stored.data.frame')
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
  models.needed <- get_kept_models(bm.em)
  if (!is.null(bm.proj)) {
    formal_pred <- get_predictions(bm.proj, full.name = models.needed)
  } else {
    # make prediction according to given environment
    tmp_dir <- paste0('Tmp', as.numeric(Sys.time()) * 100000)
    formal_pred <- BIOMOD_Projection(bm.mod = load_stored_object(bm.em@models.out),
                                     new.env = new.env,
                                     proj.name = tmp_dir,
                                     new.env.xy = NULL,
                                     models.chosen = models.needed,
                                     compress = TRUE,
                                     build.clamping.mask = FALSE,
                                     do.stack = TRUE,
                                     on_0_1000 = on_0_1000,
                                     nb.cpu = nb.cpu)
    formal_pred <- get_predictions(formal_pred, full.name = models.needed)
    
    # remove tmp directory
    unlink(file.path(bm.em@dir.name, bm.em@sp.name, paste0("proj_", tmp_dir))
           , recursive = TRUE, force = TRUE)
  }

  if (!proj_is_raster) {
    if (bm.em@data.type %in% c("multiclass", "ordinal")){
      formal_pred <- tapply(X = formal_pred$pred, INDEX = list(formal_pred$points, formal_pred$full.name), FUN = function(x){as.character(x[1])})
    } else {
      formal_pred <- tapply(X = formal_pred$pred, INDEX = list(formal_pred$points, formal_pred$full.name), FUN = mean)
    }
    formal_pred <- as.data.frame(formal_pred[, models.needed])
  }
  
  ## 4. MAKING PROJECTIONS ------------------------------------------------------------------------
  proj.em <- foreach(em.name = models.chosen) %dopar%
    {
      cat("\n\t> Projecting", em.name, "...")
      if (do.stack) {
        filename <- NULL
      } else {
        filename <- file.path(indiv_proj_dir, paste0(em.name, output.format))
      }
      
      mod <- get(BIOMOD_LoadModels(bm.out = bm.em, full.name = em.name))
      rm(list = em.name)
      ef.tmp <- predict(mod
                        , newdata = formal_pred
                        , on_0_1000 = on_0_1000
                        , data_as_formal_predictions = TRUE
                        , filename = filename
                        , mod.name = em.name
                        , na.rm = na.rm)

      ## cleaning 
      if (bm.em@data.type %in% c("count", "abundance")) {
        ef.tmp[ef.tmp < 0] <- 0
        ef.tmp <- round(ef.tmp, digits = digits)
      }
      
      if(bm.em@data.type == "ordinal" && !grepl("EMfreq|EMcv", em.name)){
        data_sp <- get_formal_data(bm.em, subinfo = "resp.var")
        if(proj_is_raster){
          ef.tmp <- terra::subst(ef.tmp, 1:length(levels(data_sp)), levels(data_sp))
        } else {
          ef.tmp <- factor(ef.tmp, levels = levels(data_sp))
        }
      }
      
      if (do.stack) {
        if (proj_is_raster) {
          return(wrap(ef.tmp)) 
        } else {
          return(ef.tmp)
        }
      } else {
        return(filename)
      }
    }
  proj_out@models.projected <- models.chosen
  
  ## putting predictions into the right format
  if (do.stack) {
    if (proj_is_raster) {
      proj.em <- rast(lapply(proj.em, rast)) # SpatRaster needs to be wrapped before saving
      names(proj.em) <- models.chosen
      proj.em <- wrap(proj.em)
      proj.trans <- proj.em
    } else {
      proj.em <- as.data.frame(proj.em)
      names(proj.em) <- models.chosen
      proj.trans <- proj.em
      proj.em <- .format_proj.df(proj.em, obj.type = "em")
    }
    
    if (keep.in.memory) {
      proj_out@proj.out@val <- proj.em
      proj_out@proj.out@inMemory <- TRUE
    }
  }
  
  ## save projections
  proj_out@type <- ifelse(is.null(new.env), bm.proj@type, .get_env_class(new.env))
  if (!do.stack) {
    saved.files = unlist(proj.em)
  } else {
    assign(x = nameProjSp, value = proj.em)
    saved.files <- file.path(namePath, paste0(nameProjSp, output.format))
    if (output.format == '.RData') {
      save(list = nameProjSp, file = saved.files, compress = compress)
    } else {
      writeRaster(x = rast(get(nameProjSp)), filename = saved.files
                  , overwrite = TRUE, NAflag = -9999
                  , datatype = ifelse(any(grepl("EMcv", models.chosen) | digits != 0), "FLT4S", "INT2S"))
    }
  }
  proj_out@proj.out@link <- saved.files
  
  ## now that proj have been saved, it can be unwrapped if it is a SpatRaster
  if (proj_is_raster && do.stack) {
    proj.trans <- rast(proj.trans) 
  }
  
  
  ## 5. Compute binary and/or filtered transformation ---------------------------------------------
  if (bm.em@data.type == "binary" && (length(metric.binary) > 0 || length(metric.filter) > 0))
  {
    cat("\n")
    saved.files.binary <- NULL
    saved.files.filtered <- NULL
    to.rm <- grepl("EMcv|EMci", models.chosen)
    if (any(to.rm)) {
      cat("\n! Binary/Filtered transformation automatically desactivated for ensemble models with coefficient of variation or confidence intervals")
      models.chosen <- models.chosen[!to.rm]
    }
    
    if (length(models.chosen) > 0) {
      thresholds <- get_evaluations(bm.em, full.name = models.chosen)
      if (!on_0_1000) { thresholds[, "cutoff"] <- thresholds[, "cutoff"] / 1000 }
      
      ## Do binary/filtering transformation -----------------------------------
      for (eval.meth in unique(c(metric.binary, metric.filter))) {
        thres.tmp <- thresholds[which(thresholds$metric.eval == eval.meth), ]
        rownames(thres.tmp) <- thres.tmp$full.name
        thres.tmp <- thres.tmp[models.chosen, "cutoff"]
        
        cat("\n\t> Building", eval.meth, "binaries / filtered")
        if (!do.stack) {
          ## NO stack + binary + filtering transformation ---------------------
          for (i in 1:length(proj_out@proj.out@link)) {
            file.tmp <- proj_out@proj.out@link[i]
            output.format.search <- paste0("\\",output.format)
            if (grepl(pattern = paste0(models.chosen, collapse = "|"), x = file.tmp)) {
              if (eval.meth %in% metric.binary) {
                file.tmp.binary <- sub(output.format.search,
                                       paste0("_", eval.meth, "bin", output.format),
                                       file.tmp)
                saved.files.binary <- c(saved.files.binary, file.tmp.binary)
                writeRaster(x = bm_BinaryTransformation(rast(file.tmp), thres.tmp[i]),
                            filename = file.tmp.binary,
                            overwrite = TRUE,
                            NAflag = -9999,
                            datatype = "INT2S")
              }
              
              if (eval.meth %in% metric.filter) {
                file.tmp.filtered <- sub(output.format.search,
                                         paste0("_", eval.meth, "filt", output.format),
                                         file.tmp)
                saved.files.filtered <- c(saved.files.filtered, file.tmp.filtered)
                writeRaster(x = bm_BinaryTransformation(rast(file.tmp), thres.tmp[i], do.filtering = TRUE),
                            filename = file.tmp.filtered,
                            overwrite = TRUE,
                            NAflag = -9999,
                            datatype = "INT2S")
              }
            }
          }
        } else {
          ## subset to remove EMcv/EMci models
          if (proj_is_raster) {
            proj.trans <- subset(proj.trans, models.chosen)
          } else {
            proj.trans <- proj.trans[, models.chosen, drop = FALSE]
          }

          ## Stack + binary transformation ------------------------------------
          if (eval.meth %in% metric.binary) {
            nameBin <- paste0(nameProjSp, "_", eval.meth, "bin")
            assign(x = nameBin, value = bm_BinaryTransformation(proj.trans, thres.tmp))
            
            file.tmp.binary <- file.path(namePath, paste0(nameBin, output.format))
            saved.files.binary <- c(saved.files.binary, file.tmp.binary)
            if (output.format == '.RData') {
              if (!proj_is_raster) {
                assign(x = nameBin, value = .format_proj.df((get(nameBin)), obj.type = "em"))
              }
              save(list = nameBin,
                   file = file.tmp.binary,
                   compress = compress)
            } else {
              writeRaster(x = get(nameBin),
                          filename = file.tmp.binary,
                          overwrite = TRUE,
                          NAflag = -9999,
                          datatype = "INT2S")
            }
          }
          
          ## Stack + filtering transformation ---------------------------------
          if (eval.meth %in% metric.filter) {
            nameFilt <- paste0(nameProjSp, "_", eval.meth, "filt")
            # assign(x = nameFilt, value = bm_BinaryTransformation(proj, thres.tmp, do.filtering = TRUE))
            assign(x = nameFilt, value = bm_BinaryTransformation(proj.trans, thres.tmp, do.filtering = TRUE))
            file.tmp.filtered <- file.path(namePath, paste0(nameFilt, output.format))
            saved.files.filtered <- c(saved.files.filtered, file.tmp.filtered)
            
            if (output.format == '.RData') {
              if (!proj_is_raster) {
                assign(x = nameFilt, value = .format_proj.df((get(nameFilt)), obj.type = "em"))
              }
              save(list = nameFilt,
                   file = file.tmp.filtered,
                   compress = compress)
            } else {
              writeRaster(x = get(nameFilt),
                          filename = file.tmp.filtered,
                          overwrite = TRUE ,
                          NAflag = -9999,
                          datatype = "INT2S")
            }
          }
        }
      }
      cat("\n")
    } else {
      metric.binary <- NULL
      metric.filter <- NULL
      cat("\n! no models left for binary/filtered transformation")
    }
  }
  
  ## Save binary/filtered files link into proj_out ----------------------------
  if (!is.null(metric.binary)) {
    proj_out@proj.out@link <- c(proj_out@proj.out@link, saved.files.binary)
  }
  if (!is.null(metric.filter)) {
    proj_out@proj.out@link <- c(proj_out@proj.out@link, saved.files.filtered)
  }
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE ----------------------------------------
  nameOut <- paste0(bm.em@sp.name, ".", proj.name, ".ensemble.projection.out")
  if (!keep.in.memory) { proj_out <- free(proj_out) }
  assign(nameOut, proj_out)
  save(list = nameOut, file = file.path(namePath, nameOut))
  
  .bm_cat("Done")
  return(proj_out)
}


###################################################################################################

.BIOMOD_EnsembleForecasting.check.args <- function(bm.em, bm.proj, proj.name
                                                   , new.env, new.env.xy
                                                   , models.chosen
                                                   , metric.binary, metric.filter
                                                   , na.rm, ...)
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
    needed_pred <- get_kept_models(bm.em)  
    missing_pred <- needed_pred[!(needed_pred %in% bm.proj@models.projected)]
    if (length(missing_pred) > 0) {
      stop("Some models predictions missing :", toString(missing_pred))
    }
  }
  
  ## 3. Check new.env ---------------------------------------------------------
  if (!is.null(new.env)) {
    .fun_testIfInherits(TRUE, "new.env", new.env, c('matrix', 'data.frame', 'SpatRaster','Raster'))
    
    if (inherits(new.env, 'matrix')) {
      if (any(sapply(get_formal_data(bm.em,"expl.var"), is.factor))) {
        stop("new.env cannot be given as matrix when model involves categorical variables")
      }
      new.env <- data.frame(new.env)
    } else if (inherits(new.env, 'data.frame')) {
      # ensure that data.table are coerced into classic data.frame
      new.env <- as.data.frame(new.env) 
    }
    
    if (inherits(new.env, 'Raster')) {
      # conversion into SpatRaster
      if(any(raster::is.factor(new.env))){
        new.env <- .categorical_stack_to_terra(raster::stack(new.env))
      } else {
        new.env <- rast(new.env)
      }
    }
    
    if (inherits(new.env, 'SpatRaster')) {
      .fun_testIfIn(TRUE, "names(new.env)", names(new.env), bm.em@expl.var.names, exact = TRUE)
      new.env <- new.env[[bm.em@expl.var.names]]
      new.env.mask <- .get_data_mask(new.env, value.out = 1)
      new.env <- mask(new.env, new.env.mask)
    } else {
      .fun_testIfIn(TRUE, "colnames(new.env)", colnames(new.env), bm.em@expl.var.names, exact = TRUE)
      new.env <- new.env[ , bm.em@expl.var.names, drop = FALSE]
    }
    
    which.factor <- which(sapply(new.env, is.factor))
    if (length(which.factor) > 0) {
      new.env <- .check_env_levels(new.env, expected_levels = head(get_formal_data(bm.em, subinfo = "expl.var")))
    }
  }
  ## 4. Check models.chosen ---------------------------------------------------
  if (models.chosen[1] == 'all') {
    models.chosen <- get_built_models(bm.em)
  } else {
    models.chosen <- intersect(models.chosen, get_built_models(bm.em))
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  ## 5. Check proj.name -------------------------------------------------------
  if (!length(proj.name) && !length(bm.proj)) {
    stop("You have to give a valid 'proj.name' if you don't work with bm.proj")
  } else if (!length(proj.name)) {
    proj.name <- bm.proj@proj.name
  }
  
  ## 6. Check metric.binary & metric.filter -----------------------------------
  if (!is.null(metric.binary) | !is.null(metric.filter)) {
    if (bm.em@data.type != "binary") {
      cat ("No metric.binary or metric.filter are needed with", bm.em@data.type, "data")
      metric.binary <- NULL
      metric.filter <- NULL
    } else {
      models.evaluation <- get_evaluations(bm.em)
      if (is.null(models.evaluation)) {
        warning("Binary and/or Filtered transformations of projection not ran because of models evaluation information missing")
      } else {
        available.evaluation <- as.character(unique(models.evaluation$metric.eval))
        if (!is.null(metric.binary) && metric.binary[1] == 'all') {
          metric.binary <- available.evaluation
        } else if (!is.null(metric.binary) && any(! metric.binary %in% available.evaluation)) {
          warning(paste0(toString(metric.binary[!(metric.binary %in% available.evaluation)]),
                         " Binary Transformation were switched off because no corresponding evaluation method found"))
          metric.binary <- metric.binary[metric.binary %in% available.evaluation]
        }
        
        if (!is.null(metric.filter) && metric.filter[1] == 'all') {
          metric.filter <- available.evaluation
        } else if (!is.null(metric.filter) && any(!(metric.filter %in% available.evaluation))) {
          warning(paste0(toString(metric.filter[!(metric.filter %in% available.evaluation)]),
                         " Filtered Transformation were switched off because no corresponding evaluation method found"))
          metric.filter <- metric.filter[metric.filter %in% available.evaluation]
        }
      }
    }
  }
  
  ## 7. Check output.format ---------------------------------------------------
  output.format <- args$output.format
  if (is.null(output.format)) {
    if (length(bm.proj) > 0) {
      output.format = ifelse(bm.proj@type != 'SpatRaster', ".RData", ".tif")
    } else {
      output.format = ifelse(!inherits(new.env, 'SpatRaster'), ".RData", ".tif")
    }
  }
  
  ## 8. Check do.stack --------------------------------------------------------
  do.stack <- args$do.stack
  if (is.null(do.stack)) {
    do.stack <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than bm.proj ones
    if (!is.null(bm.proj) && all(grepl("individual_projections", bm.proj@proj.out@link))) {
      do.stack <- FALSE
    }
  }
  
  ## 9. Check keep.in.memory --------------------------------------------------
  keep.in.memory <- args$keep.in.memory
  if (is.null(keep.in.memory)) {
    keep.in.memory <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than bm.proj ones
    if (!is.null(bm.proj)) {
      keep.in.memory <- bm.proj@proj.out@inMemory
    }
  }
  
  ## 10. Check new.env.xy ------------------------------------------------------
  if (is.null(new.env.xy)) {
    if (!is.null(bm.proj)) {
      new.env.xy <- bm.proj@coord
    } else {
      new.env.xy <- data.frame()
    }
  } else {
    new.env.xy <- as.data.frame(new.env.xy)
  }
  
  ## 11. Check na.rm ------------------------------------------------------
  stopifnot(is.logical(na.rm))
  
  ## 12.on_0_1000 --------------------------------
  on_0_1000 <- ifelse(is.null(args$on_0_1000), TRUE, args$on_0_1000)
  if (bm.em@data.type  %in% c("count", "abundance", "ordinal", "multiclass")) { on_0_1000 <- FALSE }
  
  return(list(bm.em = bm.em,
              bm.proj = bm.proj,
              proj.name = proj.name,
              new.env = new.env,
              new.env.xy = new.env.xy,
              models.chosen = models.chosen,
              metric.binary = metric.binary,
              metric.filter = metric.filter,
              digits = ifelse(is.null(args$digits), 0, args$digits),
              on_0_1000 = on_0_1000,
              keep.in.memory = keep.in.memory,
              do.stack = do.stack,
              output.format = output.format,
              compress = ifelse(is.null(args$compress), FALSE, args$compress)))
}


###################################################################################################

.BIOMOD_EnsembleForecasting.prepare.workdir <- function(dir.name, sp.name, proj.folder)
{
  cat("\nCreating suitable Workdir...\n")
  indiv_proj_dir <- file.path(dir.name, sp.name, proj.folder, "individual_projections")
  dir.create(indiv_proj_dir, showWarnings = FALSE, recursive = TRUE, mode = "777")
  return(indiv_proj_dir)
}


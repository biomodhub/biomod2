###################################################################################################
##' @name BIOMOD_RangeSize
##' @author Maya Guegen Helene Blancheteau
##' 
##' @title Analyze the range size differences between projections of species distribution models
##' 
##' @description This function allows to calculate the absolute number of locations (pixels) lost, 
##' stable and gained, as well as the corresponding relative proportions, between two (or more) 
##' binary projections of (ensemble) species distribution models (\emph{which can represent new 
##' time scales or environmental scenarios for example}).
##' 
##' 
##' @param proj.current a \code{BIOMOD.projection.out} object containing the initial projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param proj.future a \code{BIOMOD.projection.out} object containing the final binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_projected_models}} function
##' @param metric.binary (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric selected to transform predictions into binary 
##' values, must be among \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, 
##' \code{BIAS}, \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, 
##' \code{ETS}, \code{BOYCE}, \code{MPA}
##' 
##' @return
##' 
##' A \code{BIOMOD.rangesize.out} containing principaly two objects :
##' \describe{
##'   \item{Compt.By.Species}{a \code{data.frame} containing the summary of range change for each 
##'   comparison
##'   \itemize{
##'     \item \code{Loss} : number of pixels predicted to be lost
##'     \item \code{Stable_Abs} : number of pixels not currently occupied and not predicted to be
##'     \item \code{Stable_Pres} : number of pixels currently occupied and predicted to remain 
##'     occupied
##'     \item \code{Gain} : number of pixels predicted to be gained
##'     \item \code{PercLoss} : percentage of pixels currently occupied and predicted to be lost 
##'     (\code{Loss / (Loss + Stable_Pres)})
##'     \item \code{PercGain} : percentage of pixels predicted to be gained compare to the 
##'     number of pixels currently occupied (\code{Gain / (Loss + Stable_Pres)})
##'     \item \code{SpeciesRangeChange} : percentage of pixels predicted to change (loss or gain) 
##'     compare to the number of pixels currently occupied (\code{PercGain - PercLoss})
##'     \item \code{CurrentRangeSize} : number of pixels currently occupied
##'     \item \code{FutureRangeSize0Disp} : number of pixels predicted to be occupied, assuming 
##'     no migration
##'     \item \code{FutureRangeSize1Disp} : number of pixels predicted to be occupied, assuming 
##'     migration
##'   }
##'   }
##'   \item{loss.gain}{an object in the same form than the input data (\code{proj.current} and 
##'   \code{proj.future}) and containing a value for each point/pixel of each comparison among :
##'   \itemize{
##'     \item \code{-2} : predicted to be lost
##'     \item \code{-1} : predicted to remain occupied
##'     \item \code{0} : predicted to remain unoccupied
##'     \item \code{1} : predicted to be gained
##'   }
##'   }
##'   \item{Diff.By.Pixel}{an object in the same form than the input data (\code{proj.current} and 
##'   \code{proj.future}) and containing a value for each point/pixel of each comparison obtain with :
##'   \itemize{
##'   \item \code{Future - 2* Current} for binary data
##'   \item \code{Future - Current} for nonbinary after rescaling Future and Current from 0 to 1.
##'   }
##'   }
##' }
##' 
##' 
##' @keywords "species range change" projections gain loss
##' 
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}, 
##' \code{\link{bm_PlotRangeSize}}
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
##'                                       metric.eval = c('TSS', 'ROC'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' }
##' 
##' models.proj <- get_built_models(myBiomodModelOut, algo = "RF")
##'   # Project single models
##'   myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                     proj.name = 'CurrentRangeSize',
##'                                     new.env = myExpl,
##'                                     models.chosen = models.proj,
##'                                     metric.binary = 'all',
##'                                     build.clamping.mask = TRUE)
##' 
##' 
##' # --------------------------------------------------------------- #
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' data(bioclim_future)
##' myExplFuture <- terra::rast(bioclim_future)
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExplFuture <- terra::crop(myExplFuture, myExtent)
##' }
##' # Project onto future conditions
##' myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                               proj.name = 'FutureRangeSize',
##'                                               new.env = myExplFuture,
##'                                               models.chosen = models.proj,
##'                                               metric.binary = 'TSS')
##'                                               
##' # Compute differences
##' myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = myBiomodProj,
##'                                       proj.future = myBiomodProjectionFuture,
##'                                       metric.binary = "TSS")
##' 
##' 
##' # Represent main results
##' bm_PlotRangeSize(bm.range = myBiomodRangeSize)
##' 
##' 
##' 
##' @importFrom foreach foreach %do%
##' @importFrom terra rast nlyr `add<-` values
##' @importFrom dplyr count
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_RangeSize <- function(proj.current, proj.future, 
                             models.chosen = "all", metric.binary = NULL){
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future, 
                                       models.chosen, metric.binary)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  pred_current <- get_predictions(proj.current,
                                  full.name = models.chosen,
                                  metric.binary = metric.binary,
                                  model.as.col = TRUE)
  
  pred_future <- get_predictions(proj.future,
                                 full.name = models.chosen,
                                 metric.binary = metric.binary,
                                 model.as.col = TRUE)
  
  coord <- proj.future@coord
  if(is.data.frame(pred_future)){
    if(nrow(pred_current) != nrow(pred_future)){
      good <- intersect(rownames(pred_current), rownames(pred_future))
      pred_current <- as.data.frame(pred_current[rownames(pred_current) %in% good,])
      pred_future <- as.data.frame(pred_future[rownames(pred_future) %in% good,])
      #coord <- as.data.frame(proj.future@coord[rownames(proj.future@coord) %in% good,])
    } 
  }
  
  
  ordinal <- (proj.current@data.type == "ordinal")
  bm.range <- bm_RangeSize(proj.current = pred_current, proj.future = pred_future, ordinal = ordinal)
  
  link.models <- proj.current@models.out@link
  if(grepl("ensemble.models.out", link.models)){
    row.names <- c("Species", "EMmodel", "Merged_Dataset", "Merged_Run", "Merged_Algo")
  } else {
    row.names <- c("Species", "Dataset", "Run", "Algo")
  }
  
  
  out <- new("BIOMOD.rangesize.out",
             Compt.By.Models = bm.range$Compt.By.Models,
             Diff.By.Pixel = bm.range$Diff.By.Pixel,
             loss.gain = bm.range$loss.gain,
             data.type = bm.range$data.type,
             coord = coord,
             row.names = row.names)
  
  return(out)
}


.BIOMOD_RangeSize.check.args <- function(proj.current, proj.future, 
                                         models.chosen, metric.binary){
  
  if(inherits(proj.current, "SpatRaster") | inherits(proj.current, "data.frame") |inherits(proj.current, "raster")){
    stop("BIOMOD.RangeSize now accepts BIOMOD.projection.out objects directly. \n
         If you want to compare ", class(proj.current), " objects, use bm_RangeSize function. ")
  } # To Remove after a while. 
  
  .fun_testIfInherits(TRUE, "proj.current", proj.current, "BIOMOD.projection.out")
  .fun_testIfInherits(TRUE, "proj.future", proj.future, "BIOMOD.projection.out")
  
  ## Check models.chosen ---------------------------------------------------
  if (models.chosen[1] == 'all') {
    models.chosen <- proj.future@models.projected
    models.chosen <- intersect(models.chosen, proj.current@models.projected)
  } else {
    models.chosen <- intersect(models.chosen, proj.current@models.projected)
    models.chosen <- intersect(models.chosen, proj.future@models.projected)
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  if (proj.current@data.type == "binary" && is.null(metric.binary)){
    stop("Metric.binary is necessary for binary data.")
  }
  
  if (proj.current@data.type != "binary" && !is.null(metric.binary)){
    metric.binary <- NULL
    warning("metric.binary set to NULL for non binary data.")
  }
  
  
  return(list(models.chosen = models.chosen,
              metric.binary = metric.binary))
}


### Class

##' @name BIOMOD.rangesize.out
##' @aliases BIOMOD.rangesize.out-class
##' @author Helene Blancheteau
##' 
##' @title \code{BIOMOD_RangeSize()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_RangeSize}}, and used by 
##' \code{\link{bm_PlotRangeSize}}
##' 
##' 
##' @slot Compt.By.Models a \code{data.frame} containing the summary of range change for each comparison
##' @slot Diff.By.Pixel a \code{SpatRaster} or a \code{data.frame} containing a value for each point/pixel of each comparison
##' @slot loss.gain a \code{SpatRaster} or a \code{data.frame} containg for each point/pixel 
##' a value indicating a loss, a gain or a stable sitatution
##' @slot data.type a \code{character} corresponding to the data type
##' @slot coord a 2-columns \code{matrix} or \code{data.frame} containing the corresponding 
##' \code{X} and \code{Y} coordinates used to project the species distribution model(s)
##' @slot row.names A \code{vector} containing tags matching models names splitted by 
##' '_' character
##' 
##' @family Toolbox objects
##' 
NULL

##' @name BIOMOD.rangesize.out-class
##' @rdname BIOMOD.rangesize.out
##' @export
##' 

# Class Definition  -----------------------------------

setClass("BIOMOD.rangesize.out",
         representation(Compt.By.Models = 'data.frame',
                        Diff.By.Pixel = 'ANY',
                        loss.gain = 'ANY',
                        data.type = 'character',
                        coord = 'data.frame',
                        row.names = 'character'),
         prototype(data.type = "binary"),
         validity = function(object){ return(TRUE) })               
  
        

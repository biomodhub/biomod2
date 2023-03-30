###################################################################################################
##' @name BIOMOD_LoadModels
##' @author Damien Georges
##' 
##' @title Load species distribution models built with \pkg{biomod2}
##' 
##' @description This function loads individual models built with \code{\link{BIOMOD_Modeling}} 
##' or \code{\link{BIOMOD_EnsembleModeling}} functions.
##' 
##' 
##' @param bm.out a \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' object that can be obtained with the \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' 
##' @param full.name (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing model names to be kept, must be either \code{all} or a 
##' sub-selection of model names that can be obtained with the \code{\link{get_built_models}} 
##' function
##'
##' @param PA (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing pseudo-absence set to be loaded, must be among \code{PA1}, 
##' \code{PA2}, \code{...}, \code{allData}
##' @param run (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing repetition set to be loaded, must be among \code{RUN1}, 
##' \code{RUN2}, \code{...}, \code{allRun}
##' @param algo (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} containing algorithm to be loaded, must be either \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT}, \code{MAXNET}
##' 
##' @param merged.by.PA (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing merged pseudo-absence set to be loaded, must be among \code{PA1}, 
##' \code{PA2}, \code{...}, \code{mergedData}
##' @param merged.by.run (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing merged repetition set to be loaded, must be among \code{RUN1}, 
##' \code{RUN2}, \code{...}, \code{mergedRun}
##' @param merged.by.algo (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} containing merged algorithm to be loaded, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT}, \code{MAXNET}, \code{mergedAlgo}
##' @param filtered.by (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing evaluation metric selected to filter single models to build the 
##' ensemble models, must be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, 
##' \code{BIAS}, \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, 
##' \code{HK}, \code{HSS}, \code{OR}, \code{ORSS}
##' 
##' 
##' @return 
##' 
##' A \code{vector} containing the names of the loaded models.
##' 
##' 
##' @details
##' 
##' This function might be of particular use to load models and make response plot analyses. \cr \cr
##' 
##' Running the function providing only \code{bm.out} argument will load all models built by the 
##' \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} function, but a 
##' subselection of models can be done using the additional arguments (\code{full.name}, \code{PA}, 
##' \code{run}, \code{algo}, \code{merged.by.PA}, \code{merged.by.run}, \code{merged.by.algo}, 
##' \code{filtered.by}).
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}
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
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' }
##' 
##' 
##' # ---------------------------------------------------------------
##' # Loading some models built
##' BIOMOD_LoadModels(bm.out = myBiomodModelOut, algo = 'RF')
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_LoadModels <- function(bm.out, full.name = NULL, PA = NULL, run = NULL, algo = NULL
                              , merged.by.PA = NULL, merged.by.run = NULL
                              , merged.by.algo = NULL, filtered.by = NULL)
{
  # .bm_cat("Load Models")

  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_LoadModels.check.args(bm.out)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## Get names of models to load
  if (inherits(bm.out, "BIOMOD.models.out")) {
    models.to.load <- get_built_models(bm.out, full.name = full.name, PA = PA, run = run, algo = algo)
  } else if (inherits(bm.out, "BIOMOD.ensemble.models.out")) {
    models.to.load <- get_built_models(bm.out, full.name = full.name, merged.by.PA = merged.by.PA
                                       , merged.by.run = merged.by.run, merged.by.algo = merged.by.algo
                                       , filtered.by = filtered.by, algo = algo)
  }
  envir <- parent.frame()
  
  if (length(models.to.load) == 0) {
    cat("\n   ! No models computed matched, No models loaded !")
    return(NULL)
  }
  
  ## LOAD the selected models -----------------------------------------------------------
  filename = file.path(bm.out@dir.name, bm.out@sp.name, "models", bm.out@modeling.id)
  # .bm_cat("Done")
  for (mtl in models.to.load) {
    load(file = file.path(bm.out@dir.name, bm.out@sp.name, "models", bm.out@modeling.id, mtl), envir = envir)
  }
  return(models.to.load)
}


###################################################################################################

.BIOMOD_LoadModels.check.args <- function(bm.out)
{
  ## 1. Check bm.out ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
}


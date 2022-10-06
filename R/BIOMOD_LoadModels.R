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
##' @param \ldots (\emph{optional, see Details)}) 
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
##' \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} function. But a 
##' subselection of models can be done using the following additional arguments :
##' \itemize{
##'   \item{\code{models} : }{a \code{vector} containing model names to be loaded, must be among 
##'   \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, 
##'   \code{MARS}, \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}}
##'   \item{\code{run.eval} : }{a \code{vector} containing repetition set to be loaded, must be 
##'   among \code{RUN1}, \code{RUN2}, \code{...}, \code{Full}}
##'   \item{\code{data.set} : }{a \code{vector} containing pseudo-absence set to be loaded, must 
##'   be among \code{PA1}, \code{PA2}, \code{...} \cr \cr}
##'   \item{\code{path} : }{a \code{character} corresponding to the location of the species folder 
##'   (if different from the current working directory) \cr \cr}
##'   \item{\code{full.name} : }{a \code{vector} containing model names to be kept, must be either 
##'   \code{all} or a sub-selection of model names \cr \cr}
##'   \item{\code{as} : }{a \code{character} to contain the loaded models}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}
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
##' # ---------------------------------------------------------------
##' # Loading some models built
##' myLoadedModels <- BIOMOD_LoadModels(bm.out = myBiomodModelOut, models = 'RF')
##' myLoadedModels
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_LoadModels <- function(bm.out, ... )
{
  # .bm_cat("Load Models")

  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_LoadModels.check.args(bm.out, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## Get names of models to load
  models.to.load <- get_built_models(bm.out)
  envir <- parent.frame()
  
  ## Create list or a sub-list of models to load ----------------------------------------
  if (!is.null(full.name)) {
    models.to.load <- full.name
  } else { ## make a subselection
    
    ## subselection on models
    if (!is.null(models)) {
      models.to.load <- grep(paste0(models, collapse = "|"), models.to.load, value = TRUE)
    }
    
    ## subselection on run.Eval
    if (!is.null(run.eval)) {
      models.to.load <- grep(paste0(run.eval, collapse = "|"), models.to.load, value = TRUE)
    }
    
    ## subselection on data.set
    if (!is.null(data.set)) {
      models.to.load <- grep(paste0(data.set, collapse = "|"), models.to.load, value = TRUE)
    }
  }
  
  if (length(models.to.load) == 0) {
    cat("\n   ! No models computed matched, No models loaded !")
    return(NULL)
  }
  
  ## LOAD the selected models -----------------------------------------------------------
  nameFile = file.path(bm.out@dir.name, bm.out@sp.name, "models", bm.out@modeling.id)
  # .bm_cat("Done")
  if (!is.null(as) && length(models.to.load) == 1) {
    assign(x = as,
           value = get(load(file = file.path(nameFile, models.to.load))),
           envir = envir)
    invisible(TRUE)
  } else {
    for (mtl in models.to.load) {
      load(file = file.path(bm.out@dir.name, bm.out@sp.name, "models", bm.out@modeling.id, mtl), envir = envir)
    }
    return(models.to.load)
  }
}


###################################################################################################

.BIOMOD_LoadModels.check.args <- function(bm.out, ...)
{
  args <- list(...)
  
  ## 1. Check bm.out ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check args ------------------------------------------------------------
  available.args <- c("models", "run.eval", "data.set", "path", "as", "full.name")
  .fun_testIfIn(TRUE, "names(args)", names(args), available.args)
  avail_models <- get_built_models(bm.out) ## get all available model names
  
  ## 2.1 Check add.args : models ----------------------------------------------
  models <- args$models
  if (!is.null(models)) {
    infos = .extract_modelNamesInfo(model.names = avail_models, info = 'models')
    .fun_testIfIn(TRUE, "models", models, infos)
    models = paste0("_", models)
  }
  
  ## 2.2 Check add.args : run.eval --------------------------------------------
  run.eval <- args$run.eval
  if (!is.null(run.eval)) {
    infos = .extract_modelNamesInfo(model.names = avail_models, info = 'run.eval')
    .fun_testIfIn(TRUE, "run.eval", run.eval, infos)
    run.eval = paste0("_", run.eval)
  }
  
  ## 2.3 Check add.args : data.set --------------------------------------------
  data.set <- args$data.set
  if (!is.null(data.set)) {
    infos = .extract_modelNamesInfo(model.names = avail_models, info = 'data.set')
    .fun_testIfIn(TRUE, "data.set", data.set, infos)
    data.set = paste0("_", data.set, "_")
  }
  
  ## 2.4 Check path : data.set ------------------------------------------------
  path <- args$path
  if (!is.null(path) && !(file.path(bm.out@dir.name, bm.out@sp.name) %in% list.dirs(path = path))) {
    stop("invalid path given")
  } else {
    path = "."
  }
  
  ## 2.5 Check path : full.name ------------------------------------------------
  full.name <- args$full.name
  if (!is.null(full.name)) {
    .fun_testIfIn(TRUE, "full.name", full.name, avail_models)
  }
  
  return(list(models = models,
              run.eval = run.eval,
              data.set = data.set, 
              path = path, 
              as = args$as, 
              full.name = full.name))
}


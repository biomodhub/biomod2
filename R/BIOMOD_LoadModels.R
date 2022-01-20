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
##' @param bm.out a \code{BIOMOD.models.out} (or \code{BIOMOD.ensemble.models.out}) object that 
##' can be obtained from \code{\link{BIOMOD_Modeling}} (or \code{\link{BIOMOD_EnsembleModeling}}) 
##' function
##' @param ... (\emph{optional, see \href{BIOMOD_LoadModels.html#details}{Details})}) 
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
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}
##' 
##' 
##' @examples
##' 
##' # species occurrences
##' myFile <- system.file("external/species/mammals_table.csv", package="biomod2")
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0("external/bioclim/current/bio", c(3, 4, 7, 11, 12), ".grd")
##' myExpl = raster::stack(system.file(myFiles[1], package = "biomod2"),
##'                        system.file(myFiles[2], package = "biomod2"),
##'                        system.file(myFiles[3], package = "biomod2"),
##'                        system.file(myFiles[4], package = "biomod2"),
##'                        system.file(myFiles[5], package = "biomod2"))
##' 
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' # 3. Doing Modelisation
##' myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
##'                                     models = c('SRE', 'CTA', 'RF'), 
##'                                     models.options = myBiomodOption, 
##'                                     NbRunEval = 1, 
##'                                     DataSplit = 80, 
##'                                     Yweights = NULL, 
##'                                     VarImport = 0, 
##'                                     models.eval.meth = c('TSS','ROC'),
##'                                     SaveObj = TRUE,
##'                                     rescal.all.models = FALSE,
##'                                     do.full.models = FALSE)
##'                                     
##' # 4. Loading some models built
##' myLoadedModels <- BIOMOD_LoadModels(myBiomodModelOut, models = 'RF')
##' myLoadedModels
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_LoadModels <- function(bm.out, ... )
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  add.args <- list(...)
  args <- .BIOMOD_LoadModels.check.args(bm.out, add.args)
  add.args <- args$add.args
  rm(args)
  
  ## Get names of models to load
  models.to.load <- get_built_models(bm.out)
  envir <- parent.frame()
  
  ## Create list or a sub-list of models to load ----------------------------------------
  if (!is.null(add.args$full.name)) {
    models.to.load <- add.args$full.name
  } else { ## make a subselection
    
    ## subselection on models
    if (!is.null(add.args$models)) {
      model.to.load.tmp <- c()
      for (mod in add.args$models) {
        if (sum(grepl(mod, models.to.load)) > 0) {
          model.to.load.tmp <- c(model.to.load.tmp, grep(mod, models.to.load, value = TRUE))
        }
      }
      models.to.load <- model.to.load.tmp
    }
    
    ## subselection on run.Eval
    if (!is.null(add.args$run.eval)) {
      model.to.load.tmp <- c()
      for (re in add.args$run.eval) {
        if (sum(grepl(re, models.to.load)) > 0) {
          model.to.load.tmp <- c(model.to.load.tmp, grep(re, models.to.load, value = TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }
    
    ## subselection on data.set
    if (!is.null(add.args$data.set)) {
      model.to.load.tmp <- c()
      for (ds in add.args$data.set) {
        if (sum(grepl(ds,  models.to.load)) > 0) {
          model.to.load.tmp <- c(model.to.load.tmp, grep(ds, models.to.load, value = TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }
  }
  
  if (length(models.to.load) == 0) {
    cat("\n   ! No models computed matched, No models loaded !")
    return(NULL)
  }
  
  ## LOAD the selected models -----------------------------------------------------------
  nameFile = file.path(bm.out@sp.name, "models", bm.out@modeling.id)
  if (!is.null(add.args$as) && length(models.to.load) == 1) {
    assign(x = add.args$as,
           value = get(load(file = file.path(nameFile, models.to.load))),
           envir = envir)
    invisible(TRUE)
  } else {
    for (mtl in models.to.load) {
      load(file = file.path(bm.out@sp.name, "models", bm.out@modeling.id, mtl),
           envir = envir)
    }
    return(models.to.load)
  }
}


###################################################################################################

.BIOMOD_LoadModels.check.args <- function(bm.out, add.args)
{
  ## 1. Check bm.out ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check add.args --------------------------------------------------------
  available.args <- c("models", "run.eval", "data.set", "path", "as", "full.name")
  .fun_testIfIn(TRUE, "names(add.args)", names(add.args), available.args)
  avail_models <- get_built_models(bm.out) ## get all available model names
  
  ## 2.1 Check add.args : models ----------------------------------------------
  if (!is.null(add.args$models)) {
    infos = .extract_modelNamesInfo(model.names = avail_models, info = 'models')
    .fun_testIfIn(TRUE, "add.args$models", add.args$models, infos)
    add.args$models = paste0("_", add.args$models)
  }
  
  ## 2.2 Check add.args : run.eval --------------------------------------------
  if (!is.null(add.args$run.eval)) {
    infos = .extract_modelNamesInfo(model.names = avail_models, info = 'run.eval')
    .fun_testIfIn(TRUE, "add.args$run.eval", add.args$run.eval, infos)
    add.args$run.eval = paste0("_", add.args$run.eval)
  }
  
  ## 2.3 Check add.args : data.set --------------------------------------------
  if (!is.null(add.args$data.set)) {
    infos = .extract_modelNamesInfo(model.names = avail_models, info = 'data.set')
    .fun_testIfIn(TRUE, "add.args$data.set", add.args$data.set, infos)
    add.args$data.set = paste0("_", add.args$data.set, "_")
  }
  
  ## 2.4 Check path : data.set ------------------------------------------------
  if (!is.null(add.args$path) && !(bm.out@sp.name %in% list.dirs(path = add.args$path))) {
    stop("invalid path given")
  } else {
    add.args$path = "."
  }
  
  ## 2.5 Check path : full.name ------------------------------------------------
  if (!is.null(add.args$full.name)) {
    .fun_testIfIn(TRUE, "add.args$full.name", add.args$full.name, avail_models)
  }
  
  return(list(add.args = add.args))
}


###################################################################################################
##' @name BIOMOD_PresenceOnly
##' @author Frank Breiner, Maya Gueguen
##' 
##' @title Evaluate models with presence-only metrics
##' 
##' @description This function computes presence-only evaluation metrics (Boyce index and Minimal 
##' Predicted Area) for \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' objects that can be obtained with the \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions.
##' 
##' 
##' @param bm.mod a \code{\link{BIOMOD.models.out}} object returned by the 
##' \code{\link{BIOMOD_Modeling}} function
##' @param bm.em a \code{\link{BIOMOD.ensemble.models.out}} object returned by the 
##' \code{\link{BIOMOD_EnsembleModeling}} function
##' @param bg.env (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing values of 
##' environmental variables (in columns or layers) extracted from the background 
##' (\emph{if presences are to be compared to background instead of absences or 
##' pseudo-absences selected for modeling})
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' 
##' @param perc a \code{numeric} between \code{0} and \code{1} corresponding to the percentage of 
##' correctly classified presences for Minimal Predicted Area (see 
##' \code{ecospat.mpa()} in \pkg{ecospat})
## \code{\link[ecospat]{ecospat.mpa}}) # generate R CMD Check error due 
## to crossref missing ecospat package
##' 
##' 
##' @return 
##' 
##' A \code{data.frame} containing evaluation scores both for the evaluation metrics used in the 
##' \code{\link{BIOMOD_Modeling}} function and additional Boyce index and Minimal Predicted Area.
##' 
##' 
##' @details
##' 
##' \code{em.by} parameter of \code{\link{BIOMOD_EnsembleModeling}} must have been set to 
##' \code{PA+run} in order to have an ensemble for each \emph{RUN} of the 
##' \code{NbRunEval} parameter of the \code{\link{BIOMOD_Modeling}} function for evaluation. 
##' \cr \cr
##' 
##' The Boyce index returns \code{NA} values for \code{SRE} models because it can not be 
##' calculated with binary predictions. \cr This is also the reason why some \code{NA} values 
##' might appear for \code{GLM} models if they do not converge.
##' 
##' 
##' @note In order to break dependency loop between packages \pkg{biomod2} and \pkg{ecospat}, 
##' code of \code{ecospat.boyce()} and \code{ecospat.mpa()} in \pkg{ecospat})
##' functions have been copied within this file from version 3.2.2 (august 2022).
## \code{\link[ecospat]{ecospat.mpa}}) # generate R CMD Check error due 
## to crossref missing ecospat package
## code of \code{\link[ecospat]{ecospat.boyce}} and \code{\link[ecospat]{ecospat.mpa}} 
## generate R CMD Check error due to crossref missing ecospat package 
##' 
##' 
##' @references
##' 
##' \itemize{
##'   \item Engler, R., Guisan, A., and Rechsteiner L. 2004. An improved approach for predicting 
##'   the distribution of rare and endangered species from occurrence and pseudo-absence data. 
##'   \emph{Journal of Applied Ecology}, \bold{41(2)}, 263-274.
##'   \item Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., and Guisan, A. 2006. Evaluating 
##'   the ability of habitat suitability models to predict species presences. \emph{Ecological 
##'   Modelling}, \bold{199(2)}, 142-152.
##' }
##' 
##' 
##' @seealso \code{ecospat.boyce()} and \code{ecospat.mpa()} in \pkg{ecospat}, 
## @seealso \code{\link[ecospat]{ecospat.boyce}}, \code{\link[ecospat]{ecospat.mpa}}, 
##' \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{BIOMOD.ensemble.models.out}}, \code{\link{BIOMOD_EnsembleModeling}}
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
##'                                         metric.eval = c('TSS', 'ROC'),
##'                                         var.import = 3,
##'                                         seed.val = 42)
##' }
##' 
##' 
##' # --------------------------------------------------------------- #
##' # # Evaluate models with Boyce index and MPA
##' # myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
##' #                                   bm.em = myBiomodEM)
##' # myBiomodPO
##' # 
##' # # Evaluate models with Boyce index and MPA (using background data)
##' # myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
##' #                                   bm.em = myBiomodEM, 
##' #                                   bg.env = myExpl)
##' # myBiomodPO
##' 
##' 
## @importFrom ecospat ecospat.boyce ecospat.mpa
##' @importFrom PresenceAbsence presence.absence.accuracy
##' @importFrom terra rast extract
##' 
##' @export
##' 
##' 
###################################################################################################

BIOMOD_PresenceOnly <- function(bm.mod = NULL, 
                                bm.em = NULL, 
                                bg.env = NULL, 
                                perc = 0.9)
{
  .bm_cat("Do Presence-Only Evaluation")
  # if (!isNamespaceLoaded("ecospat")) { requireNamespace("ecospat", quietly = TRUE) }
  
  ## 0. Check arguments --------------------------------------------------------
  args <- .BIOMOD_PresenceOnly.check.args(bm.mod, bm.em, bg.env, perc)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Get calib.lines ------------------------------------------------------
  if (!is.null(bm.mod)) {
    calib.lines <- get_calib_lines(bm.mod)#[, 1]
    bm.data <- get_formal_data(bm.mod)
    resp.var <- get_formal_data(bm.mod, subinfo = "resp.var")
  } else {
    calib.lines <- get_calib_lines(get_formal_data(bm.em))#[, 1]
    bm.data <- get_formal_data(get_formal_data(bm.em))
    resp.var <- get_formal_data(get_formal_data(bm.em), subinfo = "resp.var")
  }
  

  calib.notNA <- which(!is.na(resp.var))
  calib.lines <- calib.lines[calib.notNA, ]
  
  ## 2. Individual models ----------------------------------
  if (!is.null(bm.mod)) {
    ## Get calibration lines and observations
    myResp <- resp.var[calib.notNA] ## keep only lines associated to sites (no pseudo-absences)
    
    ## Get evaluation scores
    myEvalMod <- get_evaluations(bm.mod)
    
    ## Get predictions on observations
    myPredMod <- .get_list_predictions(bm.out = bm.mod, evaluation = FALSE)
    if (!is.null(bg.env)) {
      myPredMod.sites <- myPredMod
      myBiomodProj.eval <- BIOMOD_Projection(bm.mod = bm.mod,
                                             proj.name = paste(bm.mod@modeling.id, "cv_EF_eval", sep = "_"), 
                                             new.env = bg.env,
                                             build.clamping.mask = FALSE)
      myPredMod <- .get_list_predictions(bm.out = myBiomodProj.eval, evaluation = FALSE)
    }
    
    ## Get predictions on evaluation data
    if (bm.mod@has.evaluation.data == TRUE) {
      myPredMod.eval <- .get_list_predictions(bm.out = bm.mod, evaluation = TRUE)
    }
  } else {
    myEvalMod = myPredMod = myPredMod.eval = NULL
  }
  
  ## 2. Ensemble models ----------------------------------
  
  if (!is.null(bm.em)) {
    
    ## Get evaluation scores
    myEvalEM <- get_evaluations(bm.em)
    if (!is.null(bm.mod)) {
      full.colnames <- c("full.name", "PA", "run", "algo",
                         "metric.eval", "cutoff", 
                         "sensitivity", "specificity", "calibration", 
                         "validation", "evaluation",
                         "merged.by.PA", "merged.by.run", "merged.by.algo", "filtered.by")
      
      myEvalMod$merged.by.PA <- NA
      myEvalMod$merged.by.run <- NA
      myEvalMod$merged.by.algo <- NA
      myEvalMod$filtered.by <- NA
      myEvalEM$PA <- NA
      myEvalEM$run <- NA
      
      myEvalMod <- rbind(myEvalMod[,full.colnames],
                         myEvalEM[,full.colnames])
    } else {
      myEvalMod <- myEvalEM
    }
    myEvalMod = as.data.frame(myEvalMod)
    for (cc in c("full.name", "PA", "run", "algo", "metric.eval")) {
      if (cc %in% colnames(myEvalMod)) {
        myEvalMod[, cc] = as.character(myEvalMod[, cc])
      }
    }
    
    ## Get predictions on observations
    myPredEM <- .get_list_predictions(bm.out = bm.em, evaluation = FALSE)
    if (!is.null(bg.env)) {
      myPredEM.sites <- myPredEM
      
      if (exists("myPredMod.sites")) {
        myPredMod.sites <- c(myPredMod.sites, myPredEM.sites)
      } else {
        myPredMod.sites <- myPredEM.sites
      }
      
      if(exists("myBiomodProj.eval")){
        myPredEM <- BIOMOD_EnsembleForecasting(
          bm.em = bm.em,
          bm.proj = myBiomodProj.eval,
          proj.name = paste(bm.em@modeling.id, "cv_EF_bg", sep = "_"))
      } else {
        myPredEM <- BIOMOD_EnsembleForecasting(
          bm.em = bm.em,
          new.env = bg.env,
          proj.name = paste(bm.em@modeling.id, "cv_EF_bg", sep = "_"))
      }
      myPredEM <- .get_list_predictions(bm.out = myPredEM, evaluation = FALSE)
    }
    if (!is.null(bm.mod)) {
      myPredMod <- c(myPredMod, myPredEM)
    } else {
      myPredMod <- myPredEM
    }
    
    ## Get predictions on evaluation data
    if (!is.null(bm.mod) && bm.mod@has.evaluation.data == TRUE) {
      myPredEM.eval <- .get_list_predictions(bm.out = bm.em, evaluation = TRUE)
      myPredMod.eval <- c(myPredMod.eval, myPredEM.eval)
    }  
  }
  
  ## CALCULATE BOYCE & MPA VALUES ---------------------------------------------
  mpa.eval <- boyce.eval <- myEvalMod[!duplicated(myEvalMod$full.name), ]
  boyce.eval$metric.eval <- "BOYCE"
  boyce.eval[, c("cutoff", "sensitivity", "specificity", "validation")] <- NA
  mpa.eval$metric.eval <- "MPA"
  mpa.eval[, c("cutoff", "sensitivity", "specificity", "validation")] <- NA
  
  for (i in 1:nrow(boyce.eval)) {
    ## Get model informations
    full.name <- boyce.eval[i, 1]
    tmp <- strsplit(as.character(full.name), split = "_")[[1]]
    length.tmp <- length(tmp)
    # run <- tmp[c(grep("RUN", tmp), grep("allRun", tmp), grep("mergedRun", tmp))]
    run <- 
      paste0(
        "_",
        paste0(tmp[c(length.tmp - 2, length.tmp - 1)], 
               collapse = "_")
      )
    # message(run)
    ## Get evaluation lines
    if (length(run) == 0) {
      ind.eval = NULL
    } else {
      if (grepl(run, pattern = "allRun") | grepl(run, pattern = "mergedRun")) {
        ind.eval = 1:nrow(calib.lines) 
      } else {
          ind.eval = which(!calib.lines[, run])
      }
    }
    
    ## Get vectors of selected observations and predictions
    if (is.null(bg.env)) {
      if (length(ind.eval) == 0) { ## No background, no evaluation : all obs and pred
        Test <- myResp
        Pred <- myPredMod[[full.name]]
      } else { ## No background : only obs and pred on eval lines
        Test <- myResp[ind.eval]
        Pred <- myPredMod[[full.name]][ind.eval]
      }
      Pred.obs <- Pred[which(Test == 1)]
    } else {
      if (length(ind.eval) == 0) { ## No evaluation : only obs and pred on presence points
        ind.obs = which(myResp == 1)
      } else { ## Only obs and pred on presence points which are also eval lines
        ind.obs = intersect(ind.eval, which(myResp == 1))
      }
      Test <- c(myResp[ind.obs], rep(0, nrow(bg.env)))
      Pred <- c(myPredMod.sites[[full.name]][ind.obs], myPredMod[[full.name]])
      Pred.obs <- Pred[1:length(ind.obs)]
    }
    
    ind.notNA = which(!is.na(Pred))
    ind.b = which(boyce.eval[, 1] == full.name)
    ind.m = which(mpa.eval[, 1] == full.name)
    
    ## Compute Boyce and MPA values -------------------------------------------
    if (length(Pred) > 0 && length(ind.notNA) > 0) {
      ## Prepare table to compute evaluation scores
      DATA <- cbind(1:length(Pred), Test, Pred / 1000)
      DATA[is.na(DATA[, 2]), 2] <- 0
      DATA <- DATA[complete.cases(DATA), ]
      
      ## Boyce index
      boy <- ecospat.boyce(fit = Pred[ind.notNA], obs = Pred.obs, PEplot = FALSE)
      boyce.eval$validation[ind.b] <- boy$cor
      if (sum(boy$F.ratio < 1, na.rm = TRUE) > 0) {
        boyce.eval$cutoff[ind.b] <- round(boy$HS[max(which(boy$F.ratio < 1))], 0)
        if (!is.na(boyce.eval$cutoff[ind.b] / 1000)) {
          EVAL <- presence.absence.accuracy(DATA, threshold = boyce.eval$cutoff[ind.b] / 1000)
          boyce.eval$sensitivity[ind.b] <- EVAL$sensitivity
          boyce.eval$specificity[ind.b] <- EVAL$specificity
        } else {
          boyce.eval[ind.b, c("sensitivity", "specificity")] < - NA
        }
      } else {
        boyce.eval[ind.b, c("cutoff", "sensitivity", "specificity")] <- NA 	
      }
      
      ## MPA index
      mpa <- ecospat.mpa(Pred.obs, perc = perc)
      mpa.eval$cutoff[ind.m] <- mpa
      EVAL <- presence.absence.accuracy(DATA, threshold = mpa / 1000)
      mpa.eval$sensitivity[ind.m] <- EVAL$sensitivity
      mpa.eval$specificity[ind.m] <- EVAL$specificity  
    }
    
    ## Compute Boyce and MPA values for evaluation data -----------------------
    if (!is.null(bm.mod) && bm.mod@has.evaluation.data == TRUE) {
      myResp.eval <- get_formal_data(bm.mod)@eval.data.species
      Pred.eval <- myPredMod.eval[[full.name]]
      
      boy <- ecospat.boyce(fit = Pred.eval,
                           obs = Pred.eval[myResp.eval == 1],
                           PEplot = FALSE)
      boyce.eval[ind.b, "evaluation"] <- boy$cor
      # mpa.eval[ind.m,"evaluation"] <- ecospat.mpa(Pred.eval[myResp.eval == 1], perc = perc)
      mpa.eval[ind.m,"evaluation"] <- NA
    }
  }
  myEvalMod[, c("sensitivity", "specificity")] <- round(myEvalMod[, c("sensitivity", "specificity")], 1)
  boyce.eval[, c("sensitivity", "specificity")] <- round(boyce.eval[, c("sensitivity", "specificity")] * 100, 1)
  mpa.eval[, c("sensitivity", "specificity")] <- round(mpa.eval[, c("sensitivity", "specificity")] * 100, 1)
  
  ## SAVE OUTPUTS ---------------------------------------------------------------------------------
  output <- rbind(myEvalMod, boyce.eval, mpa.eval)
  
  if (!is.null(bm.mod)) {
    sp <- bm.mod@sp.name
    mod.id <- bm.mod@modeling.id
  } else if (!is.null(bm.em)) {
    sp <- bm.em@sp.name
    mod.id <- bm.em@modeling.id
  }
  save(output, file = paste0(sp, "/.BIOMOD_DATA/", mod.id, "/presenceonly.evaluation_", sp))
  
  .bm_cat("Done")
  return(output)
}


.get_list_predictions <- function(bm.out, evaluation = FALSE)
{
  myPred <- get_predictions(bm.out, evaluation = evaluation)
  myPred <- tapply(X = myPred$pred, INDEX = list(myPred$points, myPred$full.name), FUN = mean)
  # myPred <- as.data.frame(myPred)
  myPred.list <- lapply(1:ncol(myPred), function(x) myPred[, x])
  names(myPred.list) <- colnames(myPred)
  return(myPred.list)
}


# Check Arguments -------------------------------------------------------------

.BIOMOD_PresenceOnly.check.args <- function(bm.mod, bm.em, bg.env, perc)
{
  
  if(is.null(bm.mod) && is.null(bm.em)){
    stop("At least one of 'bm.mod' or 'bm.em' have to be given.")
  }
  
  ## 1. Check bm.mod ----------------------------------------------------------
  if (!is.null(bm.mod)) {
    .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
    expl.var.names <- bm.mod@expl.var.names
  }
  
  ## 2. Check bm.em -----------------------------------------------------------
  if (!is.null(bm.em)) {
    .fun_testIfInherits(TRUE, "bm.em", bm.em, "BIOMOD.ensemble.models.out")
    if(is.null(bm.mod)){
      expl.var.names <- get_formal_data(bm.em)@expl.var.names
    }
  }
  
  ## 2. Check bg.env -----------------------------------------------------------
  if(!is.null(bg.env)){
    available.types.resp <- c('numeric', 'data.frame', 'matrix',
                              'SpatialPointsDataFrame', 'SpatVector',
                              'Raster','SpatRaster')
    
    .fun_testIfInherits(TRUE, "bg.env", bg.env, available.types.resp)
    
    if(is.matrix(bg.env)) {
      bg.env <- as.data.frame(bg.env)
    }
    
    if(is.numeric(bg.env)) {
      bg.env <- as.data.frame(bg.env)
      colnames(bg.env) <- expl.var.names[1]
    }
    
    if (inherits(bg.env, 'Raster')) {
      if(any(raster::is.factor(bg.env))){
        bg.env <- .categorical_stack_to_terra(bg.env)
      } else {
        bg.env <- rast(bg.env)
      }
    }
    
    if (inherits(bg.env, 'SpatRaster')) {
      bg.env <- as.data.frame(bg.env)
    }
    
    if (inherits(bg.env, 'SpatialPoints')) {
      bg.env <- as.data.frame(bg.env@data)
    }
    
    if (inherits(bg.env, 'SpatVector')) {
      bg.env <- as.data.frame(bg.env)
    }
    
    ## remove NA from background data
    if (sum(is.na(bg.env)) > 0) {
      cat("\n      ! NAs have been automatically removed from bg.env data")
      bg.env <- na.omit(bg.env)
    }
    
  }
  
  ## 4. Check perc -----------------------------------------------------------
  .fun_testIf01(TRUE, "perc", perc)
  
  
  ## 5. Check MPA & Evaluation data ---------------------------------------------
  if ((!is.null(bm.mod) && bm.mod@has.evaluation.data) ||
      (is.null(bm.mod) && get_formal_data(bm.em)@has.evaluation.data)) {
    cat("\n      ! Evaluation data will be ignored for MPA-related calculations")
  }
  
  
  return(list(bm.mod = bm.mod,
              bm.em = bm.em,
              bg.env = bg.env,
              perc = perc))
}


###################################################################################################
## FROM ECOSPAT PACKAGE VERSION 3.2.2 (august 2022)

#### Calculating Boyce index as in Hirzel et al. 2006
# fit: A vector or Raster-Layer containing the predicted suitability values 
# obs: A vector containing the predicted suitability values or xy-coordinates 
# (if fit is a Raster-Layer) of the validation points (presence records)
# nclass : number of classes or vector with classes threshold. If nclass=0, Boyce index is 
# calculated with a moving window (see next parameters)
# windows.w : width of the moving window (by default 1/10 of the suitability range)
# res : resolution of the moving window (by default 100 focals)
# PEplot : if True, plot the predicted to expected ratio along the suitability class
# rm.duplicate : if TRUE, the correlation exclude successive duplicated values
# method : correlation method used to compute the boyce index


ecospat.boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                          PEplot = TRUE, rm.duplicate = TRUE, method = 'spearman')
{
  
  #### internal function calculating predicted-to-expected ratio for each class-interval
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    return(round(pi/ei,10))
  }
  
  # here, ecospat.boyce should not receive RasterLayer
  # if (inherits(fit,"RasterLayer")) {
  #   if (is.data.frame(obs) || is.matrix(obs)) {
  #     obs <- extract(fit, obs)
  #   }
  #   fit <- getValues(fit)
  #   fit <- fit[!is.na(fit)]
  # }
  
  mini <- min(fit,obs)
  maxi <- max(fit,obs)
  
  if(length(nclass)==1){
    if (nclass == 0) { #moving window
      if (window.w == "default") {window.w <- (max(fit) - min(fit))/10}
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else{ #window based on nb of class
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else{ #user defined window
    vec.mov <- c(mini, sort(nclass[!nclass>maxi|nclass<mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  #index to remove successive duplicates
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)  # calculation of the correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  if(length(nclass)==1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  }
  HS <- HS[to.keep]  #exclude the NaN
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}


###################################################################################################
## FROM ECOSPAT PACKAGE VERSION 3.2.2 (august 2022)

## This function calculates the Minimal Predicted Area.
## 
## R. Patin, Nov. 2022 : the function does not return minimal predicted area, but
## rather the quantile of suitability corresponding to perc% of presence
##  correctly predicted
##  
## FUNCTION'S ARGUMENTS
## Pred:      numeric or RasterLayer .predicted suitabilities from a SDM prediction
## Sp.occ.xy: xy-coordinates of the species (if Pred is a RasterLayer)
## perc:      Percentage of Sp.occ.xy that should be encompassed by the binary map.

## Details:

## Value
# Returns the Minimal Predicted Area

## Author(s)
# Frank Breiner

## References
# Engler, R., A. Guisan, and L. Rechsteiner. 2004. An improved approach for predicting the 
# distribution of rare and endangered species from occurrence and pseudo-absence data. 
# Journal of Applied Ecology 41:263-274.


ecospat.mpa <- function(Pred, Sp.occ.xy = NULL, perc = 0.9)
{
  perc <- 1 - perc
  if (!is.null(Sp.occ.xy)) {
    Pred <- extract(Pred, Sp.occ.xy)
  }
  round(quantile(Pred, probs = perc,na.rm = TRUE), 3)
}

### EXAMPLE

# obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
# [which(ecospat.testData$Saxifraga_oppositifolia==1)]) ecospat.mpa(obs) ecospat.mpa(obs,perc=1) ##
# 100% of the presences encompassed


## Example script for using Ensemble Small Models ESMs according to Lomba et al. 2010 Written by
## Frank Breiner Swiss Federal Research Institute WSL, July 2014.

## References: Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J. &
## Guisan, A. (2010). Overcoming the rare species modeling paradox: A novel hierarchical framework
## applied to an Iberian endemic plant. Biological Conservation, 143, 2647-2657.

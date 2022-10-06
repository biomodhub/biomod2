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
##' A \code{matrix}, \code{data.frame} or \code{\link[raster:stack]{RasterStack}} object 
##' containing values of environmental variables extracted from the background (\emph{if 
##' presences are to be compared to background instead of absences or pseudo-absences selected 
##' for modeling})
##' @param perc a \code{numeric} between \code{0} and \code{1} corresponding to the percentage of 
##' correctly classified presences for Minimal Predicted Area (see 
##' \code{ecospat.mpa()} in \pkg{ecospat})
## \code{\link[ecospat]{ecospat.mpa}}) # generate R CMD Check error due 
## to crossref missing ecospat package
##' 
##' @param save.output (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the output is to be saved within the 
##' \code{.BIOMOD_DATA} folder or not
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
##' \code{PA_dataset+repet} in order to have an ensemble for each \emph{RUN} of the 
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
##' # Evaluate models with Boyce index and MPA
##' myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
##'                                   bm.em = myBiomodEM)
##' myBiomodPO
##' 
##' # Evaluate models with Boyce index and MPA (using background data)
##' myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
##'                                   bm.em = myBiomodEM, 
##'                                   bg.env = myExpl)
##' myBiomodPO
##' 
##' 
## @importFrom ecospat ecospat.boyce ecospat.mpa
##' @importFrom PresenceAbsence presence.absence.accuracy
##' @importFrom data.table rbindlist
##' @importFrom raster extract getValues
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_PresenceOnly <- function(bm.mod = NULL, 
                                bm.em = NULL, 
                                bg.env = NULL, 
                                perc = 0.9, 
                                save.output = TRUE)
{
  .bm_cat("Do Presence-Only Evaluation")
  # if (!isNamespaceLoaded("ecospat")) { requireNamespace("ecospat") }
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_PresenceOnly.check.args(bm.mod, bm.em, bg.env, perc, save.output)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## MODELING OUTPUT ------------------------------------------------------------------------------
  if (!is.null(bm.mod)) {
    ## Get calibration lines and observations
    calib.lines <- get_calib_lines(bm.mod)[, , 1]
    calib.notNA <- which(!is.na(calib.lines[, 1])) ## remove NA (pseudo-absences) from run1
    calib.lines <- calib.lines[calib.notNA, ] ## keep only lines associated to sites (no pseudo-absences)
    myResp <- get_formal_data(bm.mod)@data.species
    myResp <- myResp[calib.notNA] ## keep only lines associated to sites (no pseudo-absences)
    
    ## Get evaluation scores
    myModelEval <- na.omit(get_evaluations(bm.mod, as.data.frame = TRUE))
    myModelEval$Model.name <- sapply(myModelEval$Model.name, function(x) 
      paste(c(bm.mod@sp.name, strsplit(as.character(x), split = "_")[[1]][3:1]), collapse = "_"))
    
    ## Get predictions on observations
    myModelPred <- get_predictions(bm.mod, as.data.frame = TRUE)
    if (!is.null(bg.env)) {
      myModelPred.sites <- as.data.frame(myModelPred)
      myBiomodProj.eval <- BIOMOD_Projection(bm.mod = bm.mod,
                                             proj.name = paste(bm.mod@modeling.id, "cv_EF_eval", sep = "_"), 
                                             new.env = bg.env,
                                             build.clamping.mask = FALSE)
      myModelPred <- get_predictions(myBiomodProj.eval, as.data.frame = TRUE)
    }
    
    ## Get predictions on evaluation data
    if (bm.mod@has.evaluation.data == TRUE) {
      # myModelPred.eval <- as.data.frame(get(load(paste0(bm.mod@", sp.name", "/.BIOMOD_DATA/"
      #                                                   , bm.mod@modeling.id
      #                                                   , "/models.prediction.eval"))))
      # colnames(myModelPred.eval) <- sapply(colnames(myModelPred.eval), function(x) 
      #   paste(c(bm.mod@sp.name, strsplit(as.character(x), split = "[.]")[[1]][3:1]), collapse = "_"))
      myModelPred.eval <- get_predictions(bm.mod, as.data.frame = TRUE, evaluation = TRUE)
    }
  } else {
    myModelEval = myModelPred = myModelPred.eval = NULL
  }
  
  ## ENSEMBLE MODELING OUTPUT ---------------------------------------------------------------------
  if (!is.null(bm.em)) {
    
    ## Get evaluation scores
    myModelEvalEF <- na.omit(get_evaluations(bm.em, as.data.frame = TRUE))
    myModelEvalEF$Model.name <- paste(bm.em@sp.name, as.character(myModelEvalEF$Model.name), sep = "_")
    if (!is.null(bm.mod)) {
      myModelEval <- rbindlist(list(myModelEval, myModelEvalEF), fill = TRUE)
    } else {
      myModelEval <- myModelEvalEF
    }
    myModelEval = as.data.frame(myModelEval)
    for (cc in c("Model.name", "Algo", "Run", "Dataset", "Model")) {
      if (cc %in% colnames(myModelEval)) {
        myModelEval[, cc] = as.character(myModelEval[, cc])
      }
    }
    
    ## Get predictions on observations
    myBiomodProjFF <- get_predictions(bm.em, as.data.frame = TRUE)  
    if (!is.null(bg.env)) {
      myBiomodProjFF.sites <- as.data.frame(myBiomodProjFF)
      myModelPred.sites <- cbind(myModelPred.sites, myBiomodProjFF.sites)
      myBiomodProjFF <- BIOMOD_EnsembleForecasting(bm.em = bm.em,
                                                   bm.proj = myBiomodProj.eval,
                                                   proj.name = paste(bm.em@modeling.id, "cv_EF_bg", sep = "_"))
      myBiomodProjFF <- get_predictions(myBiomodProjFF, as.data.frame = TRUE)
    }
    if (!is.null(bm.mod)) {
      myModelPred <- cbind(myModelPred, myBiomodProjFF)
    } else {
      myModelPred <- myBiomodProjFF
    }
    
    ## Get predictions on evaluation data
    if (!is.null(bm.mod) && bm.mod@has.evaluation.data == TRUE) {
      myBiomodProjFF.eval <- get_predictions(bm.em, as.data.frame = TRUE, evaluation = TRUE)
      myModelPred.eval <- cbind(myModelPred.eval, myBiomodProjFF.eval)      
    }  
  }
  
  ## CALCULATE BOYCE & MPA VALUES -----------------------------------------------------------------
  mpa.eval <- boyce.eval <- myModelEval[!duplicated(myModelEval$Model.name), ]
  boyce.eval$Eval.metric <- "BOYCE"
  boyce.eval[, c("Testing.data", "Cutoff", "Sensitivity", "Specificity")] <- NA
  mpa.eval$Eval.metric <- "MPA"
  mpa.eval[, c("Testing.data", "Cutoff", "Sensitivity", "Specificity")] <- NA
  
  for (i in 1:nrow(boyce.eval))
  {
    ## Get model informations
    Model.name <- boyce.eval[i,1]
    tmp <- strsplit(as.character(Model.name), split = "_")[[1]]
    n <- length(tmp)
    tec <- paste(tmp[3:n], collapse = "_") 
    run <- tmp[c(grep("RUN", tmp), grep("Full", tmp), grep("mergedRun", tmp))]
    
    ## Get evaluation lines
    if (length(run) == 0) {
      ind.eval = NULL
    } else {
      if (run == "mergedRun") {
        ind.eval = 1:nrow(calib.lines) 
      } else {
        if (inherits(calib.lines, "matrix")) {
          ind.eval = which(calib.lines[, paste0("_", run)] == FALSE)
        } else {
          ind.eval = which(calib.lines == FALSE)
        }
      }
    }
    
    ## Get vectors of selected observations and predictions
    if (is.null(bg.env)) {
      if (length(ind.eval) == 0) { ## No background, no evaluation : all obs and pred
        Test <- myResp
        Pred <- myModelPred[, Model.name]
      } else { ## No background : only obs and pred on eval lines
        Test <- myResp[ind.eval]
        Pred <- myModelPred[ind.eval, Model.name]
      }
      Pred.obs <- Pred[which(Test == 1)]
    } else {
      if (length(ind.eval) == 0) { ## No evaluation : only obs and pred on presence points
        ind.obs = which(myResp == 1)
      } else { ## Only obs and pred on presence points which are also eval lines
        ind.obs = intersect(ind.eval, which(myResp == 1))
      }
      Test <- c(myResp[ind.obs], rep(0, nrow(bg.env)))
      Pred <- c(myModelPred.sites[ind.obs, Model.name], myModelPred[, Model.name])
      Pred.obs <- Pred[1:length(ind.obs)]
    }
    
    ind.notNA = which(!is.na(Pred))
    ind.b = which(boyce.eval[, 1] == Model.name)
    ind.m = which(mpa.eval[, 1] == Model.name)
    
    ## Compute Boyce and MPA values -------------------------------------------
    if (length(Pred) > 0 && length(ind.notNA) > 0) {
      ## Prepare table to compute evaluation scores
      DATA <- cbind(1:length(Pred), Test, Pred / 1000)
      DATA[is.na(DATA[, 2]), 2] <- 0
      DATA <- DATA[complete.cases(DATA), ]
      
      ## Boyce index
      boy <- ecospat.boyce(fit = Pred[ind.notNA], obs = Pred.obs, PEplot = FALSE)
      boyce.eval$Testing.data[ind.b] <- boy$cor
      if (sum(boy$F.ratio < 1, na.rm = TRUE) > 0) {
        boyce.eval$Cutoff[ind.b] <- round(boy$HS[max(which(boy$F.ratio < 1))], 0)
        if (!is.na(boyce.eval$Cutoff[ind.b] / 1000)) {
          EVAL <- presence.absence.accuracy(DATA, threshold = boyce.eval$Cutoff[ind.b] / 1000)
          boyce.eval$Sensitivity[ind.b] <- EVAL$sensitivity
          boyce.eval$Specificity[ind.b] <- EVAL$specificity
        } else {
          boyce.eval[ind.b, c("Sensitivity", "Specificity")] < - NA
        }
      } else {
        boyce.eval[ind.b, c("Cutoff", "Sensitivity", "Specificity")] <- NA 	
      }
      
      ## MPA index
      mpa <- ecospat.mpa(Pred.obs, perc = perc)
      mpa.eval$Cutoff[ind.m] <- mpa
      EVAL <- presence.absence.accuracy(DATA, threshold = mpa / 1000)
      mpa.eval$Sensitivity[ind.m] <- EVAL$sensitivity
      mpa.eval$Specificity[ind.m] <- EVAL$specificity  
    }
    
    ## Compute Boyce and MPA values for evaluation data -----------------------
    if (!is.null(bm.mod) && bm.mod@has.evaluation.data == TRUE) {
      myResp.eval <- get_formal_data(bm.mod)@eval.data.species
      Pred.eval <- myModelPred.eval[, Model.name]
      
      boy <- ecospat.boyce(fit = Pred.eval, obs = Pred.eval[myResp.eval == 1 & ind.1], PEplot = FALSE)
      boyce.eval[ind.b, "Evaluating.data"] <- boy$cor
      mpa.eval[ind.m,"Evaluating.data"] <- ecospat.mpa(Pred.eval[myResp.eval == 1 & ind.1], perc = perc)
    }
  }
  myModelEval[, c("Sensitivity", "Specificity")] <- round(myModelEval[, c("Sensitivity", "Specificity")], 1)
  boyce.eval[, c("Sensitivity", "Specificity")] <- round(boyce.eval[, c("Sensitivity", "Specificity")] * 100, 1)
  mpa.eval[, c("Sensitivity", "Specificity")] <- round(mpa.eval[, c("Sensitivity", "Specificity")] * 100, 1)
  
  ## SAVE OUTPUTS ---------------------------------------------------------------------------------
  output <- rbind(myModelEval, boyce.eval, mpa.eval)
  
  if (save.output) {
    if (!is.null(bm.mod)) {
      sp <- bm.mod@sp.name
      mod.id <- bm.mod@modeling.id
    } else if (!is.null(bm.em)) {
      sp <- bm.em@sp.name
      mod.id <- bm.em@modeling.id
    }
    save(output, file = paste0(sp, "/.BIOMOD_DATA/", mod.id, "/presenceonly.evaluation_", sp))
  }
  
  .bm_cat("Done")
  return(output)
}



###################################################################################################

.BIOMOD_PresenceOnly.check.args <- function(bm.mod, bm.em, bg.env, perc, save.output)
{
  ## 1. Check bm.mod ----------------------------------------------------------
  if (!is.null(bm.mod)) {
    .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
  }
  
  ## 2. Check bm.em -----------------------------------------------------------
  if (!is.null(bm.em)) {
    .fun_testIfInherits(TRUE, "bm.em", bm.em, "BIOMOD.ensemble.models.out")
  }
  
  ## 2. Check bg.env -----------------------------------------------------------
  if(is.matrix(bg.env) | is.numeric(bg.env)) {
    bg.env <- as.data.frame(bg.env)
  }
  
  if (inherits(bg.env, 'Raster')) {
    bg.env <- as.data.frame(getValues(bg.env))
  }
  
  if (inherits(bg.env, 'SpatialPoints')) {
    bg.env <- as.data.frame(bg.env@data)
  }
  
  ## remove NA from background data
  if (sum(is.na(bg.env)) > 0) {
    cat("\n      ! NAs have been automatically removed from bg.env data")
    bg.env <- na.omit(bg.env)
  }
  
  ## 4. Check perc -----------------------------------------------------------
  .fun_testIf01(TRUE, "perc", perc)
  
  
  return(list(bm.mod = bm.mod,
              bm.em = bm.em,
              bg.env = bg.env,
              perc = perc,
              save.output = save.output))
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
  
  if (inherits(fit,"RasterLayer")) {
    if (is.data.frame(obs) || is.matrix(obs)) {
      obs <- extract(fit, obs)
    }
    fit <- getValues(fit)
    fit <- fit[!is.na(fit)]
  }
  
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

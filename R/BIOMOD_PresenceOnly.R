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
##' A \code{matrix} or \code{data.frame} object containing values of environmental variables
##' extracted from the background (\emph{if presences are to be compared to background instead of
##' absences or pseudo-absences selected for modeling})
##' @param perc a \code{numeric} between \code{0} and \code{1} corresponding to the percentage of
##' correctly classified presences for Minimal Predicted Area (see
##' \code{\link[ecospat]{ecospat.mpa}})
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
##' @seealso \code{\link[ecospat]{ecospat.boyce}}, \code{\link[ecospat]{ecospat.mpa}},
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
##'
##' # ---------------------------------------------------------------
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##'
##' # Create default modeling options
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##'
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     bm.options = myBiomodOptions,
##'                                     nb.rep = 2,
##'                                     data.split.perc = 80,
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 3,
##'                                     do.full.models = FALSE)
##'
##' # Model ensemble models
##' myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
##'                                       models.chosen = 'all',
##'                                       em.by = 'all',
##'                                       metric.select = c('TSS'),
##'                                       metric.select.thresh = c(0.7),
##'                                       metric.eval = c('TSS', 'ROC'),
##'                                       var.import = 3,
##'                                       prob.mean = TRUE,
##'                                       prob.median = TRUE,
##'                                       prob.cv = TRUE,
##'                                       prob.ci = TRUE,
##'                                       prob.ci.alpha = 0.05,
##'                                       committee.averaging = TRUE,
##'                                       prob.mean.weight = TRUE,
##'                                       prob.mean.weight.decay = 'proportional')
##'
##'
##' # ---------------------------------------------------------------
##' # Evaluate models with Boyce index and MPA
##' myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
##'                                   bm.em = myBiomodEM)
##' myBiomodPO$eval
##'
##' # Evaluate models with Boyce index and MPA (using background data)
##' myBiomodPO <- BIOMOD_PresenceOnly(bm.mod = myBiomodModelOut,
##'                                   bm.em = myBiomodEM,
##'                                   bg.env = getValues(myExpl))
##' myBiomodPO$eval
##'
##'
##' @importFrom ecospat ecospat.boyce ecospat.mpa
##' @importFrom PresenceAbsence presence.absence.accuracy
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
  if (!isNamespaceLoaded("ecospat")) { requireNamespace("ecospat") }

  ## MODELING OUTPUT ------------------------------------------------------------------------------
  if (!is.null(bm.mod)) {
    ## Get calibration lines and observations
    calib.lines <- get(load(bm.mod@calib.lines@link))[, , 1]
    calib.notNA <- which(!is.na(calib.lines[, 1])) ## remove NA (pseudo-absences) from run1
    calib.lines <- calib.lines[calib.notNA, ] ## keep only lines associated to sites (no pseudo-absences)
    myResp <- get(load(bm.mod@formated.input.data@link))@data.species
    myResp <- myResp[calib.notNA] ## keep only lines associated to sites (no pseudo-absences)

    ## Get evaluation scores
    myModelEval <- get_evaluations(bm.mod, as.data.frame = TRUE)
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
      myModelPred <- as.data.frame(myBiomodProj.eval@proj.out@val)
      colnames(myModelPred) <- paste(
        bm.mod@sp.name,
        rep(dimnames(myBiomodProj.eval@proj.out@val)[[4]], prod(dim(myBiomodProj.eval@proj.out@val)[2:3])),
        rep(dimnames(myBiomodProj.eval@proj.out@val)[[3]], each = dim(myBiomodProj.eval@proj.out@val)[2]),
        rep(dimnames(myBiomodProj.eval@proj.out@val)[[2]], dim(myBiomodProj.eval@proj.out@val)[3])
        , sep = "_")
    }

    ## Get predictions on evaluation data
    if (bm.mod@has.evaluation.data == TRUE) {
      myModelPred.eval  <- as.data.frame(get(load(paste0(bm.mod@"sp.name", "/.BIOMOD_DATA/"
                                                         , bm.mod@modeling.id
                                                         , "/models.prediction.eval"))))

      colnames(myModelPred.eval) <- sapply(colnames(myModelPred.eval), function(x)
        paste(c(bm.mod@sp.name, strsplit(as.character(x), split = "[.]")[[1]][3:1]), collapse = "_"))
    }
  } else {
    myModelEval = myModelPred = myModelPred.eval = NULL
  }

  ## ENSEMBLE MODELING OUTPUT ---------------------------------------------------------------------
  if (!is.null(bm.em)) {

    ## Get evaluation scores
    myModelEvalEF <- get_evaluations(bm.em, as.data.frame = TRUE)
    myModelEvalEF$Model.name <- paste(bm.mod@sp.name, as.character(myModelEvalEF$Model.name), sep = "_")
    if (!is.null(bm.mod)) {
      myModelEval <- rbindlist(list(myModelEval, myModelEvalEF), fill = TRUE)
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
                                                   proj.name = paste(bm.mod@modeling.id, "cv_EF_bg", sep = "_"))
      myBiomodProjFF <- as.data.frame(myBiomodProjFF@proj.out@val)
    }
    if (!is.null(bm.mod)) { myModelPred <- cbind(myModelPred, myBiomodProjFF) }

    ## Get predictions on evaluation data
    if (bm.mod@has.evaluation.data == TRUE) {
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
    run <- tmp[c(grep("RUN", tmp), grep("Full", tmp))]

    ## Get evaluation lines
    if (length(run) == 0) {
      ind.eval = NULL
    } else {
      if (inherits(calib.lines, "matrix")) {
        ind.eval = which(calib.lines[, paste0("_", run)] == FALSE)
      } else {
        ind.eval = which(calib.lines == FALSE)
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
      Pred.obs <- Pred[which(Test > 0)]
    } else {
      if (length(ind.eval) == 0) { ## No evaluation : only obs and pred on presence points
        ind.obs = which(myResp > 0)
      } else { ## Only obs and pred on presence points which are also eval lines
        ind.obs = intersect(ind.eval, which(myResp > 0))
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
      boy <- ecospat::ecospat.boyce(fit = Pred[ind.notNA], obs = Pred.obs, PEplot = FALSE)
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
      mpa <- ecospat::ecospat.mpa(Pred.obs, perc = perc)
      mpa.eval$Cutoff[ind.m] <- mpa
      EVAL <- presence.absence.accuracy(DATA, threshold = mpa / 1000)
      mpa.eval$Sensitivity[ind.m] <- EVAL$sensitivity
      mpa.eval$Specificity[ind.m] <- EVAL$specificity
    }

    ## Compute Boyce and MPA values for evaluation data -----------------------
    if (bm.mod@has.evaluation.data == TRUE) {
      myResp.eval <- get(load(bm.mod@formated.input.data@link))@eval.data.species
      Pred.eval <- myModelPred.eval[, Model.name]

      boy <- ecospat::ecospat.boyce(fit = Pred.eval, obs = Pred.eval[myResp.eval > 0 & ind.1 == TRUE], PEplot = FALSE)
      boyce.eval[ind.b, "Evaluating.data"] <- boy$cor
      mpa.eval[ind.m,"Evaluating.data"] <- ecospat::ecospat.mpa(Pred.eval[myResp.eval > 0 & ind.1 == TRUE], perc = perc)
    }
  }
  myModelEval[, c("Sensitivity", "Specificity")] <- round(myModelEval[, c("Sensitivity", "Specificity")], 1)
  boyce.eval[, c("Sensitivity", "Specificity")] <- round(boyce.eval[, c("Sensitivity", "Specificity")] * 100, 1)
  mpa.eval[, c("Sensitivity", "Specificity")] <- round(mpa.eval[, c("Sensitivity", "Specificity")] * 100, 1)

  ## SAVE OUTPUTS ---------------------------------------------------------------------------------
  output <- rbind(myModelEval, boyce.eval, mpa.eval)

  if (save.output) {
    if (!is.null(bm.mod)) { sp <- bm.mod@sp.name }
    if (!is.null(bm.em)) { sp <- bm.em@sp.name }
    save(output, file = paste0(sp, "/.BIOMOD_DATA/", bm.mod@modeling.id, "/presenceonly.evaluation_", sp))
  }

  .bm_cat("Done")
  return(output)
}

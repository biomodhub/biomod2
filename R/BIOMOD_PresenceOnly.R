##' @name BIOMOD_PresenceOnly
##' @aliases BIOMOD_PresenceOnly
##' 
##' @title evaluate models with presences only metrics
##' 
##' @description This function enables to evaluate BIOMOD.models.out and 
##' BIOMOD.EnsembleModeling.out object with presence-only evaluation methods 
##' (Boyce index and Minimal Predicted Area MPA)
##' 
##' @param modeling.output  "BIOMOD.models.out" object produced by a BIOMOD_Modeling run
##' @param EM.output        a "BIOMOD.EnsembleModeling.out" returned by BIOMOD_EnsembleModeling
##' @param bg.env           a data frame or matrix of environmental variables which was extracted from the background 
##' (might be used if presences should be compared to the background instead of Absences or Pseudo-Absences selected for modelling). 
##' @param perc             Percentage of correctly classified presences for MPA (Default 90\%).  
##' @param save.output      logical. If TRUE (Default) the output is saved to the ".BIOMOD_DATA" folder
##' 
##' @details
##' 'em.by' of 'BIOMOD.EnsembleModeling' must be 'PA_dataset+repet' to have an 
##' ensemble for each RUN of the 'NbRunEval' argument (BIOMOD_Modeling funtion)
##' for evaluation.
##' The Boyce index returns NA values for 'SRE' models because it is not possible
##' to be calculated with binary predictions. This is also the reason why there 
##' are sometimes NA values for 'GLM' models if they don not converge.
##' 
##' @return
##' data.frame containing evaluation scores for the evaluation metrics used for 
##' the BIOMOD_Modeling function and additional Boyce index and MPA
##' 
##' @references
##' Engler, R., Guisan, A., and Rechsteiner L. 2004. An improved approach for predicting the distribution of rare and endangered species from occurrence and pseudo-absence data. Journal of Applied Ecology, 41(2), 263-274.
##' Hirzel, A. H., Le Lay, G., Helfer, V., Randin, C., and Guisan, A. 2006. Evaluating the ability of habitat suitability models to predict species presences. Ecological Modelling, 199(2), 142-152.
##' 
##' @author 
##' Frank Breiner \email{frank.breiner@wsl.ch}
##' 
##' @seealso 
##' \code{\link[ecospat]{ecospat.boyce}}, \code{\link[ecospat]{ecospat.mpa}}, \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{BIOMOD_EnsembleModeling}}
##' 
##' @examples
##' \dontrun{
##' requireNamespace(PresenceAbsence, 'PresenceAbsence', quietly = TRUE)
##' 
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"), row.names = 1)
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
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
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('SRE','CTA','RF'), 
##'                                      models.options = myBiomodOption, 
##'                                      NbRunEval=1, 
##'                                      DataSplit=80, 
##'                                      Yweights=NULL, 
##'                                      VarImport=3, 
##'                                      models.eval.meth = c('TSS','ROC'),
##'                                      SaveObj = TRUE,
##'                                      rescal.all.models = FALSE,
##'                                      do.full.models = FALSE)
##' 
##' # 4. Doing Ensemble Modelling
##' myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
##'                                        chosen.models = 'all',
##'                                        em.by = 'PA_dataset+repet',
##'                                        eval.metric = c('TSS'),
##'                                        eval.metric.quality.threshold = c(0.7),
##'                                        models.eval.meth = c('TSS','ROC'),
##'                                        prob.mean = TRUE,
##'                                        prob.cv = FALSE,
##'                                        prob.ci = FALSE,
##'                                        prob.ci.alpha = 0.05,
##'                                        prob.median = FALSE,
##'                                        committee.averaging = FALSE,
##'                                        prob.mean.weight = TRUE,
##'                                        prob.mean.weight.decay = 'proportional' )   
##' 
##' # evaluate Biomod models with the Boyce index and MPA
##' pres.only.eval <- BIOMOD_PresenceOnly(myBiomodModelOut, myBiomodEM)
##' pres.only.eval$eval
##'
##' # evaluate Biomod models with the Boyce index and MPA using Background data
##' bg.Values <- getValues(myExpl)
##' 
##' pres.only.eval <- BIOMOD_PresenceOnly(myBiomodModelOut, myBiomodEM, bg.env = bg.Values)
##' pres.only.eval$eval
##' }

BIOMOD_PresenceOnly <- function(modeling.output = NULL, EM.output = NULL, bg.env = NULL, perc = 0.9, save.output = T)
{
  requireNamespace('PresenceAbsence', quietly = TRUE)
  
  myModelEval <- myBiomodProjFF <- NULL
  
  ## MODELING OUTPUT ------------------------------------------------------------------------------
  if (!is.null(modeling.output)) {
    ## Get calibration lines and observations
    calib.lines <- get(load(modeling.output@calib.lines@link))[, , 1]
    calib.notNA <- which(!is.na(calib.lines[, 1])) ## remove NA (pseudo-absences) from run1
    calib.lines <- calib.lines[calib.notNA, ] ## keep only lines associated to sites (no pseudo-absences)
    myResp <- get(load(modeling.output@formated.input.data@link))@data.species
    myResp <- myResp[calib.notNA] ## keep only lines associated to sites (no pseudo-absences)
    
    ## Get evaluation scores
    myModelEval <- get_evaluations(modeling.output, as.data.frame = T)
    myModelEval[,1] <- as.character(myModelEval[, 1])
    for (i in 1:nrow(myModelEval)) {
      myModelEval[i, 1] < - paste(c(modeling.output@sp.name
                                    , strsplit(as.character(myModelEval[i,1]), split = "_")[[1]][3:1])
                                  , collapse = "_")
    }
    
    ## Get predictions on observations
    myModelPred <- get_predictions(modeling.output, as.data.frame = T)
    if (!is.null(bg.env)) {
      myModelPred.sites <- as.data.frame(myModelPred)
      myBiomodProj.eval <- BIOMOD_Projection(
        new.env = bg.env,
        proj.name = paste(modeling.output@modeling.id, "cv_EF_eval", sep = "_"), 
        modeling.output = modeling.output,
        build.clamping.mask = F)
      myModelPred <- as.data.frame(myBiomodProj.eval@proj@val)
      colnames(myModelPred) <- paste(
        modeling.output@sp.name,
        rep(dimnames(myBiomodProj.eval@proj@val)[[4]], prod(dim(myBiomodProj.eval@proj@val)[2:3])),
        rep(dimnames(myBiomodProj.eval@proj@val)[[3]], each = dim(myBiomodProj.eval@proj@val)[2]),
        rep(dimnames(myBiomodProj.eval@proj@val)[[2]], dim(myBiomodProj.eval@proj@val)[3])
        , sep = "_")
    }
    
    ## Get predictions on evaluation data
    if (modeling.output@has.evaluation.data == T) {
      myModelPred.eval  <- as.data.frame(get(load(paste0(modeling.output@"sp.name", "/.BIOMOD_DATA/"
                                                         , modeling.output@modeling.id
                                                         , "/models.prediction.eval"))))
      for (i in 1:ncol(myModelPred.eval)) {
        colnames(myModelPred.eval)[i] <- paste(c(modeling.output@sp.name
                                                 , strsplit(colnames(myModelPred.eval)[i], split="[.]")[[1]][3:1])
                                               , collapse="_")
      }       
    }
  }
  
  ## ENSEMBLE MODELING OUTPUT ---------------------------------------------------------------------
  if (!is.null(EM.output)) {
    if (EM.output@em.by != 'PA_dataset+repet') { stop("em.by of 'BIOMOD.EnsembleModeling' must be 'PA_dataset+repet'") }
    
    ## Get evaluation scores
    myModelEvalEF <- get_evaluations(EM.output, as.data.frame = T)
    myModelEvalEF[, 1] <- paste(modeling.output@sp.name, as.character(myModelEvalEF[, 1]), sep = "_")
    if (!is.null(modeling.output)) { myModelEval <- rbind(myModelEval, myModelEvalEF) }
    
    ## Get predictions on observations
    myBiomodProjFF <- get_predictions(EM.output, as.data.frame = T)  
    if (!is.null(bg.env)) {
      myBiomodProjFF.sites <- as.data.frame(myBiomodProjFF)
      myModelPred.sites <- cbind(myModelPred.sites, myBiomodProjFF.sites)
      myBiomodProjFF <- BIOMOD_EnsembleForecasting(
        proj.name = paste(modeling.output@modeling.id, "cv_EF_bg", sep = "_"), 
        projection.output = myBiomodProj.eval,
        EM.output = EM.output)    
      myBiomodProjFF <- as.data.frame(myBiomodProjFF@proj@val)     
    }
    if (!is.null(modeling.output)) { myModelPred <- cbind(myModelPred, myBiomodProjFF) }
    
    ## Get predictions on evaluation data
    if (modeling.output@has.evaluation.data == T) {
      myBiomodProjFF.eval <- get_predictions(EM.output, as.data.frame = T, evaluation = T)  
      myModelPred.eval <- cbind(myModelPred.eval, myBiomodProjFF.eval)      
    }  
  }
  
  ## CALCULATE BOYCE & MPA VALUES -----------------------------------------------------------------
  mpa.eval <- boyce.eval <- myModelEval[!duplicated(myModelEval[, 1]), ]
  boyce.eval$Eval.metric <- "boyce"
  mpa.eval$Eval.metric <- "mpa"
  boyce.eval[, 3:7] <- mpa.eval[, 3:7] <- NA
  
  for (i in 1:nrow(boyce.eval))
  {
    ## Get model informations
    Model.name <- boyce.eval[i,1]
    tmp <- strsplit(as.character(Model.name), split = "_")[[1]]
    n <- length(tmp)
    tec <- paste(tmp[3:n], collapse = "_") 
    run <- tmp[c(grep("RUN", tmp), grep("Full", tmp))]
    
    ## Get evaluation lines
    if (inherits(calib.lines, "matrix")) {
      ind.eval = which(calib.lines[, paste0("_", run)] == FALSE)
    } else {
      ind.eval = which(calib.lines == FALSE)
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
    ind.b = which(boyce.eval[,1]==Model.name)
    ind.m = which(mpa.eval[,1]==Model.name)
    
    ## Compute Boyce and MPA values -------------------------------------------
    if (length(Pred) > 0) {
      ## Prepare table to compute evaluation scores
      DATA <- cbind(1:length(Pred), Test, Pred / 1000)
      DATA[is.na(DATA[, 2]), 2] <- 0
      DATA <- DATA[stats::complete.cases(DATA), ]
      
      ## Boyce index
      boy <- ecospat::ecospat.boyce(fit = Pred[ind.notNA], obs = Pred.obs, PEplot = F)
      boyce.eval[ind.b, 3] <- boy$Spearman.cor
      if (sum(boy$F.ratio < 1, na.rm = T) > 0) {
        boyce.eval[ind.b, 5] <- round(boy$HS[max(which(boy$F.ratio < 1))], 0)
        if (!is.na(boyce.eval[ind.b, 5] / 1000)) {
          EVAL <- presence.absence.accuracy(DATA, threshold = boyce.eval[ind.b, 5] / 1000)
          boyce.eval[ind.b, 6] <- EVAL$sensitivity
          boyce.eval[ind.b, 7] <- EVAL$specificity
        } else {
          boyce.eval[ind.b, 6:7] < - NA
        }
      } else {
        boyce.eval[ind.b, 7] <- boyce.eval[ind.b, 6] <- boyce.eval[ind.b, 5] <- NA 	
      }
      
      ## MPA index
      mpa <- ecospat::ecospat.mpa(Pred.obs, perc = perc)
      mpa.eval[ind.m, 5] <- mpa
      EVAL <- presence.absence.accuracy(DATA, threshold = mpa / 1000)
      mpa.eval[ind.m, 6] <- EVAL$sensitivity
      mpa.eval[ind.m, 7] <- EVAL$specificity  
    }
    
    ## Compute Boyce and MPA values for evaluation data -----------------------
    if (modeling.output@has.evaluation.data == T) {
      myResp.eval <- get(load(modeling.output@formated.input.data@link))@eval.data.species
      Pred.eval <- myModelPred.eval[, Model.name]
      
      boy <- ecospat::ecospat.boyce(fit = Pred.eval, obs = Pred.eval[myResp.eval == 1 & ind.1], PEplot=F)
      boyce.eval[ind.b, "Evaluating.data"] <- boy$Spearman.cor
      mpa.eval[ind.m,"Evaluating.data"] <- ecospat::ecospat.mpa(Pred.eval[myResp.eval == 1 & ind.1], perc = perc)
    }
  }
  myModelEval[, 6:7] <- round(myModelEval[, 6:7], 1)
  boyce.eval[, 6:7] <- round(boyce.eval[, 6:7] * 100, 1)
  mpa.eval[, 6:7] <- round(mpa.eval[, 6:7] * 100, 1)
  
  ## SAVE OUTPUTS ---------------------------------------------------------------------------------
  if (!is.null(EM.output)) {
    if (modeling.output@has.evaluation.data == T) {
      output <- list(eval = rbind(myModelEval, boyce.eval,mpa.eval)
                     , myBiomodProjFF = myBiomodProjFF
                     , myBiomodProjEF.eval = myBiomodProjFF.eval) 
    } else {
      output <- list(eval = rbind(myModelEval, boyce.eval, mpa.eval)
                     , myBiomodProjFF = myBiomodProjFF)
    }
  } else {
    output <- list(eval = rbind(myModelEval, boyce.eval, mpa.eval))      
  }
  if (save.output) {
    if (!is.null(modeling.output)) { sp <- modeling.output@sp.name }
    if (!is.null(EM.output)) { sp <- EM.output@sp.name }
    save(output, file = paste0(sp, "/.BIOMOD_DATA/", modeling.output@modeling.id, "/presenceonly.evaluation_", sp))
  }
  return(output)
}

##' ###############################################################################################
##' @name BIOMOD_Tuning
##' @aliases BIOMOD_Tuning
##' @author Frank Breiner \email{frank.breiner@wsl.ch}
##' 
##' @title Tune models parameters
##' 
##' @description Function to tune \pkg{biomod2} single models parameters
##'
##' @param data a \code{\link{BIOMOD.formated.data} object returned by the 
##' \code{\link{BIOMOD_FormatingData}} function
##' @param models a \code{vector} containing model names to be tuned, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT.Phillips}
##' @param models.options a \code{\link{BIOMOD.models.options} object returned by the 
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param trControl global control parameters that can be obtained from 
##' \code{\link[caret]{trainControl}} function
##' @param ctrl.GAM control parameters for \code{GAM}
##' @param ctrl.GBM control parameters for \code{GBM}
##' @param ctrl.GLM control parameters for \code{GLM}
##' @param ctrl.CTA control parameters for \code{CTA}
##' @param ctrl.RF control parameters for \code{RF}
##' @param ctrl.ANN control parameters for \code{ANN}
##' @param ctrl.MARS control parameters for \code{MARS}
##' @param ctrl.FDA control parameters for \code{FDA}
##' @param metric a \code{character} corresponding to the evaluation metric used to select 
##' optimal models and tune parameters, must be either \code{ROC} or \code{TSS} 
##' (\emph{maximizing Sensitivity and Specificity)
##' @param metric.ME a \code{character} corresponding to the evaluation metric used to select 
##' optimal model and tune parameters for \code{MAXENT.Phillips}, must be either 
##' \code{auc.val.avg}, \code{auc.diff.avg}, \code{or.mtp.avg}, \code{or.10p.avg} or \code{AICc}
##' @param tuneLength (see \code{\link[caret]{train}})
##' @param method.RF a \code{character} corresponding to the classification or regression model 
##' to use for \code{RF}, must be \code{rf} (see 
##' http://topepo.github.io/caret/Random_Forest.html)
##' @param method.ANN a \code{character} corresponding to the classification or regression model 
##' to use for \code{ANN}, must be \code{avNNet} (see 
##' http://topepo.github.io/caret/Neural_Network.html)
##' @param method.MARS a \code{character} corresponding to the classification or regression model 
##' to use for \code{MARS}, must be \code{earth} (see 
##' http://topepo.github.io/caret/Multivariate_Adaptive_Regression_Splines.html)
##' @param method.GAM a \code{character} corresponding to the classification or regression model 
##' to use for \code{GAM}, must be \code{gam} (see 
##' http://topepo.github.io/caret/Generalized_Additive_Model.html)
##' @param method.GLM a \code{character} corresponding to the classification or regression model 
##' to use for \code{GLM}, must be \code{glmStepAIC} (see 
##' http://topepo.github.io/caret/Generalized_Linear_Model.html)
##' @param type.GLM a \code{vector} of \code{character} corresponding to modeling types for 
##' \code{GLM}, must be among \code{simple}, \code{quadratic}, \code{polynomial}, 
##' \code{s_smoother}
##' @param interaction.GLM a \code{vector} of interaction types, must be amond \code{0}, \code{1}
##' @param cvmethod.ME a \code{character} corresponding to the method used to partition data for 
##' \code{MAXENT.Phillips}, must be \code{randomkfold}
##' @param kfolds.ME an \code{integer} corresponding to the number of bins for k-fold 
##' cross-validation for \code{MAXENT.Phillips}
##' @param overlap.ME (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether to calculate pairwise metric of niche overlap or not 
##' (see \code{\link{calc.niche.overlap}})
##' @param clamp.ME (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether \emph{Features are constrained to remain within the 
##' range of values in the training data} (Elith et al. 2011) or not
##' @param n.bg.ME an \code{integer} corresponding to the number of background points used to run 
##' \code{MAXENT.Phillips}
##' @param env.ME a \code{RasterStack} object containing model predictor variables
##' @param size.tune.ANN a \code{vector} of size parameters (number of units in the hidden layer) 
##' for \code{ANN}
##' @param decay.tune.ANN a \code{vector} of weight decay parameters for \code{ANN}
##' @param maxit.ANN an \code{integer} corresponding to the maximum number of iterations for 
##' \code{ANN}
##' @param MaxNWts.ANN an \code{integer} corresponding to the maximum allowable number of weights 
##' for \code{ANN}
##' @param parallel.ME (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether to enable parallel computing for 
##' \code{MAXENT.Phillips} or not
##' @param numCores.ME an \code{integer} corresponding to the number of cores to be used to 
##' train \code{MAXENT.Phillips}
##' @param Yweights a \code{vector} of response points weights for models allowing case weights
##' 
##' 
##' @value
##' 
##' A \code{BIOMOD.models.options} object (see \code{\link{BIOMOD_ModelingOptions}}) with 
##' optimized parameters
##' 
##' 
##' @details 
##' 
##' \code{trControl} parameter is set by default to :
##' \code{caret::trainControl(method = 'cv', summaryFunction = caret::twoClassSummary, 
##' classProbs = TRUE, returnData = FALSE)}.
##' 
##' All control parameters for other models are set to \code{trControl} if unspecified.
##' 
##' For more details on \code{MAXENT.Phillips} tuning, please refer to  
##' \code{\link[ENMeval]{ENMevaluate}}.
##' For more details on other models tuning, please refer to \code{\link[caret]{train}}.
##' 
##' 
##' @references
##' 
##' \itemize{
##'   \item Kuhn, Max. 2008. Building predictive models in R using the caret package. 
##'   \emph{Journal of Statistical Software} \bold{28}, 1-26.
##'   \item Kuhn, Max, and Kjell Johnson. 2013. Applied predictive modeling. New York: Springer.
##'   \item Muscarella, Robert, et al. 2014. ENMeval: An R package for conducting spatially 
##'   independent evaluations and estimating optimal model complexity for Maxent ecological 
##'   niche models. \emph{Methods in Ecology and Evolution}, \bold{5}, 1198-1205.
##' }
##' 
##' 
##' @seealso \code{\link[caret]{train}}, \code{\link[ENMeval]{ENMevaluate}}, 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions}}
##' 
##' 
##' @examples
##' 
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv", package="biomod2"), row.names = 1)
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
##' # 1. Formating Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options
##' ### Duration for turing all models sequential with default settings 
##' ### on 3.4 GHz processor: approx. 45 min tuning all models in parallel
##' ### (on 8 cores) using foreach loops runs much faster: approx. 14 min
##' 
##' # library(doParallel)
##' # cl<-makeCluster(8)
##' # doParallel::registerDoParallel(cl) 
##' 
##' time.seq <- system.time(
##'   bm.tuning <- BIOMOD_Tuning(myBiomodData,
##'                              env.ME = myExpl,
##'                              n.bg.ME = ncell(myExpl))
##' )
##' 
##' # stopCluster(cl)
##' 
##' par(mfrow = c(1,3))
##' plot(bm.tuning$tune.CTA.rpart)
##' plot(bm.tuning$tune.CTA.rpart2)
##' plot(bm.tuning$tune.RF)
##' 
##' myBiomodModelOut <- BIOMOD_Modeling(myBiomodData, 
##'                                     models = c('RF', 'CTA'), 
##'                                     models.options = bm.tuning$models.options, 
##'                                     NbRunEval = 1, 
##'                                     DataSplit = 100, 
##'                                     VarImport = 0,
##'                                     models.eval.meth = c('ROC'),
##'                                     do.full.models = FALSE,
##'                                     modeling.id = 'test')
##' 
##' 
##' 
##' ###############################################################################################


BIOMOD_tuning <- function(data,
                          models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips'),
                          models.options = BIOMOD_ModelingOptions(),
                          method.ANN = 'avNNet',
                          method.RF = 'rf',
                          method.MARS = 'earth',
                          method.GAM = 'gam',
                          method.GLM = 'glmStepAIC',
                          trControl = NULL,
                          metric = 'ROC',
                          ctrl.CTA = NULL,
                          ctrl.RF = NULL,
                          ctrl.ANN = NULL,
                          ctrl.MARS = NULL,                                  
                          ctrl.FDA = NULL,
                          ctrl.GAM = NULL,
                          ctrl.GBM = NULL,
                          ctrl.GLM = NULL,
                          tuneLength = 30,
                          decay.tune.ANN = c(0.001, 0.01, 0.05, 0.1),
                          size.tune.ANN = c(2, 4, 6, 8),
                          maxit.ANN = 500,
                          MaxNWts.ANN = 10 * (ncol(data@data.env.var) + 1) + 10 + 1,
                          type.GLM = c('simple', 'quadratic', 'polynomial', 's_smoother'), 
                          interaction.GLM = c(0, 1),
                          cvmethod.ME = 'randomkfold',
                          overlap.ME = FALSE,
                          kfolds.ME = 10,
                          n.bg.ME = 10000,
                          env.ME = NULL,
                          metric.ME = 'ROC',
                          clamp.ME = TRUE,
                          parallel.ME = FALSE,
                          numCores.ME = NULL,
                          Yweights = NULL)
{
  
  ## MAXENT: http://cran.r-project.org/web/packages/ENMeval/ENMeval.pdf --> ENMevaluate()
  ## or:    http://cran.r-project.org/web/packages/maxent/maxent.pdf -->  tune.maxent()
  
  ## 0. Check namespaces --------------------------------------------------------------------------
  mod.names = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips')
  
  if (sum(mod.names %in% models) > 0) {
    if (!isNamespaceLoaded("caret")) { requireNamespace("caret", quietly = TRUE) }
    if (!isNamespaceLoaded('dplyr')) { requireNamespace("dplyr", quietly = TRUE) }
    if (is.null(trControl)) {
      trControl <- caret::trainControl(method = "cv",
                                       summaryFunction = caret::twoClassSummary,
                                       classProbs = TRUE,
                                       returnData = FALSE)
    }
    if ("MAXENT.Phillips" %in% models && !isNamespaceLoaded('ENMeval')) { requireNamespace("ENMeval", quietly = TRUE) }
    # if ("MAXENT.Tsuruoka" %in% models && !isNamespaceLoaded('maxent')) { requireNamespace("maxent", quietly = TRUE) }
  }
  
  tune.SRE <- tune.GLM <- tune.MAXENT.Phillips <- tune.GAM <- tune.GBM <- 
    tune.CTA.rpart <- tune.CTA.rpart2 <- tune.RF <- tune.ANN <- tune.MARS <- tune.FDA <- NULL
  # tune.MAXENT.Tsuruoka <- NULL
  
  
  
  ## 1.1 SRE --------------------------------------------------------------------------------------
  resp <- data@data.species
  
  if ('SRE' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning SRE\n"))
    
    tune.SRE = foreach(rep = 1:trControl$repeats, .combine = "rbind") %do%
      {
        fold <- dismo::kfold(resp, by = resp, k = trControl$number)
        RES = foreach (quant = c(0, 0.0125, 0.025, 0.05, 0.1), .combine = "rbind") %:%
          foreach (i = 1:trControl$number, .combine = "rbind") %do%
          {
            DATA <- cbind(1:sum(fold == i)
                          , resp[fold == i]
                          , sre(Response = resp[fold != i],
                                Explanatory = data@data.env.var[fold != i, ],
                                NewData = data@data.env.var[fold == i, ],
                                Quant = quant,
                                return_extremcond = FALSE))
            
            RES = presence.absence.accuracy(DATA, threshold = as.vector(
              PresenceAbsence::optimal.thresholds(DATA, opt.methods = 3)[2], mode = "numeric"))
            return(data.frame(RES, quant = quant))
          }
        return(RES)
      }
    
    t <- aggregate(tune.SRE, by = list(quant = tune.SRE$quant), mean)
    if (metric == 'ROC') {
      models.options@SRE$quant <- t[which.max(t$AUC), "quant"]
    } else if (metric == 'TSS') {
      models.options@SRE$quant <- t[which.max(t$sensitivity + t$specificity - 1), "quant"]
    }
    cat(paste("Finished tuning SRE\n", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if(metric == 'ROC' | metric == 'TSS'){ resp <- as.factor(ifelse(resp == 1 & !is.na(resp), "Presence", "Absence")) }
  
  ## 1.2 GBM --------------------------------------------------------------------------------------
  
  if ('GBM' %in% models)
  {  
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning GBM. Start coarse tuning\n"))
    if (is.null(ctrl.GBM)) { ctrl.GBM <- trControl }
    
    tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
                             .n.trees = c(500, 1000, 2500),
                             .shrinkage = c(0.001, 0.01, 0.1),
                             .n.minobsinnode = 10)
    try(tune.GBM <- caret::train(data@data.env.var, 
                                 resp,
                                 method = "gbm",
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.GBM,
                                 verbose = FALSE,
                                 weights = Yweights))
    
    cat("Best optimization of coarse tuning:\n")
    cat(paste(tune.GBM$bestTune,"\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if (!is.null(tune.GBM)) {
      cat("Start fine tuning\n")
      
      if (tune.GBM$bestTune$n.trees == 2500) {
        cat("Best optimization with large trees! Tuning GBM will take a while.\n")
        n.trees <- seq(2500, 10000, by = 2500)
      } else if(tune.GBM$bestTune$n.trees == 1000) {
        n.trees <- seq(750, 2000, by = 250)
      } if(tune.GBM$bestTune$n.trees == 500) {
        n.trees <- seq(100, 1000, by = 50)
      }
      
      tune.grid <- expand.grid(.interaction.depth = c(tune.GBM$bestTune$interaction.depth - 1
                                                      , tune.GBM$bestTune$interaction.depth
                                                      , tune.GBM$bestTune$interaction.depth + 1),
                               .n.trees = n.trees,
                               .shrinkage = c(tune.GBM$bestTune$shrinkage / 2,
                                              tune.GBM$bestTune$shrinkage,
                                              tune.GBM$bestTune$shrinkage * 5), 
                               .n.minobsinnode = 10)
      tune.GBM <- NULL
      try(tune.GBM <- caret::train(data@data.env.var, 
                                   resp,
                                   method = "gbm",
                                   tuneGrid = tune.grid,
                                   trControl = ctrl.GBM,
                                   verbose = FALSE,
                                   weights = Yweights))  
    }
    
    if (!is.null(tune.GBM)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.GBM$results[, c("Sens", "Spec")], 1, sum) - 1)
        models.options@GBM$n.trees <- tune.GBM$results[tmp, "n.trees"]
        models.options@GBM$interaction.depth <- tune.GBM$results[tmp, "interaction.depth"]
        models.options@GBM$shrinkage <- tune.GBM$results[tmp, "shrinkage"]
      } else {
        models.options@GBM$n.trees <- tune.GBM$bestTune$n.trees
        models.options@GBM$interaction.depth <- tune.GBM$bestTune$interaction.depth
        models.options@GBM$shrinkage <- tune.GBM$bestTune$shrinkage
      }
    } else {
      cat("Tuning GBM failed!")
      tune.GBM <- "FAILED"
    }
    cat(paste("\n Finished tuning GBM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.3 RF ---------------------------------------------------------------------------------------
  
  if ('RF' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning RF\n"))
    if (is.null(ctrl.RF)) { ctrl.RF <- trControl }
    tuneLength.rf <- min(tuneLength, ncol(data@data.env.var))
    
    try(tune.RF <- caret::train(data@data.env.var, 
                                resp,
                                method = method.RF,
                                tuneLength = tuneLength.rf,
                                trControl = ctrl.RF,
                                metric = metric,
                                weights = Yweights))
    
    if (!is.null(tune.RF)) { ## give both mtry as bestTune
      if (metric == 'TSS') {
        models.options@RF$mtry <- tune.RF$results[which.max(apply(tune.RF$results[, c("Sens", "Spec")], 1, sum) - 1), "mtry"]
      } else {
        models.options@RF$mtry <- tune.RF$bestTune$mtry
      }
    } else {
      cat("Tuning RF failed!")
      tune.RF <- "FAILED"
    }
    cat(paste("Finished tuning RF\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.4 ANN --------------------------------------------------------------------------------------
  
  if ('ANN' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning ANN\n"))
    if (is.null(ctrl.ANN)) { ctrl.ANN <- trControl }
    ## already tuning: 
    # size: optimised by cross validation based on model AUC (NbCv cross validation; tested size will be the following c(2,4,6, 8))
    # decay: optimised by cross validation on model AUC (NbCv cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1)).
    # could increase maxit from 200 to 500
    # a nice option would be to use model averaging for ann: avNNet in package(caret)
    
    ## Create a specific candidate set of models to evaluate:
    tune.grid <- expand.grid(.decay = decay.tune.ANN,
                             .size = size.tune.ANN,
                             .bag = FALSE)
    try(tune.ANN <- caret::train(data@data.env.var, 
                                 resp, 
                                 method = method.ANN,
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.ANN,
                                 ## Automatically standardize data prior to modeling and prediction
                                 preProc = c("center", "scale"),
                                 linout = TRUE,
                                 trace = FALSE,
                                 MaxNWts.ANN = MaxNWts.ANN,
                                 maxit = maxit.ANN,
                                 metric = metric,
                                 weights = Yweights))
    if (!is.null(tune.ANN)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.ANN$results[, c("Sens", "Spec")], 1, sum) - 1)
        models.options@ANN$size <- tune.ANN$results[tmp, "size"]
        models.options@ANN$decay <- tune.ANN$results[tmp, "decay"]
        models.options@ANN$maxit <- maxit.ANN
      } else{
        models.options@ANN$size <- tune.ANN$bestTune$size
        models.options@ANN$decay <- tune.ANN$bestTune$decay
        models.options@ANN$maxit <- maxit.ANN
      }
    } else {
      cat("Tuning ANN failed!")
      tune.ANN <- "FAILED"
    }
    cat(paste("Finished tuning ANN\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.5 GAM --------------------------------------------------------------------------------------
  
  if ('GAM' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning GAM\n"))
    if (is.null(ctrl.GAM)) { ctrl.GAM <- trControl }
    
    try(tune.GAM <-   caret::train(data@data.env.var, 
                                   resp, 
                                   method = method.GAM,
                                   trControl = ctrl.GAM,
                                   weights = Yweights))
    if (!is.null(tune.GAM)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.GAM$results[, c("Sens", "Spec")], 1, sum) - 1)
        models.options@GAM$select <- tune.GAM$results[tmp, "select"]
        models.options@GAM$method <- as.character(tune.GAM$results[tmp, "method"])
      } else {
        models.options@GAM$select <- tune.GAM$bestTune$select
        models.options@GAM$method <- as.character(tune.GAM$bestTune$method)
      }
    } else {
      cat("Tuning GAM failed!")
      tune.GAM <- "FAILED"
    }
    cat(paste("Finished tuning GAM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.6 MARS -------------------------------------------------------------------------------------
  
  if ('MARS' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning MARS\n"))
    if (is.null(ctrl.MARS)) { ctrl.MARS <- trControl }
    if (is.null(models.options@MARS$nk)) {
      nprune <- 2:max(21, 2 * ncol(data@data.env.var) + 1)
    } else {
      nprune <- 2:min(models.options@MARS$nk, 38)
    }
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = nprune)
    try(tune.MARS <- caret::train(data@data.env.var, 
                                  resp, 
                                  method = method.MARS,
                                  tuneGrid = tune.grid,
                                  trControl = ctrl.MARS,
                                  weights = Yweights))
    
    if (!is.null(tune.MARS)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.MARS$results[, c("Sens", "Spec")], 1, sum) - 1)
        if ("degree" %in% names(models.options@MARS)) {
          models.options@MARS$degree <- tune.MARS$results[tmp, "degree"]
        } else {
          models.options@MARS$interaction.level <- tune.MARS$results[tmp, "degree"] - 1
        }
        models.options@MARS$nprune <- tune.MARS$results[tmp, "nprune"]
      } else {
        if ("degree" %in% names(models.options@MARS)) {
          models.options@MARS$degree <- tune.MARS$bestTune$degree
        } else {
          models.options@MARS$interaction.level <- tune.MARS$bestTune$degree - 1
        }
        models.options@MARS$nprune <- tune.MARS$bestTune$nprune
      }
    } else {
      cat("Tuning MARS failed!")
      tune.MARS <- "FAILED"
    }
    cat(paste("Finished tuning MARS\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.7 GLM --------------------------------------------------------------------------------------
  
  if ('GLM' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning GLM\n"))
    if (is.null(ctrl.GLM)) { ctrl.GLM <- trControl }
    if ("s_smoother" %in% type.GLM) { requireNamespace("gam", quietly = TRUE) }
    
    fm <- list()
    GLM.results = foreach (type = type.GLM, .combine = "rbind") %:%
      foreach (IA = interaction.GLM, .combine = "rbind") %do%
      {
        try(tune.GLM <- caret::train(bm_MakeFormula("resp",
                                                    data@data.env.var,
                                                    type = type,
                                                    interaction.level = IA),
                                     data = cbind(data@data.env.var, resp = resp),
                                     method = method.GLM,
                                     trControl = ctrl.GLM,
                                     weights = Yweights))
        try(fm[[length(fm) + 1]] <- formula(tune.GLM$finalModel))
        try(RES <- cbind(tune.GLM$results, il = IA, type = type))
        return(RES)
      }
    
    glm.best <- which.max(GLM.results$ROC)
    models.options@GLM$interaction.level <- GLM.results[glm.best, "il"]
    models.options@GLM$type <- as.character(GLM.results[glm.best, "type"])
    models.options@GLM$myFormula <- formula(paste(data@sp.name, "~", gsub("`", "", as.character(fm[[glm.best]])[3])))
    models.options@GLM$test <- "none" 
    
    cat(paste("Finished tuning GLM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }      
  
  ## 1.8 FDA --------------------------------------------------------------------------------------
  
  if ('FDA' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning FDA\n"))
    if (is.null(ctrl.FDA)) { ctrl.FDA <- trControl }
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    try(tune.FDA <- caret::train(data@data.env.var, 
                                 factor(resp), 
                                 method = "fda",
                                 tuneGrid = tune.grid,                  
                                 trControl = ctrl.FDA,
                                 weights = Yweights))
    
    if (!is.null(tune.FDA)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.FDA$results[, c("Sens", "Spec")], 1, sum) - 1)
        models.options@FDA$add_args <- list(degree = tune.FDA$results[tmp, "degree"],
                                            nprune = tune.FDA$results[tmp, "nprune"])
      } else {
        models.options@FDA$add_args <- list(degree = tune.FDA$bestTune$degree,
                                            nprune = tune.FDA$bestTune$nprune)
      }
    } else {
      cat("Tuning FDA failed!")
      tune.FDA <- "FAILED"
    }
    cat(paste("Finished tuning FDA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.9 CTA --------------------------------------------------------------------------------------
  
  if ('CTA' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning CTA\n"))
    if (is.null(ctrl.CTA)) { ctrl.CTA <- trControl }
    
    cat("Tuning Complexity Parameter")    
    try(tune.CTA.rpart <- caret::train(data@data.env.var, 
                                       resp, 
                                       method = "rpart",
                                       tuneLength = tuneLength,
                                       trControl = ctrl.CTA,
                                       metric = metric,
                                       weights = Yweights))
    
    cat("Tuning Max Tree Depth")
    try(tune.CTA.rpart2 <-  caret::train(data@data.env.var, 
                                         resp,
                                         method = "rpart2",
                                         tuneLength = tuneLength,
                                         trControl = ctrl.CTA,
                                         metric = metric,
                                         weights = Yweights))
    
    if (!is.null(tune.CTA.rpart)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.CTA.rpart$results[, c("Sens", "Spec")], 1, sum) - 1)
        models.options@CTA$control$cp <- tune.CTA.rpart$results[tmp, "cp"]
      } else {
        models.options@CTA$control$cp <- tune.CTA.rpart$bestTune
      }
    } else {
      cat("Tuning CTA cp failed!")
      tune.CTA.rpart <- "FAILED"
    }
    
    if (!is.null(tune.CTA.rpart2)) {
      if (metric == 'TSS') {
        tmp = which.max(apply(tune.CTA.rpart2$results[, c("Sens", "Spec")], 1, sum) - 1)
        models.options@CTA$control$maxdepth <- tune.CTA.rpart2$results[tmp, "maxdepth"]
      } else {
        models.options@CTA$control$maxdepth <- tune.CTA.rpart2$bestTune
      }
    } else {
      cat("Tuning CTA maxdepth failed!")
      tune.CTA.rpart2 <- "FAILED"
    }
    cat(paste("Finished tuning CTA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.10 MAXENT.Phillips -------------------------------------------------------------------------
  
  if ('MAXENT.Phillips' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning MAXENT.Phillips\n"))
    if (cvmethod.ME != 'randomkfold') { kfolds.ME <- NA }
    
    try(tune.MAXENT.Phillips <- tuning.maxent(pres = data@data.env.var[data@data.species == 1 & !is.na(data@data.species), ],
                                              bg = data@data.env.var[data@data.species == 0 | is.na(data@data.species), ],
                                              method = cvmethod.ME, 
                                              kfolds = kfolds.ME, #env.ME,
                                              clamp = clamp.ME, 
                                              parallel = parallel.ME, 
                                              numCores = numCores.ME,
                                              categoricals = NULL))
    
    if (!is.null(tune.MAXENT.Phillips)) {
      if (!metric.ME %in% c("auc.val.avg", "auc.diff.avg", "or.mtp.avg", "or.10p.avg", "AICc")) {
        metric.ME <- "auc.val.avg"
        cat("Invalid metric.ME argument! metric.ME was set to auc.val.avg")
      }
      if (metric.ME == 'auc.val.avg') {
        tmp = which.max(tune.MAXENT.Phillips@results[, metric.ME])
      } else {
        tmp = which.min(tune.MAXENT.Phillips@results[, metric.ME])
      }
      models.options@MAXENT.Phillips$linear <- grepl("L", tune.MAXENT.Phillips@results[tmp, "fc"])
      models.options@MAXENT.Phillips$quadratic <- grepl("Q", tune.MAXENT.Phillips@results[tmp, "fc"])
      models.options@MAXENT.Phillips$hinge <- grepl("H", tune.MAXENT.Phillips@results[tmp, "fc"])
      models.options@MAXENT.Phillips$product <- grepl("P", tune.MAXENT.Phillips@results[tmp, "fc"])
      models.options@MAXENT.Phillips$threshold <- grepl("T", tune.MAXENT.Phillips@results[tmp, "fc"])
      models.options@MAXENT.Phillips$betamultiplier <- tune.MAXENT.Phillips@results[tmp, "rm"]
    } else {
      cat("Tuning MAXENT.Phillips failed!")
      tune.MAXENT.Phillips <- "FAILED"
    }
    cat(paste("Finished tuning MAXENT.Phillips\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  
  # if ('MAXENT.Tsuruoka' %in% models) {
  #   cat("Start tuning MAXENT.Tsuruoka\n")
  #   try(tune.MAXENT.Tsuruoka <- as.data.frame(tune.maxent(data@data.env.var,data@data.species,nfold=kfolds.ME,showall=T)))
  #   cat(paste("Finished tuning MAXENT.Tsuruoka\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  #   
  #   if(!is.null(tune.MAXENT.Tsuruoka)){
  #     models.options@MAXENT.Tsuruoka$l1_regularizer <- tune.MAXENT.Tsuruoka$l1_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     models.options@MAXENT.Tsuruoka$l2_regularizer <- tune.MAXENT.Tsuruoka$l2_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     models.options@MAXENT.Tsuruoka$use_sgd <- ifelse(tune.MAXENT.Tsuruoka[which.max(tune.MAXENT.Tsuruoka$accuracy),]$use_sgd==0,F,T)
  #     models.options@MAXENT.Tsuruoka$set_heldout <- tune.MAXENT.Tsuruoka$set_heldout[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #   } else { if('MAXENT.Tsuruoka' %in% models){cat("Tuning MAXENT.Tsuruoka failed!"); tune.MAXENT.Tsuruoka <- "FAILED"}}
  # }  
  
  
  return(list(models.options=models.options, tune.SRE =tune.SRE,  tune.CTA.rpart = tune.CTA.rpart, tune.CTA.rpart2 = tune.CTA.rpart2,
              tune.RF = tune.RF, tune.ANN = tune.ANN,  tune.MARS = tune.MARS, tune.FDA = tune.FDA, tune.GBM=tune.GBM,
              tune.GAM = tune.GAM, tune.MAXENT.Phillips = tune.MAXENT.Phillips, 
              # tune.MAXENT.Tsuruoka = tune.MAXENT.Tsuruoka, 
              tune.GLM=tune.GLM))
}


###################################################################################################
#### Modified tuning function from the ENMeval package to tune MAXENT.Phillips (internal function for BIOMOD_tuning)

tuning.maxent <- function(pres,
                          bg,
                          method,
                          kfolds,
                          clamp,
                          parallel,
                          numCores,
                          categoricals,
                          tune.args = list(rm = seq(0.5, 1, 0.5), fc = c("L")))
{
  requireNamespace("ENMeval", quietly = TRUE)
  results <- ENMeval::ENMevaluate(occs = pres,
                                  bg = bg,
                                  tune.args = tune.args,
                                  partitions = method,
                                  algorithm = "maxent.jar",
                                  partition.settings = list(kfolds = kfolds),
                                  doClamp = clamp,
                                  parallel = parallel,
                                  numCores = numCores)
  return(results)
}




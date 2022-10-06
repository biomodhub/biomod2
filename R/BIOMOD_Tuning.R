###################################################################################################
##' @name BIOMOD_Tuning
##' @author Frank Breiner
##' 
##' @title Tune models parameters
##' 
##' @description Function to tune \pkg{biomod2} single models parameters
##'
##'
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param bm.options a \code{\link{BIOMOD.models.options}} object returned by the  
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param models a \code{vector} containing model names to be tuned, \cr 
##' must be among \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, 
##' \code{FDA}, \code{MARS}, \code{RF}, \code{MAXENT.Phillips}
##' @param metric.eval a \code{character} corresponding to the evaluation metric used to select 
##' optimal models and tune parameters, must be either \code{ROC} or \code{TSS} 
##' (\emph{maximizing Sensitivity and Specificity})
##' @param ctrl.train global control parameters that can be obtained from the 
##' \code{\link[caret]{trainControl}} function
##' @param ctrl.train.tuneLength (see \code{tuneLength} parameter in \code{\link[caret]{train}})
##' @param ctrl.ANN control parameters for \code{ANN}
##' @param ctrl.CTA control parameters for \code{CTA}
##' @param ctrl.FDA control parameters for \code{FDA}
##' @param ctrl.GAM control parameters for \code{GAM}
##' @param ctrl.GBM control parameters for \code{GBM}
##' @param ctrl.GLM control parameters for \code{GLM}
##' @param ctrl.MARS control parameters for \code{MARS}
##' @param ctrl.RF control parameters for \code{RF}
##' @param ANN.method a \code{character} corresponding to the classification or regression model 
##' to use for \code{ANN}, \cr 
##' must be \code{avNNet} (see \url{http://topepo.github.io/caret/train-models-by-tag.html#Neural_Network})
##' @param ANN.size.tune a \code{vector} of size parameters (number of units in the hidden layer) 
##' for \code{ANN}
##' @param ANN.decay.tune a \code{vector} of weight decay parameters for \code{ANN}
##' @param ANN.maxit an \code{integer} corresponding to the maximum number of iterations for 
##' \code{ANN}
##' @param ANN.MaxNWts an \code{integer} corresponding to the maximum allowable number of weights 
##' for \code{ANN}
##' @param GAM.method a \code{character} corresponding to the classification or regression model 
##' to use for \code{GAM}, \cr 
##' must be \code{gam} (see \url{http://topepo.github.io/caret/train-models-by-tag.html#generalized-additive-model})
##' @param GLM.method a \code{character} corresponding to the classification or regression model 
##' to use for \code{GLM}, \cr 
##' must be \code{glmStepAIC} (see 
##' \url{http://topepo.github.io/caret/train-models-by-tag.html#Generalized_Linear_Model})
##' @param GLM.type a \code{vector} of \code{character} corresponding to modeling types for 
##' \code{GLM}, \cr must be among \code{simple}, \code{quadratic}, \code{polynomial}, 
##' \code{s_smoother}
##' @param GLM.interaction a \code{vector} of interaction types, must be among \code{0}, \code{1}
##' @param MARS.method a \code{character} corresponding to the classification or regression model 
##' to use for \code{MARS}, \cr 
##' must be \code{earth} (see 
##' \url{http://topepo.github.io/caret/train-models-by-tag.html#Multivariate_Adaptive_Regression_Splines})
##' @param ME.metric a \code{character} corresponding to the evaluation metric used to select 
##' optimal model and tune parameters for \code{MAXENT.Phillips}, must be either 
##' \code{auc.val.avg}, \code{auc.diff.avg}, \code{or.mtp.avg}, \code{or.10p.avg} or \code{AICc}
##' @param ME.cvmethod a \code{character} corresponding to the method used to partition data for 
##' \code{MAXENT.Phillips}, \cr must be \code{randomkfold}
##' @param ME.kfolds an \code{integer} corresponding to the number of bins for k-fold 
##' cross-validation for \code{MAXENT.Phillips}
##' @param ME.overlap (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether to calculate pairwise metric of niche overlap or not 
##' (see \code{\link[ENMeval]{calc.niche.overlap}})
##' @param ME.clamp (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether \emph{Features are constrained to remain within the 
##' range of values in the training data} (Elith et al. 2011) or not
##' @param ME.n.bg an \code{integer} corresponding to the number of background points used to run 
##' \code{MAXENT.Phillips}
##' @param ME.env a \code{RasterStack} object containing model predictor variables
##' @param ME.parallel (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether to enable parallel computing for 
##' \code{MAXENT.Phillips} or not
##' @param ME.numCores an \code{integer} corresponding to the number of cores to be used to 
##' train \code{MAXENT.Phillips}
##' @param RF.method a \code{character} corresponding to the classification or regression model 
##' to use for \code{RF}, \cr 
##' must be \code{rf} (see \url{http://topepo.github.io/caret/train-models-by-tag.html#random-forest})
##' @param weights a \code{vector} of \code{numeric} values corresponding to observation weights
##' 
##' 
##' @return 
##' 
##' A \code{\link{BIOMOD.models.options}} object (see \code{\link{BIOMOD_ModelingOptions}}) with 
##' optimized parameters
##' 
##' 
##' @details 
##' 
##' \itemize{
##'   \item \code{ctrl.train} parameter is set by default to : \cr
##'   \code{caret::trainControl(method = 'cv', summaryFunction = caret::twoClassSummary,} \cr
##'   \code{classProbs = TRUE, returnData = FALSE)}.
##'   \item All control parameters for other models are set to \code{ctrl.train} if unspecified.
##'   \item For more details on \code{MAXENT.Phillips} tuning, please refer to 
##'   \code{\link[ENMeval]{ENMevaluate}}.
##'   \item For more details on other models tuning, please refer to \code{\link[caret]{train}}.
##' }
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
##' @seealso \code{\link[caret]{trainControl}}, \code{\link[caret]{train}}, 
##' \code{\link[ENMeval]{calc.niche.overlap}}, \code{\link[ENMeval]{ENMevaluate}}, 
##' \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Modeling}}
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
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' 
##' # ---------------------------------------------------------------
##' ### Duration for turing all models sequential with default settings 
##' ### on 3.4 GHz processor: approx. 45 min tuning all models in parallel
##' ### (on 8 cores) using foreach loops runs much faster: approx. 14 min
##' 
##' \dontrun{
##' # library(doParallel)
##' # cl <- makeCluster(8)
##' # doParallel::registerDoParallel(cl) 
##' 
##' time.seq <- system.time(
##'   bm.tuning <- BIOMOD_Tuning(bm.format = myBiomodData, ME.env = myExpl, ME.n.bg = ncell(myExpl))
##' )
##' 
##' # stopCluster(cl)
##' 
##' plot(bm.tuning$tune.CTA.rpart)
##' plot(bm.tuning$tune.CTA.rpart2)
##' plot(bm.tuning$tune.RF)
##' plot(bm.tuning$tune.ANN)
##' plot(bm.tuning$tune.MARS)
##' plot(bm.tuning$tune.FDA)
##' plot(bm.tuning$tune.GBM)
##' plot(bm.tuning$tune.GAM)
##' 
##' # Get tuned modeling options
##' myBiomodOptions <- bm.tuning$models.options
##' }
##' 
##' 
##' @importFrom foreach foreach
## @importFrom caret trainControl train twoClassSummary
## @importFrom dismo kfold
##' @importFrom PresenceAbsence optimal.thresholds presence.absence.accuracy
## @importFrom ENMeval ENMevaluate
##' @importFrom stats aggregate as.formula binomial complete.cases cor formula glm 
##' median na.exclude na.omit qt quantile sd
##' 
##' 
##' @export
##' 
###################################################################################################


BIOMOD_Tuning <- function(bm.format,
                          bm.options = BIOMOD_ModelingOptions(),
                          models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips'),
                          metric.eval = 'ROC',
                          ctrl.train = NULL,
                          ctrl.train.tuneLength = 30,
                          ctrl.ANN = NULL,
                          ctrl.CTA = NULL,
                          ctrl.FDA = NULL,
                          ctrl.GAM = NULL,
                          ctrl.GBM = NULL,
                          ctrl.GLM = NULL,
                          ctrl.MARS = NULL,                                  
                          ctrl.RF = NULL,
                          ANN.method = 'avNNet',
                          ANN.decay.tune = c(0.001, 0.01, 0.05, 0.1),
                          ANN.size.tune = c(2, 4, 6, 8),
                          ANN.maxit = 500,
                          ANN.MaxNWts = 10 * (ncol(bm.format@data.env.var) + 1) + 10 + 1,
                          MARS.method = 'earth',
                          GAM.method = 'gam',
                          GLM.method = 'glmStepAIC',
                          GLM.type = c('simple', 'quadratic', 'polynomial', 's_smoother'), 
                          GLM.interaction = c(0, 1),
                          ME.cvmethod = 'randomkfold',
                          ME.overlap = FALSE,
                          ME.kfolds = 10,
                          ME.n.bg = 10000,
                          ME.env = NULL,
                          ME.metric = 'ROC',
                          ME.clamp = TRUE,
                          ME.parallel = FALSE,
                          ME.numCores = NULL,
                          RF.method = 'rf',
                          weights = NULL)
{
  .bm_cat("Tune Modeling Options")
  
  ## MAXENT: http://cran.r-project.org/web/packages/ENMeval/ENMeval.pdf --> ENMevaluate()
  ## or:    http://cran.r-project.org/web/packages/maxent/maxent.pdf -->  tune.maxent()
  
  ## 0. Check namespaces --------------------------------------------------------------------------
  mod.names = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips')
  
  if (sum(mod.names %in% models) > 0) {
    if (!isNamespaceLoaded("caret")) { requireNamespace("caret", quietly = TRUE) }
    if (!isNamespaceLoaded('dplyr')) { requireNamespace("dplyr", quietly = TRUE) }
    if (is.null(ctrl.train)) {
      ctrl.train <- caret::trainControl(method = "cv",
                                        repeats = 3,
                                        summaryFunction = caret::twoClassSummary,
                                        classProbs = TRUE,
                                        returnData = FALSE)
    }
    if ("MAXENT.Phillips" %in% models && !isNamespaceLoaded('ENMeval')) { requireNamespace("ENMeval", quietly = TRUE) }
    # if ("MAXENT.Tsuruoka" %in% models && !isNamespaceLoaded('maxent')) { requireNamespace("maxent", quietly = TRUE) }
    if ("SRE" %in% models && !isNamespaceLoaded('dismo')) { requireNamespace("dismo", quietly = TRUE) }
  }
  
  tune.SRE <- tune.GLM <- tune.MAXENT.Phillips <- tune.GAM <- tune.GBM <- 
    tune.CTA.rpart <- tune.CTA.rpart2 <- tune.RF <- tune.ANN <- tune.MARS <- tune.FDA <- NULL
  # tune.MAXENT.Tsuruoka <- NULL
  
  resp <- bm.format@data.species
  # if (is.null(weights)) { weights = rep(1, length(bm.format@data.species))}
  
  
  ## 1.1 SRE --------------------------------------------------------------------------------------
  
  if ('SRE' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning SRE\n"))
    
    tune.SRE = foreach(rep = 1:ctrl.train$repeats, .combine = "rbind") %do%
      {
        fold <- dismo::kfold(resp, by = resp, k = ctrl.train$number)
        RES = foreach (quant = c(0, 0.0125, 0.025, 0.05, 0.1), .combine = "rbind") %:%
          foreach (i = 1:ctrl.train$number, .combine = "rbind") %do%
          {
            DATA <- cbind(1:sum(fold == i)
                          , resp[fold == i]
                          , bm_SRE(resp.var = resp[fold != i],
                                   expl.var = bm.format@data.env.var[fold != i, ],
                                   new.env = bm.format@data.env.var[fold == i, ],
                                   quant = quant,
                                   do.extrem = FALSE))
            
            RES = presence.absence.accuracy(DATA, threshold = as.vector(
              optimal.thresholds(DATA, opt.methods = 3)[2], mode = "numeric"))
            return(data.frame(RES, quant = quant))
          }
        return(RES)
      }
    
    t <- aggregate(tune.SRE, by = list(quant = tune.SRE$quant), mean)
    if (metric.eval == 'ROC') {
      bm.options@SRE$quant <- t[which.max(t$AUC), "quant"]
    } else if (metric.eval == 'TSS') {
      bm.options@SRE$quant <- t[which.max(t$sensitivity + t$specificity - 1), "quant"]
    }
    cat(paste("Finished tuning SRE", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  if(metric.eval == 'ROC' | metric.eval == 'TSS'){ resp <- as.factor(ifelse(resp == 1 & !is.na(resp), "Presence", "Absence")) }
  
  ## 1.2 GBM --------------------------------------------------------------------------------------
  
  if ('GBM' %in% models)
  {  
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start coarse tuning GBM\n"))
    if (is.null(ctrl.GBM)) { ctrl.GBM <- ctrl.train }
    
    tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
                             .n.trees = c(500, 1000, 2500),
                             .shrinkage = c(0.001, 0.01, 0.1),
                             .n.minobsinnode = 10)
    try(tune.GBM <- caret::train(bm.format@data.env.var, 
                                 resp,
                                 method = "gbm",
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.GBM,
                                 verbose = FALSE,
                                 weights = weights))
    
    cat("Best optimization of coarse tuning:\n")
    cat(paste(tune.GBM$bestTune, collapse = ' / '), "\n")
    
    if (!is.null(tune.GBM)) {
      cat(" Start fine tuning GBM\n")
      
      if (tune.GBM$bestTune$n.trees == 2500) {
        cat("Best optimization with large trees! Tuning GBM will take a while.\n")
        n.trees <- seq(2500, 10000, by = 2500)
      } else if(tune.GBM$bestTune$n.trees == 1000) {
        n.trees <- seq(750, 2000, by = 250)
      } else if(tune.GBM$bestTune$n.trees == 500) {
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
      try(tune.GBM <- caret::train(bm.format@data.env.var, 
                                   resp,
                                   method = "gbm",
                                   tuneGrid = tune.grid,
                                   trControl = ctrl.GBM,
                                   verbose = FALSE,
                                   weights = weights))  
    }
    
    if (!is.null(tune.GBM)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.GBM$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@GBM$n.trees <- tune.GBM$results[tmp, "n.trees"]
        bm.options@GBM$interaction.depth <- tune.GBM$results[tmp, "interaction.depth"]
        bm.options@GBM$shrinkage <- tune.GBM$results[tmp, "shrinkage"]
      } else {
        bm.options@GBM$n.trees <- tune.GBM$bestTune$n.trees
        bm.options@GBM$interaction.depth <- tune.GBM$bestTune$interaction.depth
        bm.options@GBM$shrinkage <- tune.GBM$bestTune$shrinkage
      }
    } else {
      cat("Tuning GBM failed!")
      tune.GBM <- "FAILED"
    }
    cat(paste("Finished tuning GBM", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.3 RF ---------------------------------------------------------------------------------------
  
  if ('RF' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning RF\n"))
    if (is.null(ctrl.RF)) { ctrl.RF <- ctrl.train }
    tuneLength.rf <- min(ctrl.train.tuneLength, ncol(bm.format@data.env.var))
    
    try(tune.RF <- caret::train(bm.format@data.env.var, 
                                resp,
                                method = RF.method,
                                tuneLength = tuneLength.rf,
                                trControl = ctrl.RF,
                                metric = metric.eval,
                                weights = weights))
    
    if (!is.null(tune.RF)) { ## give both mtry as bestTune
      if (metric.eval == 'TSS') {
        bm.options@RF$mtry <- tune.RF$results[which.max(apply(tune.RF$results[, c("Sens", "Spec")], 1, sum) - 1), "mtry"]
      } else {
        bm.options@RF$mtry <- tune.RF$bestTune$mtry
      }
    } else {
      cat("Tuning RF failed!")
      tune.RF <- "FAILED"
    }
    cat(paste("Finished tuning RF", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.4 ANN --------------------------------------------------------------------------------------
  
  if ('ANN' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning ANN\n"))
    if (is.null(ctrl.ANN)) { ctrl.ANN <- ctrl.train }
    ## already tuning: 
    # size: optimised by cross validation based on model AUC (NbCv cross validation; tested size will be the following c(2,4,6, 8))
    # decay: optimised by cross validation on model AUC (NbCv cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1)).
    # could increase maxit from 200 to 500
    # a nice option would be to use model averaging for ann: avNNet in package(caret)
    
    ## Create a specific candidate set of models to evaluate:
    tune.grid <- expand.grid(.decay = ANN.decay.tune,
                             .size = ANN.size.tune,
                             .bag = FALSE)
    try(tune.ANN <- caret::train(bm.format@data.env.var, 
                                 resp, 
                                 method = ANN.method,
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.ANN,
                                 ## Automatically standardize data prior to modeling and prediction
                                 preProc = c("center", "scale"),
                                 linout = TRUE,
                                 trace = FALSE,
                                 MaxNWts.ANN = ANN.MaxNWts,
                                 maxit = ANN.maxit,
                                 metric = metric.eval,
                                 weights = weights))
    if (!is.null(tune.ANN)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.ANN$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@ANN$size <- tune.ANN$results[tmp, "size"]
        bm.options@ANN$decay <- tune.ANN$results[tmp, "decay"]
        bm.options@ANN$maxit <- ANN.maxit
      } else{
        bm.options@ANN$size <- tune.ANN$bestTune$size
        bm.options@ANN$decay <- tune.ANN$bestTune$decay
        bm.options@ANN$maxit <- ANN.maxit
      }
    } else {
      cat("Tuning ANN failed!")
      tune.ANN <- "FAILED"
    }
    cat(paste("Finished tuning ANN", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.5 GAM --------------------------------------------------------------------------------------
  
  if ('GAM' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning GAM\n"))
    if (is.null(ctrl.GAM)) { ctrl.GAM <- ctrl.train }
    
    try(tune.GAM <- caret::train(bm.format@data.env.var, 
                                 resp, 
                                 method = GAM.method,
                                 trControl = ctrl.GAM,
                                 weights = weights))
    if (!is.null(tune.GAM)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.GAM$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@GAM$select <- tune.GAM$results[tmp, "select"]
        bm.options@GAM$method <- as.character(tune.GAM$results[tmp, "method"])
      } else {
        bm.options@GAM$select <- tune.GAM$bestTune$select
        bm.options@GAM$method <- as.character(tune.GAM$bestTune$method)
      }
    } else {
      cat("Tuning GAM failed!")
      tune.GAM <- "FAILED"
    }
    cat(paste("Finished tuning GAM", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.6 MARS -------------------------------------------------------------------------------------
  
  if ('MARS' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning MARS\n"))
    if (is.null(ctrl.MARS)) { ctrl.MARS <- ctrl.train }
    if (is.null(bm.options@MARS$nk)) {
      nprune <- 2:max(21, 2 * ncol(bm.format@data.env.var) + 1)
    } else {
      nprune <- 2:min(bm.options@MARS$nk, 38)
    }
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = nprune)
    try(tune.MARS <- caret::train(bm.format@data.env.var, 
                                  resp, 
                                  method = MARS.method,
                                  tuneGrid = tune.grid,
                                  trControl = ctrl.MARS,
                                  weights = weights))
    
    if (!is.null(tune.MARS)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.MARS$results[, c("Sens", "Spec")], 1, sum) - 1)
        if ("degree" %in% names(bm.options@MARS)) {
          bm.options@MARS$degree <- tune.MARS$results[tmp, "degree"]
        } else {
          bm.options@MARS$interaction.level <- tune.MARS$results[tmp, "degree"] - 1
        }
        bm.options@MARS$nprune <- tune.MARS$results[tmp, "nprune"]
      } else {
        if ("degree" %in% names(bm.options@MARS)) {
          bm.options@MARS$degree <- tune.MARS$bestTune$degree
        } else {
          bm.options@MARS$interaction.level <- tune.MARS$bestTune$degree - 1
        }
        bm.options@MARS$nprune <- tune.MARS$bestTune$nprune
      }
    } else {
      cat("Tuning MARS failed!")
      tune.MARS <- "FAILED"
    }
    cat(paste("Finished tuning MARS", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.7 GLM --------------------------------------------------------------------------------------
  
  if ('GLM' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning GLM\n"))
    if (is.null(ctrl.GLM)) { ctrl.GLM <- ctrl.train }
    if ("s_smoother" %in% GLM.type) { requireNamespace("gam", quietly = TRUE) }
    
    fm <- list()
    GLM.results = foreach (type = GLM.type, .combine = "rbind") %:%
      foreach (IA = GLM.interaction, .combine = "rbind") %do%
      {
        try(tune.GLM <- caret::train(bm_MakeFormula(resp.name = "resp",
                                                    expl.var = bm.format@data.env.var,
                                                    type = type,
                                                    interaction.level = IA),
                                     data = cbind(bm.format@data.env.var, resp = resp),
                                     method = GLM.method,
                                     trControl = ctrl.GLM))
        # weights = weights))
        try(fm[[length(fm) + 1]] <- formula(tune.GLM$finalModel))
        try(RES <- cbind(tune.GLM$results, il = IA, type = type))
        return(RES)
      }
    tune.GLM <- GLM.results
    
    glm.best <- which.max(GLM.results$ROC)
    bm.options@GLM$interaction.level <- GLM.results[glm.best, "il"]
    bm.options@GLM$type <- as.character(GLM.results[glm.best, "type"])
    bm.options@GLM$myFormula <- formula(paste(bm.format@sp.name, "~", gsub("`", "", as.character(fm[[glm.best]])[3])))
    bm.options@GLM$test <- "none" 
    
    cat(paste("Finished tuning GLM", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }      
  
  ## 1.8 FDA --------------------------------------------------------------------------------------
  
  if ('FDA' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning FDA\n"))
    if (is.null(ctrl.FDA)) { ctrl.FDA <- ctrl.train }
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    try(tune.FDA <- caret::train(bm.format@data.env.var, 
                                 factor(resp), 
                                 method = "fda",
                                 tuneGrid = tune.grid,                  
                                 trControl = ctrl.FDA,
                                 weights = weights))
    
    if (!is.null(tune.FDA)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.FDA$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@FDA$add_args <- list(degree = tune.FDA$results[tmp, "degree"],
                                        nprune = tune.FDA$results[tmp, "nprune"])
      } else {
        bm.options@FDA$add_args <- list(degree = tune.FDA$bestTune$degree,
                                        nprune = tune.FDA$bestTune$nprune)
      }
    } else {
      cat("Tuning FDA failed!")
      tune.FDA <- "FAILED"
    }
    cat(paste("Finished tuning FDA", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.9 CTA --------------------------------------------------------------------------------------
  
  if ('CTA' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning CTA\n"))
    if (is.null(ctrl.CTA)) { ctrl.CTA <- ctrl.train }
    
    cat("Tuning Complexity Parameter")    
    try(tune.CTA.rpart <- caret::train(bm.format@data.env.var, 
                                       resp, 
                                       method = "rpart",
                                       tuneLength = ctrl.train.tuneLength,
                                       trControl = ctrl.CTA,
                                       metric = metric.eval,
                                       weights = weights))
    
    cat("Tuning Max Tree Depth")
    try(tune.CTA.rpart2 <-  caret::train(bm.format@data.env.var, 
                                         resp,
                                         method = "rpart2",
                                         tuneLength = ctrl.train.tuneLength,
                                         trControl = ctrl.CTA,
                                         metric = metric.eval,
                                         weights = weights))
    
    if (!is.null(tune.CTA.rpart)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.CTA.rpart$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@CTA$control$cp <- tune.CTA.rpart$results[tmp, "cp"]
      } else {
        bm.options@CTA$control$cp <- tune.CTA.rpart$bestTune
      }
    } else {
      cat("Tuning CTA cp failed!")
      tune.CTA.rpart <- "FAILED"
    }
    
    if (!is.null(tune.CTA.rpart2)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.CTA.rpart2$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@CTA$control$maxdepth <- tune.CTA.rpart2$results[tmp, "maxdepth"]
      } else {
        bm.options@CTA$control$maxdepth <- tune.CTA.rpart2$bestTune
      }
    } else {
      cat("Tuning CTA maxdepth failed!")
      tune.CTA.rpart2 <- "FAILED"
    }
    cat(paste("Finished tuning CTA", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  ## 1.10 MAXENT.Phillips -------------------------------------------------------------------------
  
  if ('MAXENT.Phillips' %in% models)
  {
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n", "Start tuning MAXENT.Phillips\n"))
    if (ME.cvmethod != 'randomkfold') { ME.kfolds <- NA }
    
    try(tune.MAXENT.Phillips <- .maxent_tuning(pres = bm.format@data.env.var[bm.format@data.species == 1 & !is.na(bm.format@data.species), ],
                                               bg = bm.format@data.env.var[bm.format@data.species == 0 | is.na(bm.format@data.species), ],
                                               method = ME.cvmethod, 
                                               kfolds = ME.kfolds, #ME.env,
                                               clamp = ME.clamp, 
                                               parallel = ME.parallel, 
                                               numCores = ME.numCores,
                                               categoricals = NULL))
    
    if (!is.null(tune.MAXENT.Phillips)) {
      if (!ME.metric %in% c("auc.val.avg", "auc.diff.avg", "or.mtp.avg", "or.10p.avg", "AICc")) {
        ME.metric <- "auc.val.avg"
        cat("Invalid ME.metric argument! ME.metric was set to auc.val.avg")
      }
      if (ME.metric == 'auc.val.avg') {
        tmp = which.max(tune.MAXENT.Phillips@results[, ME.metric])
      } else {
        tmp = which.min(tune.MAXENT.Phillips@results[, ME.metric])
      }
      bm.options@MAXENT.Phillips$linear <- grepl("L", tune.MAXENT.Phillips@results[tmp, "fc"])
      bm.options@MAXENT.Phillips$quadratic <- grepl("Q", tune.MAXENT.Phillips@results[tmp, "fc"])
      bm.options@MAXENT.Phillips$hinge <- grepl("H", tune.MAXENT.Phillips@results[tmp, "fc"])
      bm.options@MAXENT.Phillips$product <- grepl("P", tune.MAXENT.Phillips@results[tmp, "fc"])
      bm.options@MAXENT.Phillips$threshold <- grepl("T", tune.MAXENT.Phillips@results[tmp, "fc"])
      bm.options@MAXENT.Phillips$betamultiplier <- tune.MAXENT.Phillips@results[tmp, "rm"]
    } else {
      cat("Tuning MAXENT.Phillips failed!")
      tune.MAXENT.Phillips <- "FAILED"
    }
    cat(paste("Finished tuning MAXENT.Phillips", "\n-=-=-=-=-=-=-=-=-=-=\n"))
  }
  
  
  # if ('MAXENT.Tsuruoka' %in% models) {
  #   cat("Start tuning MAXENT.Tsuruoka\n")
  #   try(tune.MAXENT.Tsuruoka <- as.data.frame(tune.maxent(bm.format@data.env.var,bm.format@data.species,nfold=ME.kfolds,showall=T)))
  #   cat(paste("Finished tuning MAXENT.Tsuruoka\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  #   
  #   if(!is.null(tune.MAXENT.Tsuruoka)){
  #     bm.options@MAXENT.Tsuruoka$l1_regularizer <- tune.MAXENT.Tsuruoka$l1_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     bm.options@MAXENT.Tsuruoka$l2_regularizer <- tune.MAXENT.Tsuruoka$l2_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     bm.options@MAXENT.Tsuruoka$use_sgd <- ifelse(tune.MAXENT.Tsuruoka[which.max(tune.MAXENT.Tsuruoka$accuracy),]$use_sgd==0,F,T)
  #     bm.options@MAXENT.Tsuruoka$set_heldout <- tune.MAXENT.Tsuruoka$set_heldout[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #   } else { if('MAXENT.Tsuruoka' %in% models){cat("Tuning MAXENT.Tsuruoka failed!"); tune.MAXENT.Tsuruoka <- "FAILED"}}
  # }  
  
  .bm_cat("Done")
  return(list(models.options = bm.options, tune.SRE = tune.SRE,  tune.CTA.rpart = tune.CTA.rpart, tune.CTA.rpart2 = tune.CTA.rpart2,
              tune.RF = tune.RF, tune.ANN = tune.ANN,  tune.MARS = tune.MARS, tune.FDA = tune.FDA, tune.GBM = tune.GBM,
              tune.GAM = tune.GAM, tune.GLM = tune.GLM, tune.MAXENT.Phillips = tune.MAXENT.Phillips))
  # tune.MAXENT.Tsuruoka = tune.MAXENT.Tsuruoka, 
}


###################################################################################################
#### Modified tuning function from the ENMeval package to tune MAXENT.Phillips (internal function for BIOMOD_tuning)

.maxent_tuning <- function(pres,
                           bg,
                           method,
                           kfolds,
                           clamp,
                           parallel,
                           numCores,
                           categoricals,
                           tune.args = list(rm = seq(0.5, 1, 0.5), fc = c("L")))
{
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


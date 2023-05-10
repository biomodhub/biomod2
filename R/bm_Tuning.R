# BIOMOD_Tuning Documentation -------------------------------------------------
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
##' \code{FDA}, \code{MARS}, \code{RF}, \code{MAXENT}
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
##' optimal model and tune parameters for \code{MAXENT}, must be either 
##' \code{auc.val.avg}, \code{auc.diff.avg}, \code{or.mtp.avg}, \code{or.10p.avg} or \code{AICc}
##' @param ME.cvmethod a \code{character} corresponding to the method used to partition data for 
##' \code{MAXENT}, \cr must be \code{randomkfold}
##' @param ME.kfolds an \code{integer} corresponding to the number of bins for k-fold 
##' cross-validation for \code{MAXENT}
##' @param ME.overlap (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether to calculate pairwise metric of niche overlap or not 
##' (see \code{\link[ENMeval]{calc.niche.overlap}})
##' @param ME.clamp (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether \emph{Features are constrained to remain within the 
##' range of values in the training data} (Elith et al. 2011) or not
##' @param ME.n.bg an \code{integer} corresponding to the number of background points used to run 
##' \code{MAXENT}
##' @param ME.env a \code{\link[terra:rast]{SpatRaster}} object 
##' containing model predictor variables
##' @param ME.parallel (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether to enable parallel computing for 
##' \code{MAXENT} or not
##' @param ME.numCores an \code{integer} corresponding to the number of cores to be used to 
##' train \code{MAXENT}
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
##'   \item For more details on \code{MAXENT} tuning, please refer to 
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
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' 
##' # --------------------------------------------------------------- #
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
#------------------------------------------------------------------------------#


BIOMOD_Tuning <- function(bm.format,
                          bm.options = BIOMOD_ModelingOptions(),
                          models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT'),
                          metric.eval = 'ROC', ## to correct
                          ctrl.train = NULL,
                          tuning.length = 30, ## to correct
                          ME.cvmethod = 'randomkfold',
                          ME.kfolds = 10,
                          
                          ANN.method = 'avNNet',
                          ANN.decay.tune = c(0.001, 0.01, 0.05, 0.1),
                          ANN.size.tune = c(2, 4, 6, 8),
                          ANN.maxit = 500,
                          ANN.MaxNWts = 10 * (ncol(bm.format@data.env.var) + 1) + 10 + 1,
                          GLM.method = 'glmStepAIC',
                          GLM.type = c('simple', 'quadratic', 'polynomial', 's_smoother'), 
                          GLM.interaction = c(0, 1),
                          weights = NULL)
{
  .bm_cat("Tune Modeling Options")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  # args <- .bm_Tuning.check.args(model = model, bm.format = bm.format, weights = weights)
  # for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  # rm(args)
  .bm_Tuning.check.args(model = model, bm.format = bm.format, weights = weights)
  
  if (!(model %in% c("MAXENT", "SRE"))) {
    ctrl.train <- caret::trainControl(method = "repeatedcv",
                                      repeats = 3,
                                      number = 10,
                                      summaryFunction = caret::twoClassSummary,
                                      classProbs = TRUE,
                                      # savePredictions = "all", ## for GLM ??
                                      returnData = FALSE)
  }
  
  
  ## create dataset -----------------------------------------------------------
  # myResp = bm.format@data.species
  # myExpl = bm.format@data.env.var
  mySpExpl <- get_species_data(bm.format)
  myResp = mySpExpl[, 1]
  myExpl = mySpExpl[, 4:ncol(mySpExpl)]
  if (model == "FDA") myResp = factor(myResp)
  if (model != "SRE" && metric.eval %in% c("ROC", "TSS")) myResp <- as.factor(ifelse(myResp == 1 & !is.na(myResp), "Presence", "Absence"))
  
  ## check weights ------------------------------------------------------------
  ## FOR FDA and CTA only ??
  if (is.null(weights)) { weights = rep(1, length(myResp)) }
  
  ## check tuning.length ------------------------------------------------------
  # ctrl.mod <- ctrl.train ## if no control given for the specific model
  tuning.length <- 1
  if (model == "CTA") tuning.length <- 30
  if (model == "RF") tuning.length <- min(30, ncol(myExpl))
  
  
  
  ## 2. Run tuning --------------------------------------------------------------------------------
  if ("MAXENT" %in% models) { 
    try(tune.MAXENT <- ENMeval::ENMevaluate(occs = mySpExpl[mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), ],
                                            bg = mySpExpl[mySpExpl[, 1] == 0 | is.na(mySpExpl[, 1]), ],
                                            tune.args = list(rm = seq(0.5, 1, 0.5), fc = c("L")),
                                            algorithm = "maxent.jar",
                                            partitions = ME.cvmethod,
                                            partition.settings = list(kfolds = ME.kfolds),
                                            doClamp = TRUE, ## allow to change or not ?
                                            parallel = TRUE,
                                            numCores = nb.cpu, ## default to 1 or NULL (all available cores used then) ?
                                            categoricals = NULL))
    
    if (!is.null(tune.MAXENT)) {
      if (ME.metric == 'auc.val.avg') {
        tmp = which.max(tune.MAXENT@results[, ME.metric])
      } else {
        tmp = which.min(tune.MAXENT@results[, ME.metric])
      }
      bm.options@MAXENT$linear <- grepl("L", tune.MAXENT@results[tmp, "fc"])
      bm.options@MAXENT$quadratic <- grepl("Q", tune.MAXENT@results[tmp, "fc"])
      bm.options@MAXENT$hinge <- grepl("H", tune.MAXENT@results[tmp, "fc"])
      bm.options@MAXENT$product <- grepl("P", tune.MAXENT@results[tmp, "fc"])
      bm.options@MAXENT$threshold <- grepl("T", tune.MAXENT@results[tmp, "fc"])
      bm.options@MAXENT$betamultiplier <- tune.MAXENT@results[tmp, "rm"]
    }
  } else if ("SRE" %in% models) {
    SRE.quant = c(0, 0.0125, 0.025, 0.05, 0.1)
    tune.SRE = foreach(rep = 1:ctrl.train$repeats, .combine = "rbind") %do%
      {
        fold <- dismo::kfold(myResp, by = myResp, k = ctrl.train$number) ## by = to keep prevalence
        RES = foreach (quant = SRE.quant, .combine = "rbind") %:%
          foreach (i = 1:ctrl.train$number, .combine = "rbind") %do%
          {
            DATA <- cbind(1:sum(fold == i)
                          , myResp[fold == i]
                          , bm_SRE(resp.var = myResp[fold != i],
                                   expl.var = myExpl[fold != i, ],
                                   new.env = myExpl[fold == i, ],
                                   quant = quant,
                                   do.extrem = FALSE))
            RES = presence.absence.accuracy(DATA, threshold = as.vector(
              optimal.thresholds(DATA, opt.methods = 3)[2], mode = "numeric"))
            return(data.frame(RES, quant = quant))
          }
        return(data.frame(RES, rep = rep))
      }
    tune.SRE$TSS <- tune.SRE$sensitivity + tune.SRE$specificity - 1
    tmp <- aggregate(tune.SRE[, c("sensitivity", "specificity", "Kappa", "AUC", "TSS")], by = list(quant = tune.SRE$quant), mean)
    bm.options@SRE$quant <- tmp[which.max(tmp[, metric.eval]), "quant"]
  } else {
    
    ## Parameters to be tested
    all.fun <- c('avNNet', 'rpart', 'rpart2', 'fda', 'gamSpline', 'bam', 'gam', 'gbm', 'glm', 'earth', 'rf', 'xgbTree')
    all.params <- foreach (fi = all.fun) %do% {
      params <- getModelInfo(model = fi)
      return(list(pkg = params[[fi]]$library, params = params[[fi]]$parameters$parameter))
    }
    names(all.params) <- all.fun
    train.params <- all.params[[tuning.fun]]
    
    
    tuning.grid = NULL
    ## expand.grid with different parameter values to be tested TODO
    # ## GBM
    # tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
    #                          .n.trees = c(500, 1000, 2500),
    #                          .shrinkage = c(0.001, 0.01, 0.1),
    #                          .n.minobsinnode = 10)
    # tune.grid <- expand.grid(.interaction.depth = c(tune.GBM$bestTune$interaction.depth - 1
    #                                                 , tune.GBM$bestTune$interaction.depth
    #                                                 , tune.GBM$bestTune$interaction.depth + 1),
    #                          .n.trees = n.trees,
    #                          .shrinkage = c(tune.GBM$bestTune$shrinkage / 2,
    #                                         tune.GBM$bestTune$shrinkage,
    #                                         tune.GBM$bestTune$shrinkage * 5), 
    #                          .n.minobsinnode = 10)
    # ## ANN
    ANN.decay.tune = c(0.001, 0.01, 0.05, 0.1)
    ANN.size.tune = c(2, 4, 6, 8)
    ANN.maxit = 500
    ANN.MaxNWts = 10 * (ncol(myExpl) + 1) + 10 + 1
    tuning.grid <- expand.grid(.decay = c(0.001, 0.01, 0.05, 0.1),
                               .size = c(2, 4, 6, 8),
                               .bag = FALSE)
    sizetmp <- c(2, 4, 6, 8)
    decaytmp <- c(0.001, 0.01, 0.05, 0.1)
    maxittmp <- 200
    nbCVtmp <- 5
    ## do cross validation test to find the optimal values of size and decay parameters (prevent from overfitting)
    CV_nnet <- bm_CVnnet(Input = myExpl, 
                         Target = myResp,
                         size = sizetmp,
                         decay = decaytmp,
                         maxit = maxittmp,
                         nbCV = nbCVtmp,
                         weights = weights)
                         # seedval = seed.val)
    ## MARS
    tuning.grid <- expand.grid(.degree = 1:2, .nprune = 2:max(38, 2 * ncol(myExpl) + 1))
    ## FDA
    tuning.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    ## GLM
    GLM.type = c('simple', 'quadratic', 'polynomial', 's_smoother')
    GLM.interaction = c(0, 1)
    
    
    if (model == "GLM") {
      RES = foreach (typ = GLM.type, .combine = "rbind") %:%
        foreach (intlev = GLM.interaction.level, .combine = "rbind") %do%
        {
          tuned.mod <- caret::train(formula = bm_MakeFormula(resp.name = bm.format@sp.name,
                                                             expl.var = myExpl,
                                                             type = typ,
                                                             interaction.level = intlev),
                                    # data = cbind(myExpl, resp = factor(myResp, c("Presence", "Absence"))),
                                    data = tot, #[, -c(2, 3)],
                                    method = tuning.fun,
                                    # preProcess = c('center', 'scale'),
                                    metric = "ROC",
                                    trControl = ctrl.train)
          return(tuned.mod$results)
        }
      
      
      # fm <- list()
      # GLM.results = foreach (type = GLM.type, .combine = "rbind") %:%
      #   foreach (IA = GLM.interaction, .combine = "rbind") %do%
      #   {
      #     try(tune.GLM <- caret::train(bm_MakeFormula(resp.name = "resp",
      #                                                 expl.var = bm.format@data.env.var,
      #                                                 type = type,
      #                                                 interaction.level = IA),
      #                                  data = cbind(bm.format@data.env.var, resp = resp),
      #                                  method = GLM.method,
      #                                  trControl = ctrl.GLM))
      #     # weights = weights))
      #     try(fm[[length(fm) + 1]] <- formula(tune.GLM$finalModel))
      #     try(RES <- cbind(tune.GLM$results, il = IA, type = type))
      #     return(RES)
      #   }
      # tune.GLM <- GLM.results
    } else {
      try(tuned.mod <- caret::train(x = myExpl, 
                                    y = myResp,
                                    method = tuning.fun,
                                    tuneGrid = tuning.grid, ## null for RF, GAM
                                    tuneLength = tuning.length,
                                    trControl = ctrl.train,
                                    metric = metric.eval, ## RF
                                    verbose = FALSE, ## RF (rm for earth, bam, fda)
                                    weights = weights))
      
      try(tuned.mod <- caret::train(x = myExpl, 
                                    y = myResp,
                                    method = tuning.fun,
                                    tuneGrid = tuning.grid, ## null for RF, GAM
                                    # tuneLength = tuning.length,
                                    trControl = ctrl.train,
                                    metric = metric.eval, ## RF
                                    # verbose = FALSE, ## RF (rm for earth, bam, fda)
                                    weights = weights,
                                    ## Automatically standardize data prior to modeling and prediction
                                    preProc = c("center", "scale"),
                                    linout = TRUE,
                                    trace = FALSE,
                                    MaxNWts.ANN = ANN.MaxNWts,
                                    maxit = ANN.maxit))
    }
    
    ## RF, MARS
    if (!is.null(tuned.mod)) {
      tmp = tuned.mod$results
      tmp$TSS = tmp$Sens + tmp$Spec - 1
      
      # bm.options@RF$mtry <- tmp[which.max(tmp[, metric.eval]), "mtry"]
      for (param in all.params[[tuning.fun]]$params) {
        bm.options@RF[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
      }
    }
  }
  
  
  if ('GBM' %in% models)
  {  
    # tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
    #                          .n.trees = c(500, 1000, 2500),
    #                          .shrinkage = c(0.001, 0.01, 0.1),
    #                          .n.minobsinnode = 10)
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
      
      # tune.grid <- expand.grid(.interaction.depth = c(tune.GBM$bestTune$interaction.depth - 1
      #                                                 , tune.GBM$bestTune$interaction.depth
      #                                                 , tune.GBM$bestTune$interaction.depth + 1),
      #                          .n.trees = n.trees,
      #                          .shrinkage = c(tune.GBM$bestTune$shrinkage / 2,
      #                                         tune.GBM$bestTune$shrinkage,
      #                                         tune.GBM$bestTune$shrinkage * 5), 
      #                          .n.minobsinnode = 10)
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
    }
  }
  
  if ('GLM' %in% models)
  {
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
  }      
  
  if ('CTA' %in% models)
  {
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
    }
    
    if (!is.null(tune.CTA.rpart2)) {
      if (metric.eval == 'TSS') {
        tmp = which.max(apply(tune.CTA.rpart2$results[, c("Sens", "Spec")], 1, sum) - 1)
        bm.options@CTA$control$maxdepth <- tune.CTA.rpart2$results[tmp, "maxdepth"]
      } else {
        bm.options@CTA$control$maxdepth <- tune.CTA.rpart2$bestTune
      }
    }
  }
  
  .bm_cat("Done")
}



# ---------------------------------------------------------------------------- #

.bm_Tuning.check.args <- function(model, bm.format, weights = NULL)
{
  ## check namespace ----------------------------------------------------------
  if (!(model %in% c("MAXENT", "SRE"))) {
    if (!isNamespaceLoaded("caret")) { 
      if(!requireNamespace('caret', quietly = TRUE)) stop("Package 'caret' not found")
    }
  } else if (model == "MAXENT" && !isNamespaceLoaded('ENMeval')) { 
    if(!requireNamespace('ENMeval', quietly = TRUE)) stop("Package 'ENMeval' not found")
  } else if (model == "SRE" && !isNamespaceLoaded('dismo')) { 
    if(!requireNamespace('dismo', quietly = TRUE)) stop("Package 'dismo' not found")
  }
  if (model == "GLM" && "s_smoother" %in% GLM.type) { 
    if(!requireNamespace('gam', quietly = TRUE)) stop("Package 'gam' not found")
  }
  
  ## check bm.format ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  
  ## check evaluation metric --------------------------------------------------
  if (model == "MAXENT") {
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("auc.val.avg", "auc.diff.avg", "or.mtp.avg", "or.10p.avg", "AICc"))
  } else if (model == "SRE") {
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("AUC", "Kappa", "TSS"))
  } else {
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("ROC", "TSS")) ## TO CHECK !!!!!
  }
  
  # ## create dataset -----------------------------------------------------------
  # myResp = bm.format@data.species
  # myExpl = bm.format@data.env.var
  # if (model == "FDA") myResp = factor(myResp)
  # if (metric.eval %in% c("ROC", "TSS")) myResp <- as.factor(ifelse(myResp == 1 & !is.na(myResp), "Presence", "Absence"))
  # 
  # ## check weights ------------------------------------------------------------
  # if (is.null(weights)) { weights = rep(1, length(myResp))}
  # 
  # ## check tuning.length ------------------------------------------------------
  # # ctrl.mod <- ctrl.train ## if no control given for the specific model
  # tuning.length <- 1
  # if (model == "CTA") tuning.length <- 30
  # if (model == "RF") tuning.length <- min(30, ncol(myExpl))
  
  
  # return(list(myResp = myResp
  #             , myExpl = myExpl
  #             , weights = weights
  #             , tuning.length = tuning.length))
}


# bm_Tuning Documentation -------------------------------------------------
##' @name bm_Tuning
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
##' @param metric.eval a \code{character} corresponding to the evaluation metric used to select 
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
##'   bm.tuning <- bm_Tuning(bm.format = myBiomodData, ME.env = myExpl, ME.n.bg = ncell(myExpl))
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


bm_Tuning <- function(bm.format,
                      bm.opt.def,
                      model,
                      metric.eval = 'TSS',
                      weights = NULL,
                      ctrl.train = NULL,
                      params.train = list(ANN.size = c(2, 4, 6, 8),
                                          ANN.decay = c(0.001, 0.01, 0.05, 0.1),
                                          ANN.bag = FALSE, 
                                          FDA.degree = 1:2, 
                                          FDA.nprune = 2:38,
                                          GBM.n.trees = c(500, 1000, 2500),
                                          GBM.interaction.depth = seq(2, 8, by = 3),
                                          GBM.shrinkage = c(0.001, 0.01, 0.1),
                                          GBM.n.minobsinnode = 10,
                                          MARS.degree = 1:2, 
                                          MARS.nprune = 2:max(38, 2 * ncol(bm.format@data.env.var) + 1),
                                          SRE.quant = c(0, 0.0125, 0.025, 0.05, 0.1),
                                          XGBOOST.nrounds = 50,
                                          XGBOOST.max_depth = 1,
                                          XGBOOST.eta = c(0.3, 0.4),
                                          XGBOOST.gamma = 0,
                                          XGBOOST.colsample_bytree = c(0.6, 0.8),
                                          XGBOOST.min_child_weight = 1,
                                          XGBOOST.subsample = 0.5),
                      ## ADD : mtry for rf ? select + method for bam / gam ?
                      ANN.maxit = 500,
                      ANN.MaxNWts = 10 * (ncol(bm.format@data.env.var) + 1) + 10 + 1,
                      ME.cvmethod = 'randomkfold',
                      ME.kfolds = 10)
{
  .bm_cat("Tune Modeling Options")
  
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_Tuning.check.args(model = model, bm.format = bm.format, metric.eval = metric.eval, weights = weights)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## LOOP OVER CALIB LINES 
  if (is.null(calib.lines)) {
    calib.lines = data.frame("_allData_allRun" = rep(TRUE, length(bm.format@data.species)))
  }
  
  argsval <- foreach(calib.i = 1:ncol(calib.lines)) %do%
    {
      argstmp <- BOM@args.default
      
      ## 1. SPECIFIC CASE OF MAXENT OR SRE ------------------------------------------------------------
      if (model %in% c("MAXENT", "SRE")) {
        
        ## create dataset ---------------------------------------------------------
        mySpExpl <- get_species_data(bm.format)
        mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
        myResp = mySpExpl[, 1]
        myExpl = mySpExpl[, 4:ncol(mySpExpl)]
        
        
        if (model == "MAXENT") { # --------------------------------------------------------------------
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
            if (metric.eval == 'auc.val.avg') {
              tmp = which.max(tune.MAXENT@results[, metric.eval])
            } else {
              tmp = which.min(tune.MAXENT@results[, metric.eval])
            }
            argstmp$linear <- grepl("L", tune.MAXENT@results[tmp, "fc"])
            argstmp$quadratic <- grepl("Q", tune.MAXENT@results[tmp, "fc"])
            argstmp$hinge <- grepl("H", tune.MAXENT@results[tmp, "fc"])
            argstmp$product <- grepl("P", tune.MAXENT@results[tmp, "fc"])
            argstmp$threshold <- grepl("T", tune.MAXENT@results[tmp, "fc"])
            argstmp$betamultiplier <- tune.MAXENT@results[tmp, "rm"]
          }
        } else if (model == "SRE") { # ----------------------------------------------------------------
          tune.SRE = foreach(rep = 1:ctrl.train$repeats, .combine = "rbind") %do%
            {
              fold <- dismo::kfold(myResp, by = myResp, k = ctrl.train$number) ## by = to keep prevalence
              RES = foreach (quant = params.train$SRE.quant, .combine = "rbind") %:%
                foreach (i = 1:ctrl.train$number, .combine = "rbind") %do%
                {
                  DATA <- cbind(1:sum(fold == i)
                                , myResp[fold == i]
                                , bm_SRE(resp.var = myResp[fold != i],
                                         expl.var = myExpl[fold != i, ],
                                         new.env = myExpl[fold == i, ],
                                         quant = quant,
                                         do.extrem = FALSE))
                  RES = presence.absence.accuracy(DATA, threshold = as.vector(optimal.thresholds(DATA, opt.methods = 3)[2], mode = "numeric"))
                  return(data.frame(RES, quant = quant))
                }
              return(data.frame(RES, rep = rep))
            }
          tune.SRE$TSS <- tune.SRE$sensitivity + tune.SRE$specificity - 1
          tmp <- aggregate(tune.SRE[, c("sensitivity", "specificity", "Kappa", "AUC", "TSS")], by = list(quant = tune.SRE$quant), mean)
          argstmp$quant <- tmp[which.max(tmp[, metric.eval]), "quant"]
        }
        
      } else {
        ## 2. ALL OTHER MODELS ------------------------------------------------------------------------
        
        ## create dataset ---------------------------------------------------------
        mySpExpl <- get_species_data(bm.format)
        mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
        myResp = mySpExpl[, 1]
        myExpl = mySpExpl[, 4:ncol(mySpExpl)]
        myResp <- as.factor(ifelse(myResp == 1 & !is.na(myResp), "Presence", "Absence"))
        
        
        ## check control ----------------------------------------------------------
        ctrl.train <- caret::trainControl(method = "repeatedcv",
                                          repeats = 3,
                                          number = 10,
                                          summaryFunction = caret::twoClassSummary,
                                          classProbs = TRUE,
                                          returnData = FALSE)
        
        
        ## run tuning ---------------------------------------------------------------------------------
        if (model == "GLM") {
          # RES = foreach (typ = GLM.type, .combine = "rbind") %:%
          #   foreach (intlev = GLM.interaction.level, .combine = "rbind") %do%
          #   {
          #     tuned.mod <- caret::train(form = bm_MakeFormula(resp.name = "resp",
          #                                                     expl.var = myExpl,
          #                                                     type = typ,
          #                                                     interaction.level = intlev),
          #                               data = cbind(myExpl, resp = myResp),
          #                               method = tuning.fun,
          #                               metric = metric.eval,
          #                               trControl = ctrl.train)
          #     if (!is.null(tuned.mod)) {
          #       tmp = tuned.mod$results
          #       tmp$TSS = tmp$Sens + tmp$Spec - 1
          #       # formu = as.character(tuned.mod$finalModel$formula)
          #       # formu = paste0(formu[c(2,1,3)], collapse = " ")
          #       formu = tuned.mod$coefnames
          #       formu = paste0(bm.format@sp.name, " ~ 1 + ", paste0(formu, collapse = " + "))
          #       return(data.frame(tmp, type = typ, interaction.level = intlev, formula = formu))
          #     }
          #   }
          # for (param in train.params$params) { ## must be type, interaction.level, formula
          #   argstmp[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
          # }
        } else {
          if (model == "ANN") {
            try(tuned.mod <- caret::train(x = myExpl,
                                          y = myResp,
                                          method = tuning.fun,
                                          tuneGrid = tuning.grid,
                                          trControl = ctrl.train,
                                          metric = metric.eval,
                                          weights = weights,
                                          ## Automatically standardize data prior to modeling and prediction
                                          preProc = c("center", "scale"),
                                          linout = TRUE,
                                          trace = FALSE,
                                          MaxNWts.ANN = ANN.MaxNWts,
                                          maxit = ANN.maxit))
          } else if (tuning.fun %in% c("earth", "bam", "fda", "rpart")) { ## remove verbose
            try(tuned.mod <- caret::train(x = myExpl, 
                                          y = myResp,
                                          method = tuning.fun,
                                          tuneGrid = tuning.grid,
                                          tuneLength = tuning.length,
                                          trControl = ctrl.train,
                                          metric = metric.eval,
                                          weights = weights))
            
            if (model == "CTA" && !is.null(tuned.mod)) {
              tmp = tuned.mod$results
              tmp$TSS = tmp$Sens + tmp$Spec - 1
              argstmp[["cp"]] <- tmp[which.max(tmp[, metric.eval]), "cp"]
              
              try(tuned.mod <- caret::train(x = myExpl, 
                                            y = myResp,
                                            method = "rpart2",
                                            tuneGrid = tuning.grid,
                                            tuneLength = tuning.length,
                                            trControl = ctrl.train,
                                            metric = metric.eval,
                                            weights = weights))
              if (!is.null(tuned.mod)) {
                tmp = tuned.mod$results
                tmp$TSS = tmp$Sens + tmp$Spec - 1
                argstmp[["maxdepth"]] <- tmp[which.max(tmp[, metric.eval]), "maxdepth"]
              }
            }
          } else {
            try(tuned.mod <- caret::train(x = myExpl, 
                                          y = myResp,
                                          method = tuning.fun,
                                          tuneGrid = tuning.grid,
                                          tuneLength = tuning.length,
                                          trControl = ctrl.train,
                                          metric = metric.eval,
                                          verbose = FALSE,
                                          weights = weights))
          }
          
          ## GET tuned parameter values ---------------------------------------------------------------
          if (!is.null(tuned.mod)) {
            tmp = tuned.mod$results
            tmp$TSS = tmp$Sens + tmp$Spec - 1
            for (param in train.params$params) {
              argstmp[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
            }
          }
        }
      }
      return(argstmp)
    }
  names(argsval) <- colnames(calib.lines)
  
  .bm_cat("Done")
  return(argsval)
}



# ---------------------------------------------------------------------------- #

TABLE_MODELS <- data.frame(model = c('ANN', 'CTA', 'FDA', 'GAM', 'GAM', 'GAM', 'GBM', 'GLM'
                                     , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
                           , type = 'binary'
                           , package = c('nnet', 'rpart', 'mda', 'gam', 'mgcv', 'mgcv', 'gbm', 'stats'
                                         , 'earth', 'MAXENT', 'maxnet', 'randomForest', 'biomod2', 'xgboost')
                           , func = c('nnet', 'rpart', 'fda', 'gam', 'bam', 'gam', 'gbm', 'glm'
                                      , 'earth', 'MAXENT', 'maxnet', 'randomForest', 'bm_SRE', 'xgboost')
                           , train = c('avNNet', 'rpart', 'fda', 'gamSpline', 'bam', 'gam', 'gbm', 'glm'
                                       , 'earth', 'ENMevaluate', '', 'rf', 'bm_SRE', 'xgbTree'))

.bm_Tuning.check.args <- function(model, bm.format, metric.eval, weights = NULL)
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
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("ROC", "TSS"))
  }
  
  ## check weights ------------------------------------------------------------
  if (model %in% c("CTA", "FDA") && is.null(weights)) { weights = rep(1, length(bm.format@data.species)) }
  
  
  
  ## get tuning function and parameters ---------------------------------------
  all.fun <- c('avNNet', 'rpart', 'rpart2', 'fda', 'gamSpline', 'bam', 'gam', 'gbm', 'glm', 'earth', 'rf', 'xgbTree')
  all.params <- foreach (fi = all.fun) %do% {
    params <- getModelInfo(model = fi)
    return(list(pkg = params[[fi]]$library, params = params[[fi]]$parameters$parameter))
  }
  names(all.params) <- all.fun
  # tuning.fun <- TABLE_MODELS$train[which(TABLE_MODELS$model == model)]
  train.params <- all.params[[tuning.fun]]
  
  ## get tuning grid through params.train -------------------------------------
  tuning.grid <- NULL
  if (model %in% c("ANN", "FDA", "GBM", "MARS", "XGBOOST")) {
    params.train = params.train[grep(model, names(params.train))]
    .fun_testIfIn(TRUE, "names(params.train)", names(params.train), paste0(model, ".", train.params$params))
    names(params.train) = sub(model, "", names(params.train))
    tuning.grid <- do.call(expand.grid, params.train)
  }
  
  ## get tuning length --------------------------------------------------------
  tuning.length <- 1
  if (model == "CTA") tuning.length <- 30
  if (model == "RF") tuning.length <- min(30, ncol(bm.format@data.env.var))
  
  
  return(list(tuning.fun = tuning.fun
              , train.params = train.params
              , tuning.length = tuning.length
              , tuning.grid = tuning.grid))
}


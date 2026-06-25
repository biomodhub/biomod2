###################################################################################################
##' @name bm_Tuning
##' @author Frank Breiner, Maya Guéguen, Hélène Blancheteau
##' 
##' @title Tune models parameters
##' 
##' @description This internal \pkg{biomod2} function allows to tune single model parameters and 
##' select more efficient ones based on an evaluation metric.
##' 
##'
##' @param model a \code{character} corresponding to the  algorithm to be tuned, must be either 
##' \code{ANN}, \code{CTA}, \code{DNN}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{RFd}, \code{SRE}, \code{XGBOOST}
##' @param tuning.fun a \code{character} corresponding to the model function name to be called 
##' through \code{\link[caret]{train}} function for tuning parameters (see \code{\link{ModelsTable}} 
##' dataset)
##' @param do.formula (\emph{optional, default} \code{FALSE}) \cr  
##' A \code{logical} value defining whether formula is to be optimized or not
##' @param do.stepAIC (\emph{optional, default} \code{FALSE}) \cr  
##' A \code{logical} value defining whether variables selection is to be performed for 
##' \code{GLM} and \code{GAM} models or not
##' @param bm.options a \code{\link{BIOMOD.options.default}} or \code{\link{BIOMOD.options.dataset}} 
##' object returned by the \code{\link{bm_ModelingOptions}} function
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' A \code{data.frame} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{bm_CrossValidation}} functions
##' @param metric.eval a \code{character} corresponding to the evaluation metric to be used, must 
##' be either \code{AUC}, \code{Kappa} or \code{TSS} for \code{SRE} only ; \code{auc.val.avg}, 
##' \code{auc.diff.avg}, \code{or.mtp.avg}, \code{or.10p.avg}, \code{AICc} for \code{MAXENT} only ; 
##' \code{ROC} or \code{TSS} for all other models
##' @param metric.AIC a \code{character} corresponding to the AIC metric to be used, must 
##' be either \code{AIC} or \code{BIC}
##' @param weights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to observation weights (one per 
##' observation, see Details)
##' @param ctrl.train (\emph{optional, default} \code{NULL}) \cr 
##' A \code{\link[caret]{trainControl}} object
##' @param params.train a \code{list} containing values of model parameters to be tested 
##' (see Details)
##' 
##' 
##' @return 
##' 
##' A \code{\link{BIOMOD.models.options}} object (see \code{\link{bm_ModelingOptions}}) with 
##' optimized parameters
##' 
##' 
##' @details 
##' 
##' \bold{Concerning \code{ctrl.train} parameter :}
##' 
##' Set by default to : \cr
##' 
##' \code{ctrl.train <- caret::trainControl(method = "repeatedcv", repeats = 3, number = 10,} \cr
##' \code{summaryFunction = caret::twoClassSummary,} \cr
##' \code{classProbs = TRUE, returnData = FALSE)} \cr \cr
##' 
##' 
##' \bold{Concerning \code{params.train} parameter :}
##' 
##' All elements of the \code{list} must have names matching \code{model.parameter_name} format, 
##' \code{parameter_name} being one of the parameter of the \code{tuning.fun} function called by 
##' \code{caret} package and that can be found through the \code{\link[caret]{getModelInfo}} 
##' function.
##' 
##' Currently, the available parameters to be tuned are the following :
##' \describe{
##'   \item{ANN}{\code{size}, \code{decay}, \code{bag}}
##'   \item{CTA}{\code{maxdepth}}
##'   \item{DNN}{\code{hidden} , \code{bias}, \code{lambda}, \code{alpha}, \code{lr}, \code{batchsize}, \code{150}}
##'   \item{FDA}{\code{degree}, \code{nprune}}
##'   \item{GAM.gam}{\code{span}, \code{degree}}
##'   \item{GAM.mgcv}{\code{select}, \code{method}}
##'   \item{GBM}{\code{n.trees}, \code{interaction.depth}, \code{shrinkage}, \code{n.minobsinnode}}
##'   \item{MARS}{\code{degree}, \code{nprune}}
##'   \item{MAXENT}{\code{algorithm},\code{tune.args}, \code{parallel}, \code{partitions}, \code{kfolds}, 
##'   \code{user.grp}}
##'   \item{RF}{\code{mtry}}
##'   \item{RFd}{\code{mtry}}
##'   \item{SRE}{\code{quant}}
##'   \item{XGBOOST}{\code{nrounds}, \code{max_depth}, \code{eta}, \code{gamma}, 
##'   \code{colsampl_bytree}, \code{min_child_weight}, \code{subsample}}
##' }
##' 
##' 
##' The \code{\link{expand.grid}} function is used to build a \code{matrix} containing all 
##' combinations of parameters to be tested.
##' 
##' @note 
##' \itemize{
##'   \item No tuning for \code{GLM} and \code{MAXNET}
##'   \item \code{MAXENT} is tuned through \code{\link[ENMeval]{ENMevaluate}} function which is
##'   calling either :
##'   \itemize{
##'     \item maxnet (by defining \code{MAXENT.algorithm = 'maxnet'}) (\emph{default})
##'     \item Java version of Maxent defined in \pkg{dismo} package (by defining 
##'     \code{MAXENT.algorithm = 'maxent.jar'})
##'   }
##'   \item \code{DNN} is tuned through \code{\link[cito]{tune}} function. 
##'   The values include in \code{params.train} are the lower or upper range which hyperparameters are sampled.
##'   If there is only one value, the hyperparameter is fixed by biomod2 (inclunding the width and depth of \code{hidden} parameters.) 
##'   \item \code{SRE} is tuned through \code{\link{bm_SRE}} function
##'   \item All other models are tuned through \code{\link[caret]{train}} function
##'   \item No optimization of formula for \code{DNN}, \code{MAXENT}, \code{MAXNET}, \code{SRE} and 
##'   \code{XGBOOST}
##'   \item No interaction included in formula for \code{CTA}
##'   \item Variables selection only for \code{GAM.gam} and \code{GLM}
##' }
##' 
##' 
##' @seealso \code{\link[caret]{trainControl}}, \code{\link[caret]{train}}, 
##' \code{\link[ENMeval]{ENMevaluate}}, 
##' \code{\link{ModelsTable}}, \code{\link{BIOMOD.models.options}}, 
##' \code{\link{bm_ModelingOptions}}, \code{\link{BIOMOD_Modeling}}
##' @family Secondary functions
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
##' myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
##'                                      resp.var = myResp,
##'                                      resp.xy = myRespXY,
##'                                      expl.var = myExpl)
##' 
##' 
##' # --------------------------------------------------------------- #
##' # List of all models currently available in `biomod2` (and their related package and function)
##' # Some of them can be tuned through the `train` function of the `caret` package 
##' # (and corresponding training function to be used is indicated)
##' data(ModelsTable)
##' ModelsTable
##' 
##' allModels <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM'
##'                , 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
##' 
##' # default parameters
##' opt.d <- bm_ModelingOptions(data.type = 'binary',
##'                             models = allModels,
##'                             strategy = 'default')
##'                             
##' # tune parameters for Random Forest model
##' tuned.rf <- bm_Tuning(model = 'RF',
##'                       tuning.fun = 'rf', ## see in ModelsTable
##'                       do.formula = FALSE,
##'                       bm.options = opt.d@options$RF.binary.randomForest.randomForest,
##'                       bm.format = myBiomodData)
##' tuned.rf
##' 
##' \dontrun{
##' # tune parameters for GAM (from mgcv package) model
##' tuned.gam <- bm_Tuning(model = 'GAM',
##'                        tuning.fun = 'gam', ## see in ModelsTable
##'                        do.formula = TRUE,
##'                        do.stepAIC = TRUE,
##'                        bm.options = opt.d@options$GAM.binary.mgcv.gam,
##'                        bm.format = myBiomodData)
##' tuned.gam
##' }                  
##' 
##' 
##' 
##' @importFrom foreach foreach %do% %:%
##' @importFrom stats aggregate formula 
##' @importFrom PresenceAbsence optimal.thresholds presence.absence.accuracy
##' @importFrom dplyr mutate_at
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


bm_Tuning <- function(model,
                      workflow,
                      bm.options,
                      bm.format,
                      calib.lines = NULL,
                      # tuning.fun,
                      do.tuning = TRUE,
                      do.formula = FALSE,
                      do.stepAIC = FALSE,
                      metric.eval = 'AUCroc',
                      metric.AIC = 'AIC',
                      weights = NULL,
                      ctrl.train = NULL,
                      params.train = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_Tuning.check.args(model = model, workflow = workflow, bm.options = bm.options, bm.format = bm.format
                                , do.tuning = do.tuning, do.formula = do.formula, do.stepAIC = do.stepAIC
                                , metric.eval = metric.eval, metric.AIC = metric.AIC
                                , weights = weights, params.train = params.train)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## LOOP OVER CALIB LINES 
  if (is.null(calib.lines)) {
    calib.lines <- data.frame(rep(TRUE, length(bm.format@data.species)))
    colnames(calib.lines) <- "_allData_allRun"
  }
  ## LOOP OVER PA DATASETS
  if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
    PA.lines <- colnames(bm.format@PA.table)
  } else {
    PA.lines <- "_allData_allRun"
  }
  ## LOOP OVER ALL COMBINED
  test_forPA <- sapply(PA.lines, function(xx) any(grepl(xx, colnames(calib.lines))))
  if (inherits(bm.format, "BIOMOD.formated.data.PA") && 
      sum(test_forPA) == length(PA.lines)) {
    combi <- expand.grid(PA = NA, calib = colnames(calib.lines), stringsAsFactors = FALSE)
    combi$PA <- sapply(combi$calib, function(ii) strsplit(ii, "_")[[1]][2])
    if (length(which(combi$calib == "_allData_allRun")) > 0) {
      combi$PA[which(combi$calib == "_allData_allRun")] <- "_allData_allRun"
    }
    combi$name_dataset <- combi$calib
  } else {
    combi <- expand.grid(PA = PA.lines, calib = colnames(calib.lines), stringsAsFactors = FALSE)
    combi$name_dataset <- sapply(1:nrow(combi), function(ii) {
      tmp1 <- combi$PA[ii]
      tmp2 <- combi$calib[ii]
      if (!grepl("PA", tmp2) || (grepl("PA", tmp2) && grepl(tmp1, tmp2))) {
        if (tmp1 == "_allData_allRun") tmp1 <- "allData"
        tmp2 <- strsplit(tmp2, "_")[[1]][3]
        return(paste0("_", tmp1, "_", tmp2))
      } else {
        return(NA)
      }
    })
    combi <- na.exclude(combi)
  }
  
  
  ## ----------------------------------------------------------------------------------------------
  
  argsval <- foreach(PA.i = combi$PA, calib.i = combi$calib, dataset.i = combi$name_dataset) %do%
    {
      cat("\n\t\t> Dataset", dataset.i)
      
      ## If some values already available (defined by user), keep them
      if (inherits(bm.options, "BIOMOD.options.dataset") && !is.null(bm.options@args.values[[dataset.i]])) {
        argstmp <- bm.options@args.values[[dataset.i]]
      } else {
        argstmp <- bm.options@args.default
      }
      
      ## create dataset -------------------------------------------------------
      mySpExpl <- get_species_data(bm.format)
      mySpExpl[["_allData_allRun"]] <- ifelse(model %in% c("DNN", "MAXENT", "SRE"), TRUE, 1) ## CONCAT
      current.weights <- weights[which(calib.lines[, calib.i] == TRUE & mySpExpl[, PA.i] == TRUE)] ## !(DNN, MAXENT, SRE)
      mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
      mySpExpl <- mySpExpl[which(mySpExpl[, PA.i] == TRUE), ]
      
      if (model == "MAXENT") { # ----------------------------------------------
        if (params.train$MAXENT.algorithm == "maxnet") {
          mySpExpl[["_allData_allRun"]] <- NULL
          mySpExpl[, 1] <- ifelse(mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), 1, 0)
          mySpExpl <- mySpExpl[, 1:(3 + ncol(bm.format@data.env.var))]
        }
      } else { # --------------------------------------------------------------
        if (!(model %in% c("DNN", "SRE")) && bm.format@data.type == "binary") { ## CONCAT
          mySpExpl[, 1] <- as.factor(ifelse(mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), "presence", "absence"))
        }
        myResp <- mySpExpl[, 1]
        myExpl <- mySpExpl[, 4:(3 + ncol(bm.format@data.env.var))]
        if (model == "SRE") {
          myResp <- sapply(myResp, function(xx) ifelse(xx == 0 || is.na(xx), 0, 1))
        }
      }
      
      ## run tuning -------------------------------------------------------------------------------
      
      if (do.tuning) {
        cat("\n\t\t\t> Tuning parameters...")
        
        if (model == "DNN") { ## DNN case  # -------------------------------------- #
          
          ## Preparation of data
          scale_data <- scale(myExpl)
          argstmp$data <- cbind(myResp, as.data.frame(scale_data))
          argstmp$formula <- as.formula("myResp ~.")
          
          ## Preparation of parameters to be tuned
          names_DNN <- grep(paste0(model,"\\."), names(params.train), value = TRUE)
          for (param.n in names_DNN) {
            real.name <- unlist(strsplit(param.n, split = "\\."))[2]
            argstmp[[real.name]] <- params.train[[param.n]]
          }
          argstmp$tuning <- cito::config_tuning(steps = 5)
          
          ## Tune model
          tune.DNN <- do.call(cito::dnn, argstmp)
          
          ## Keep the tuned parameters
          argstmp$hidden = tune.DNN$model_properties$hidden
          argstmp$bias = tune.DNN$model_properties$bias
          argstmp$lambda = tune.DNN$training_properties$lambda
          argstmp$alpha = tune.DNN$training_properties$alpha
          argstmp$lr = as.numeric(tune.DNN$training_properties$lr)
          argstmp$batchsize = tune.DNN$training_properties$batchsize
          argstmp$epochs = tune.DNN$training_properties$epochs
          
        } else if (model == "MAXENT") { # ------------------------------------------#
          
          try(tune.MAXENT <- ENMeval::ENMevaluate(occs = mySpExpl[mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), ],
                                                  bg = mySpExpl[mySpExpl[, 1] == 0 | is.na(mySpExpl[, 1]), ],
                                                  tune.args = params.train$MAXENT.tune.args,
                                                  algorithm = params.train$MAXENT.algorithm, ## ??
                                                  partitions = params.train$MAXENT.partitions,
                                                  partition.settings = list(kfolds = params.train$MAXENT.kfolds),
                                                  user.grp = params.train$MAXENT.user.grp,
                                                  doClamp = TRUE, ## allow to change or not ?
                                                  parallel = params.train$MAXENT.parallel,
                                                  numCores = NULL, ## default to 1 or NULL (all available cores used then) ?
                                                  categoricals = NULL))
          
          if (!is.null(tune.MAXENT)) {
            if (metric.eval == 'auc.val.avg') {
              tmp <- which.max(tune.MAXENT@results[, metric.eval])
            } else {
              tmp <- which.min(tune.MAXENT@results[, metric.eval])
            }
            argstmp$linear <- grepl("L", tune.MAXENT@results[tmp, "fc"])
            argstmp$quadratic <- grepl("Q", tune.MAXENT@results[tmp, "fc"])
            argstmp$hinge <- grepl("H", tune.MAXENT@results[tmp, "fc"])
            argstmp$product <- grepl("P", tune.MAXENT@results[tmp, "fc"])
            argstmp$threshold <- grepl("T", tune.MAXENT@results[tmp, "fc"])
            argstmp$betamultiplier <- tune.MAXENT@results[tmp, "rm"]
          }
          
        } else if (model == "SRE") { # -------------------------------------------- #
          
          tune.SRE <- foreach(rep = 1:ctrl.train$repeats, .combine = "rbind") %do%
            {
              fold <- dismo::kfold(myResp, by = myResp, k = ctrl.train$number) ## by = to keep prevalence
              RES <- foreach (quant = params.train$SRE.quant, .combine = "rbind") %:%
                foreach (i = 1:ctrl.train$number, .combine = "rbind") %do%
                {
                  DATA <- cbind(1:sum(fold == i)
                                , myResp[fold == i]
                                , bm_SRE(resp.var = myResp[fold != i],
                                         expl.var = myExpl[fold != i, ],
                                         new.env = myExpl[fold == i, ],
                                         quant = quant,
                                         do.extrem = FALSE))
                  thresh <- as.vector(optimal.thresholds(DATA, opt.methods = 3)[2], mode = "numeric")
                  RES <- presence.absence.accuracy(DATA, threshold = thresh)
                  return(data.frame(RES, quant = quant))
                }
              return(data.frame(RES, rep = rep))
            }
          tune.SRE$TSS <- tune.SRE$sensitivity + tune.SRE$specificity - 1
          tmp <- aggregate(tune.SRE[, c("sensitivity", "specificity", "Kappa", "AUC", "TSS")]
                           , by = list(quant = tune.SRE$quant), mean)
          argstmp$quant <- tmp[which.max(tmp[, metric.eval]), "quant"]
          
        } else if (workflow == "caret") {
          
          ## prepare tuning command -------------------------------------------
          cmd.tuning <- "caret::train(x = myExpl, y = myResp, method = tuning.fun, tuneGrid = tuning.grid,"
          if (bm.format@data.type == "binary") {
            cmd.tuning <- paste0(cmd.tuning, " trControl = ctrl.train, metric = 'ROC',")
          } else {
            cmd.tuning <- paste0(cmd.tuning, " trControl = ctrl.train, metric ='", metric.eval, "',")
          }
          
          if (tuning.fun %in% c("fda", "rpart")) { ## add weights
            cmd.tuning <- paste0(cmd.tuning, " weights = current.weights,")
          }
          if (tuning.fun == "avNNet") {
            maxit = 500
            maxnwts = 10 * (ncol(myExpl) + 1) + 10 + 1
            ## Automatically standardize data prior to modeling and prediction
            cmd.tuning <- paste0(cmd.tuning, " preProc = c('center', 'scale'), linout = TRUE, trace = FALSE,")
            cmd.tuning <- paste0(cmd.tuning, " MaxNWts.ANN = maxnwts, maxit = maxit))")
          } else if (tuning.fun %in% c("earth", "bam", "fda", "rpart", "glm")) { ## remove verbose
            cmd.tuning <- paste0(cmd.tuning, " tuneLength = tuning.length))")
          } else {
            cmd.tuning <- paste0(cmd.tuning, " tuneLength = tuning.length, verbose = FALSE))")
          }
          
          ## run tuning -------------------------------------------------------
          eval(parse(text = paste0("try(tuned.mod <- ", cmd.tuning)))
          
          ## GET tuned parameter values -------------------------------------------------------------
          if (!is.null(tuned.mod)) {
            tmp <- tuned.mod$results
            if (bm.format@data.type == 'binary') { tmp$TSS <- tmp$Sens + tmp$Spec - 1 }
            
            if (metric.eval == "RMSE") {
              selected <- which.min(tmp[, metric.eval])
            } else {
              selected <- which.max(tmp[, metric.eval])
            }
            
            if (model == "XGBOOST") {
              ## TO CHECK
              # for (param in train.params$params) {
              #   if (is.null(argstmp[[param]])){
              #     argstmp$params[[param]] <- tmp[selected, param]
              #   } else {
              #     argstmp[[param]] <- tmp[selected, param]}
              # }
            } else {
              for (param in train.params$params) {
                argstmp[[param]] <- tmp[selected, param]
              }
            }
            
            if (model == "CTA") {
              tuning.fun = "rpart2"
              eval(parse(text = paste0("try(tuned.mod <- ", cmd.tuning)))
              tuning.fun = "rpart" # needed to reset the tuning function in non parallel mode
              if (!is.null(tuned.mod)) {
                tmp = tuned.mod$results
                if (bm.format@data.type == 'binary') { tmp$TSS = tmp$Sens + tmp$Spec - 1 }
                if(metric.eval == "RMSE"){
                  argstmp[["maxdepth"]] <- tmp[which.min(tmp[, metric.eval]), "maxdepth"]
                } else {
                  argstmp[["maxdepth"]] <- tmp[which.max(tmp[, metric.eval]), "maxdepth"]
                }
              }
            }
          }
        } else if (workflow == "tidymodels") {
          
          tab.train = mySpExpl[, c(names(mySpExpl)[1], colnames(myExpl))]
          
          eval(parse(text = paste0("tm_recipe <- recipes::recipe(", names(mySpExpl)[1], " ~ ., data = tab.train)")))
          eval(parse(text = paste0("tm_model <- ", tuning.fun, "("
                                   , toString(paste0(colnames(tuning.grid), ' = tune()'))
                                   , ")")))
          tm_model = set_mode(tm_model, tuning.mod)
          tm_model = set_engine(tm_model, tuning.eng)
          tm_workflow <- workflows::workflow(preprocessor = tm_recipe, spec = tm_model)
          tm_folds <- rsample::vfold_cv(tab.train, v = 10, repeats = 3)
          eval(parse(text = paste0("tm_metrics <- yardstick::metric_set(", toString(metric.eval), ")")))
          
          ## run tuning -------------------------------------------------------
          tuned.mod <- tune::tune_grid(tm_workflow, resamples = tm_folds, grid = tuning.grid, metrics = tm_metrics)
          
          ## GET tuned parameter values -------------------------------------------------------------
          if (!is.null(tuned.mod)) {
            tmp = as.data.frame(tune::collect_metrics(tuned.mod))
            # select_best(tuned.mod, metric = met)
            
            if (metric.eval == "RMSE") {
              selected <- which.min(tmp[, "mean"])
            } else {
              selected <- which.max(tmp[, "mean"])
            }
            
            for (param in colnames(tuning.grid)) {
              tmp1 <- TuningParam[which(TuningParam$model == model), ]
              tmp2 <- tmp1$caret[which(tmp1$tidymodels == param)]
              if (is.na(tmp2)) {
                argstmp[[param]] <- tmp[selected, param]
              } else {
                argstmp[[tmp2]] <- tmp[selected, param]
              }
            }
          }
        }
      }
      
      ## run formula selection --------------------------------------------------------------------
      if (do.formula) {
        cat("\n\t\t\t> Tuning formula...")
        
        typ.vec = c('simple', 'quadratic', 'polynomial', 's_smoother')
        max.intlev <- min(ncol(myExpl) - 1, 3)
        if (model == "CTA") { max.intlev <- 0}
        if (model == "FDA") {
          typ.vec = c('simple', 's_smoother')
          max.intlev <- 0
          argstmp$method <- eval(parse(text = paste0("quote(", argstmp$method, ")")))
        }
        if (model %in% c("RF", "RFd", "GBM")) { typ.vec = c('simple', 'quadratic', 'polynomial') }
        
        if (bm.format@data.type == "binary") {
          myObs <- as.numeric(factor(myResp, labels = c("absence" = 0, "presence" = 1))) - 1
        } else {
          myObs <- myResp
        }
        
        data <- cbind(myExpl, resp = myObs)
        if (bm.format@data.type == "binary" && model %in% c("RF", "RFd")) {
          data <- data %>% mutate_at("resp", factor)
          argstmp$strata <- data[, "resp"]
          argstmp$sampsize <- unlist(ifelse(!is.null(argstmp$sampsize), list(argstmp$sampsize), nrow(data)))
        }
        argstmp$data <- data
        
        
        TMP <- foreach (typ = typ.vec, .combine = "rbind") %:%
          foreach (intlev = 0:max.intlev, .combine = "rbind") %do%
          {
            tuned.form <- NULL
            model.call <- paste0(bm.options@package, "::", bm.options@func)
            formu <- bm_MakeFormula(resp.name = "resp", expl.var = myExpl, type = typ, interaction.level = intlev)
            argstmp$formula <- formu
            argstmp <- argstmp[c("formula", "data", names(argstmp)[which(!(names(argstmp) %in% c("formula", "data")))])]
            tuned.form <- try(do.call(eval(parse(text = model.call)), argstmp), silent = TRUE)
            
            if (!inherits(tuned.form, "try-error") && !inherits(tuned.form, "data.frame")) {
              fit <- predict(tuned.form)
              if (bm.format@data.type %in% c("multiclass","ordinal")) {
                if (model %in% c("GAM", "GLM")) {
                  fit <- .threshold_ordinal(myObs, fit, metric.bm)$fit_factor
                } else if (model %in% c("CTA", "MARS")) {
                  fit <- predict(tuned.form, type = "class")
                }
              }
              if (bm.format@data.type == "binary" && model %in% c("CTA", "FDA", "RF", "RFd")) {
                fit <- predict(tuned.form, type = "class")
                fit <- as.numeric(fit) - 1
              }
              if (model == "MARS") {
                if (bm.format@data.type %in% c("multiclass","ordinal")) {
                  fit <- as.factor(fit[, 1])
                } else if (bm.format@data.type == "binary") {
                  fit <- as.numeric(fit[, 1])
                }
              }
              tmp <- bm_FindOptimStat(metric.bm, obs = myObs, fit = fit)[, "best.stat"]
              formu <- paste0(bm.format@sp.name, "~" , as.character(formu)[3])
              return(data.frame(stat = tmp, type = typ, interaction.level = intlev, formula = formu))
            }
          }
        
        if (metric.eval == "RMSE") {
          argstmp$formula <- formula(TMP[which.min(TMP[, 'stat']), "formula"])
        } else {
          argstmp$formula <- formula(TMP[which.max(TMP[, 'stat']), "formula"])
        }
      } else {
        if (length(argstmp$formula) <= 1 && model %in% c("CTA", "FDA", "GAM", "GBM", "GLM")) {
          argstmp$formula <- bm_MakeFormula(resp.name = bm.format@sp.name
                                            , expl.var = myExpl
                                            , type = 'simple'
                                            , interaction.level = 0)
        }
      }
      
      ## run variable selection -------------------------------------------------------------------
      if (do.stepAIC) {
        cat("\n\t\t\t> Tuning variables (AIC)...")
        
        if (model == "GLM") {
          glmStart <- glm(as.formula(paste0(bm.format@sp.name, " ~ 1")),
                          data = mySpExpl,
                          family = argstmp$family,
                          control = argstmp$control,
                          weights = current.weights,
                          mustart = rep(ifelse(!is.null(argstmp$mustart) & nchar(argstmp$mustart) > 0
                                               , argstmp$mustart, 0.5), length(myResp)),
                          model = TRUE)
          
          try(tuned.AIC <- MASS::stepAIC(glmStart,
                                         scope = list(upper = argstmp$formula, lower = ~1), ##upper = (sub(".*~", "~", argstmp$formula))
                                         k = criteria.AIC,
                                         direction = "both",
                                         trace = FALSE,
                                         steps = 10000))
          
          if (!is.null(tuned.AIC)) { argstmp$formula <- deparse(tuned.AIC$formula) }
          
        } else if (model == "GAM") {
          gamStart <- do.call(gam::gam, list(formula = as.formula(paste0(bm.format@sp.name, " ~ 1")),
                                             data = mySpExpl,
                                             family = argstmp$family,
                                             control = argstmp$control,
                                             weights = current.weights))
          tuned.AIC <- NULL
          try(tuned.AIC <-gam::step.Gam(gamStart,
                                        scope = .scope(head(myExpl), "gam::s", 6),
                                        direction = "both",
                                        trace = FALSE))
          if (!is.null(tuned.AIC)) { argstmp$formula <- formula(deparse(tuned.AIC$formula)) }
        }
      }
      
      cat("\n")
      return(argstmp)
    }
  names(argsval) <- combi$name_dataset
  return(argsval)
}


###################################################################################################

.bm_Tuning.check.args <- function(model, workflow, bm.options, bm.format, do.tuning, do.formula, do.stepAIC
                                  , metric.eval, metric.AIC, weights = NULL, ctrl.train = NULL, params.train)
{
  tuning.fun <- tuning.eng <- tuning.mod <- tuning.grid <- tuning.length <- train.params <- NULL
  metric.eval2 <- metric.eval
  
  ## check model --------------------------------------------------------------
  .fun_testIfInOnlyOne("model", model, c("ANN", "CTA", "DNN", "FDA", "GAM", "GBM", "GLM"
                                         , "MARS", "MAXENT", "MAXNET", "RF", "RFd", "SRE", "XGBOOST"))
  
  ## check bm.options ----------------------------------------------------------
  .fun_testIfInherits("bm.options", bm.options, c("BIOMOD.options.default", "BIOMOD.options.dataset"))
  ## check bm.format ----------------------------------------------------------
  .fun_testIfInherits("bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  
  ## check actions requested --------------------------------------------------
  if (do.tuning) {
    # if (model %in% c("GLM", "MAXNET")) {
    if (model %in% c("MAXNET")) {
      do.tuning <- FALSE
      .message("No tuning available for that model. Sorry.")
    } else if (model == "XGBOOST") {
      do.tuning <- FALSE
      .message("Due to upgrade of xgboost package, currently, no tuning available for that model. Sorry.")
    }
  }
  if (do.formula && model %in% c("DNN", "MAXENT", "MAXNET", "SRE", "XGBOOST")) {
    do.formula <- FALSE
    .message("No optimization of formula available for that model. Sorry.")
  }
  if (do.stepAIC && model %in% c("GAM", "GLM")) {
    if ((model == "GAM" && bm.options@package == "gam") ||
        bm.format@data.type %in% c("multiclass", "ordinal", "relative")) {
      do.stepAIC <- FALSE
      .message("No variables selection available for that combination of model/datatype. Sorry.")
    }
  } else {
    do.stepAIC <- FALSE
    .message("No variables selection available for that model. Sorry.")
  }
  
  ## check workflow -----------------------------------------------------------
  if (do.tuning) {
    .fun_testIfIn("workflow", workflow, c("caret", "tidymodels", "cito", "ENMeval", "biomod2"))
    if (model == "DNN" && workflow != "cito") {
      workflow <- "cito"
      .message("workflow set to cito (only one tuning package available for DNN)")
    } else if (model == "MAXENT" && workflow != "ENMeval") {
      workflow <- "ENMeval"
      .message("workflow set to ENMeval (only one tuning package available for MAXENT)")
    } else if (model == "SRE" && workflow != "biomod2") {
      workflow <- "biomod2"
      .message("workflow set to biomod2 (only one tuning package available for SRE)")
    } else if (!(model %in% c("DNN", "MAXENT", "SRE"))) {
      .fun_testIfInOnlyOne("workflow", workflow, c("caret", "tidymodels"))
    }
  }
  
  ## check namespace ----------------------------------------------------------
  if (do.tuning) {
    if (workflow == "caret") {
      .fun_testIfNamespace("caret")
    } else if (workflow == "tidymodels") {
      # if (isNamespaceLoaded("caret")) { unloadNamespace("caret") } ## need to unload caret before recipes
      # if (isNamespaceLoaded("tune")) { unloadNamespace("tune") } ## need to unload tune before recipes
      # .fun_testIfNamespace("tidymodels")
      .fun_testIfNamespace("parsnip")
      .fun_testIfNamespace("recipes")
      .fun_testIfNamespace("workflows")
      .fun_testIfNamespace("dials")
      .fun_testIfNamespace("rsample")
      .fun_testIfNamespace("yardstick")
      .fun_testIfNamespace("tune")
    } else if (model == "MAXENT") {
      .fun_testIfNamespace("ENMeval")
    } else if (model == "SRE") {
      .fun_testIfNamespace("dismo")
    }
  }
  if (do.formula == TRUE) { .fun_testIfNamespace("gam") }
  if (do.stepAIC == TRUE) { .fun_testIfNamespace("MASS") }
  
  ## check weights ------------------------------------------------------------
  if (model %in% c("CTA", "FDA", "GAM", "GLM") && is.null(weights)) { ## tuning or stepAIC
    weights <- rep(1, length(bm.format@data.species))
  }
  
  
  ## Do tuning ----------------------------------------------------------------
  if (do.tuning) {
    
    ## check combination [model + datatype + workflow]
    tmp <- unlist(ifelse(bm.format@data.type == "binary", "binary"
                         , list(c("nonbinary", bm.format@data.type))))
    TuningTable_ <- TuningTable[which(TuningTable$model == model &
                                        TuningTable$type %in% tmp &
                                        TuningTable$tuning_package == workflow), , drop = FALSE]
    if (model == "GAM") {
      tmp <- ifelse(bm.options@package == "gam", "gamLoess", bm.options@func)
      TuningTable_ <- TuningTable_[which(TuningTable_$tuning_function == tmp), , drop = FALSE]
    } else if (model == "MAXENT") {
      TuningTable_ <- TuningTable_[which(TuningTable_$tuning_engine == "maxnet"), , drop = FALSE]
      .message("Tuning with 'maxent.jar' currently not working. tuning_engine has been set to 'maxnet'.")
    }
    
    if (nrow(TuningTable_) == 0) {
      do.tuning <- FALSE
      .message("No tuning available for that combination of model/datatype. Sorry.")
    } else if (nrow(TuningTable_) > 1) {
      stop("Petiprobleum")
    }
    tuning.fun <- TuningTable_$tuning_function
    tuning.eng <- TuningTable_$tuning_engine
    tuning.mod <- TuningTable_$tuning_mode
  }
  
  ## Do tuning ----------------------------------------------------------------
  if (do.tuning) {
    
    ## check params.train ---------------------------------
    tmp <- TuningParamTrain
    if (!is.null(params.train)) {
      if (any(!names(params.train) %in% names(TuningParamTrain))) {
        tmp <- setdiff(names(params.train), names(TuningParamTrain)) ## ou l'inverse ? to test !!
        .message(tmp, " provided in params.train will be ignored (not contained within TuningParamTrain)")
      }
      ## complete TuningParamTrain
      for (i in names(params.train)) {
        if (i %in% names(tmp)) {
          tmp[[i]] = params.train[[i]]
        }
      }
    }
    params.train = tmp
    print("yi")
    print(names(params.train))
    
    
    if (workflow == "caret") { ## A quel point c'est utile ?
      params <- caret::getModelInfo(model = tuning.fun)
      train.params <- list(pkg = params[[tuning.fun]]$library
                           , params = params[[tuning.fun]]$parameters$parameter)
      ## ATTENTION CTA rpart + rpart2 ??
    }
    
    ## get tuning grid through params.train ---------------
    # browser()
    # if (model %in% c("CTA")) { ##, "GLM")) {
    if (model == "CTA" && workflow == "caret") { ##, "GLM")) {
      tuning.grid <- NULL
    } else if (workflow %in% c("caret", "tidymodels")) {
      # browser()
      params.train <- params.train[grep(paste0(model, "\\."), names(params.train))]
      
      # tmp1 <- grep(paste0(model, "\\."), names(params.train), value = TRUE)
      # tmp2 <- sub(paste0(model, "\\."), "", tmp1)
      # tmp3 <- TuningParam[which(TuningParam$model == model), c("caret", "tidymodels")]
      # tmp3 <- as.vector(na.exclude(unlist(tmp3)))
      # tmp <- tmp1[which(tmp2 %in% tmp3)]
      # # browser()
      # params.train <- params.train[tmp]
      
      # tmp1 <- grep(paste0(model, "\\."), names(params.train), value = TRUE)
      # tmp2 <- sub(paste0(model, "\\."), "", tmp1)
      # tmp3 <- which(TuningParam$model == model & !is.na(TuningParam[, workflow]))
      # tmp4 <- TuningParam$caret[tmp3]
      # if (length(which(is.na(tmp4))) > 0) tmp4[which(is.na(tmp4))] <- TuningParam$tidymodels[tmp3[which(is.na(tmp4))]]
      # tmp <- tmp1[which(tmp2 %in% tmp4)] ## PROBLEUM : param dans tidymodels mais NA dans caret
      # # tmp <- tmp1[which(tmp2 %in% TuningParam$caret[tmp3])] ## PROBLEUM : param dans tidymodels mais NA dans caret
      # # tmp <- paste0(model, ".", TuningParam[tmp3[which(TuningParam[tmp3, workflow] %in% tmp2)], "caret"])
      # # browser()
      # # print(tmp)
      # # stop()
      # params.train <- params.train[tmp]
      print("yo")
      print(names(params.train))
      
      
      if (workflow == "caret") {
        toCompare <- paste0(model, ".", train.params$params)
      } else if (workflow == "tidymodels") {
        tmp <- TuningParam[which(TuningParam$model == model), c("caret", "tidymodels")]
        tmp <- tmp[which(!is.na(tmp$tidymodels)), ]
        ind.NA <- is.na(tmp$caret)
        toCompare <- c(tmp$caret[!ind.NA], tmp$tidymodels[ind.NA])
        toCompare <- paste0(model, ".", toCompare)
      }
      # print(toCompare)
      ## remove parameters not found in function definition
      if (any(!names(params.train) %in% toCompare)) {
        tmp <- setdiff(names(params.train), toCompare)
        params.train <- params.train[which(names(params.train) %in% toCompare)]
        .message("params.train set to ", toString(names(params.train)), " (", toString(tmp), " removed)")
      }
      
      if (model %in% c("RF", "RFd")) {
        ## reduce mtry not to be superior to number of variables
        ind <- which(params.train[[paste0(model, ".mtry")]] <= ncol(bm.format@data.env.var))
        params.train[[paste0(model, ".mtry")]] <- params.train[[paste0(model, ".mtry")]][ind]
      }
      
      if (workflow == "caret") {
        ## subselect some parameters for some models
        if (model == "GAM") {
          if (tuning.fun == "gamLoess"){
            params.train <- params.train[c("GAM.span", "GAM.degree")]
          } else { ## bam, gam
            params.train <- params.train[c("GAM.select", "GAM.method")]
          }
        } else if (model %in% c("RF", "RFd")) {
          params.train <- params.train[paste0(model, ".mtry")]
        }
        # ## remove parameters not found in function definition
        # if (any(!names(params.train) %in% paste0(model, ".", train.params$params))) {
        #   tmp <- setdiff(names(params.train), paste0(model, ".", train.params$params))
        #   params.train <- params.train[which(names(params.train) %in% paste0(model, ".", train.params$params))]
        #   .message("params.train set to ", toString(names(params.train)), " (", toString(tmp), " removed)")
        # }
        names(params.train) <- sub(model, "", names(params.train)) ## caret
      } else if (workflow == "tidymodels") {
        names(params.train) <- sub(paste0(model, "."), "", names(params.train)) ## tidymodels
        for (i in 1:length(params.train)) {
          ind <- which(TuningParam$model == model & TuningParam$caret == names(params.train)[i])
          if (length(ind) > 0) {
            names(params.train)[i] <- TuningParam$tidymodels[ind] ## tidymodels
          }
        }
      }
      print("ye")
      print(names(params.train))
      tuning.grid <- do.call(expand.grid, params.train)
      print(tuning.grid)
      
    } else if (model == "DNN") {
      names_DNN <- grep(paste0(model,"\\."), names(params.train), value = TRUE)
      # browser()
      for (param.n in names_DNN) {
        print("--------------------------------------------------------")
        print(param.n)
        print(params.train[[param.n]])
        real.name <- unlist(strsplit(param.n, split = "\\."))[2]
        if (real.name == "hidden") {
          hidden_d <- params.train$DNN.hidden$depth
          hidden_w <- params.train$DNN.hidden$width
          if (length(hidden_d) == 1 && length(hidden_w) == 1) {
            DNN_hidden <- rep(hidden_w, hidden_d)
          } else if (length(hidden_d) == 1) {
            DNN_hidden <- cito::tune(lower = hidden_w[1], upper = hidden_w[2], additional = hidden_d, fixed = 'depth')
          } else if (length(hidden_w) == 1) {
            DNN_hidden <- cito::tune(lower = hidden_d[1], upper = hidden_d[2], additional = hidden_w, fixed = 'width')
          } else {
            DNN_hidden <- cito::tune(lower = c(hidden_w[1], hidden_d[1]), upper = c(hidden_w[2], hidden_d[2]))
          }
          params.train$DNN.hidden <- DNN_hidden
        } else if (length(params.train[[param.n]]) == 2) {
          params.train[[param.n]] <- cito::tune(lower = params.train[[param.n]][1], upper = params.train[[param.n]][2])
        # } else if (length(params.train[[param.n]]) == 1) {
        } else {
          params.train[[param.n]] <- cito::tune(values = params.train[[param.n]])
        }
        print(params.train[[param.n]])
      }
    } else if (model == "MAXENT") {
      .fun_testIfIn("params.train$MAXENT.algorithm", params.train$MAXENT.algorithm, c("maxent.jar", "maxnet"))
    } else if (model == "SRE") {
      sapply(params.train$SRE.quant, FUN = .fun_testIf0X, objName =  "params.train$SRE.quant")
    }
    
    
    if (workflow == "caret") {
      ## get tuning length --------------------------------
      tuning.length <- 1
      if (model == "CTA") tuning.length <- 30
      if (model == "RF" || model == "RFd") tuning.length <- min(30, ncol(bm.format@data.env.var))
      
      ## check control ------------------------------------
      if (!(model %in% c("DNN", "MAXENT", "SRE")) && is.null(ctrl.train)) {
        if (bm.format@data.type == "binary") {
          ctrl.train <- caret::trainControl(method = "repeatedcv",
                                            repeats = 3,
                                            number = 10,
                                            summaryFunction = caret::twoClassSummary,
                                            classProbs = TRUE,
                                            returnData = FALSE)
        } else {
          ctrl.train <- caret::trainControl(method = "repeatedcv",
                                            repeats = 3,
                                            number = 10,
                                            summaryFunction = caret::defaultSummary,
                                            classProbs = FALSE,
                                            returnData = FALSE)
        }
      }
    } else {
      tuning.length <- NULL
      ctrl.train <- NULL
    }
    
    ## check metric.eval ----------------------------------
    if (model == "DNN") {
      metric.eval <- NULL
      .message("metric.eval set to NULL (DNN)")
    } else {
      tmp1 = ifelse(bm.format@data.type == "binary", "binary"
                    , ifelse(bm.format@data.type %in% c("multiclass", "ordinal"), "class", "nonbinary"))
      tmp2 = ifelse(model %in% c("MAXENT", "SRE"), model, workflow)
      tmp = na.exclude(TuningMetric[which(TuningMetric$type == tmp1), c("biomod2", tmp2)])
      rownames(tmp) <- tmp[, "biomod2"]
      .fun_testIfIn("metric.eval", metric.eval, tmp[, "biomod2"])
      metric.eval2 <- tmp[metric.eval, tmp2]
    }
  }
  
  ## Do formula ---------------------------------------------------------------
  if (do.formula) {
    if (bm.format@data.type == "binary") {
      metric.bm <- metric.eval
    } else if (bm.format@data.type %in% c("multiclass", "ordinal")) {
      metric.bm <- "Accuracy"
    } else {
      metric.bm <- ifelse(metric.eval == "RMSE", "RMSE", "Rsquared")
    }
  } else { metric.bm <- NULL }
  
  ## Do step AIC --------------------------------------------------------------
  if (do.stepAIC) {
    .fun_testIfIn("metric.AIC", metric.AIC, c("AIC", "BIC"))
    if (metric.AIC == "AIC") { criteria.AIC <- 2 }
    if (metric.AIC == "BIC")  { criteria.AIC <- log(ncol(bm.format@data.env.var)) }
  } else { criteria.AIC <- NA }
  
  return(list(weights = weights
              , do.tuning = do.tuning
              , do.formula = do.formula
              , do.stepAIC = do.stepAIC
              , criteria.AIC = criteria.AIC
              , tuning.fun = tuning.fun
              , tuning.eng = tuning.eng
              , tuning.mod = tuning.mod
              , tuning.length = tuning.length
              , tuning.grid = tuning.grid
              , train.params = train.params
              , params.train = params.train
              , ctrl.train = ctrl.train
              , metric.eval = metric.eval2
              , metric.bm = metric.bm))
}


###################################################################################################

.scope <- function(enviroTrain, Smoother, degree)
{
  XXX <- enviroTrain
  deg <- degree
  vnames <- names(XXX[])
  step.list <- as.list(vnames)
  names(step.list) <- vnames
  NbVar <- dim(enviroTrain)[2]
  i <- 1
  while (i <= NbVar)
  {
    vname <- names(XXX)[i]
    # loops through independent variable names
    junk <- paste0("1 + ", vname)
    # minimum scope
    if (is.numeric(XXX[, i])) {
      junk <- c(junk, paste0(Smoother, "(", vname, ",", deg, ")"))
      junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
    } else if (is.factor(XXX[, i])) {
      junk <- c(junk, vname)
      junk <- eval(parse(text = paste("~", paste(junk, collapse = "+"))))
    }
    step.list[[vname]] <- junk
    i <- i + 1
  }
  
  return(step.list)
}

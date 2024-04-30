## --------------------------------------------------------------------------- ##
##' @name bm_Tuning
##' @author Frank Breiner, Maya Gueguen, Helene Blancheteau
##' 
##' @title Tune models parameters
##' 
##' @description This internal \pkg{biomod2} function allows to tune single model parameters and 
##' select more efficient ones based on an evaluation metric.
##' 
##'
##' @param model a \code{character} corresponding to the  algorithm to be tuned, must be either 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}
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
##' \code{                                  summaryFunction = caret::twoClassSummary,} \cr
##' \code{                                  classProbs = TRUE, returnData = FALSE)} \cr \cr
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
##'   \item{FDA}{\code{degree}, \code{nprune}}
##'   \item{GAM.gam}{\code{span}, \code{degree}}
##'   \item{GAM.mgcv}{\code{select}, \code{method}}
##'   \item{GBM}{\code{n.trees}, \code{interaction.depth}, \code{shrinkage}, \code{n.minobsinnode}}
##'   \item{MARS}{\code{degree}, \code{nprune}}
##'   \item{MAXENT}{\code{algorithm}, \code{parallel}}
##'   \item{RF}{\code{mtry}}
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
##'   \item \code{SRE} is tuned through \code{\link{bm_SRE}} function
##'   \item All other models are tuned through \code{\link[caret]{train}} function
##'   \item No optimization of formula for \code{MAXENT}, \code{MAXNET}, \code{SRE} and 
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
##' @family Secundary functions
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
##' 
##' 
##' @export
##' 
#------------------------------------------------------------------------------#


bm_Tuning <- function(model,
                      tuning.fun,
                      do.formula = FALSE,
                      do.stepAIC = FALSE,
                      bm.options,
                      bm.format,
                      calib.lines = NULL,
                      metric.eval = 'TSS',
                      metric.AIC = 'AIC',
                      weights = NULL,
                      ctrl.train = NULL,
                      params.train = list(ANN.size = c(2, 4, 6, 8),
                                          ANN.decay = c(0.001, 0.01, 0.05, 0.1),
                                          ANN.bag = FALSE, 
                                          FDA.degree = 1:2, 
                                          FDA.nprune = 2:38,
                                          GAM.select = c(TRUE, FALSE),
                                          GAM.method = c('GCV.Cp', 'GACV.Cp', 'REML', 'P-REML', 'ML', 'P-ML'),
                                          GAM.span = c(0.3, 0.5, 0.7),
                                          GAM.degree = 1,
                                          GBM.n.trees = c(500, 1000, 2500),
                                          GBM.interaction.depth = seq(2, 8, by = 3),
                                          GBM.shrinkage = c(0.001, 0.01, 0.1),
                                          GBM.n.minobsinnode = 10,
                                          MARS.degree = 1:2, 
                                          MARS.nprune = 2:max(38, 2 * ncol(bm.format@data.env.var) + 1),
                                          MAXENT.algorithm = 'maxnet',
                                          MAXENT.parallel = TRUE,
                                          RF.mtry = 1:min(10, ncol(bm.format@data.env.var)),
                                          SRE.quant = c(0, 0.0125, 0.025, 0.05, 0.1),
                                          XGBOOST.nrounds = 50,
                                          XGBOOST.max_depth = 1,
                                          XGBOOST.eta = c(0.3, 0.4),
                                          XGBOOST.gamma = 0,
                                          XGBOOST.colsample_bytree = c(0.6, 0.8),
                                          XGBOOST.min_child_weight = 1,
                                          XGBOOST.subsample = 0.5))
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_Tuning.check.args(model = model, tuning.fun = tuning.fun
                                , do.formula = do.formula, do.stepAIC = do.stepAIC
                                , bm.options = bm.options, bm.format = bm.format
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
  
  
  if (model != "MAXENT" && is.null(ctrl.train)) {
    ## check control
    ctrl.train <- caret::trainControl(method = "repeatedcv",
                                      repeats = 3,
                                      number = 10,
                                      summaryFunction = caret::twoClassSummary,
                                      classProbs = TRUE,
                                      returnData = FALSE)
  }
  
  
  argsval <- foreach(PA.i = combi$PA, calib.i = combi$calib, dataset.i = combi$name_dataset) %do%
    {
      cat(paste0("\n\t\t> Dataset ", dataset.i))
      argstmp <- bm.options@args.default
      
      if (model == "MAXNET") {
        warning("No tuning available for that model. Sorry.")
      } else {
        ## 1. SPECIFIC CASE OF MAXENT OR SRE ------------------------------------------------------------
        if (model %in% c("MAXENT", "SRE")) {
          cat("\n\t\t\t> Tuning parameters...")
          
          ## create dataset ---------------------------------------------------------
          mySpExpl <- get_species_data(bm.format)
          mySpExpl[["_allData_allRun"]] <- TRUE
          mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
          mySpExpl <- mySpExpl[which(mySpExpl[, PA.i] == TRUE), ]
          
          ## SRE case
          myResp <- mySpExpl[, 1]
          myExpl <- mySpExpl[, 4:(3 + ncol(bm.format@data.env.var))]
          
          ## MAXENT case
          if (params.train$MAXENT.algorithm == "maxnet") {
            mySpExpl[["_allData_allRun"]] <- NULL
            mySpExpl[, 1] <- ifelse(mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), 1, 0)
            mySpExpl <- mySpExpl[, 1:(3 + ncol(bm.format@data.env.var))]
          }
          
          if (model == "MAXENT") { # ------------------------------------------#
            try(tune.MAXENT <- ENMeval::ENMevaluate(occs = mySpExpl[mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), ],
                                                    bg = mySpExpl[mySpExpl[, 1] == 0 | is.na(mySpExpl[, 1]), ],
                                                    tune.args = list(rm = seq(0.5, 1, 0.5), fc = c("L")),
                                                    algorithm = params.train$MAXENT.algorithm,
                                                    partitions = "randomkfold",
                                                    partition.settings = list(kfolds = 10),
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
          } else if (model == "SRE") { # -------------------------------------- #
            myResp <- sapply(myResp, function(xx) ifelse(xx == 0 || is.na(xx), 0, 1))
            
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
          }
          
        } else {
          ## 2. ALL OTHER MODELS ------------------------------------------------------------------------
          
          ## create dataset
          mySpExpl <- get_species_data(bm.format)
          mySpExpl[["_allData_allRun"]] <- 1
          current.weights <- weights[which(calib.lines[, calib.i] == TRUE &
                                             mySpExpl[, PA.i] == TRUE)]
          mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
          mySpExpl <- mySpExpl[which(mySpExpl[, PA.i] == TRUE), ]
          mySpExpl[, 1] <- as.factor(ifelse(mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), "presence", "absence"))
          myResp <- mySpExpl[, 1]
          myExpl <- mySpExpl[, 4:(3 + ncol(bm.format@data.env.var))]
          
          ## run tuning -------------------------------------------------------
          
          cmd.tuning <- "caret::train(x = myExpl, y = myResp, method = tuning.fun, tuneGrid = tuning.grid,"
          cmd.tuning <- paste0(cmd.tuning, " trControl = ctrl.train, metric = 'ROC',")
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
          
          if (model != "GLM") {
            cat("\n\t\t\t> Tuning parameters...")
            eval(parse(text = paste0("try(tuned.mod <- ", cmd.tuning)))
          
            ## GET tuned parameter values -------------------------------------------------------------
            if (!is.null(tuned.mod)) {
              tmp <- tuned.mod$results
              tmp$TSS <- tmp$Sens + tmp$Spec - 1
              
              if (model == "XGBOOST") {
                for (param in train.params$params) {
                  if (is.null(argstmp[[param]])){
                    argstmp$params[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
                  } else {
                    argstmp[[param]] <- tmp[which.max(tmp[, metric.eval]), param]}
                }
              } else {
                for (param in train.params$params) {
                  argstmp[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
                }
              }
              
              tuning.form <- tuning.grid[which.max(tmp[, metric.eval]), ]
              
              if (model == "RF") {
                tuning.form <- data.frame(mtry = tuning.grid[which.max(tmp[, metric.eval]), ])
              }
              
              if (model == "CTA") {
                tuning.fun = "rpart2"
                eval(parse(text = paste0("try(tuned.mod <- ", cmd.tuning)))
                tuning.fun = "rpart" # needed to reset the tuning function in non parallel mode
                if (!is.null(tuned.mod)) {
                  tmp = tuned.mod$results
                  tmp$TSS = tmp$Sens + tmp$Spec - 1
                  argstmp[["maxdepth"]] <- tmp[which.max(tmp[, metric.eval]), "maxdepth"]
                }
              }
            }
          } else { tuning.form <- tuning.grid }
          
          ## run formula selection ------------------------------------------------------------------
          if (do.formula) {
            cat("\n\t\t\t> Tuning formula...")
            
            cmd.form <- sub("tuneGrid = tuning.grid", "tuneGrid = tuning.form", cmd.tuning)
            cmd.form <- sub("weights = current.weights,", "", cmd.form)
            cmd.init <- "form = bm_MakeFormula(resp.name = 'resp', expl.var = myExpl, type = typ, interaction.level = intlev),"
            cmd.init <- paste0(cmd.init, " data = cbind(myExpl, resp = myResp),")
            cmd.form <- sub("x = myExpl, y = myResp,", cmd.init, cmd.form)
            
            max.intlev <- min(ncol(myExpl) - 1, 3)
            typ.vec = c('simple', 'quadratic', 'polynomial', 's_smoother')
            
            if (model %in% c("CTA", "FDA")) {
              if (model == "CTA") { typ.vec = c('simple', 'quadratic', 'polynomial', 's_smoother') }
              if (model == "FDA") { typ.vec = c('simple', 's_smoother') }
              
              TMP <- foreach (typ = typ.vec, .combine = "rbind") %do%
                {
                  tuned.form <- NULL
                  intlev <- 0
                  eval(parse(text = paste0("capture.output("
                                           , "try(tuned.form <- ", sub(")$", ", silent = TRUE)", cmd.form)
                                           , ")")))
                  if (!is.null(tuned.form)) {
                    tmp <- tuned.form$results
                    tmp$TSS <- tmp$Sens + tmp$Spec - 1
                    formu <- tuned.form$coefnames
                    formu <- paste0(bm.format@sp.name, " ~ 1 + ", paste0(formu, collapse = " + "))
                    return(data.frame(tmp, type = typ, interaction.level = intlev, formula = formu))
                  }
                }
            } else {
              if (model == "RF") { typ.vec = c('simple','quadratic', 'polynomial') }
              
              TMP <- foreach (typ = typ.vec, .combine = "rbind") %:%
                foreach (intlev = 0:max.intlev, .combine = "rbind") %do%
                {
                  tuned.form <- NULL
                  eval(parse(text = paste0("capture.output("
                                           , "try(tuned.form <- ", sub(")$", ", silent = TRUE)", cmd.form)
                                           , ")")))
                  if (!is.null(tuned.form)) {
                    tmp <- tuned.form$results
                    tmp$TSS <- tmp$Sens + tmp$Spec - 1
                    formu <- tuned.form$coefnames
                    formu <- paste0(bm.format@sp.name, " ~ 1 + ", paste0(formu, collapse = " + "))
                    return(data.frame(tmp, type = typ, interaction.level = intlev, formula = formu))
                  }
                }
            }
            argstmp$formula <- TMP[which.max(TMP[, metric.eval]), "formula"]
            if (model %in% c("ANN", "GAM", "GBM", "MARS", "RF")) {
              argstmp$formula <- formula(argstmp$formula)
            }
          } else {
            if (model %in% c("CTA", "FDA", "GAM", "GBM", "GLM")) {
              argstmp$formula <- bm_MakeFormula(resp.name = bm.format@sp.name
                                                , expl.var = myExpl
                                                , type = 'simple'
                                                , interaction.level = 0)
            }
          }
          
          ## run variable selection -----------------------------------------------------------------
          if (do.stepAIC &&
              (model == "GLM" || 
               (model == "GAM" && bm.options@package == "gam")) ) {
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
                                             scope = list(upper = (sub(".*~", "~", argstmp$formula)), lower = ~1),
                                             k = criteria.AIC,
                                             direction = "both",
                                             trace = FALSE,
                                             steps = 10000))
              if (!is.null(tuned.AIC)) { argstmp$formula <- deparse(tuned.AIC$formula) }
              
            } else if (model == "GAM") { # if (bm.options@GAM$algo == 'GAM_gam') { ## gam package
              gamStart <- do.call(gam::gam, list(formula = as.formula(paste0(bm.format@sp.name, " ~ 1")),
                                                 data = mySpExpl, 
                                                 family = argstmp$family,
                                                 control = argstmp$control,
                                                 weights = current.weights))          
              tuned.AIC <- NULL
              try(tuned.AIC <- 
                    gam::step.Gam(
                      gamStart,
                      scope = .scope(head(myExpl), "gam::s", 6),
                      direction = "both",
                      trace = FALSE))
              if (!is.null(tuned.AIC)) { argstmp$formula <- formula(deparse(tuned.AIC$formula)) }
            }
          }
        }
      }
      return(argstmp)
    }
  names(argsval) <- combi$name_dataset
  return(argsval)
}



# Check arguments ------------------------------------------------------------

.bm_Tuning.check.args <- function(model, tuning.fun, do.formula, do.stepAIC
                                  , bm.options, bm.format, metric.eval, metric.AIC
                                  , weights = NULL, params.train)
{
  ## check model --------------------------------------------------------------
  .fun_testIfIn(TRUE, "model", model, c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM"
                                        , "MARS", "MAXENT", "MAXNET", "RF", "SRE", "XGBOOST"))
  
  ## check namespace ----------------------------------------------------------
  if (!isNamespaceLoaded("caret")) { 
    if (!requireNamespace('caret', quietly = TRUE)) stop("Package 'caret' not found")
  }
  if (model == "MAXENT" && !isNamespaceLoaded('ENMeval')) { 
    if (!requireNamespace('ENMeval', quietly = TRUE)) stop("Package 'ENMeval' not found")
  } else if (model == "SRE" && !isNamespaceLoaded('dismo')) { 
    if (!requireNamespace('dismo', quietly = TRUE)) stop("Package 'dismo' not found")
  }
  if (do.formula == TRUE) {
    if (!requireNamespace('gam', quietly = TRUE)) stop("Package 'gam' not found")
  }
  if (do.stepAIC == TRUE) {
    if (!requireNamespace('MASS', quietly = TRUE)) stop("Package 'MASS' not found")
  }
  
  ## check bm.options ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.options", bm.options, c("BIOMOD.options.default", "BIOMOD.options.dataset"))
  ## check bm.format ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  ## check params.train -------------------------------------------------------
  params.train_init = list(ANN.size = c(2, 4, 6, 8),
                           ANN.decay = c(0.001, 0.01, 0.05, 0.1),
                           ANN.bag = FALSE, 
                           FDA.degree = 1:2, 
                           FDA.nprune = 2:38,
                           GAM.select = c(TRUE, FALSE),
                           GAM.method = c('GCV.Cp', 'GACV.Cp', 'REML', 'P-REML', 'ML', 'P-ML'),
                           GAM.span = c(0.3, 0.5, 0.7),
                           GAM.degree = 1,
                           GBM.n.trees = c(500, 1000, 2500),
                           GBM.interaction.depth = seq(2, 8, by = 3),
                           GBM.shrinkage = c(0.001, 0.01, 0.1),
                           GBM.n.minobsinnode = 10,
                           MARS.degree = 1:2, 
                           MARS.nprune = 2:max(38, 2 * ncol(bm.format@data.env.var) + 1),
                           MAXENT.algorithm = 'maxnet',
                           MAXENT.parallel = TRUE,
                           RF.mtry = 1:min(10, ncol(bm.format@data.env.var)),
                           SRE.quant = c(0, 0.0125, 0.025, 0.05, 0.1),
                           XGBOOST.nrounds = 50,
                           XGBOOST.max_depth = 1,
                           XGBOOST.eta = c(0.3, 0.4),
                           XGBOOST.gamma = 0,
                           XGBOOST.colsample_bytree = c(0.6, 0.8),
                           XGBOOST.min_child_weight = 1,
                           XGBOOST.subsample = 0.5)
  
  for (i in names(params.train)) {
    if (i %in% names(params.train_init)) {
      params.train_init[[i]] = params.train[[i]]
    }
  }
  params.train = params.train_init
  ## check evaluation metric --------------------------------------------------
  if (model == "MAXENT") {
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("auc.val.avg", "auc.diff.avg", "or.mtp.avg", "or.10p.avg", "AICc"))
    .fun_testIfIn(TRUE, "params.train$MAXENT.algorithm", params.train$MAXENT.algorithm, c("maxent.jar", "maxnet"))
  } else if (model == "SRE") {
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("AUC", "Kappa", "TSS"))
    sapply(params.train$SRE.quant,FUN=.fun_testIf01,test = TRUE,objName =  "params.train$SRE.quant")
  } else {
    .fun_testIfIn(TRUE, "metric.eval", metric.eval, c("ROC", "TSS"))
  }
  ## check weights ------------------------------------------------------------
  if (model %in% c("CTA", "FDA", "GAM", "GLM") && is.null(weights)) { 
    weights = rep(1, length(bm.format@data.species))
  }
  
  
  ## get tuning function and parameters ---------------------------------------
  all.fun <- c('avNNet', 'rpart', 'rpart2', 'fda', 'gamLoess', 'bam', 'gam', 'gbm', 'glm', 'earth', 'rf', 'xgbTree')
  all.params <- foreach (fi = all.fun) %do% {
    params <- caret::getModelInfo(model = fi)
    return(list(pkg = params[[fi]]$library, params = params[[fi]]$parameters$parameter))
  }
  names(all.params) <- all.fun
  
  .fun_testIfIn(TRUE, "tuning.fun", tuning.fun, c(all.fun, "bm_SRE", "ENMevaluate", "maxnet"))
  train.params <- all.params[[tuning.fun]]
  ## get tuning grid through params.train -------------------------------------
  tuning.grid <- NULL
  if (model %in% c("ANN", "FDA", "GAM", "GBM", "MARS", "RF", "XGBOOST")) {
    if (!(model == "GAM")) {
      params.train = params.train[grep(model, names(params.train))]
      .fun_testIfIn(TRUE, "names(params.train)", names(params.train), paste0(model, ".", train.params$params))
    } else if (tuning.fun == "gamLoess"){
      params.train = params.train[c('GAM.span', "GAM.degree")]
    } else {
      params.train = params.train[c('GAM.select', 'GAM.method')]
    }
    names(params.train) = sub(model, "", names(params.train))
    tuning.grid <- do.call(expand.grid, params.train)
  }
  
  ## get tuning length --------------------------------------------------------
  tuning.length <- 1
  if (model == "CTA") tuning.length <- 30
  if (model == "RF") tuning.length <- min(30, ncol(bm.format@data.env.var))
  
  ## Do formula ---------------------------------------------------------------
  if (model %in% c("MAXENT", "MAXNET", "SRE", "XGBOOST") && do.formula == TRUE) {
    do.formula <- FALSE 
    cat("\n No optimization of formula for", model)
  }
  
  ## get criteria -------------------------------------------------------------
  if (do.stepAIC && (model == "GLM" || 
                     (model == "GAM" && bm.options@package == "gam"))) {
    .fun_testIfIn(TRUE, "metric.AIC", metric.AIC, c("AIC", "BIC"))
    if (metric.AIC == "AIC") criteria.AIC <- 2
    if (metric.AIC == "BIC") criteria.AIC <- log(ncol(bm.format@data.env.var))
  } else {
    do.stepAIC <- FALSE
    criteria.AIC <- NA
  }
  
  return(list(weights = weights
              , do.formula = do.formula
              , criteria.AIC = criteria.AIC
              , tuning.fun = tuning.fun
              , train.params = train.params
              , tuning.length = tuning.length
              , tuning.grid = tuning.grid
              , params.train = params.train))
}


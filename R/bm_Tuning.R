## --------------------------------------------------------------------------- ##
##' @name bm_Tuning
##' @author Frank Breiner, Maya Gueguen
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
##' be either \code{ROC}, \code{TSS}, \code{KAPPA} for \code{SRE} only, and \code{auc.val.avg}, 
##' \code{auc.diff.avg}, \code{or.mtp.avg}, \code{or.10p.avg}, \code{AICc} for \code{MAXENT} only
##' @param metric.AIC a \code{character} corresponding to the AIC metric to be used, must 
##' be either \code{AIC} or \code{BIC}
##' @param weights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to observation weights (one per 
##' observation, see Details)
##' 
##' @param params.train a \code{list} containing values of model parameters to be tested 
##' (see Details)
## @param ctrl.train (\emph{optional, default} \code{NULL}) \cr to be added ?
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
##'   \item{FDA}{\code{degree}, \code{nprune}}
##'   \item{GAM}{\code{select}, \code{method}}
##'   \item{GBM}{\code{n.trees}, \code{interaction.depth}, \code{shrinkage}, \code{n.minobsinnode}}
##'   \item{MARS}{\code{degree}, \code{nprune}}
##'   \item{RF}{\code{mtry}}
##'   \item{SRE}{\code{quant}}
##'   \item{XGBOOST}{\code{nrounds}, \code{max_depth}, \code{eta}, \code{gamma}, 
##'   \code{colsampl_bytree}, \code{min_child_weight}, \code{subsample}}
##' }
##' 
##' The \code{\link{expand.grid}} function is used to build a \code{matrix} containing all 
##' combinations of parameters to be tested.
##' 
##' @note 
##' \itemize{
##'   \item No tuning for \code{GLM} and \code{MAXNET}
##'   \item \code{MAXENT} is tuned through \code{\link[ENMeval]{ENMevaluate}} function
##'   \item \code{SRE} is tuned through \code{\link{bm_SRE}} function
##'   \item All other models are tuned through \code{\link[caret]{train}} function
##'   \item No optimization of formula for \code{MAXENT}, \code{MAXNET} and \code{SRE}
##'   \item Selection variables only for \code{GAM} and \code{GLM}
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
##' @importFrom foreach foreach %do%
##' @importFrom stats aggregate  
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
                                          GBM.n.trees = c(500, 1000, 2500),
                                          GBM.interaction.depth = seq(2, 8, by = 3),
                                          GBM.shrinkage = c(0.001, 0.01, 0.1),
                                          GBM.n.minobsinnode = 10,
                                          MARS.degree = 1:2, 
                                          MARS.nprune = 2:max(38, 2 * ncol(bm.format@data.env.var) + 1),
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
  combi <- expand.grid(PA = PA.lines, calib = colnames(calib.lines), stringsAsFactors = FALSE)
  combi$name_dataset <- sapply(1:nrow(combi), function(ii) {
    tmp1 <- combi$PA[ii]
    tmp2 <- combi$calib[ii]
    if (tmp2 == "_allData_allRun") tmp2 <- "_allRun"
    if (tmp1 == "_allData_allRun") tmp1 <- "allData"
    paste0("_", tmp1, tmp2)
  })
    
  #   expected_CVnames <- paste0("_allData_RUN", seq_len(ncol(calib.lines)))
  # 
  # if (!is.null(bm.format) && inherits(bm.format, "BIOMOD.formated.data.PA")) {
  #   expected_CVnames <- c(expected_CVnames
  #                         , sapply(1:ncol(bm.format@PA.table)
  #                                  , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
  #                                                        , paste0("_PA", this_PA, "_allRun"))))
  # } 
  
  if (model != "MAXENT") {
    ## check control
    ctrl.train <- caret::trainControl(method = "repeatedcv",
                                      repeats = 3,
                                      number = 10,
                                      summaryFunction = caret::twoClassSummary,
                                      classProbs = TRUE,
                                      returnData = FALSE)
  }
  
  argsval <- foreach(PA.i = combi$PA, calib.i = combi$calib) %do%
    {
      # cat("\n\t\t> Dataset..")
      
      argstmp <- bm.options@args.default
      
      if (model == "MAXNET") {
        warning("No tuning available for that model. Sorry.")
      } else {
        ## 1. SPECIFIC CASE OF MAXENT OR SRE ------------------------------------------------------------
        if (model %in% c("MAXENT", "SRE")) {
          cat("\n\t\t> Tuning parameters...")
          
          ## create dataset ---------------------------------------------------------
          mySpExpl <- get_species_data(bm.format)
          mySpExpl[["_allData_allRun"]] <- TRUE
          mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
          mySpExpl <- mySpExpl[which(mySpExpl[, PA.i] == TRUE), ]
          myResp <- mySpExpl[, 1]
          myExpl <- mySpExpl[, 4:ncol(mySpExpl)]
          mySpExpl[["_allData_allRun"]] <- NULL
          
          
          if (model == "MAXENT") { # ------------------------------------------#
            try(tune.MAXENT <- ENMeval::ENMevaluate(occs = mySpExpl[mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), ],
                                                    bg = mySpExpl[mySpExpl[, 1] == 0 | is.na(mySpExpl[, 1]), ],
                                                    tune.args = list(rm = seq(0.5, 1, 0.5), fc = c("L")),
                                                    algorithm = "maxent.jar",
                                                    partitions = "randomkfold",
                                                    partition.settings = list(kfolds = 10),
                                                    doClamp = TRUE, ## allow to change or not ?
                                                    parallel = TRUE,
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
                    RES <- presence.absence.accuracy(DATA, threshold = as.vector(optimal.thresholds(DATA, opt.methods = 3)[2], mode = "numeric"))
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
          
          ## create dataset
          mySpExpl <- get_species_data(bm.format)
          mySpExpl[["_allData_allRun"]] <- 1
          mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
          mySpExpl <- mySpExpl[which(mySpExpl[, PA.i] == TRUE), ]
          mySpExpl[, 1] <- as.factor(ifelse(mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), "Presence", "Absence"))
          myResp <- mySpExpl[, 1]
          myExpl <- mySpExpl[, 4:ncol(mySpExpl)]
          myExpl[["_allData_allRun"]] <- NULL
          
          ## run tuning -----------------------------------------------------------------------------
          cmd.tuning <- "caret::train(x = myExpl, y = myResp, method = tuning.fun, tuneGrid = tuning.grid,"
          cmd.tuning <- paste0(cmd.tuning, " trControl = ctrl.train, metric = 'ROC',")
          if (tuning.fun %in% c("fda", "rpart")) { ## add weights
            cmd.tuning <- paste0(cmd.tuning, " weights = weights,")
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
            cat("\n\t\t> Tuning parameters...")
            eval(parse(text = paste0("try(tuned.mod <- ", cmd.tuning)))
            
            ## GET tuned parameter values -------------------------------------------------------------
            if (!is.null(tuned.mod)) {
              tmp <- tuned.mod$results
              tmp$TSS <- tmp$Sens + tmp$Spec - 1
              if (model == "XGBOOST") {
                for (param in train.params$params) {
                  argstmp$params[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
                }
              } else {
                for (param in train.params$params) {
                  argstmp[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
                }
              }
              tuning.form <- tuning.grid[which.max(tmp[, metric.eval]), ]
              
              if (model == "CTA") {
                tuning.fun = "rpart2"
                eval(parse(text = paste0("try(tuned.mod <- ", cmd.tuning)))
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
            cat("\n\t\t> Tuning formula...")
            
            cmd.form <- sub("tuneGrid = tuning.grid", "tuneGrid = tuning.form", cmd.tuning)
            cmd.init <- "form = bm_MakeFormula(resp.name = 'resp', expl.var = myExpl, type = typ, interaction.level = intlev),"
            cmd.init <- paste0(cmd.init, " data = cbind(myExpl, resp = myResp),")
            cmd.form <- sub("x = myExpl, y = myResp,", cmd.init, cmd.form)
            
            TMP <- foreach (typ = c('simple', 'quadratic', 'polynomial', 's_smoother'), .combine = "rbind") %:%
              foreach (intlev = 0:(ncol(myExpl) - 1), .combine = "rbind") %do%
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
            argstmp$formula <- TMP[which.max(TMP[, metric.eval]), "formula"]
          }
          
          ## run variable selection -----------------------------------------------------------------
          if (do.stepAIC && (model == "GLM" || 
                             (model == "GAM" && bm.options@package == "gam")) ) {
            cat("\n\t\t> Tuning variables (AIC)...")
            
            if (model == "GLM") {
              glmStart <- glm(as.formula(paste0(bm.format@sp.name, " ~ 1")),
                              data = mySpExpl, 
                              family = argstmp$family, 
                              control = argstmp$control, 
                              weights = weights,
                              mustart = rep(ifelse(!is.null(argstmp$mustart) & nchar(argstmp$mustart) > 0, argstmp$mustart, 0.5), length(myResp)),
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
                                                 weights = weights))          
              tuned.AIC <- NULL
              try(tuned.AIC <- 
                    gam::step.Gam(
                      gamStart,
                      # scope = list(upper = (sub(".*~", "~", argstmp$formula)), lower = ~1),
                      scope = .scope(head(myExpl),
                                     "gam::s", 6),
                      direction = "both",
                      trace = FALSE))
              if (!is.null(tuned.AIC)) { argstmp$formula <- deparse(tuned.AIC$formula) }
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
  
  ##check bm.options ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.options", bm.options, c("BIOMOD.options.default", "BIOMOD.options.dataset"))
  
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
  if (model %in% c("CTA", "FDA", "GLM", "GAM") && is.null(weights)) { 
    weights = rep(1, length(bm.format@data.species))
  }
  
  
  ## get tuning function and parameters ---------------------------------------
  all.fun <- c('avNNet', 'rpart', 'rpart2', 'fda', 'gamSpline', 'bam', 'gam', 'gbm', 'glm', 'earth', 'rf', 'xgbTree')
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
    if (!(model == "GAM" && tuning.fun == "gamSpline")) {
      params.train = params.train[grep(model, names(params.train))]
      .fun_testIfIn(TRUE, "names(params.train)", names(params.train), paste0(model, ".", train.params$params))
      names(params.train) = sub(model, "", names(params.train))
      tuning.grid <- do.call(expand.grid, params.train)
    }
  }
  
  ## get tuning length --------------------------------------------------------
  tuning.length <- 1
  if (model == "CTA") tuning.length <- 30
  if (model == "RF") tuning.length <- min(30, ncol(bm.format@data.env.var))
  
  ## get criteria -------------------------------------------------------------
  if (do.stepAIC) {
    .fun_testIfIn(TRUE, "metric.AIC", metric.AIC, c("AIC", "BIC"))
    if (metric.AIC == "AIC") criteria.AIC <- 2
    if (metric.AIC == "BIC") criteria.AIC <- log(ncol(bm.format@data.env.var))
  }
  
  
  return(list(weights = weights
              , criteria.AIC = criteria.AIC
              , tuning.fun = tuning.fun
              , train.params = train.params
              , tuning.length = tuning.length
              , tuning.grid = tuning.grid))
}


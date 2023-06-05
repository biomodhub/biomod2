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
##' \code{\link{bm_ModelingOptions}} function
##' @param models a \code{vector} containing model names to be tuned, must be among 
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}
##' @param metric.eval a \code{character} corresponding to the evaluation metric used to select 
##' optimal models and tune parameters, must be either \code{ROC} or \code{TSS} 
##' (\emph{maximizing Sensitivity and Specificity})
##' \code{auc.val.avg}, \code{auc.diff.avg}, \code{or.mtp.avg}, \code{or.10p.avg} or \code{AICc}
##' @param ctrl.train global control parameters that can be obtained from the 
##' \code{\link[caret]{trainControl}} function
##' @param ctrl.train.tuneLength (see \code{tuneLength} parameter in \code{\link[caret]{train}})
##' for \code{ANN}
##' @param ANN.MaxNWts an \code{integer} corresponding to the maximum allowable number of weights 
##' for \code{ANN}
##' @param GLM.method a \code{character} corresponding to the classification or regression model 
##' to use for \code{GLM}, \cr 
##' must be \code{glmStepAIC} (see 
##' \url{http://topepo.github.io/caret/train-models-by-tag.html#Generalized_Linear_Model})
##' @param ME.cvmethod a \code{character} corresponding to the method used to partition data for 
##' \code{MAXENT}, \cr must be \code{randomkfold}
##' @param ME.kfolds an \code{integer} corresponding to the number of bins for k-fold 
##' cross-validation for \code{MAXENT}
##' @param weights a \code{vector} of \code{numeric} values corresponding to observation weights
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
##' \code{\link{bm_ModelingOptions}}, \code{\link{BIOMOD_Modeling}}
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


bm_Tuning <- function(model,
                      tuning.fun,
                      do.formula = FALSE,
                      bm.opt.def,
                      bm.format,
                      calib.lines = NULL,
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
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_Tuning.check.args(model = model, tuning.fun = tuning.fun, do.formula = do.formula, bm.format = bm.format
                                , metric.eval = metric.eval, weights = weights, params.train = params.train)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## LOOP OVER CALIB LINES 
  if (is.null(calib.lines)) {
    calib.lines <- data.frame(rep(TRUE, length(bm.format@data.species)))
    colnames(calib.lines) <- "_allData_allRun"
  }
  
  if (model != "MAXENT") {
    ## check control
    ctrl.train <- caret::trainControl(method = "repeatedcv",
                                      repeats = 3,
                                      number = 10,
                                      summaryFunction = caret::twoClassSummary,
                                      classProbs = TRUE,
                                      returnData = FALSE)
  }
  
  argsval <- foreach(calib.i = 1:ncol(calib.lines)) %do%
    {
      argstmp <- bm.opt.def@args.default
      
      ## 1. SPECIFIC CASE OF MAXENT OR SRE ------------------------------------------------------------
      if (model %in% c("MAXENT", "SRE")) {
        cat("\n\t\t> Tuning parameters...")
        
        ## create dataset ---------------------------------------------------------
        mySpExpl <- get_species_data(bm.format)
        mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
        myResp <- mySpExpl[, 1]
        myExpl <- mySpExpl[, 4:ncol(mySpExpl)]
        
        
        if (model == "MAXENT") { # --------------------------------------------------------------------
          try(tune.MAXENT <- ENMeval::ENMevaluate(occs = mySpExpl[mySpExpl[, 1] == 1 & !is.na(mySpExpl[, 1]), ],
                                                  bg = mySpExpl[mySpExpl[, 1] == 0 | is.na(mySpExpl[, 1]), ],
                                                  tune.args = list(rm = seq(0.5, 1, 0.5), fc = c("L")),
                                                  algorithm = "maxent.jar",
                                                  partitions = ME.cvmethod,
                                                  partition.settings = list(kfolds = ME.kfolds),
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
        } else if (model == "SRE") { # ----------------------------------------------------------------
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
        mySpExpl <- mySpExpl[which(calib.lines[, calib.i] == TRUE), ]
        myResp <- mySpExpl[, 1]
        myExpl <- mySpExpl[, 4:ncol(mySpExpl)]
        myResp <- as.factor(ifelse(myResp == 1 & !is.na(myResp), "Presence", "Absence"))
        
        ## run tuning -----------------------------------------------------------------------------
        cmd.tuning <- "caret::train(x = myExpl, y = myResp, method = tuning.fun, tuneGrid = tuning.grid,"
        cmd.tuning <- paste0(cmd.tuning, " trControl = ctrl.train, metric = 'ROC', weights = weights,")
        
        if (model == "ANN") {
          ## Automatically standardize data prior to modeling and prediction
          cmd.tuning <- paste0(cmd.tuning, " preProc = c('center', 'scale'), linout = TRUE, trace = FALSE,")
          cmd.tuning <- paste0(cmd.tuning, " MaxNWts.ANN = ANN.MaxNWts, maxit = ANN.maxit))")
          
        } else if (tuning.fun %in% c("earth", "bam", "fda", "rpart")) { ## remove verbose
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
            for (param in train.params$params) {
              argstmp[[param]] <- tmp[which.max(tmp[, metric.eval]), param]
            }
            tuning.form <- tuning.grid[which.max(tmp[, metric.eval]), ]
          }
        } else { tuning.form <- tuning.grid }
        
        ## run formula selection ------------------------------------------------------------------
        if (do.formula) {
          cat("\n\t\t> Tuning formula...")
          
          cmd.form <- sub("tuneGrid = tuning.grid", "tuneGrid = tuning.form", cmd.tuning)
          cmd.init <- "form = bm_MakeFormula(resp.name = 'resp', expl.var = myExpl, type = typ, interaction.level = intlev),"
          cmd.init <- paste0(cmd.init, " data = cbind(myExpl, resp = myResp),")
          cmd.form <- sub("x = myExpl, y = myResp,", cmd.init, cmd.form)
          
          TMP <- foreach (typ = c("simple", "quadratic", "polynomial", "s_smoother"), .combine = "rbind") %:%
            foreach (intlev = 0:3, .combine = "rbind") %do%
            {
              eval(parse(text = paste0("try(tuned.form <- ", cmd.form)))
              
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
        
        # if (model == "CTA" && !is.null(tuned.mod)) {
        #   tmp = tuned.mod$results
        #   tmp$TSS = tmp$Sens + tmp$Spec - 1
        #   argstmp[["cp"]] <- tmp[which.max(tmp[, metric.eval]), "cp"]
        #   
        #   try(tuned.mod <- caret::train(x = myExpl, 
        #                                 y = myResp,
        #                                 method = "rpart2",
        #                                 tuneGrid = tuning.grid,
        #                                 tuneLength = tuning.length,
        #                                 trControl = ctrl.train,
        #                                 metric = metric.eval,
        #                                 weights = weights))
        #   if (!is.null(tuned.mod)) {
        #     tmp = tuned.mod$results
        #     tmp$TSS = tmp$Sens + tmp$Spec - 1
        #     argstmp[["maxdepth"]] <- tmp[which.max(tmp[, metric.eval]), "maxdepth"]
        #   }
        # }
      }
      return(argstmp)
    }
  names(argsval) <- colnames(calib.lines)
  
  return(argsval)
}



# ---------------------------------------------------------------------------- #

.bm_Tuning.check.args <- function(model, tuning.fun, do.formula, bm.format, metric.eval, weights = NULL, params.train)
{
  ## check model --------------------------------------------------------------
  .fun_testIfIn(TRUE, "model", model, c("ANN", "CTA", "FDA", "GAM", "GBM", "GLM"
                                        , "MARS", "MAXENT", "MAXNET", "RF", "SRE", "XGBOOST"))
  
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
  if (do.formula == TRUE) {
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
  
  .fun_testIfIn(TRUE, "tuning.fun", tuning.fun, c(all.fun, "bm_SRE", "ENMevaluate"))
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


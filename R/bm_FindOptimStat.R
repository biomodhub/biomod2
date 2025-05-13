###################################################################################################
##' @name bm_FindOptimStat
##' @aliases bm_FindOptimStat
##' @aliases bm_CalculateStatBin
##' @aliases bm_CalculateStatAbund
##' @aliases get_optim_value
##' @author Damien Georges
##' 
##' @title Calculate the best score according to a given evaluation method
##'
##' @description This internal \pkg{biomod2} function allows the user to find the threshold to 
##' convert continuous values into binary ones leading to the best score for a given evaluation 
##' metric.
##'
##'
##' @param metric.eval a \code{character} corresponding to the evaluation metric to be used, must 
##' be either \code{AUCroc}, \code{AUCprg}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{OR}, \code{ORSS}, 
##' \code{BOYCE}, \code{MPA} (\emph{binary data}), 
##' \code{RMSE}, \code{MAE}, \code{MSE}, \code{Rsquared}, \code{Rsquared_aj}, \code{Max_error} 
##' (\emph{abundance / count / relative data}), 
##' \code{Accuracy}, \code{Recall}, \code{Precision}, \code{F1} (\emph{multiclass/ordinal data})
##' @param obs a \code{vector} of observed values (binary, \code{0} or \code{1})
##' @param fit a \code{vector} of fitted values (continuous)
##' @param nb.thresh an \code{integer} corresponding to the number of thresholds to be 
##' tested over the range of fitted values
##' @param threshold (\emph{optional, default} \code{NULL}) \cr 
##' A \code{numeric} corresponding to the threshold used to convert the given data
##' @param boyce.bg.env (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing values of 
##' environmental variables (in columns or layers) extracted from the background 
##' (\emph{if presences are to be compared to background instead of absences or 
##' pseudo-absences selected for modeling})
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' @param mpa.perc a \code{numeric} between \code{0} and \code{1} corresponding to the percentage 
##' of correctly classified presences for Minimal Predicted Area (see \code{ecospat.mpa()} in 
##' \pkg{ecospat})
##' @param k an \code{integer} corresponding to the number of environmental variables used in the 
##' model
##' @param misc a \code{matrix} corresponding to a contingency table
##'
##'
##' @return 
##' 
##' A \code{1} row x \code{5} columns \code{data.frame} containing :
##' \itemize{
##'   \item \code{metric.eval} : the chosen evaluation metric
##'   \item \code{cutoff} : the associated cut-off used to transform the continuous values into 
##'   binary
##'   \item \code{sensitivity} : the sensibility obtained on fitted values with this threshold
##'   \item \code{specificity} : the specificity obtained on fitted values with this threshold
##'   \item \code{best.stat} : the best score obtained for the chosen evaluation metric
##' }
##' 
##'
##' @details
##' 
##' \emph{Please refer to \href{https://www.cawcr.gov.au/projects/verification/}{CAWRC website 
##' ("Methods for dichotomous forecasts")} to get detailed description (simple/complex metrics).} \cr
##' Optimal value of each method can be obtained with the \code{\link{get_optim_value}} function.
##' 
##' \describe{
##'   \item{simple}{
##'   \itemize{
##'     \item \code{POD} : Probability of detection (hit rate)
##'     \item \code{FAR} : False alarm ratio
##'     \item \code{POFD} : Probability of false detection (fall-out)
##'     \item \code{SR} : Success ratio
##'     \item \code{ACCURACY} : Accuracy (fraction correct)
##'     \item \code{BIAS} : Bias score (frequency bias)
##'   }
##'   }
##'   \item{complex}{
##'   \itemize{
##'     \item \code{AUCroc} : Area Under Curve of Relative operating characteristic
##'     \item \code{AUCprg} : Area Under Curve of Precision-Recall-Gain curve
##'     \item \code{TSS} : True skill statistic (Hanssen and Kuipers discriminant, Peirce's 
##'     skill score)
##'     \item \code{KAPPA} : Cohen's Kappa (Heidke skill score)
##'     \item \code{OR} : Odds Ratio
##'     \item \code{ORSS} : Odds ratio skill score (Yule's Q)
##'     \item \code{CSI} : Critical success index (threat score)
##'     \item \code{ETS} : Equitable threat score (Gilbert skill score)
##'   }
##'   }
##'   \item{presence-only}{
##'   \itemize{
##'     \item \code{BOYCE} : Boyce index
##'     \item \code{MPA} : Minimal predicted area (cutoff optimizing MPA to predict 90\% of 
##'     presences)
##'   }
##'   }
##'   \item{abundance / count / relative data}{
##'   \itemize{
##'     \item \code{RMSE} : Root Mean Square Error
##'     \item \code{MSE} : Mean Square Error
##'     \item \code{MAE} : Mean Absolute Error
##'     \item \code{Rsquared} : R squared
##'     \item \code{Rsquared_aj} : R squared adjusted
##'     \item \code{Max_error} : Maximum error
##'   }
##'   }
##'   \item{multiclass/ordinal data}{
##'   \itemize{
##'     \item \code{Accuracy} : Accuracy
##'     \item \code{Recall} : Macro average Recall
##'     \item \code{Precision} : Macro average Precision
##'     \item \code{F1} : Macro F1 score
##'   }
##'   }
##' }
##'   
##' Note that if a value is given to \code{threshold}, no optimization will be done, and 
##' only the score for this threshold will be returned.
##' 
##' The Boyce index returns \code{NA} values for \code{SRE} models because it can not be 
##' calculated with binary predictions. \cr This is also the reason why some \code{NA} values 
##' might appear for \code{GLM} models if they did not converge.
##'
##' 
##' @note In order to break dependency loop between packages \pkg{biomod2} and \pkg{ecospat}, 
##' code of \code{ecospat.boyce()} and \code{ecospat.mpa()} in \pkg{ecospat})
##' functions have been copied within this file from version 3.2.2 (august 2022).
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
##' @keywords models options evaluation auc tss boyce mpa
##' 
##'
##' @seealso \code{ecospat.boyce()} and \code{ecospat.mpa()} in \pkg{ecospat}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{bm_RunModelsLoop}}, 
##' \code{\link{BIOMOD_EnsembleModeling}}
##' @family Secondary functions
##' 
##' 
##' @examples
##' ## Generate a binary vector
##' vec.a <- sample(c(0, 1), 100, replace = TRUE)
##' 
##' ## Generate a 0-1000 vector (random drawing)
##' vec.b <- runif(100, min = 0, max = 1000)
##' 
##' ## Generate a 0-1000 vector (biased drawing)
##' BiasedDrawing <- function(x, m1 = 300, sd1 = 200, m2 = 700, sd2 = 200) {
##'   return(ifelse(x < 0.5, rnorm(1, m1, sd1), rnorm(1, m2, sd2)))
##' }
##' vec.c <- sapply(vec.a, BiasedDrawing)
##' vec.c[which(vec.c < 0)] <- 0
##' vec.c[which(vec.c > 1000)] <- 1000
##' 
##' ## Find optimal threshold for a specific evaluation metric
##' bm_FindOptimStat(metric.eval = 'TSS', fit = vec.b, obs = vec.a)
##' bm_FindOptimStat(metric.eval = 'TSS', fit = vec.c, obs = vec.a, nb.thresh = 100)
##' bm_FindOptimStat(metric.eval = 'TSS', fit = vec.c, obs = vec.a, threshold = 280)
##'
##' 
##' @importFrom stats median complete.cases
##' @importFrom pROC roc coords auc
##' @importFrom PresenceAbsence presence.absence.accuracy
##' 
##' @export
##' 
##' 
###################################################################################################


bm_FindOptimStat <- function(metric.eval = 'TSS',
                             obs,
                             fit,
                             nb.thresh = 100,
                             threshold = NULL,
                             boyce.bg.env = NULL,
                             mpa.perc = 0.9,
                             k = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_FindOptimStat.check.args(boyce.bg.env, mpa.perc)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## Check obs and fit 
  if (length(obs) != length(fit)) {
    cat("obs and fit are not the same length  => model evaluation skipped !")
    eval.out <- matrix(NA, 1, 4, dimnames = list(metric.eval, c("best.stat", "cutoff", "sensitivity", "specificity")))
    return(eval.out)
  }
  
  ## remove all unfinite values
  to_keep <- (is.finite(fit) & is.finite(obs))
  obs <- obs[to_keep]
  fit <- fit[to_keep]
  
  ## check some data is still here
  if (!length(obs) | !length(fit)) {
    cat("\n Non finite obs or fit available => model evaluation skipped !")
    eval.out <- matrix(NA, 1, 4, dimnames = list(metric.eval, c("best.stat", "cutoff", "sensitivity", "specificity")))
    return(eval.out)
  }
  
  if (length(unique(obs)) == 1 | length(unique(fit)) == 1) {
    warning("\nObserved or fitted data contains a unique value... Be careful with these model predictions\n", immediate. = TRUE)
  }
  
  
  ## 1. Compute metrics ---------------------------------------------------------------------------
  cutoff <- sensitivity <- specificity <- best.stat <- NA
  abundance_metrics <- c("RMSE", "MSE", "MAE", "AIC","Rsquared", "Rsquared_aj"
                         , "Max_error", "Accuracy", "Recall", "Precision", "F1")
  
  if (!(metric.eval %in% c('AUCroc', "AUCprg", abundance_metrics))) ## BINARY METRICS OTHER THAN AUC -----------
  {
    if (!(metric.eval %in% c('BOYCE', 'MPA'))) ## BINARY METRICS OTHER THAN BOYCE, MPA ------------
    {
      ## A. get threshold values to be tested -----------------------------------
      if (is.null(threshold)) { # test a range of threshold to get the one giving the best score
        
        ## guess fit value scale (e.g. 0-1 for a classic fit or 0-1000 for a biomod2 model fit)
        fit.scale <- .guess_scale(fit)
        
        if (length(unique(fit)) == 1) { ## ONE VALUE predicted
          valToTest <- unique(fit)
          ## add 2 values to test based on mean with 0 and the guessed max of fit (1 or 1000)
          valToTest <- round(c(mean(c(fit.scale["min"], valToTest)),
                               mean(c(fit.scale["max"], valToTest))))
        } else { ## SEVERAL VALUES predicted
          mini <- max(min(fit, na.rm = TRUE), fit.scale["min"], na.rm = TRUE)
          maxi <- min(max(fit, na.rm = TRUE), fit.scale["max"], na.rm = TRUE)
          valToTest <- try(unique(round(c(seq(mini, maxi, length.out = nb.thresh), mini, maxi))))
          if (inherits(valToTest, "try-error")) {
            valToTest <- seq(fit.scale["min"], fit.scale["max"], length.out = nb.thresh)
          }
          # deal with unique value to test case
          if (length(valToTest) < 3) {
            valToTest <- round(c(mean(fit.scale["min"], mini), valToTest, mean(fit.scale["max"], maxi)))
          }
        }
      } else { # test only one value
        valToTest <- threshold
      }
      
      ## B. apply the bm_CalculateStatBin function ------------------------------
      calcStat <- sapply(lapply(valToTest, function(x) {
        return(table(fit >= x, obs))
      }), bm_CalculateStatBin, metric.eval = metric.eval)
      
      ## C. scale obtained scores and find best value ---------------------------
      StatOptimum <- get_optim_value(metric.eval)
      calcStat <- 1 - abs(StatOptimum - calcStat)
      best.stat <- max(calcStat, na.rm = TRUE)
      cutoff <- median(valToTest[which(calcStat == best.stat)]) # if several values are selected
      
    } else if (metric.eval == "BOYCE") ## BOYCE ---------------------------------------------------
    {
      boy <- ecospat.boyce(obs = fit[which(obs == 1)], fit = fit, PEplot = FALSE)
      best.stat <- boy$cor
      if (sum(boy$F.ratio < 1, na.rm = TRUE) > 0) {
        cutoff <- round(boy$HS[max(which(boy$F.ratio < 1))], 0)
      }
    } else if (metric.eval == "MPA") ## MPA -------------------------------------------------------
    {
      cutoff <- ecospat.mpa(fit[which(obs == 1)], perc = mpa.perc)
    }
    
    if (!is.na(cutoff / 1000) & cutoff <= 1000) {
      DATA <- cbind(1:length(fit), obs, fit / 1000)
      DATA[is.na(DATA[, 2]), 2] <- 0
      DATA <- DATA[complete.cases(DATA), ]
      EVAL <- presence.absence.accuracy(DATA, threshold = cutoff / 1000, find.auc = FALSE)
      sensitivity <- EVAL$sensitivity * 100
      specificity <- EVAL$specificity * 100
    }
    
  } else if (metric.eval == 'AUCroc') ## ROC ---------------------------------------------------------
  {
    roc1 <- roc(obs, fit, percent = TRUE, direction = "<", levels = c(0, 1))
    roc1.out <- coords(roc1, "best", ret = c("threshold", "sens", "spec")
                       , transpose = TRUE, best.method = "closest.topleft")
    ## if two optimal values are returned keep only the first one
    if (!is.null(ncol(roc1.out))) { roc1.out <- roc1.out[, 1] }
    
    best.stat <- as.numeric(auc(roc1)) / 100
    cutoff <- as.numeric(roc1.out["threshold"])
    specificity <- as.numeric(roc1.out["specificity"])
    sensitivity <- as.numeric(roc1.out["sensitivity"])
  
  } else if (metric.eval == 'AUCprg') ## PRG ---------------------------------------------------------
  {
    prg1 <- create_prg_curve(labels = obs, pos_scores = fit)
    prg1.out <- best_point_prg(prg1)
    
    best.stat <- calc_auprg(prg1) 
    cutoff <- as.numeric(prg1.out["pos_score"])
    specificity <- as.numeric(prg1.out["TN"]) / (as.numeric(prg1.out["TN"]) + as.numeric(prg1.out["FP"])) * 100
    sensitivity <- as.numeric(prg1.out["TP"]) / (as.numeric(prg1.out["TP"]) + as.numeric(prg1.out["FN"])) * 100
    
  } else if (metric.eval %in% abundance_metrics) ## ABUNDANCE METRICS -----------------------------
  {
    if (metric.eval %in% c("Accuracy", "Recall", "Precision", "F1") &&
        length(levels(obs)) != length(levels(fit))) {
      fit <- factor(fit, levels = levels(obs), ordered = TRUE)
      cat("\n \t\t Careful : some categories are not predicted! ")
    }
    best.stat <- bm_CalculateStatAbun(metric.eval, obs, fit, k)
  }
  
  eval.out <- data.frame(metric.eval, cutoff, sensitivity, specificity, best.stat)
  return(eval.out)
}


###################################################################################################

.bm_FindOptimStat.check.args <- function(boyce.bg.env = NULL, mpa.perc = 0.9)
{
  ## 1. Check boyce.bg.env ----------------------------------------------------
  if (!is.null(boyce.bg.env)) {
    available.types.resp <- c('numeric', 'data.frame', 'matrix', 'Raster',
                              'SpatialPointsDataFrame', 'SpatVector', 'SpatRaster')
    .fun_testIfInherits(TRUE, "boyce.bg.env", boyce.bg.env, available.types.resp)
    
    if (is.numeric(boyce.bg.env)) {
      boyce.bg.env <- as.data.frame(boyce.bg.env)
      colnames(boyce.bg.env) <- expl.var.names[1]
    } else if (inherits(boyce.bg.env, 'Raster')) {
      if (any(raster::is.factor(boyce.bg.env))) {
        boyce.bg.env <- .categorical_stack_to_terra(boyce.bg.env)
      } else {
        boyce.bg.env <- rast(boyce.bg.env)
      }
    }
    if (is.matrix(boyce.bg.env) || inherits(boyce.bg.env, 'SpatRaster') || inherits(boyce.bg.env, 'SpatVector')) {
      boyce.bg.env <- as.data.frame(boyce.bg.env)
    } else if (inherits(boyce.bg.env, 'SpatialPoints')) {
      boyce.bg.env <- as.data.frame(boyce.bg.env@data)
    }
    
    ## remove NA from background data
    if (sum(is.na(boyce.bg.env)) > 0) {
      cat("\n      ! NAs have been automatically removed from boyce.bg.env data")
      boyce.bg.env <- na.omit(boyce.bg.env)
    }
  }
  
  ## 2. Check mpa.perc --------------------------------------------------------
  .fun_testIf01(TRUE, "mpa.perc", mpa.perc)
  
  return(list(boyce.bg.env = boyce.bg.env, mpa.perc = mpa.perc))
}


###################################################################################################

##'
##' @rdname bm_FindOptimStat
##' @export
##'

get_optim_value <- function(metric.eval)
{
  switch(metric.eval
         , 'POD' = 1
         , 'FAR' = 0
         , 'POFD' = 0
         , 'SR' = 1
         , 'ACCURACY' = 1
         , 'BIAS' = 1
         , 'AUCroc' = 1
         , 'AUCprg' = 1
         , 'TSS' = 1
         , 'KAPPA' = 1
         , 'OR' = 1000000
         , 'ORSS' = 1
         , 'CSI' = 1
         , 'ETS' = 1
         , 'BOYCE' = 1
         , 'MPA' = 1
  )
}


###################################################################################################

.guess_scale <- function(fit)
{
  min <- 0
  max <- ifelse(max(fit, na.rm = TRUE) <= 1, 1, 1000)
  out <- c(min, max)
  names(out) <- c("min", "max")
  return(out)
}

.contingency_table_check <- function(misc)
{
  # Contingency table checking
  if (dim(misc)[1] == 1) {
    if (row.names(misc)[1] == "FALSE") {
      misc <- rbind(misc, c(0, 0))
      rownames(misc) <- c('FALSE', 'TRUE')
    } else {
      a <- misc
      misc <- c(0, 0)
      misc <- rbind(misc, a)
      rownames(misc) <- c('FALSE', 'TRUE')
    }
  }
  
  if (ncol(misc) != 2 | nrow(misc) != 2) {
    misc = matrix(0, ncol = 2, nrow = 2, dimnames = list(c('FALSE', 'TRUE'), c('0', '1')))
  }
  
  if ((sum(colnames(misc) %in% c('FALSE', 'TRUE', '0', '1')) < 2) | 
      (sum(rownames(misc) %in% c('FALSE', 'TRUE', '0', '1')) < 2) ) {
    stop("Unavailable contingency table given")
  }
  
  if ('0' %in% rownames(misc)) { rownames(misc)[which(rownames(misc) == '0')] <- 'FALSE' }
  if ('1' %in% rownames(misc)) { rownames(misc)[which(rownames(misc) == '1')] <- 'TRUE' }
  
  return(misc)
}

.contingency_table_ordinal <- function(obs, fit)
{
  nblevels <- length(levels(obs))
  m <- matrix(0,nrow = nblevels, ncol = nblevels, dimnames = list(levels(obs), levels(obs)))
  for (i in 1:length(obs)) {
    m[obs[i], fit[i]] <- m[obs[i], fit[i]] + 1
  }
  return(m)
}


##'
##' @rdname bm_FindOptimStat
##' @export
##'

bm_CalculateStatBin <- function(misc, metric.eval = 'TSS')
{
  ## check contingency table
  misc <- .contingency_table_check(misc)
  
  ## calculate basic classification information -------------------------------
  hits <- misc['TRUE', '1'] ## true positives
  misses <- misc['FALSE', '1'] ##  false positives
  false_alarms <- misc['TRUE', '0'] ## false negatives 
  correct_negatives <- misc['FALSE', '0'] ## true negatives
  
  total <- sum(misc)
  forecast_1 <- sum(misc['TRUE', ]) ## hits + false_alarms
  forecast_0 <- sum(misc['FALSE', ]) ## misses + correct_negatives
  observed_1 <- sum(misc[, '1']) ## hits + misses
  observed_0 <- sum(misc[, '0']) ## false_alarms + correct_negatives
  
  ## calculate chosen evaluation metric ---------------------------------------
  out = switch(metric.eval
               , 'POD' =  hits / observed_1
               , 'POFD' = false_alarms / observed_0
               , 'FAR' = false_alarms / forecast_1
               , 'SR' = hits / forecast_1
               , 'ACCURACY' = (hits + correct_negatives) / total
               , 'BIAS' = forecast_1 / observed_1
               
               , 'TSS' = (hits / observed_1) + (correct_negatives / observed_0) - 1
               , 'KAPPA' = { ## PAREIL ?
                 Po <- (1 / total) * (hits + correct_negatives)
                 Pe <- ((1 / total) ^ 2) * ((forecast_1 * observed_1) + (forecast_0 * observed_0))
                 return((Po - Pe) / (1 - Pe))
               }
               , 'OR' = (hits * correct_negatives) / (misses * false_alarms)
               , 'ORSS' = (hits * correct_negatives - misses * false_alarms) / (hits * correct_negatives + misses * false_alarms)
               , 'CSI' = hits / (observed_1 + false_alarms)
               , 'ETS' = {
                 hits_rand <- (observed_1 * forecast_1) / total
                 return((hits - hits_rand) / (observed_1 + false_alarms - hits_rand))
               }
  )
}

##'
##' @rdname bm_FindOptimStat
##' @export
##'

bm_CalculateStatAbun <- function(metric.eval, obs, fit, k)
{
  ## get contingency table
  if (metric.eval %in% c("Accuracy", "Recall", "Precision", "F1")) {
    m <- .contingency_table_ordinal(obs, fit)
  }
  
  ## calculate chosen evaluation metric ---------------------------------------
  out = switch(metric.eval
               , 'RMSE' = sqrt(mean((obs - fit) ^ 2))
               , 'MSE' = mean((obs - fit) ^ 2)
               , 'MAE' = mean(abs(obs - fit))
               , 'Rsquared' = cor(obs,fit) ^ 2
               , 'Rsquared_aj' = 1 - (1 - cor(obs, fit) ^ 2) * (length(obs) - 1) / (length(obs) - k - 1) 
               , 'Max_error' = max(abs(obs - fit), na.rm = TRUE)
               , 'Accuracy' = sum(diag(m)) / sum(m)
               , 'Recall' = sum(diag(m) / colSums(m), na.rm = TRUE) / nrow(m)
               , 'Precision' = sum(diag(m) / rowSums(m), na.rm = TRUE) / nrow(m)
               , 'F1' = {
                 p <- sum(diag(m) / rowSums(m), na.rm = TRUE) / nrow(m)
                 r <- sum(diag(m) / colSums(m), na.rm = TRUE) / nrow(m)
                 return(2 * p * r / (p + r))
               }
  )
}


###################################################################################################
## FROM ECOSPAT PACKAGE VERSION 3.2.2 (august 2022)

#### Calculating Boyce index as in Hirzel et al. 2006
# fit: A vector or Raster-Layer containing the predicted suitability values 
# obs: A vector containing the predicted suitability values or xy-coordinates 
# (if fit is a Raster-Layer) of the validation points (presence records)
# nclass : number of classes or vector with classes threshold. If nclass=0, Boyce index is 
# calculated with a moving window (see next parameters)
# windows.w : width of the moving window (by default 1/10 of the suitability range)
# res : resolution of the moving window (by default 100 focals)
# PEplot : if True, plot the predicted to expected ratio along the suitability class
# rm.duplicate : if TRUE, the correlation exclude successive duplicated values
# method : correlation method used to compute the boyce index


ecospat.boyce <- function(fit, obs, nclass = 0, window.w = "default", res = 100, 
                          PEplot = TRUE, rm.duplicate = TRUE, method = 'spearman')
{
  
  #### internal function calculating predicted-to-expected ratio for each class-interval
  boycei <- function(interval, obs, fit) {
    pi <- sum(as.numeric(obs >= interval[1] & obs <= interval[2])) / length(obs)
    ei <- sum(as.numeric(fit >= interval[1] & fit <= interval[2])) / length(fit)
    return(round(pi/ei,10))
  }
  
  # here, ecospat.boyce should not receive RasterLayer
  # if (inherits(fit,"RasterLayer")) {
  #   if (is.data.frame(obs) || is.matrix(obs)) {
  #     obs <- extract(fit, obs)
  #   }
  #   fit <- getValues(fit)
  #   fit <- fit[!is.na(fit)]
  # }
  
  mini <- min(fit,obs)
  maxi <- max(fit,obs)
  
  if(length(nclass)==1){
    if (nclass == 0) { #moving window
      if (window.w == "default") {window.w <- (max(fit) - min(fit))/10}
      vec.mov <- seq(from = mini, to = maxi - window.w, by = (maxi - mini - window.w)/res)
      vec.mov[res + 1] <- vec.mov[res + 1] + 1  #Trick to avoid error with closed interval in R
      interval <- cbind(vec.mov, vec.mov + window.w)
    } else{ #window based on nb of class
      vec.mov <- seq(from = mini, to = maxi, by = (maxi - mini)/nclass)
      interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
    }
  } else{ #user defined window
    vec.mov <- c(mini, sort(nclass[!nclass>maxi|nclass<mini]))
    interval <- cbind(vec.mov, c(vec.mov[-1], maxi))
  }
  
  f <- apply(interval, 1, boycei, obs, fit)
  to.keep <- which(f != "NaN")  # index to keep no NaN data
  f <- f[to.keep]
  if (length(f) < 2) {
    b <- NA  #at least two points are necessary to draw a correlation
  } else {
    r<-1:length(f)
    if(rm.duplicate == TRUE){
      r <- c(1:length(f))[f != c( f[-1],TRUE)]  #index to remove successive duplicates
    }
    b <- cor(f[r], vec.mov[to.keep][r], method = method)  # calculation of the correlation (i.e. Boyce index) after removing successive duplicated values
  }
  HS <- apply(interval, 1, sum)/2  # mean habitat suitability in the moving window
  if(length(nclass)==1 & nclass == 0) {
    HS[length(HS)] <- HS[length(HS)] - 1  #Correction of the 'trick' to deal with closed interval
  }
  HS <- HS[to.keep]  #exclude the NaN
  if (PEplot == TRUE) {
    plot(HS, f, xlab = "Habitat suitability", ylab = "Predicted/Expected ratio", col = "grey", cex = 0.75)
    points(HS[r], f[r], pch = 19, cex = 0.75)
  }
  return(list(F.ratio = f, cor = round(b, 3), HS = HS))
}


###################################################################################################
## FROM ECOSPAT PACKAGE VERSION 3.2.2 (august 2022)

## This function calculates the Minimal Predicted Area.
## 
## R. Patin, Nov. 2022 : the function does not return minimal predicted area, but
## rather the quantile of suitability corresponding to perc% of presence
##  correctly predicted
##  
## FUNCTION'S ARGUMENTS
## Pred:      numeric or RasterLayer .predicted suitabilities from a SDM prediction
## Sp.occ.xy: xy-coordinates of the species (if Pred is a RasterLayer)
## perc:      Percentage of Sp.occ.xy that should be encompassed by the binary map.

## Details:

## Value
# Returns the Minimal Predicted Area

## Author(s)
# Frank Breiner

## References
# Engler, R., A. Guisan, and L. Rechsteiner. 2004. An improved approach for predicting the 
# distribution of rare and endangered species from occurrence and pseudo-absence data. 
# Journal of Applied Ecology 41:263-274.


ecospat.mpa <- function(Pred, Sp.occ.xy = NULL, perc = 0.9)
{
  perc <- 1 - perc
  if (!is.null(Sp.occ.xy)) {
    Pred <- extract(Pred, Sp.occ.xy)
  }
  round(quantile(Pred, probs = perc,na.rm = TRUE), 3)
}

### EXAMPLE

# obs <- (ecospat.testData$glm_Saxifraga_oppositifolia
# [which(ecospat.testData$Saxifraga_oppositifolia==1)]) ecospat.mpa(obs) ecospat.mpa(obs,perc=1) ##
# 100% of the presences encompassed


## Example script for using Ensemble Small Models ESMs according to Lomba et al. 2010 Written by
## Frank Breiner Swiss Federal Research Institute WSL, July 2014.

## References: Lomba, A., Pellissier, L., Randin, C.F., Vicente, J., Moreira, F., Honrado, J. &
## Guisan, A. (2010). Overcoming the rare species modeling paradox: A novel hierarchical framework
## applied to an Iberian endemic plant. Biological Conservation, 143, 2647-2657.

###################################################################################################
## FROM PRG PACKAGE VERSION 0.5.1 (April 2025)

# Precision Gain
#
# This function calculates Precision Gain from the entries of the contingency table: number of true positives (TP), false negatives (FN), false positives (FP), and true negatives (TN). More information on Precision-Recall-Gain curves and how to cite this work is available at http://www.cs.bris.ac.uk/~flach/PRGcurves/.
# @param TP number of true positives, can be a vector
# @param FN number of false negatives, can be a vector
# @param FP number of false positives, can be a vector
# @param TN number of true negatives, can be a vector
# @return Precision Gain (a numeric value less than or equal to 1; or -Inf or NaN, see the details below)
# @details Precision Gain (PrecGain) quantifies by how much precision is improved over the default precision of the always positive predictor, equal to the proportion of positives (pi). PrecGain=1 stands for maximal improvement (Prec=1) and PrecGain=0 stands for no improvement (Prec=pi). If Prec=0, then PrecGain=-Inf. It can happen that PrecGain=NaN, for instance if there are no positives (TP=0 and FN=0) and TN>0.
# @examples 
# precision_gain(3,0,1,2)
# # [1] 0.6666667
# TP = c(0,2,3)
# FN = 3-TP
# FP = c(0,1,2)
# TN = 2-FP
# precision_gain(TP,FN,FP,TN)
# # [1]  NaN 0.25 0.00
precision_gain = function(TP,FN,FP,TN) {
  n_pos = TP+FN
  n_neg = FP+TN
  prec_gain = 1-(n_pos*FP)/(n_neg*TP)
  prec_gain[TN+FN==0] = 0
  return(prec_gain)
}

# Recall Gain
#
# This function calculates Recall Gain from the entries of the contingency table: number of true positives (TP), false negatives (FN), false positives (FP), and true negatives (TN). More information on Precision-Recall-Gain curves and how to cite this work is available at http://www.cs.bris.ac.uk/~flach/PRGcurves/.
# @param TP number of true positives, can be a vector
# @param FN number of false negatives, can be a vector
# @param FP number of false positives, can be a vector
# @param TN number of true negatives, can be a vector
# @return Recall Gain (a numeric value less than or equal to 1; or -Inf or NaN, see the details below)
# @details Recall Gain (RecGain) quantifies by how much recall is improved over the recall equal to the proportion of positives (pi). RecGain=1 stands for maximal improvement (Rec=1) and RecGain=0 stands for no improvement (Rec=pi). If Rec=0, then RecGain=-Inf. It can happen that RecGain=NaN, for instance if there are no negatives (FP=0 and TN=0) and FN>0 and TP=0.
# @examples 
# recall_gain(3,0,1,2)
# # [1] 1
# TP = c(0,2,3)
# FN = 3-TP
# FP = c(0,1,2)
# TN = 2-FP
# recall_gain(TP,FN,FP,TN)
# # [1]  -Inf 0.25 1.00
recall_gain = function(TP,FN,FP,TN) {
  n_pos = TP+FN
  n_neg = FP+TN
  rg = 1-(n_pos*FN)/(n_neg*TP)
  rg[TN+FN==0] = 1
  return(rg)
}

# create a table of segments
.create.segments = function(labels,pos_scores,neg_scores) {
  # reorder labels and pos_scores by decreasing pos_scores, using increasing neg_scores in breaking ties
  new_order = order(pos_scores,-neg_scores,decreasing=TRUE)
  labels = labels[new_order]
  pos_scores = pos_scores[new_order]
  neg_scores = neg_scores[new_order]
  # create a table of segments
  segments = data.frame(pos_score=NA,neg_score=NA,pos_count=0,neg_count=rep(0,length(labels)))
  j = 0
  for (i in seq_along(labels)) {
    if ((i==1)||(pos_scores[i-1]!=pos_scores[i])||(neg_scores[i-1]!=neg_scores[i])) {
      j = j + 1
      segments$pos_score[j] = pos_scores[i]
      segments$neg_score[j] = neg_scores[i]
    }
    segments[j,4-labels[i]] = segments[j,4-labels[i]] + 1
  }
  segments = segments[1:j,]
  return(segments)
}

.create.crossing.points = function(points,n_pos,n_neg) {
  n = n_pos+n_neg
  points$is_crossing = 0
  # introduce a crossing point at the crossing through the y-axis
  j = min(which(points$recall_gain>=0))
  if (points$recall_gain[j]>0) { # otherwise there is a point on the boundary and no need for a crossing point
    delta = points[j,,drop=FALSE]-points[j-1,,drop=FALSE]
    if (delta$TP>0) {
      alpha = (n_pos*n_pos/n-points$TP[j-1])/delta$TP
    } else {
      alpha = 0.5
    }
    new_point = points[j-1,,drop=FALSE] + alpha*delta
    new_point$precision_gain = precision_gain(new_point$TP,new_point$FN,new_point$FP,new_point$TN)
    new_point$recall_gain = 0
    new_point$is_crossing = 1
    points = rbind(points,new_point)
    points = points[order(points$index),,drop=FALSE]
  }   
  # now introduce crossing points at the crossings through the non-negative part of the x-axis
  crossings = data.frame()
  x = points$recall_gain
  y = points$precision_gain
  for (i in which((c(y,0)*c(0,y)<0)&(c(1,x)>=0))) {
    cross_x = x[i-1]+(-y[i-1])/(y[i]-y[i-1])*(x[i]-x[i-1])
    delta = points[i,,drop=FALSE]-points[i-1,,drop=FALSE]
    if (delta$TP>0) {
      alpha = (n_pos*n_pos/(n-n_neg*cross_x)-points$TP[i-1])/delta$TP
    } else {
      alpha = (n_neg/n_pos*points$TP[i-1]-points$FP[i-1])/delta$FP
    }
    new_point = points[i-1,,drop=FALSE] + alpha*delta
    new_point$precision_gain = 0
    new_point$recall_gain = recall_gain(new_point$TP,new_point$FN,new_point$FP,new_point$TN)
    new_point$is_crossing = 1
    crossings = rbind(crossings,new_point)
  }
  # add crossing points to the 'points' data frame
  points = rbind(points,crossings)
  points = points[order(points$index,points$recall_gain),2:ncol(points),drop=FALSE]
  rownames(points) = NULL
  points$in_unit_square = 1
  points$in_unit_square[points$recall_gain<0] = 0
  points$in_unit_square[points$precision_gain<0] = 0
  return(points)
}


# Precision-Recall-Gain curve
#
# This function creates the Precision-Recall-Gain curve from the vector of labels and vector of scores where higher score indicates a higher probability to be positive. More information on Precision-Recall-Gain curves and how to cite this work is available at http://www.cs.bris.ac.uk/~flach/PRGcurves/.
# @param labels a vector of labels, where 1 marks positives and 0 or -1 marks negatives
# @param pos_scores vector of scores for the positive class, where a higher score indicates a higher probability to be a positive
# @param neg_scores vector of scores for the negative class, where a higher score indicates a higher probability to be a negative (by default, equal to -pos_scores)
# @return A data.frame which lists the points on the PRG curve with the following columns: pos_score, neg_score, TP, FP, FN, TN, precision_gain, recall_gain, is_crossing and in_unit_square. All the points are listed in the order of increasing thresholds on the score to be positive (the ties are broken by decreasing thresholds on the score to be negative).
# @details The PRG-curve is built by considering all possible score thresholds, starting from -Inf and then using all scores that are present in the given data. The results are presented as a data.frame which includes the following columns: pos_score, neg_score, TP, FP, FN, TN, precision_gain, recall_gain, is_crossing and in_unit_square. The resulting points include the points where the PRG curve crosses the y-axis and the positive half of the x-axis. The added points have is_crossing=1 whereas the actual PRG points have is_crossing=0. To help in visualisation and calculation of the area under the curve the value in_unit_square=1 marks that the point is within the unit square [0,1]x[0,1], and otherwise, in_unit_square=0.
# @examples
# create_prg_curve(c(1,1,0,0,1,0),c(0.8,0.8,0.6,0.4,0.4,0.2))
# #   pos_score neg_score  TP FP  FN TN precision_gain recall_gain is_crossing in_unit_square
# # 1      -Inf       Inf 0.0  0 3.0  3            NaN        -Inf           0              0
# # 2       NaN       NaN 1.5  0 1.5  3      1.0000000         0.0           1              1
# # 3       0.8      -0.8 2.0  0 1.0  3      1.0000000         0.5           0              1
# # 4       0.6      -0.6 2.0  1 1.0  2      0.5000000         0.5           0              1
# # 5       0.4      -0.4 3.0  2 0.0  1      0.3333333         1.0           0              1
# # 6       0.2      -0.2 3.0  3 0.0  0      0.0000000         1.0           0              1
#
# create_prg_curve(c(1,1,0,0,1,0),c(1,1,1,1e-20,1e-40,1e-60),c(0,1e-20,1e-20,1,1,1))
# #   pos_score neg_score  TP  FP  FN  TN precision_gain recall_gain is_crossing in_unit_square
# # 1      -Inf       Inf 0.0 0.0 3.0 3.0            NaN        -Inf           0              0
# # 2     1e+00     0e+00 1.0 0.0 2.0 3.0      1.0000000        -1.0           0              0
# # 3     1e+00     5e-21 1.5 0.5 1.5 2.5      0.6666667         0.0           1              1
# # 4     1e+00     1e-20 2.0 1.0 1.0 2.0      0.5000000         0.5           0              1
# # 5     1e-20     1e+00 2.0 2.0 1.0 1.0      0.0000000         0.5           0              1
# # 6     1e-40     1e+00 3.0 2.0 0.0 1.0      0.3333333         1.0           0              1
# # 7     1e-60     1e+00 3.0 3.0 0.0 0.0      0.0000000         1.0           0              1
create_prg_curve = function(labels,pos_scores,neg_scores=-pos_scores) { 
  create_crossing_points = TRUE
  n = length(labels)
  n_pos = sum(labels)
  n_neg = n - n_pos
  # convert negative labels into 0s
  labels = 1*(labels==1) 
  segments = .create.segments(labels,pos_scores,neg_scores)
  # calculate recall gains and precision gains for all thresholds
  points = data.frame(index=1:(1+nrow(segments)))
  points$pos_score=c(Inf,segments$pos_score)
  points$neg_score=c(-Inf,segments$neg_score)
  points$TP = c(0,cumsum(segments$pos_count))
  points$FP = c(0,cumsum(segments$neg_count))
  points$FN = n_pos-points$TP
  points$TN = n_neg-points$FP
  points$precision_gain = precision_gain(points$TP,points$FN,points$FP,points$TN)
  points$recall_gain = recall_gain(points$TP,points$FN,points$FP,points$TN)
  if (create_crossing_points) {
    points = .create.crossing.points(points,n_pos,n_neg)
  } else {
    points = points[,2:ncol(points)]
  }
  return(points)
}

# Calculate area under the Precision-Recall-Gain curve
#
# This function calculates the area under the Precision-Recall-Gain curve from the results of the function create_prg_curve. More information on Precision-Recall-Gain curves and how to cite this work is available at http://www.cs.bris.ac.uk/~flach/PRGcurves/. 
# @param prg_curve the data structure resulting from the function create_prg_curve
# @return A numeric value representing the area under the Precision-Recall-Gain curve.
# @details This function calculates the area under the Precision-Recall-Gain curve, taking into account only the part of the curve with non-negative recall gain. The regions with negative precision gain (PRG-curve under the x-axis) contribute as negative area.
# @examples
# calc_auprg(create_prg_curve(c(1,1,0,0,1,0),c(0.8,0.8,0.6,0.4,0.4,0.2)))
# # [1] 0.7083333
calc_auprg = function(prg_curve) {
  area = 0
  for (i in 2:nrow(prg_curve)) {
    if (!is.na(prg_curve$recall_gain[i-1]) && (prg_curve$recall_gain[i-1]>=0)) {
      width = prg_curve$recall_gain[i]-prg_curve$recall_gain[i-1]
      height = (prg_curve$precision_gain[i]+prg_curve$precision_gain[i-1])/2
      area = area + width*height
    }
  }
  return(area)
}


#######################################################################################################

#Best stat for AUCprg

# This functions calculates the optimized point of the PRG curves 
# Author : Helene Blancheteau

best_point_prg <- function(prg_curve) {
  
  #Distance of the optimal point coord(1,1)
  distance <- (prg_curve$recall_gain - 1)^2 + (prg_curve$precision_gain - 1)^2
  distance <- sqrt(distance)
  
  optimal_point <- as.vector(prg_curve[which.min(distance),])
  return(optimal_point)
}

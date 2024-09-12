###################################################################################################
##' @name bm_FindOptimStat
##' @aliases bm_FindOptimStat
##' @aliases bm_CalculateStat
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
##' be either \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, 
##' \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, 
##' \code{BOYCE}, \code{MPA}
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
##'     \item \code{ROC} : Relative operating characteristic
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
##' }
##'   
##' Optimal value of each method can be obtained with the \code{\link{get_optim_value}} function. \cr
##' \emph{Please refer to the \href{https://www.cawcr.gov.au/projects/verification/}{CAWRC website 
##' (section "Methods for dichotomous forecasts")} to get detailed description of each metric.}
##' 
##' Note that if a value is given to \code{threshold}, no optimization will be done., and 
##' only the score for this threshold will be returned.
##' 
##' The Boyce index returns \code{NA} values for \code{SRE} models because it can not be 
##' calculated with binary predictions. \cr This is also the reason why some \code{NA} values 
##' might appear for \code{GLM} models if they do not converge.
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
                             mpa.perc = 0.9)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_FindOptimStat.check.args(threshold, boyce.bg.env, mpa.perc)
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
    cat("Non finite obs or fit available => model evaluation skipped !")
    eval.out <- matrix(NA, 1, 4, dimnames = list(metric.eval, c("best.stat", "cutoff", "sensitivity", "specificity")))
    return(eval.out)
  }
  
  if (length(unique(obs)) == 1 | length(unique(fit)) == 1) {
    warning("\nObserved or fitted data contains a unique value... Be careful with this models predictions\n", immediate. = TRUE)
  }
  
  
  if (!(metric.eval %in% c('ROC', 'BOYCE', 'MPA')))
  { ## for all evaluation metrics other than ROC, BOYCE, MPA --------------------------------------
    
    ## 1. get threshold values to be tested -----------------------------------
    if (is.null(threshold)) { # test a range of threshold to get the one giving the best score
      
      ## guess fit value scale (e.g. 0-1 for a classic fit or 0-1000 for a biomod2 model fit)
      fit.scale <- .guess_scale(fit)
      
      if (length(unique(fit)) == 1) {
        valToTest <- unique(fit)
        ## add 2 values to test based on mean with 0 and the guessed max of fit (1 or 1000)
        valToTest <- round(c(mean(c(fit.scale["min"], valToTest)),
                             mean(c(fit.scale["max"], valToTest))))
      } else {
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
    
    ## 2. apply the bm_CalculateStat function ---------------------------------
    calcStat <- sapply(lapply(valToTest, function(x) {
      return(table(fit > x, obs))
    }), bm_CalculateStat, metric.eval = metric.eval)
    
    ## 3. scale obtained scores and find best value ---------------------------
    StatOptimum <- get_optim_value(metric.eval)
    calcStat <- 1 - abs(StatOptimum - calcStat)
    best.stat <- max(calcStat, na.rm = TRUE)
    cutoff <- median(valToTest[which(calcStat == best.stat)]) # if several values are selected
    
    misc <- table(fit >= cutoff, obs)
    misc <- .contingency_table_check(misc)
    true.pos <- misc['TRUE', '1']
    true.neg <- misc['FALSE', '0']
    specificity <- (true.neg * 100) / sum(misc[, '0'])
    sensitivity <- (true.pos * 100) / sum(misc[, '1'])
    
  } else if (metric.eval == 'ROC') { ## specific procedure for ROC value --------------------------
    roc1 <- roc(obs, fit, percent = TRUE, direction = "<", levels = c(0, 1))
    roc1.out <- coords(roc1, "best", ret = c("threshold", "sens", "spec"), transpose = TRUE)
    ## if two optimal values are returned keep only the first one
    if (!is.null(ncol(roc1.out))) { roc1.out <- roc1.out[, 1] }
    best.stat <- as.numeric(auc(roc1)) / 100
    cutoff <- as.numeric(roc1.out["threshold"])
    specificity <- as.numeric(roc1.out["specificity"])
    sensitivity <- as.numeric(roc1.out["sensitivity"])
    
  } else { ## specific procedure for BOYCE, MPA values --------------------------------------------
    ## Prepare table to compute evaluation scores
    DATA <- cbind(1:length(fit), obs, fit / 1000)
    DATA[is.na(DATA[, 2]), 2] <- 0
    DATA <- DATA[complete.cases(DATA), ]
    cutoff <- sensitivity <- specificity <- best.stat <- NA
    
    if (metric.eval == "BOYCE") { ## Boyce index
      boy <- ecospat.boyce(obs = fit[which(obs == 1)], fit = fit, PEplot = FALSE)
      best.stat <- boy$cor
      if (sum(boy$F.ratio < 1, na.rm = TRUE) > 0) {
        cutoff <- round(boy$HS[max(which(boy$F.ratio < 1))], 0)
      }
    } else if (metric.eval == "MPA") {
      cutoff <- ecospat.mpa(fit[which(obs == 1)], perc = mpa.perc)
    }
    if (!is.na(cutoff / 1000) & cutoff <= 1000) {
      EVAL <- presence.absence.accuracy(DATA, threshold = cutoff / 1000)
      sensitivity <- EVAL$sensitivity * 100
      specificity <- EVAL$specificity * 100
    }
  }
  
  eval.out <- data.frame(metric.eval, cutoff, sensitivity, specificity, best.stat)
  return(eval.out)
}


# Check Arguments ---------------------------------------------------------------------------------
.bm_FindOptimStat.check.args <- function(threshold, boyce.bg.env = NULL, mpa.perc = 0.9)
{
  ## 1. Check boyce.bg.env -----------------------------------------------------------
  if (!is.null(boyce.bg.env)) {
    available.types.resp <- c('numeric', 'data.frame', 'matrix', 'Raster',
                              'SpatialPointsDataFrame', 'SpatVector', 'SpatRaster')
    .fun_testIfInherits(TRUE, "boyce.bg.env", boyce.bg.env, available.types.resp)
    
    if(is.numeric(boyce.bg.env)) {
      boyce.bg.env <- as.data.frame(boyce.bg.env)
      colnames(boyce.bg.env) <- expl.var.names[1]
    } else if (inherits(boyce.bg.env, 'Raster')) {
      if(any(raster::is.factor(boyce.bg.env))){
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
  
  ## 2. Check perc -----------------------------------------------------------
  .fun_testIf01(TRUE, "mpa.perc", mpa.perc)
  
  return(list(threshold = threshold
              , boyce.bg.env = boyce.bg.env
              , mpa.perc = mpa.perc))
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
         , 'ROC' = 1
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

##'
##' @rdname bm_FindOptimStat
##' @export
##'

bm_CalculateStat <- function(misc, metric.eval = 'TSS')
{
  ## check contagency table
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


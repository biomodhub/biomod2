###################################################################################################
##' @name bm_FindOptimStat
##' @aliases bm_FindOptimStat
##' @aliases bm_CalculateStat
##' @aliases get_optim_value
##' @author Damien Georges
##' 
##' @title Calculate the best score according to a given evaluation method
##'
##' @description 
##' 
##' This internal \pkg{biomod2} function allows the user to find the threshold to convert 
##' continuous values into binary ones leading to the best score for a given evaluation metric.
##'
##' @param Stat a \code{character} corresponding to the evaluation metric to be used, must be 
##' either \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR} 
##' or \code{ORSS}
##' @param Fit a \code{vector} of fitted values (continuous)
##' @param Obs a \code{vector} of observed values (binary, \code{0} or \code{1})
##' @param Nb.thresh.test an \code{integer} corresponding to the number of thresholds to be 
##' tested over the range of fitted values
##' @param Fixed.thresh (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the only threshold value to be tested
##'
##'
##' @return 
##' 
##' A \code{1} row x \code{4} columns \code{matrix} containing :
##' \itemize{
##'   \item{\code{best.iter}}{ : the best score obtained for the chosen evaluation metric}
##'   \item{\code{cutoff}}{ : the associated cut-off used to transform the continuous values into 
##'   binary}
##'   \item{\code{sensibility}}{ : the sensibility obtained on fitted values with this threshold}
##'   \item{\code{specificity}}{ : the specificity obtained on fitted values with this threshold}
##' }
##' 
##'
##' @details
##' 
##' \emph{Please refer to \code{\link{BIOMOD_Modeling}} to get more information about these 
##' evaluation metrics.}
##' 
##' Note that if a value is given to \code{Fixed.thresh}, no optimisation will be done., and 
##' only the score for this threshold will be returned.
##'
##'
##' @keywords models, options, evaluation
##' 
##'
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}
##'
##'
##' @examples
##' 
##' ## generate a binary vector
##' a <- sample(c(0, 1), 100, replace = TRUE)
##'
##' ## random drawing
##' b <- runif(100, min = 0, max = 1000)
##' bm_FindOptimStat(Stat = 'TSS', Fit = b, Obs = a)
##'
##' ## biased drawing
##' BiasedDrawing <- function(x, m1 = 300, sd1 = 200, m2 = 700, sd2 = 200) {
##'   return(ifelse(x < 0.5, rnorm(1, m1, sd1), rnorm(1, m2, sd2)))
##' }
##' c <- sapply(a, BiasedDrawing)
##' bm_FindOptimStat(Stat = 'TSS', Fit = c, Obs = a, Nb.thresh.test = 100)
##'
##' 
##' @importFrom pROC roc coords auc
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


bm_FindOptimStat <- function(Stat = 'TSS',
                             Fit,
                             Obs,
                             Nb.thresh.test = 100,
                             Fixed.thresh = NULL)
{
  ## remove all unfinite values
  to_keep <- (is.finite(Fit) & is.finite(Obs))
  Fit <- Fit[to_keep]
  Obs <- Obs[to_keep]
  
  ## check some data is still here
  if (!length(Obs) | !length(Fit)) {
    cat("Non finite obs or fit available => model evaluation skipped !")
    eval.out <- matrix(NA, 1, 4, dimnames = list(Stat, c("best.stat", "cutoff", "sensitivity", "specificity")))
    return(eval.out)
  }
  
  if (length(unique(Obs)) == 1 | length(unique(Fit)) == 1) {
    warning("\nObserved or fitted data contains a unique value... Be careful with this models predictions\n", immediate. = TRUE)
  }
  
  
  if (Stat != 'ROC') { ## for all evaluation metrics other than ROC -------------------------------
    StatOptimum <- get_optim_value(Stat)
    
    ## 1. get threshold values to be tested -----------------------------------
    if (is.null(Fixed.thresh)) { # test a range of threshold to get the one giving the best score
      
      ## guess fit value scale (e.g. 0-1 for a classic fit or 0-1000 for a biomod2 model fit)
      fit.scale <- .guess_scale(Fit)
      
      if (length(unique(Fit)) == 1) {
        valToTest <- unique(Fit)
        ## add 2 values to test based on mean with 0 and the guessed max of Fit (1 or 1000)
        valToTest <- round(c(mean(c(fit.scale["min"], valToTest)),
                             mean(c(fit.scale["max"], valToTest))))
      } else {
        mini <- max(min(Fit, na.rm = TRUE), fit.scale["min"], na.rm = TRUE)
        maxi <- min(max(Fit, na.rm = TRUE), fit.scale["max"], na.rm = TRUE)
        valToTest <- try(unique(round(c(seq(mini, maxi, length.out = Nb.thresh.test), mini, maxi))))
        if (inherits(valToTest, "try-error")) {
          valToTest <- seq(fit.scale["min"], fit.scale["max"], length.out = Nb.thresh.test)
        }
        # deal with unique value to test case
        if (length(valToTest) < 3) {
          valToTest <- round(c(mean(fit.scale["min"], mini), valToTest, mean(fit.scale["max"], maxi)))
        }
      }
    } else { # test only one value
      valToTest <- Fixed.thresh
    }
    
    ## 2. apply the bm_CalculateStat function ---------------------------------
    calcStat <- sapply(lapply(valToTest, function(x) {
      return(table(Fit > x, Obs))
    }), bm_CalculateStat, stat = Stat)
    
    ## 3. scale obtained scores and find best value ---------------------------
    calcStat <- 1 - abs(StatOptimum - calcStat)
    best.stat <- max(calcStat, na.rm = TRUE)
    cutoff <- median(valToTest[which(calcStat == best.stat)]) # if several values are selected
    
    misc <- table(Fit >= cutoff, Obs)
    misc <- .contingency_table_check(misc)
    true.pos <- misc['TRUE', '1']
    true.neg <- misc['FALSE', '0']
    specificity <- (true.neg * 100) / sum(misc[, '0'])
    sensitivity <- (true.pos * 100) / sum(misc[, '1'])
    
  } else { ## specific procedure for ROC value ----------------------------------------------------
    roc1 <- roc(Obs, Fit, percent = TRUE, direction = "<", levels = c(0, 1))
    roc1.out <- coords(roc1, "best", ret = c("threshold", "sens", "spec"), transpose = TRUE)
    ## if two optimal values are returned keep only the first one
    if (!is.null(ncol(roc1.out))) { roc1.out <- roc1.out[, 1] }
    best.stat <- as.numeric(auc(roc1)) / 100
    cutoff <- as.numeric(roc1.out["threshold"])
    sensitivity <- as.numeric(roc1.out["sensitivity"])
    specificity <- as.numeric(roc1.out["specificity"])
  }

  eval.out <- cbind(best.stat, cutoff, sensitivity, specificity)
  rownames(eval.out) <- Stat
  
  return(eval.out)
}

###################################################################################################

##'
##' @rdname bm_FindOptimStat
##' @export
##'

get_optim_value <- function(stat)
{
  switch(stat
         , 'TSS' = 1
         , 'KAPPA' = 1
         , 'ACCURACY' = 1
         , 'BIAS' = 1
         , 'POD' = 1
         , 'FAR' = 0
         , 'POFD' = 0
         , 'SR' = 1
         , 'CSI' = 1
         , 'ETS' = 1
         , 'HK' = 1
         , 'HSS' = 1
         , 'OR' = 1000000
         , 'ORSS' = 1
         , 'ROC' = 1
  )
}

###################################################################################################

.guess_scale <- function(Fit)
{
  min <- 0
  max <- ifelse(max(Fit, na.rm = TRUE) <= 1, 1, 1000)
  out <- c(min, max)
  names(out) <- c("min", "max")
  return(out)
}

.contingency_table_check <- function(Misc)
{
  # Contingency table checking
  if (dim(Misc)[1] == 1) {
    if (row.names(Misc)[1] == "FALSE") {
      Misc <- rbind(Misc, c(0, 0))
      rownames(Misc) <- c('FALSE', 'TRUE')
    } else {
      a <- Misc
      Misc <- c(0, 0)
      Misc <- rbind(Misc, a)
      rownames(Misc) <- c('FALSE', 'TRUE')
    }
  }
  
  if (ncol(Misc) != 2 | nrow(Misc) != 2) {
    Misc = matrix(0, ncol = 2, nrow = 2, dimnames = list(c('FALSE', 'TRUE'), c('0', '1')))
  }
  
  if ((sum(colnames(Misc) %in% c('FALSE', 'TRUE', '0', '1')) < 2) | 
      (sum(rownames(Misc) %in% c('FALSE', 'TRUE', '0', '1')) < 2) ) {
    stop("Unavailable contingency table given")
  }
  
  if ('0' %in% rownames(Misc)) { rownames(Misc)[which(rownames(Misc) == '0')] <- 'FALSE' }
  if ('1' %in% rownames(Misc)) { rownames(Misc)[which(rownames(Misc) == '1')] <- 'TRUE' }
  
  return(Misc)
}

##'
##' @rdname bm_FindOptimStat
##' @export
##'

bm_CalculateStat <- function(Misc, stat = 'TSS')
{
  ## check contagency table
  Misc <- .contingency_table_check(Misc)
  
  ## calculate basic classification information -------------------------------
  hits <- Misc['TRUE', '1'] ## true positives
  misses <- Misc['FALSE', '1'] ##  false positives
  false_alarms <- Misc['TRUE', '0'] ## false negatives 
  correct_negatives <- Misc['FALSE', '0'] ## true negatives
  
  total <- sum(Misc)
  forecast_1 <- sum(Misc['TRUE', ])
  forecast_0 <- sum(Misc['FALSE', ])
  observed_1 <- sum(Misc[, '1'])
  observed_0 <- sum(Misc[, '0'])
  
  ## calculate chosen evaluation metric ---------------------------------------
  out = switch(stat
               , 'TSS' = (hits / (hits + misses)) + (correct_negatives / (false_alarms + correct_negatives)) - 1
               , 'KAPPA' = {
                 Po <- (1 / total) * (hits + correct_negatives)
                 Pe <- ((1 / total) ^ 2) * ((forecast_1 * observed_1) + (forecast_0 * observed_0))
                 return((Po - Pe) / (1 - Pe))
               }
               , 'ACCURACY' = (hits + correct_negatives) / total
               , 'BIAS' = (hits + false_alarms) / (hits + misses)
               , 'POD' =  hits / (hits + misses)
               , 'FAR' = false_alarms / (hits + false_alarms)
               , 'POFD' = false_alarms / (correct_negatives + false_alarms)
               , 'SR' = hits / (hits + false_alarms)
               , 'CSI' = hits / (hits + misses + false_alarms)
               , 'ETS' = {
                 hits_rand <- ((hits + misses) * (hits + false_alarms)) / total
                 return((hits - hits_rand) / (hits + misses + false_alarms - hits_rand))
               }
               , 'HK' = (hits / (hits + misses)) - (false_alarms / (false_alarms + correct_negatives))
               , 'HSS' = {
                 expected_correct_rand <- (1 / total) * 
                   (((hits + misses) * (hits + false_alarms)) + ((correct_negatives + misses) * (correct_negatives + false_alarms)))
                 return((hits + correct_negatives - expected_correct_rand) / (total - expected_correct_rand))
               }
               , 'OR' = (hits * correct_negatives) / (misses * false_alarms)
               , 'ORSS' = (hits * correct_negatives - misses * false_alarms) / (hits * correct_negatives + misses * false_alarms)
  )
}


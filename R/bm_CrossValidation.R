###################################################################################################
##' @name BIOMOD_CrossValidation
##' @author Frank Breiner
##' 
##' @title Custom models cross-validation procedure
##' 
##' @description This function creates a \code{matrix} or \code{data.frame} that can be given to 
##' \code{data.split.table} parameter of \code{\link{BIOMOD_Modeling}} function to evaluate 
##' models with repeated k-fold or stratified cross-validation (CV) instead of repeated split samples.
##' 
##' 
##' @param bm.format a \code{\link{BIOMOD.formated.data-class}} or \code{\link{BIOMOD.formated.data.PA-class}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param k an \code{integer} corresponding to the number of bins/partitions for k-fold CV
##' @param nb.rep an \code{integer} corresponding to the number of repetitions of k-fold CV 
##' (\emph{set to \code{1} if \code{do.stratification = TRUE}})
##' @param do.stratification a \code{logical} defining whether stratified CV should be run 
##' @param method a \code{character} corresponding to the CV stratification method (\emph{if 
##' \code{do.stratification = TRUE}}), must be \code{x}, \code{y}, \code{both}, \code{block} 
##' or the name of a predictor for environmental stratified CV
##' @param balance a \code{character} defining whether partitions should be balanced for 
##' \code{presences} or \code{absences} (resp. pseudo-absences or background)
##' @param do.full.models (\emph{optional, default} \code{TRUE}) \cr  
##' A \code{logical} value defining whether models should be also calibrated and validated over 
##' the whole dataset or not
##' 
##' 
##' @return 
##' 
##' A \code{matrix} or \code{data.frame} with \code{k * nb.rep} (\emph{+ 1 if 
##' \code{do.full.models = TRUE}}) columns that can be given to \code{data.split.table} 
##' parameter of \code{\link{BIOMOD_Modeling}} function.
##' 
##' 
##' @details
##' 
##' \bold{Stratified cross-validation} may be used to test for model overfitting and to assess 
##' transferability in geographic and environmental space : 
##' \itemize{
##'   \item \code{x} and \code{y} stratification was described in \emph{Wenger and Olden 2012} 
##'   (see  \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_CrossValidation.html#References}{References}). While \code{y} 
##'   stratification uses \code{k} partitions along the y-gradient, \code{x} stratification does 
##'   the same for the x-gradient, and \code{both} combines them.
##'   \item \code{block} stratification was described in \emph{Muscarella et al. 2014} (see 
##'   \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_CrossValidation.html#References}{References}). Four bins of equal size are 
##'   partitioned (bottom-left, bottom-right, top-left and top-right).
##' }
##' 
##' If \code{balance = 'presences'}, presences are divided (balanced) equally over the 
##' partitions (e.g. \emph{Fig. 1b in Muscarelly et al. 2014}). Pseudo-absences will however be 
##' unbalanced over the partitions especially if the presences are clumped on an edge of the 
##' study area.
##' 
##' If \code{balance = 'absences'}, absences (resp. pseudo-absences or background) are divided 
##' (balanced) as equally as possible between the partitions (geographical balanced bins given 
##' that absences are spread over the study area equally, approach similar to \emph{Fig. 1 in 
##' Wenger et Olden 2012}). Presences will however be unbalanced over the partitions especially 
##' if the presences are clumped on an edge of the study area.
##' 
##' 
##' @references
##' 
##' \itemize{
##'   \item Muscarella, R., Galante, P.J., Soley-Guardia, M., Boria, R.A., Kass, J.M., Uriarte, M. 
##'   & Anderson, R.P. (2014). ENMeval: An R package for conducting spatially independent 
##'   evaluations and estimating optimal model complexity for Maxent ecological niche models. 
##'   \emph{Methods in Ecology and Evolution}, \bold{5}, 1198-1205.
##'   \item Wenger, S.J. & Olden, J.D. (2012). Assessing transferability of ecological models: an 
##'   underappreciated aspect of statistical validation. \emph{Methods in Ecology and Evolution}, 
##'   \bold{3}, 260-267.
##' }
##' 
##' 
##' @seealso \code{\link[ENMeval]{get.block}}, \code{\link[dismo]{kfold}}, 
##' \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_Modeling}}
##' @family Main functions
##'
##'
##' @examples
##' 
##' library(terra)
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
##' # ---------------------------------------------------------------
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # Create default modeling options
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##'  
##' # ---------------------------------------------------------------
##' # Create the different validation datasets
##' myBiomodCV <- BIOMOD_CrossValidation(bm.format = myBiomodData)
##' head(myBiomodCV)
##' 
##' # Several validation strategies can be combined
##' DataSplitTable.b <- BIOMOD_CrossValidation(bm.format = myBiomodData,
##'                                            k = 5,
##'                                            nb.rep = 2,
##'                                            do.full.models = FALSE)
##' DataSplitTable.y <- BIOMOD_CrossValidation(bm.format = myBiomodData,
##'                                            k = 2,
##'                                            do.stratification = TRUE,
##'                                            method = "y")
##' colnames(DataSplitTable.y)[1:2] <- c("RUN11", "RUN12")
##' myBiomodCV <- cbind(DataSplitTable.b, DataSplitTable.y)
##' head(myBiomodCV)
##' 
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'mod.CV',
##'                                     models = c('RF'),
##'                                     bm.options = myBiomodOptions,
##'                                     nb.rep = 2,
##'                                     data.split.table = myBiomodCV,
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 0,
##'                                     do.full.models = FALSE,
##'                                     seed.val = 42)
##' 
##' # Get evaluation scores & variables importance
##' myEval <- get_evaluations(myBiomodModelOut)
##' myEval$CV.strategy <- "Random"
##' myEval$CV.strategy[grepl("13", myEval$full.name)] <- "Full"
##' myEval$CV.strategy[grepl("11|12", myEval$full.name)] <- "Stratified"
##' head(myEval)
##' 
##' boxplot(myEval$calibration ~ interaction(myEval$algo, myEval$CV.strategy),
##'         xlab = "", ylab = "ROC AUC", col = rep(c("brown", "cadetblue"), 3))
##' boxplot(myEval$validation ~ interaction(myEval$algo, myEval$CV.strategy),
##'         xlab = "", ylab = "ROC AUC", col = rep(c("brown", "cadetblue"), 3))
##'          
##' 
## @importFrom ENMeval get.block
## @importFrom dismo kfold
##' 
##' @export
##' 
##' 
###################################################################################################


bm_CrossValidation <- function(bm.format, strategy = 'random', nb.rep = 1, perc, k, balance, strat
                               , do.full.models = TRUE, user.table = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_CrossValidation.check.args(bm.format, strategy, nb.rep, perc, k, balance, strat
                                         , do.full.models, seed.val, user.table)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  if ((nb.rep == 0 || k == 0) & strategy != 'user.defined') {
    out <- NULL
  } else {
    out <- switch(strategy,
                  user.defined = bm_CrossValidation_user.defined(user.table),
                  random = bm_CrossValidation_random(bm.format, nb.rep, perc),
                  kfold = bm_CrossValidation_kfold(bm.format, nb.rep, k),
                  block = bm_CrossValidation_block(bm.format),
                  strat = bm_CrossValidation_strat(bm.format, balance, strat, k),
                  env = bm_CrossValidation_env(bm.format, balance))
  }
  
  ## 2. CLEAN FINAL TABLE ----------------------------------------------------------------------------
  # colnames(DataSplitTable) <- paste0("RUN", 1:ncol(DataSplitTable))
  # 
  # if (isTRUE(do.full.models)) {
  #   DataSplitTable <- cbind(DataSplitTable, TRUE)
  #   colnames(DataSplitTable)[ncol(DataSplitTable)] <- "allRun"
  # }
  
  cat("\n")
  return(out)
}


# Argument Check ----------------------------------------------------------------------------------

.bm_CrossValidation.check.args <- function(bm.format, strategy, nb.rep, perc, k, balance, strat
                                           , do.full.models, seed.val, user.table)
{
  cat('\n\nChecking Cross-Validation arguments...\n')
  
  ## 0. Check bm.format argument ------------------------------------
  .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  
  ## 1. Check strategy argument -------------------------------------
  .fun_testIfIn(TRUE, "strategy", strategy, c("random", "kfold", "block", "strat", "env", "user.defined"))
  
  ## 2.a Check nb.rep / perc argument -------------------------------
  if (strategy %in% c("random", "kfold")) {
    .fun_testIfPosInt(TRUE, "nb.rep", nb.rep)
    if (nb.rep < 1) { stop("nb.rep must be an integer >= 1") }
    
    if (strategy == "random") {
      if (perc < 0 || perc > 100) {
        stop("perc must be a value between 0 and 100")
      } else if (perc < 50) {
        warning("You chose to allocate more data to evaluation than to calibration of your model
                (perc<50)\nMake sure you really wanted to do that. \n", immediate. = TRUE)
      } else if (perc == 100) {
        nb.rep <- 0
        warning(paste0("The models will be evaluated on the calibration data only "
                       , "(nb.rep=0 and no independent data) \n\t "
                       , "It could lead to over-optimistic predictive performances.\n")
                , immediate. = TRUE)
      }
    }
  }
  
  ## 2.b Check k argument -------------------------------------------
  if (strategy %in% c("kfold", "strat")) {
    .fun_testIfPosInt(TRUE, "k", k)
    if (k < 2) { stop("k must be an integer >= 2") }
  }
  
  ## 3. Check balance / strat argument ------------------------------
  if (strategy %in% c("strat", "env")) {
    .fun_testIfIn(TRUE, "balance", balance, c("presences","absences"))
    ind.NA  <- which(is.na(bm.format@data.species))
    tmp  <- bm.format@data.species
    tmp[ind.NA] <- 0 # was 2 before
    if (balance == "absences") {
      balance <- (tmp == 1 | tmp == 0)
    } else {
      balance <- (tmp == 1)
    }
    
    if (strategy == "strat") {
      .fun_testIfIn(TRUE, "strat", strat, c("x", "y", "both"))
    }
  }
  
  ## 4. Check user.table argument -----------------------------------
  if (strategy == "user.defined") {
    if (is.null(user.table)) {
      stop("user.table must be a matrix or a data.frame") 
    } else {
      .fun_testIfInherits(TRUE, "user.table", user.table, c("matrix", "data.frame"))
      if (inherits(user.table, 'data.frame')) {
        user.table <- as.matrix(user.table)
      }
      if (dim(user.table)[1] != length(bm.format@data.species)) { 
        stop("user.table must have as many rows (dim1) than your species as data")
      }
      nb.rep <- dim(user.table)[2]
    }
  }
  
  return(list(bm.format = bm.format,
              strategy = strategy,
              nb.rep = nb.rep,
              perc = perc,
              k = k,
              balance = balance,
              strat = strat,
              do.full.models = do.full.models,
              seed.val = seed.val,
              user.table = user.table))
}

# ---------------------------------------------------------------------------- #

.sample_mat <- function(data.sp, data.split, nb.rep = 1, data.env = NULL, seed.val = NULL)
{
  # data.sp is a 0, 1 vector
  # return a matrix with nb.rep columns of boolean (T: calib, F= eval)
  
  pres <- which(data.sp == 1)
  abs <- (1:length(data.sp))[-pres]
  
  nbPresEval <- round(length(pres) * data.split / 100)
  nbAbsEval <- round(length(abs) * data.split / 100)
  
  mat.out <- matrix(FALSE, nrow = length(data.sp), ncol = nb.rep)
  colnames(mat.out) <- paste0('_RUN', 1:nb.rep)
  
  set.seed(seed.val)
  for (i in 1:ncol(mat.out)) {
    ## force to sample at least one level of each factorial variable for calibration
    fact.cell.samp <- NULL
    if (!is.null(data.env)) {
      fact.cell.samp <- bm_SampleFactorLevels(expl.var = data.env)
      mat.out[fact.cell.samp, i] <- TRUE ## in fact.cell.samp
    }
    mat.out[sample(setdiff(pres, fact.cell.samp), ## in pres, not in fact.cell.samp
                   max(nbPresEval - length(fact.cell.samp), 0)), i] <- TRUE
    mat.out[sample(setdiff(abs, fact.cell.samp), ## in abs, not in fact.cell.samp
                   max(nbAbsEval - length(fact.cell.samp), 0)), i] <- TRUE
  }
  return(mat.out)
}



# bm_CrossValidation user-defined methods ---------------------------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_user.defined",
           def = function(bm.format, ...) {
             standardGeneric( "bm_CrossValidation_user.defined")
           })

## bm_CrossValidation user-defined BIOMOD.formated.data methods -------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_user.defined', signature(bm.format = "BIOMOD.formated.data"),
          function(bm.format, user.table) {
            cat("\n   > User defined cross-validation selection")
            calib.lines <- user.table
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })

## bm_CrossValidation user-defined BIOMOD.formated.data.PA methods ----------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_user.defined', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format, user.table) {
            cat("\n   > User defined cross-validation selection")
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                calib.pa <- user.table
                calib.pa[which(bm.format@PA.table[, pa] == FALSE | is.na(bm.format@PA.table[, pa])), ] <- NA
                return(calib.pa)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })


# bm_CrossValidation random methods ---------------------------------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_random",
           def = function(bm.format, ...) {
             standardGeneric( "bm_CrossValidation_random")
           })

## bm_CrossValidation random BIOMOD.formated.data methods -------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_random', signature(bm.format = "BIOMOD.formated.data"),
          function(bm.format, nb.rep, perc) {
            cat("\n   > Random cross-validation selection")
            if (nb.rep == 0) { # take all available data
              calib.lines <- matrix(rep(TRUE, length(bm.format@data.species)), ncol = 1)
              colnames(calib.lines) <- '_allRun'
            } else {
              calib.lines <- .sample_mat(data.sp = bm.format@data.species,
                                         data.split = perc,
                                         nb.rep = nb.rep,
                                         data.env = bm.format@data.env.var,
                                         seed.val = seed.val)
              # if (do.full.models) {
              #   calib.lines <- cbind(calib.lines, rep(TRUE, length(bm.format@data.species)))
              #   colnames(calib.lines)[nb.rep + 1] <- '_allRun'
              # }
            }
            return(calib.lines)
          })

## bm_CrossValidation random BIOMOD.formated.data.PA methods ----------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_random', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format, nb.rep, perc) {
            cat("\n   > Random cross-validation selection")
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                
                if (nb.rep == 0) { # take all available data
                  calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = 1)
                  calib.pa[ind.PA, ] <- TRUE
                  colnames(calib.pa) <- '_allRun'
                } else {
                  sampled.mat <- .sample_mat(data.sp = bm.format@data.species[ind.pa],
                                             data.split = perc,
                                             nb.rep = nb.rep,
                                             data.env = bm.format@data.env.var[ind.PA, , drop = FALSE],
                                             seed.val = seed.val)
                  calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = nb.rep)
                  calib.pa[ind.PA, ] <- sampled.mat
                  # colnames(calib.lines) <- colnames(sampled.mat)
                  # if (do.full.models) {
                  #   calib.lines <- cbind(calib.lines, rep(NA, length(bm.format@data.species)))
                  #   calib.lines[ind.PA, nb.rep + 1] <- TRUE
                  #   colnames(calib.lines)[nb.rep + 1] <- '_allRun'
                  # }
                }
                return(calib.pa)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })


# bm_CrossValidation kfold methods ----------------------------------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_kfold",
           def = function(bm.format, ...) {
             standardGeneric( "bm_CrossValidation_kfold")
           })

## bm_CrossValidation kfold BIOMOD.formated.data methods --------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_kfold', signature(bm.format = "BIOMOD.formated.data"),
          function(bm.format, nb.rep, k) {
            cat("\n   > k-fold cross-validation selection")
            if (!isNamespaceLoaded("dismo")) { 
              if(!requireNamespace('dismo', quietly = TRUE)) stop("Package 'dismo' not found")
            }
            
            ind.NA  <- which(is.na(bm.format@data.species))
            tmp  <- bm.format@data.species
            tmp[ind.NA] <- 0 # was 2 before
            
            calib.lines <- foreach(rep = 1:nb.rep) %do%
              {
                fold <- dismo::kfold(tmp, by = tmp, k = k)
                calib.rep <- NULL
                for (i in 1:k) {
                  calib.rep <- cbind(calib.rep, fold != i)
                }
                return(calib.rep)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })

## bm_CrossValidation kfold BIOMOD.formated.data.PA methods -----------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_kfold', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format, nb.rep, k) {
            cat("\n   > k-fold cross-validation selection")
            if (!isNamespaceLoaded("dismo")) {
              if(!requireNamespace('dismo', quietly = TRUE)) stop("Package 'dismo' not found")
            }
            
            ind.NA  <- which(is.na(bm.format@data.species))
            tmp  <- bm.format@data.species
            tmp[ind.NA] <- 0 # was 2 before
            
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                
                if (nb.rep == 0) { # take all available data
                  calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = 1)
                  calib.pa[ind.PA, ] <- TRUE
                  colnames(calib.pa) <- '_allRun'
                } else {
                  tmp.pa <- tmp[ind.PA]
                  calib.pa_rep <- foreach(rep = 1:nb.rep) %do%
                    {
                      fold <- dismo::kfold(tmp.pa, by = tmp.pa, k = k)
                      calib.rep <- NULL
                      for (i in 1:k) {
                        calib.rep <- cbind(calib.rep, fold != i)
                      }
                      return(calib.rep)
                    }
                  calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = nb.rep * k)
                  calib.pa[ind.PA, ] <- calib.pa_rep
                }
                return(calib.pa)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })


# bm_CrossValidation block methods ----------------------------------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_block",
           def = function(bm.format, ...) {
             standardGeneric( "bm_CrossValidation_block")
           })

## bm_CrossValidation block BIOMOD.formated.data methods --------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_block', signature(bm.format = "BIOMOD.formated.data"),
          function(bm.format) {
            cat("\n   > Block cross-validation selection")
            if (!isNamespaceLoaded("ENMeval")) { 
              if(!requireNamespace('ENMeval', quietly = TRUE)) stop("Package 'ENMeval' not found")
            }
            
            ind.NA  <- which(is.na(bm.format@data.species))
            tmp  <- bm.format@data.species
            tmp[ind.NA] <- 0 # was 2 before
            
            blocks <- ENMeval::get.block(bm.format@coord[tmp == 1, ], bm.format@coord[tmp == 0, ])
            calib.lines <- matrix(NA, nrow = length(tmp), ncol = 4)
            for (i in 1:4) {
              calib.lines[tmp == 1, i] <- blocks[[1]] != i
              calib.lines[tmp == 0, i] <- blocks[[2]] != i     
            }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })

## bm_CrossValidation block BIOMOD.formated.data.PA methods -----------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_block', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format) {
            cat("\n   > Block cross-validation selection")
            if (!isNamespaceLoaded("ENMeval")) { 
              if(!requireNamespace('ENMeval', quietly = TRUE)) stop("Package 'ENMeval' not found")
            }
            
            ind.NA  <- which(is.na(bm.format@data.species))
            tmp  <- bm.format@data.species
            tmp[ind.NA] <- 0 # was 2 before
            
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                tmp.pa <- tmp[ind.PA]
                
                blocks <- ENMeval::get.block(bm.format@coord[tmp.pa == 1, ], bm.format@coord[tmp.pa == 0, ])
                calib.pa_rep <- matrix(NA, nrow = length(tmp.pa), ncol = 4)
                for (i in 1:4) {
                  calib.pa_rep[tmp == 1, i] <- blocks[[1]] != i
                  calib.pa_rep[tmp == 0, i] <- blocks[[2]] != i     
                }
                
                calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = 4)
                calib.pa[ind.PA, ] <- calib.pa_rep
                return(calib.pa)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })


# bm_CrossValidation strat methods ----------------------------------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_strat",
           def = function(bm.format, ...) {
             standardGeneric( "bm_CrossValidation_strat")
           })

## bm_CrossValidation strat BIOMOD.formated.data methods --------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_strat', signature(bm.format = "BIOMOD.formated.data"),
          function(bm.format, balance, strat, k) {
            cat("\n   > Stratified cross-validation selection")
            tmp.coord <- bm.format@coord
              
            if (strat == "x" || strat == "both") {
              bands <- quantile(tmp.coord[balance, 1], probs = seq(0, 100, 100 / k) / 100)
              bands[1] <- bands[1] - 1
              bands[k + 1] <- bands[k + 1] + 1
              calib.x <- matrix(NA, nrow(tmp.coord), k)
              for (i in 1:k) {
                calib.x[, i] <- (tmp.coord[, 1] >= bands[i] & tmp.coord[, 1] < bands[i + 1])
              }
              if (strat == "x") { calib.lines <- calib.x }
            }
            
            if (strat == "y" || strat == "both") {
              bands <- quantile(tmp.coord[balance, 2], probs = seq(0, 100, 100 / k) / 100)
              bands[1] <- bands[1] - 1
              bands[k + 1] <- bands[k + 1] + 1
              calib.y <- matrix(NA, nrow(tmp.coord), k)
              for (i in 1:k) {
                calib.y[, i] <- (tmp.coord[, 2] >= bands[i] & tmp.coord[, 2] < bands[i + 1])
              }
              if (strat == "y") { calib.lines <- calib.y }
            }
            
            if (strat == "both") { ## Merge X and Y tables
              calib.lines <- cbind(calib.x, calib.y)
            }
            
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })

## bm_CrossValidation block BIOMOD.formated.data.PA methods -----------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_strat', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format, balance, strat, k) {
            cat("\n   > Stratified cross-validation selection")
            tmp.coord <- bm.format@coord
            
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                tmp.pa <- tmp.coord[ind.PA, ]
                
                if (strat == "x" || strat == "both") {
                  bands <- quantile(tmp.pa[balance, 1], probs = seq(0, 100, 100 / k) / 100)
                  bands[1] <- bands[1] - 1
                  bands[k + 1] <- bands[k + 1] + 1
                  calib.x <- matrix(NA, nrow(tmp.pa), k)
                  for (i in 1:k) {
                    calib.x[, i] <- (tmp.pa[, 1] >= bands[i] & tmp.pa[, 1] < bands[i + 1])
                  }
                  if (strat == "x") { calib.pa_rep <- calib.x }
                }
                if (strat == "y" || strat == "both") {
                  bands <- quantile(tmp.pa[balance, 2], probs = seq(0, 100, 100 / k) / 100)
                  bands[1] <- bands[1] - 1
                  bands[k + 1] <- bands[k + 1] + 1
                  calib.y <- matrix(NA, nrow(tmp.pa), k)
                  for (i in 1:k) {
                    calib.y[, i] <- (tmp.pa[, 2] >= bands[i] & tmp.pa[, 2] < bands[i + 1])
                  }
                  if (strat == "y") { calib.pa_rep <- calib.y }
                }
                if (strat == "both") { ## Merge X and Y tables
                  calib.pa_rep <- cbind(calib.x, calib.y)
                }
                
                calib.pa <- matrix(NA, nrow = nrow(tmp.coord), ncol = 4)
                calib.pa[ind.PA, ] <- calib.pa_rep
                return(calib.pa)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })


# bm_CrossValidation env methods ------------------------------------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_env",
           def = function(bm.format, ...) {
             standardGeneric( "bm_CrossValidation_env")
           })

## bm_CrossValidation env BIOMOD.formated.data methods ----------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_env', signature(bm.format = "BIOMOD.formated.data"),
          function(bm.format, balance) {
            cat("\n   > Environmental cross-validation selection")
            calib.lines <- foreach(env = colnames(bm.format@data.env.var), .combine = "cbind") %do%
              {
                tmp.env <- bm.format@data.env.var[balance, env]
                
                bands <- quantile(tmp.env, probs = seq(0, 100, 100 / k) / 100)
                bands[1] <- bands[1] - 1
                bands[k + 1] <- bands[k + 1] + 1
                calib.env <- matrix(NA, nrow = length(tmp.env), ncol = k)
                for (i in 1:k) {
                  calib.env[, i] <- (tmp.env <= bands[i] | tmp.env > bands[i + 1])
                }
                return(calib.env)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })

## bm_CrossValidation env BIOMOD.formated.data.PA methods -------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_env', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format, balance) {
            cat("\n   > Environmental cross-validation selection")
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                calib.env <- foreach(env = colnames(bm.format@data.env.var), .combine = "cbind") %do%
                  {
                    tmp.env <- bm.format@data.env.var[balance, env]
                    tmp.pa <- tmp.env[ind.PA]
                    
                    bands <- quantile(tmp.pa, probs = seq(0, 100, 100 / k) / 100)
                    bands[1] <- bands[1] - 1
                    bands[k + 1] <- bands[k + 1] + 1
                    calib.pa_env <- matrix(NA, nrow = length(tmp.pa), ncol = k)
                    for (i in 1:k) {
                      calib.pa_env[, i] <- (tmp.pa <= bands[i] | tmp.pa > bands[i + 1])
                    }
                    
                    calib.pa <- matrix(NA, nrow = nrow(bm.format@data.env.var), ncol = k)
                    calib.pa[ind.PA, ] <- calib.pa_env
                    return(calib.pa)
                  }
                return(calib.env)
              }
            colnames(calib.lines) <- paste('_RUN', 1:ncol(calib.lines), sep = '')
            return(calib.lines)
          })

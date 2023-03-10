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
                               , do.full.models = TRUE, seed.val = NULL, user.table = NULL)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_CrossValidation.check.args(bm.format, strategy, nb.rep, perc, k, balance, strat
                                         , do.full.models, seed.val, user.table)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  if ((nb.rep == 0 || k == 0) & strategy != 'user.defined') { ##TODO
    out <- NULL
  } else {
    out <- switch(strategy, ##TODO
                  user.defined = bm_CrossValidation_user.defined(resp.var, expl.var, user.table),
                  random = bm_CrossValidation_random(resp.var, expl.var, nb.absences, nb.rep),
                  kfold = bm_CrossValidation_kfold(),
                  block = bm_CrossValidation_block(),
                  strat = bm_CrossValidation_strat(),
                  env = bm_CrossValidation_env())
  }
  
  cat("\n")
  return(out)
}

  
# Argument Check --------------------------------------------------------------

.bm_CrossValidation.check.args <- function(bm.format, strategy, nb.rep, perc, k, balance, strat
                                           , do.full.models, seed.val, user.table)
{
  cat('\n\nChecking Cross-Validation arguments...\n')
  
  ## 0. Check bm.format argument ------------------------------------
  .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  
  ## 1. Check strategy argument -------------------------------------
  .fun_testIfIn(TRUE, "strategy", strategy, c("random", "kfold", "block", "strat", "env", "user.defined"))
  
  ## 2. Check nb.rep / perc / k argument ----------------------------
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
    
    if (strategy == "kfold") {
      .fun_testIfPosInt(TRUE, "k", k)
      if (k < 2) { stop("k must be an integer >= 2") }
    }
  }
  
  ## 3. Check balance / strat argument ------------------------------
  if (strategy %in% c("block", "strat", "env")) {
    .fun_testIfIn(TRUE, "balance", balance, c("presences","absences"))
    
    if (strategy == "strat") {
      .fun_testIfIn(TRUE, "strat", strat, c("x", "y", "both"))
    }
  }
  
  ## 4. Check user.table argument -----------------------------------
  if (strategy == "user.defined") {
    if (is.null(user.table)) {
      stop("user.table must be a matrix, a data.frame or a 3D array") 
    } else {
      if (!(length(dim(user.table) %in% c(2, 3)))) {
        stop("user.table must be a matrix, a data.frame or a 3D array") 
      }
      if (inherits(user.table,'data.frame')) {
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



# bm_CrossValidation user-defined methods --------------------------------------

##'
##' @rdname bm_CrossValidation
##' @export
##'

setGeneric("bm_CrossValidation_user.defined",
           def = function(user.table, ...) {
             standardGeneric( "bm_CrossValidation_user.defined")
           })

## bm_CrossValidation user-defined data.frame methods --------------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_user.defined', signature(user.table = "data.frame"),
          function(user.table) {
            cat("\n   > User defined cross-validation selection")
            # return(list(xy = crds(resp.var),
            #             sp = as.numeric(unlist(values(resp.var), use.names = FALSE)),
            #             env = as.data.frame(expl.var),
            #             pa.tab = user.table))
          })





setGeneric("bm_CrossValidation", def = function(bm.format, ...) { standardGeneric("bm_CrossValidation") })

setMethod("bm_CrossValidation", signature("BIOMOD.formated.data"),
          function(bm.format, nb.rep, data.split.perc
                   , weights = NULL, prevalence = NULL
                   , do.full.models = TRUE, data.split.table = NULL, seed.val = NULL)
          {
            
          }
)

setMethod("bm_CrossValidation", signature("BIOMOD.formated.data.PA"),
          function(bm.format, nb.rep, data.split.perc, weights = NULL, prevalence = NULL
                   , do.full.models = TRUE, data.split.table = NULL, seed.val = NULL)
          {
            
          }
)

## method / strategy ? : kfold / random / stratified / block / environment
## nb.rep
## k
## perc
## balance
## strat






# BIOMOD_CrossValidation <- function(bm.format,
                                   k = 5,
                                   # nb.rep = 5,
                                   do.stratification = FALSE,
                                   method = "both",
                                   balance = "presences",
                                   do.full.models = TRUE) {
  
  args <- .BIOMOD_CrossValidation.check.args(bm.format,
                                             k = k,
                                             nb.rep = nb.rep,
                                             do.stratification = do.stratification,
                                             method = method,
                                             balance = balance,
                                             do.full.models = do.full.models)
  
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  .bm_cat("Build Cross-Validation Table")
  DataSplitTable.y <- DataSplitTable.x <- DataSplitTable <- NULL
  ind.NA  <- which(is.na(bm.format@data.species))
  tmp  <- bm.format@data.species
  tmp[ind.NA] <- 0 # was 2 before
  
  ## STRATIFIED (X, Y, BOTH) / BLOCK / ENVIRONMENTAL CROSS VALIDATION -----------------------------
  if (do.stratification) {
    if (balance == "absences") {
      balance <- (tmp == 1 | tmp == 0)
    } else {
      balance <- (tmp == 1)
    }
    
    ## (X, Y, BOTH) STRATIFIED CROSS VALIDATION  -------------------------------
    if (method == "x" | method == "both") {
      DataSplitTable.x <- matrix(NA, nrow(bm.format@coord), k)
      bands <- quantile(bm.format@coord[balance, 1], probs = seq(0, 100, 100 / k) / 100)
      bands[1] <- bands[1] - 1
      bands[k + 1] <- bands[k + 1] + 1
      for (i in 1:k) {
        DataSplitTable.x[, i] <- (bm.format@coord[, 1] >= bands[i] & bm.format@coord[, 1] < bands[i + 1])
      }
      if (method == "x") { DataSplitTable <- DataSplitTable.x }
    }
    if (method == "y" | method == "both") {
      DataSplitTable.y <- matrix(NA, nrow(bm.format@coord), k)
      bands <- quantile(bm.format@coord[balance, 2], probs = seq(0, 100, 100 / k) / 100)
      bands[1] <- bands[1] - 1
      bands[k + 1] <- bands[k + 1] + 1
      for (i in 1:k) {
        DataSplitTable.y[, i] <- (bm.format@coord[, 2] >= bands[i] & bm.format@coord[, 2] < bands[i + 1])
      }
      if (method == "y") { DataSplitTable <- DataSplitTable.y }
    }
    if (method == "both") { ## Merge X and Y tables
      DataSplitTable <- cbind(DataSplitTable.x, DataSplitTable.y)
    }
    
    ## BLOCK STRATIFIED CROSS VALIDATION --------------------------------------
    if (method == "block") {
      if (!isNamespaceLoaded("ENMeval")) { 
        if(!requireNamespace('ENMeval', quietly = TRUE)) stop("Package 'ENMeval' not found")
      }
      DataSplitTable <- as.data.frame(matrix(NA, nrow(bm.format@coord), 4))
      blocks <- ENMeval::get.block(bm.format@coord[tmp == 1, ]
                                   , bm.format@coord[tmp == 0, ])
      for (i in 1:4) {
        DataSplitTable[tmp == 1, i] <- blocks[[1]] != i
        DataSplitTable[tmp == 0, i] <- blocks[[2]] != i     
      }
      DataSplitTable <- as.matrix(DataSplitTable)
    }
    
    ## ENVIRONMENTAL STRATIFIED CROSS VALIDATION ------------------------------
    if (method != "block" & method != "x" & method != "y" & method != "both") {
      DataSplitTable2 <- as.data.frame(matrix(NA, nrow(bm.format@coord), k))
      bands <- quantile(bm.format@data.env.var[balance, method], probs = seq(0, 100, 100 / k) / 100)
      bands[1] <- bands[1] - 1
      bands[k + 1] <- bands[k + 1] + 1
      for (i in 1:k) {
        DataSplitTable2[, i] <- (bm.format@data.env.var[balance, method] <= bands[i] | 
                                   bm.format@data.env.var[balance, method] > bands[i + 1])
      }
    }
  } else {
    ## K-FOLD CROSS VALIDATION --------------------------------------------------------------------
    if (!isNamespaceLoaded("dismo")) { 
      if(!requireNamespace('dismo', quietly = TRUE)) stop("Package 'dismo' not found")
    }
    for (rep in 1:nb.rep) {
      fold <- dismo::kfold(tmp, by = tmp, k = k)
      for (i in 1:k) {
        DataSplitTable <- cbind(DataSplitTable, fold != i)
      }
    }
  }
  
  ## CLEAN FINAL TABLE ----------------------------------------------------------------------------
  colnames(DataSplitTable) <- paste0("RUN", 1:ncol(DataSplitTable))
  
  if (isTRUE(do.full.models)) {
    DataSplitTable <- cbind(DataSplitTable, TRUE)
    colnames(DataSplitTable)[ncol(DataSplitTable)] <- "allRun"
  }
  
  .bm_cat("Done")
  return(DataSplitTable)
}



# Argument check ----------------------------------------------------------

.BIOMOD_CrossValidation.check.args <- 
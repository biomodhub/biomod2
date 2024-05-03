###################################################################################################
##' @name bm_CrossValidation
##' @author Frank Breiner, Maya Gueguen
##' 
##' @title Build cross-validation table
##' 
##' @description This internal \pkg{biomod2} function allows to build a cross-validation table 
##' according to 6 different methods : \code{random}, \code{kfold}, \code{block}, \code{strat}, 
##' \code{env} or \code{user.defined} (see Details).
##' 
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param strategy a \code{character} corresponding to the cross-validation selection strategy, 
##' must be among \code{random}, \code{kfold}, \code{block}, \code{strat}, \code{env} or 
##' \code{user.defined}
##' @param \ldots (\emph{optional, one or several of the following arguments depending on the 
##' selected method}) 
##' 
##' @param nb.rep (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'random'} or \code{strategy = 'kfold'}, an \code{integer} corresponding 
##' to the number of sets (repetitions) of cross-validation points that will be drawn
##' @param perc (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'random'}, a \code{numeric} between \code{0} and \code{1} defining the 
##' percentage of data that will be kept for calibration
##' @param k (\emph{optional, default} \code{0}) \cr
##' If \code{strategy = 'kfold'} or \code{strategy = 'strat'} or \code{strategy = 'env'}, an 
##' \code{integer} corresponding to the number of partitions 
##' @param balance (\emph{optional, default} \code{'presences'}) \cr 
##' If \code{strategy = 'strat'} or \code{strategy = 'env'}, a \code{character} corresponding 
##' to how data will be balanced between partitions, must be either \code{presences} or 
##' \code{absence}
##' @param env.var (\emph{optional}) \cr 
##' If \code{strategy = 'env'}, a \code{character} corresponding to the environmental variables 
##' used to build the partition. \code{k} partitions will be built for each environmental 
##' variables. \emph{By default the function uses all environmental variables available.}
##' @param strat (\emph{optional, default} \code{'both'}) \cr 
##' If \code{strategy = 'env'}, a \code{character} corresponding to how data will partitioned 
##' along gradient, must be among \code{x}, \code{y}, \code{both}
##' @param user.table (\emph{optional, default} \code{NULL}) \cr
##' If \code{strategy = 'user.defined'}, a \code{matrix} or \code{data.frame} defining for each 
##' repetition (in columns) which observation lines should be used for models calibration 
##' (\code{TRUE}) and validation (\code{FALSE})
##' @param do.full.models (\emph{optional, default} \code{TRUE}) \cr  
##' A \code{logical} value defining whether models should be also calibrated and validated over 
##' the whole dataset (and pseudo-absence datasets) or not
##' 
##' 
##' @return 
##' 
##' A \code{matrix} or \code{data.frame} defining for each repetition (in columns) which 
##' observation lines should be used for models calibration (\code{TRUE}) and validation 
##' (\code{FALSE}).
##' 
##' 
##' @details
##' 
##' Several parameters are available within the function and some of them can be used with 
##' different cross-validation strategies :
##' 
##' \code{| ....... | random | kfold | block | strat | env |} \cr
##' __________________________________________________ \cr
##' \code{| nb.rep. | x..... | x.... | ..... | ..... | ... |} \cr
##' \code{| perc... | x..... | ..... | ..... | ..... | ... |} \cr
##' \code{| k...... | ...... | x.... | ..... | x.... | x.. |} \cr
##' \code{| balance | ...... | ..... | ..... | x.... | x.. |} \cr
##' \code{| strat.. | ...... | ..... | ..... | x.... | ... |} \cr \cr \cr
##' 
##' 
##' \bold{Concerning column names of \code{matrix} output :}
##' 
##' The number of columns depends on the strategy selected. 
##' The column names are given \emph{a posteriori} of the selection, ranging from 1 to the 
##' number of columns. 
##' If \code{do.full.models = TRUE}, columns merging runs (and/or pseudo-absence datasets) 
##' are added at the end. \cr \cr
##' 
##' 
##' 
##' \bold{Concerning cross-validation strategies :}
##' 
##' \describe{
##'   \item{random}{Most simple method to calibrate and validate a model is to split the original 
##'   dataset in two datasets : one to calibrate the model and the other one to validate it. The 
##'   splitting can be repeated \code{nb.rep} times.}
##'   \item{k-fold}{The k-fold method splits the original dataset in \code{k} datasets of equal 
##'   sizes : each part is used successively as the validation dataset while the other \code{k-1} 
##'   parts are used for the calibration, leading to \code{k} calibration/validation ensembles. 
##'   This multiple splitting can be repeated \code{nb.rep} times.}
##'   \item{block}{It may be used to test for model overfitting and to assess transferability in 
##'   geographic space. \code{block} stratification was described in \emph{Muscarella et al. 2014} 
##'   (see References). Four bins of equal size are partitioned (bottom-left, bottom-right, 
##'   top-left and top-right).}
##'   \item{stratified}{It may be used to test for model overfitting and to assess transferability 
##'   in geographic space. \code{x} and \code{y} stratification was described in \emph{Wenger and 
##'   Olden 2012} (see References). \code{y} stratification uses \code{k} partitions along the 
##'   y-gradient, \code{x} stratification does the same for the x-gradient. \code{both} returns 
##'   \code{2k} partitions: \code{k} partitions stratified along the x-gradient and \code{k} 
##'   partitions stratified along the y-gradient.}
##'   \item{environmental}{It may be used to test for model overfitting and to assess 
##'   transferability in environmental space. It returns \code{k} partitions for each variable 
##'   given in \code{env.var}.}
##'   \item{user-defined}{Allow the user to give its own crossvalidation table. For a 
##'   presence-absence dataset, column names must be formatted as: \code{_allData_RUNx} with 
##'   \code{x} an integer. For a presence-only dataset for which several pseudo-absence dataset 
##'   were generated, column names must be formatted as: \code{_PAx_RUNy} with \code{x} an 
##'   integer and \code{PAx} an existing pseudo-absence dataset and \code{y} an integer \cr \cr}
##' }
##' 
##' 
##' \bold{Concerning balance parameter :}
##' 
##' If \code{balance = 'presences'}, presences are divided (balanced) equally over the partitions 
##' (e.g. \emph{Fig. 1b in Muscarelly et al. 2014}). 
##' Absences or pseudo-absences will however be unbalanced over the partitions especially if the 
##' presences are clumped on an edge of the study area.
##' 
##' If \code{balance = 'absences'}, absences (resp. pseudo-absences or background) are divided 
##' (balanced) as equally as possible between the partitions (geographical balanced bins given 
##' that absences are spread over the study area equally, approach similar to \emph{Fig. 1 in 
##' Wenger et Olden 2012}). Presences will however be unbalanced over the partitions especially
##' if the presences are clumped on an edge of the study area.
##' 
##' 
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
##' @family Secundary functions
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
##' # --------------------------------------------------------------- #
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # --------------------------------------------------------------- #
##' # Create the different validation datasets
##' 
##' # random selection
##' cv.r <- bm_CrossValidation(bm.format = myBiomodData,
##'                            strategy = "random",
##'                            nb.rep = 3,
##'                            k = 0.8)
##' 
##' # k-fold selection
##' cv.k <- bm_CrossValidation(bm.format = myBiomodData,
##'                            strategy = "kfold",
##'                            nb.rep = 2,
##'                            k = 3)
##' 
##' # block selection
##' cv.b <- bm_CrossValidation(bm.format = myBiomodData,
##'                            strategy = "block")
##' 
##' # stratified selection (geographic)
##' cv.s <- bm_CrossValidation(bm.format = myBiomodData,
##'                            strategy = "strat",
##'                            k = 2,
##'                            balance = "presences",
##'                            strat = "x")
##' 
##' # stratified selection (environmental)
##' cv.e <- bm_CrossValidation(bm.format = myBiomodData,
##'                            strategy = "env",
##'                            k = 2,
##'                            balance = "presences")
##' 
##' head(cv.r)
##' apply(cv.r, 2, table)
##' head(cv.k)
##' apply(cv.k, 2, table)
##' head(cv.b)
##' apply(cv.b, 2, table)
##' head(cv.s)
##' apply(cv.s, 2, table)
##' head(cv.e)
##' apply(cv.e, 2, table)
##' 
##' 
##' @importFrom foreach foreach %do%
##' 
##' @export
##' 
##' 



bm_CrossValidation <- function(bm.format,
                               strategy = 'random',
                               nb.rep = 0, 
                               perc = 0.8,
                               k = 0, 
                               balance = 'presences', 
                               env.var = NULL,
                               strat = 'both',
                               user.table = NULL,
                               do.full.models = FALSE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_CrossValidation.check.args(bm.format = bm.format,
                                         strategy = strategy,
                                         nb.rep = nb.rep,
                                         perc = perc, 
                                         k = k, 
                                         balance = balance, 
                                         env.var = env.var,
                                         strat = strat,
                                         user.table = user.table,
                                         do.full.models = do.full.models)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  if ((strategy %in% c('random', 'kfold') && nb.rep == 0) ||
      (strategy %in% c('kfold', 'strat', 'env') && k == 0)) {
    # take all available data
    out <- matrix(rep(TRUE, length(bm.format@data.species)), ncol = 1)
    colnames(out) <- '_allRun'
  } else {
    out <- switch(strategy,
                  user.defined = bm_CrossValidation_user.defined(bm.format, user.table),
                  random = bm_CrossValidation_random(bm.format, nb.rep, perc),
                  kfold = bm_CrossValidation_kfold(bm.format, nb.rep, k),
                  block = bm_CrossValidation_block(bm.format),
                  strat = bm_CrossValidation_strat(bm.format, balance, strat, k),
                  env = bm_CrossValidation_env(bm.format, balance, k, env.var))
  }
  
  ## 2. CLEAN FINAL TABLE ----------------------------------------------------------------------------
  if (!inherits(bm.format, "BIOMOD.formated.data.PA") & strategy != "user.defined") {
    colnames(out) <- paste0("_allData", colnames(out))
  }
  if (do.full.models) {
    out <- cbind(out, TRUE)
    colnames(out)[ncol(out)] <- "_allData_allRun"
    if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
      for (pa in 1:ncol(bm.format@PA.table)) {
        ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
        out <- cbind(out, NA)
        out[ind.PA, ncol(out)] <- TRUE
        colnames(out)[ncol(out)] <- paste0("_PA", pa, "_allRun")
      }
    }
  }
  
  # check for unbalanced dataset (dataset missing presences or absences)
  which.calibration.unbalanced <-
    which(
      apply(out, 2, 
            function(x) {
              length(unique(bm.format@data.species[which(x)]))
            }
      ) != 2)
  
  if(length(which.calibration.unbalanced) > 0) {
    cat("\n   !!! Some calibration dataset do not have both presences and absences: ", 
        paste0(colnames(out)[which.calibration.unbalanced], collapse = ", "))
    warning("Some calibration repetion do not have both presences and absences")
  }
  
  which.validation.unbalanced <-
    which(
      apply(out, 2, 
            function(x) {
              length(unique(bm.format@data.species[which(!x)]))
            }
      ) != 2)
  
  # Models with allRun have no validation
  which.validation.unbalanced <- 
    which.validation.unbalanced[which(!grepl(names(which.validation.unbalanced), pattern = "allRun"))]
  
  if (length(which.validation.unbalanced) > 0) {
    cat("\n   !!! Some validation dataset do not have both presences and absences: ", 
        paste0(colnames(out)[which.validation.unbalanced], collapse = ", "))
    warning("Some validation repetion do not have both presences and absences")
  }
  
  cat("\n")
  return(out)
}


# Argument Check ----------------------------------------------------------------------------------

.bm_CrossValidation.check.args <- function(bm.format, strategy, nb.rep, perc, k, balance,
                                           env.var, strat, user.table, do.full.models) {
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
      if (is.null(perc)) {
        stop("perc (or CV.perc) is required when strategy = 'random'")
      }
      .fun_testIf01(TRUE, "perc", perc)
      if (perc < 0.5) {
        warning("You chose to allocate more data to validation than to calibration of your model
                (perc<0.5)\nMake sure you really wanted to do that. \n", immediate. = TRUE)
      } else if (perc == 1) {
        nb.rep <- 0
        warning(paste0("The models will be evaluated on the calibration data only "
                       , "(nb.rep=0 and no independent data) \n\t "
                       , "It could lead to over-optimistic predictive performances.\n")
                , immediate. = TRUE)
      }
    }
  }
  
  ## 2.b Check k argument -------------------------------------------
  if (strategy %in% c("kfold", "strat", "env")) {
    .fun_testIfPosInt(TRUE, "k", k)
    if (k < 2) { stop("k must be an integer >= 2") }
  }
  
  ## 2.c Check env.var argument -------------------------------------------
  if (strategy %in% c("env")) {
    if (is.null(env.var)) {
      env.var <- colnames(bm.format@data.env.var)
    } else {
      .fun_testIfIn(TRUE, "env.var", env.var, colnames(bm.format@data.env.var))
    }
  }
  ## 3. Check balance / strat argument ------------------------------
  if (strategy %in% c("strat", "env")) {
    .fun_testIfIn(TRUE, "balance", balance, c("presences","absences"))
    ind.NA  <- which(is.na(bm.format@data.species))
    tmp  <- bm.format@data.species
    tmp[ind.NA] <- 0
    if (balance == "absences") {
      balance <- (tmp == 0) 
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
      # nb.rep <- dim(user.table)[2]
    }
  }
  
  return(list(bm.format = bm.format,
              strategy = strategy,
              nb.rep = nb.rep,
              perc = perc,
              k = k,
              balance = balance,
              env.var = env.var,
              strat = strat,
              do.full.models = do.full.models,
              # seed.val = seed.val,
              user.table = user.table))
}

# ---------------------------------------------------------------------------- #

.sample_mat <- function(data.sp, data.split, nb.rep = 1, data.env = NULL, seed.val = NULL)
{
  # data.sp is a 0, 1 vector
  # return a matrix with nb.rep columns of boolean (T: calib, F= eval)
  
  pres <- which(data.sp == 1)
  abs <- (1:length(data.sp))[-pres]
  
  nbPresEval <- round(length(pres) * data.split)
  nbAbsEval <- round(length(abs) * data.split)
  
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
            .check_calib.lines_names(calib.lines)
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
            nb_PA <- ncol(bm.format@PA.table)
            names_to_remove <- c("_allData_allRun", paste0("_PA", 1:nb_PA, "_allRun")) #We suppose that if they are in the format it should be fine
            table_to_test <- user.table[,!(colnames(user.table) %in% names_to_remove)]
            
            expected_PA.names <- colnames(bm.format@PA.table)
            .check_calib.lines_names(table_to_test, expected_PA.names = expected_PA.names)
            
            PA.table <- cbind(bm.format@PA.table, "allData" = TRUE) #just to pass the check
            calib.lines <-
              foreach(this.colnames = colnames(user.table), .combine = "cbind") %do% {
                        this.pa = strsplit(this.colnames, split = "_")[[1]][2]
                        calib.pa <- user.table[,this.colnames, drop = FALSE]
                        which.not.pa <- which(PA.table[, this.pa] == FALSE |
                                                is.na(PA.table[, this.pa]))
                        if(length(which.not.pa) > 0){
                          calib.pa[which.not.pa, ] <- NA
                        }
                        return(calib.pa)
                      }
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
                                         data.env = bm.format@data.env.var)
              # seed.val = seed.val)
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
                  sampled.mat <- .sample_mat(data.sp = bm.format@data.species[ind.PA],
                                             data.split = perc,
                                             nb.rep = nb.rep,
                                             data.env = bm.format@data.env.var[ind.PA, , drop = FALSE])
                  # seed.val = seed.val)
                  calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = nb.rep)
                  calib.pa[ind.PA, ] <- sampled.mat
                  colnames(calib.pa) <- paste0('_PA', pa, '_RUN', 1:ncol(calib.pa))
                }
                return(calib.pa)
              }
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
            tmp[ind.NA] <- 0
            
            calib.lines <- foreach(rep = 1:nb.rep, .combine = "cbind") %do%
              {
                fold <- dismo::kfold(tmp, by = tmp, k = k)
                calib.rep <- NULL
                for (i in 1:k) {
                  calib.rep <- cbind(calib.rep, fold != i)
                }
                return(calib.rep)
              }
            colnames(calib.lines) <- paste0('_RUN', 1:ncol(calib.lines))
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
            tmp[ind.NA] <- 0
            
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                
                if (nb.rep == 0) { # take all available data
                  calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = 1)
                  calib.pa[ind.PA, ] <- TRUE
                  colnames(calib.pa) <- '_allRun'
                } else {
                  tmp.pa <- tmp[ind.PA]
                  calib.pa_rep <- foreach(rep = 1:nb.rep, .combine = "cbind") %do%
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
                  colnames(calib.pa) <- paste0('_PA', pa, '_RUN', 1:ncol(calib.pa))
                }
                return(calib.pa)
              }
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
            tmp[ind.NA] <- 0
            
            tab.coord <- bm.format@coord
            blocks <- ENMeval::get.block(tab.coord[tmp == 1, ], tab.coord[tmp == 0, ])
            calib.lines <- matrix(NA, nrow = length(tmp), ncol = 4)
            for (i in 1:4) {
              calib.lines[tmp == 1, i] <- blocks[[1]] != i
              calib.lines[tmp == 0, i] <- blocks[[2]] != i     
            }
            colnames(calib.lines) <- paste0('_RUN', 1:ncol(calib.lines))
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
            tmp[ind.NA] <- 0
            
            calib.lines <- foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do%
              {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                tmp.pa <- tmp[ind.PA]
                
                tab.coord <- bm.format@coord[ind.PA, ]
                blocks <- ENMeval::get.block(tab.coord[tmp.pa == 1, ], tab.coord[tmp.pa == 0, ])
                calib.pa_rep <- matrix(NA, nrow = length(tmp.pa), ncol = 4)
                for (i in 1:4) {
                  calib.pa_rep[tmp.pa == 1, i] <- blocks[[1]] != i
                  calib.pa_rep[tmp.pa == 0, i] <- blocks[[2]] != i     
                }
                
                calib.pa <- matrix(NA, nrow = length(bm.format@data.species), ncol = 4)
                calib.pa[ind.PA, ] <- calib.pa_rep
                colnames(calib.pa) <- paste0('_PA', pa, '_RUN', 1:ncol(calib.pa))
                return(calib.pa)
              }
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
              bands[1] <- -Inf
              bands[k + 1] <- Inf
              calib.x <- matrix(NA, nrow(tmp.coord), k)
              for (i in 1:k) {
                calib.x[, i] <- !(tmp.coord[, 1] >= bands[i] & tmp.coord[, 1] < bands[i + 1])
              }
              if (strat == "x") { calib.lines <- calib.x }
            }
            
            if (strat == "y" || strat == "both") {
              bands <- quantile(tmp.coord[balance, 2], probs = seq(0, 100, 100 / k) / 100)
              bands[1] <- -Inf
              bands[k + 1] <- Inf
              calib.y <- matrix(NA, nrow(tmp.coord), k)
              for (i in 1:k) {
                calib.y[, i] <- !(tmp.coord[, 2] >= bands[i] & tmp.coord[, 2] < bands[i + 1])
              }
              if (strat == "y") { calib.lines <- calib.y }
            }
            
            if (strat == "both") { ## Merge X and Y tables
              calib.lines <- cbind(calib.x, calib.y)
            }

            colnames(calib.lines) <- paste0('_RUN', 1:ncol(calib.lines))
            return(calib.lines)
          })

## bm_CrossValidation strat BIOMOD.formated.data.PA methods --------------------------------
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
                tmp.pa.balance <- tmp.coord[intersect(ind.PA, which(balance)), ]
                if (strat == "x" || strat == "both") {
                  bands <- quantile(tmp.pa.balance[ , 1], probs = seq(0, 100, 100 / k) / 100)
                  bands[1] <- -Inf
                  bands[k + 1] <- Inf
                  calib.x <- matrix(NA, nrow(tmp.pa), k)
                  for (i in 1:k) {
                    calib.x[, i] <- !(tmp.pa[, 1] >= bands[i] & tmp.pa[, 1] < bands[i + 1])
                  }
                  if (strat == "x") { calib.pa_rep <- calib.x }
                }
                if (strat == "y" || strat == "both") {
                  bands <- quantile(tmp.pa[ , 2], probs = seq(0, 100, 100 / k) / 100)
                  bands[1] <- -Inf
                  bands[k + 1] <- Inf
                  calib.y <- matrix(NA, nrow(tmp.pa), k)
                  for (i in 1:k) {
                    calib.y[, i] <- !(tmp.pa[, 2] >= bands[i] & tmp.pa[, 2] < bands[i + 1])
                  }
                  if (strat == "y") { calib.pa_rep <- calib.y }
                }
                if (strat == "both") { ## Merge X and Y tables
                  calib.pa_rep <- cbind(calib.x, calib.y)
                }
                calib.pa <- matrix(NA, nrow = nrow(tmp.coord), 
                                   ncol = ifelse(strat =="both", 2*k, k))
                calib.pa[ind.PA, ] <- calib.pa_rep
                colnames(calib.pa) <- paste0('_PA', pa, '_RUN', 1:ncol(calib.pa))
                return(calib.pa)
              }
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
          function(bm.format, balance, k, env.var) {
            cat("\n   > Environmental cross-validation selection")
            calib.lines <- foreach(env = env.var, .combine = "cbind") %do%
              {
                tmp.env <- bm.format@data.env.var[balance, env]
                full.env <- bm.format@data.env.var[ , env]
                bands <- quantile(tmp.env, probs = seq(0, 100, 100 / k) / 100)
                bands[1] <- bands[1] - 1
                bands[k + 1] <- bands[k + 1] + 1
                calib.env <- matrix(NA, nrow = length(full.env), ncol = k)
                for (i in 1:k) {
                  calib.env[, i] <- (full.env <= bands[i] | full.env > bands[i + 1])
                }
                return(calib.env)
              }
            colnames(calib.lines) <- paste0('_RUN', 1:ncol(calib.lines))
            return(calib.lines)
          })

## bm_CrossValidation env BIOMOD.formated.data.PA methods -------------------------------
##'
##' @rdname bm_CrossValidation
##' @export
##'

setMethod('bm_CrossValidation_env', signature(bm.format = "BIOMOD.formated.data.PA"),
          function(bm.format, balance, k, env.var) {
            cat("\n   > Environmental cross-validation selection")
            calib.lines <- 
              foreach(pa = 1:ncol(bm.format@PA.table), .combine = "cbind") %do% {
                ind.PA <- which(bm.format@PA.table[, pa] == TRUE)
                calib.env <- foreach(env = env.var, .combine = "cbind") %do% {
                  tmp.pa <- bm.format@data.env.var[intersect(which(balance), ind.PA), env]
                  full.pa <- bm.format@data.env.var[ind.PA, env]
                  bands <- quantile(tmp.pa, probs = seq(0, 100, 100 / k) / 100)
                  bands[1] <- bands[1] - 1
                  bands[k + 1] <- bands[k + 1] + 1
                  calib.tmp.env <- matrix(NA, nrow = length(full.pa), ncol = k)
                  for (i in 1:k) {
                    calib.tmp.env[, i] <- (full.pa <= bands[i] | full.pa > bands[i + 1])
                  }
                  calib.pa.env <- matrix(NA, nrow = nrow(bm.format@data.env.var), ncol = k)
                  calib.pa.env[ind.PA, ] <- calib.tmp.env
                  return(calib.pa.env)
                }
                colnames(calib.env) <- paste0('_PA', pa, '_RUN', 1:ncol(calib.env))
                
                return(calib.env)
              }
            return(calib.lines)
          })

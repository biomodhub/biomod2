###################################################################################################
##' @name BIOMOD_CrossValidation
##' @author Frank Breiner
##' 
##' @title Custom models cross-validation procedure
##' 
##' @description This function creates a \code{DataSplitTable} that can be given as parameter to 
##' the \code{\link{BIOMOD_Modeling}} function to evaluate models with repeated k-fold or 
##' stratified cross-validation (CV) instead of repeated split samples.
##' 
##' 
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
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
##' A \code{DataSplitTable} {matrix} with \code{k * nb.rep} (\emph{+ 1 if 
##' \code{do.full.models = TRUE}}) columns that can be given as parameter to the 
##' \code{\link{BIOMOD_Modeling}} function.
##' 
##' 
##' @details
##' 
##' \bold{Stratified cross-validation} may be used to test for model overfitting and to assess 
##' transferability in geographic and environmental space : 
##' \itemize{
##'   \item \code{x} and \code{y} stratification was described in \emph{Wenger and Olden 2012} 
##'   (see \href{BIOMOD_CrossValidation.html#References}{References}). While \code{y} 
##'   stratification uses \code{k} partitions along the y-gradient, \code{x} stratification does 
##'   the same for the x-gradient, and \code{both} combines them.
##'   \item \code{block} stratification was described in \emph{Muscarella et al. 2014} (see 
##'   \href{BIOMOD_CrossValidation.html#References}{References}). Four bins of equal size are 
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
##' # Load species occurrences (6 species available)
##' myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
##' DataSpecies <- read.csv(myFile, row.names = 1)
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
##' myFiles = paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl = raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
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
##' myBiomodCV <- BIOMOD_CrossValidation(myBiomodData)
##' head(myBiomodCV)
##' 
##' # Several validation strategies can be combined
##' DataSplitTable.b <- BIOMOD_CrossValidation(myBiomodData,
##'                                           k = 5,
##'                                           rep = 2,
##'                                           do.full.models = FALSE)
##' DataSplitTable.y <- BIOMOD_CrossValidation(myBiomodData,
##'                                           k = 2,
##'                                           do.stratification = TRUE,
##'                                           method = "y")
##' colnames(DataSplitTable.y)[1:2] <- c("RUN11", "RUN12")
##' myBiomodCV <- cbind(DataSplitTable.b, DataSplitTable.y)
##' head(myBiomodCV)
##' 
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
##'                                     models = c('RF', 'GLM'),
##'                                     models.options = myBiomodOptions,
##'                                     NbRunEval = 2,
##'                                     DataSplitTable = myBiomodCV,
##'                                     VarImport = 3,
##'                                     models.eval.meth = c('TSS','ROC'),
##'                                     do.full.models = FALSE,
##'                                     modeling.id = 'mod.CV')
##' 
##' # Get evaluation scores & variables importance
##' myEval <- get_evaluations(myBiomodModelOut.CV, as.data.frame = TRUE)
##' myEval$CV.strategy <- "Random"
##' myEval$CV.strategy[grepl("13", myEval$Model.name)] <- "Full"
##' myEval$CV.strategy[grepl("11|12", myEval$Model.name)] <- "Stratified"
##' head(myEval)
##' 
##' boxplot(myEval$Testing.data ~ interaction(myEval$Algo, myEval$CV.strategy),
##'         xlab = "", ylab = "ROC AUC", col = rep(c("brown", "cadetblue"), 3))
##' 
##' 
##' @importFrom ENMeval get.block
##' @importFrom dismo kfold
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_CrossValidation <- function(bm.format,
                                   k = 5,
                                   nb.rep = 5,
                                   do.stratification = FALSE,
                                   method = "both",
                                   balance = "presences",
                                   do.full.models = TRUE)
{
  .bm_cat("Build Cross-Validation Table")
  DataSplitTable.y <- DataSplitTable.x <- DataSplitTable <- NULL
  
  ## STRATIFIED (X, Y, BOTH) / BLOCK / ENVIRONMENTAL CROSS VALIDATION -----------------------------
  if (do.stratification)
  {
    if (balance == "absences") {
      balance <- (bm.format@data.species == 1 | bm.format@data.species == 0)
    } else {
      balance <- (bm.format@data.species == 1)
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
      DataSplitTable <- as.data.frame(matrix(NA, nrow(bm.format@coord), 4))
      blocks <- get.block(bm.format@coord[bm.format@data.species == 1, ]
                          , bm.format@coord[bm.format@data.species == 0, ])
      for (i in 1:4) {
        DataSplitTable[bm.format@data.species == 1, i] <- blocks[[1]] != i
        DataSplitTable[bm.format@data.species == 0, i] <- blocks[[2]] != i     
      }
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
    for (rep in 1:nb.rep) {
      fold <- kfold(bm.format@data.species, by = bm.format@data.species, k = k)
      for (i in 1:k) {
        DataSplitTable <- cbind(DataSplitTable, fold != i)
      }
    }
  }
  
  ## CLEAN FINAL TABLE ----------------------------------------------------------------------------
  colnames(DataSplitTable) <- paste0("RUN", 1:ncol(DataSplitTable))
  if (do.full.models == TRUE) {
    DataSplitTable <- cbind(DataSplitTable, T)
    colnames(DataSplitTable)[ncol(DataSplitTable)] <- "Full"
  }
  
  .bm_cat("Done")
  return(DataSplitTable)
}

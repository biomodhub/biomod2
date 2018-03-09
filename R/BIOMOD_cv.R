##' @name BIOMOD_cv
##' @aliases BIOMOD_cv
##' 
##' @title Custom models cross-validation procedure
##' 
##' @description This function creates a DataSplitTable which could be used to evaluate models in Biomod with repeated
##'   k-fold cross-validation (cv) or stratified cv instead of repeated split sample runs
##' 
##' @param data            BIOMOD.formated.data object returned by BIOMOD_FormatingData
##' @param k               number of bins/partitions for k-fold cv
##' @param stratified.cv   logical. run a stratified cv 
##' @param stratify        stratification method of the cv. Could be "x", "y", "both" (default), "block" or the name of a predictor for environmental stratified cv.
##' @param balance         make balanced particions for "presences" (default) or "absences" (resp. pseudo-absences or background).
##' @param repetition      number of repetitions of k-fold cv (1 if stratified.cv=TRUE)
##' @param do.full.models  if true, models calibrated and evaluated with the whole dataset are done
##' 
##' @details
##'   Stratified cv could be used to test for model overfitting and for assessing transferability in geographic and environmental space. 
##'   If balance = "presences" presences are divided (balanced) equally over the particions (e.g. Fig. 1b in Muscarelly et al. 2014).
##'   Pseudo-Absences will however be unbalanced over the particions especially if the presences are clumped on an edge of the study area.
##'   If balance = "absences" absences (resp. Pseudo-Absences or background) are divided (balanced) as equally as possible for the particions
##'   (geographical balanced bins given that absences are spread over the study area equally, approach similar to Fig. 1 in Wenger et Olden 2012).
##'   Presences will however be unbalanced over the particians. Be careful: If the presences are clumped on an edge of the study area it is possible that all presences are in one bin.
##' 
##' @return
##' DataSplitTable matrix with k*repetition (+ 1 for Full models if  do.full.models = TRUE) columns for BIOMOD_Modeling function.
##' Stratification "x" and "y" was described in Wenger and Olden 2012. While Stratification "y" uses k partitions along the y-gradient, "x" does the same for the x-gradient and "both" combines them.
##' Stratification "block" was described in Muscarella et al. 2014. For bins of equal number are partitioned (bottom-left, bottom-right, top-left and top-right).
##' 
##' @references
##' Muscarella, R., Galante, P.J., Soley-Guardia, M., Boria, R.A., Kass, J.M., Uriarte, M. & Anderson, R.P. (2014). ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}, 1198-1205.
##' Wenger, S.J. & Olden, J.D. (2012). Assessing transferability of ecological models: an underappreciated aspect of statistical validation. \emph{Methods in Ecology and Evolution}, \bold{3}, 260-267.
##' 
##' @author Frank Breiner \email{frank.breiner@wsl.ch}
##' 
##' @seealso
##' \code{\link[ENMeval]{get.block}}
##'
##' @examples
##' \dontrun{
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"))
##' head(DataSpecies)
##' 
##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
##' 
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' 
##' # 3. Creating DataSplitTable
##' 
##' DataSplitTable <- BIOMOD_cv(myBiomodData, k=5, rep=2, do.full.models=F)
##' DataSplitTable.y <- BIOMOD_cv(myBiomodData,stratified.cv=T, stratify="y", k=2)
##' colnames(DataSplitTable.y)[1:2] <- c("RUN11","RUN12")
##' DataSplitTable <- cbind(DataSplitTable,DataSplitTable.y)
##' head(DataSplitTable)
##' 
##' # 4. Doing Modelisation
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('RF'), 
##'                                      models.options = myBiomodOption, 
##'                                      DataSplitTable = DataSplitTable,
##'                                      VarImport=0, 
##'                                      models.eval.meth = c('ROC'),
##'                                      do.full.models=FALSE,
##'                                      modeling.id="test")
##' 
##' ## get cv evaluations
##' eval <- get_evaluations(myBiomodModelOut,as.data.frame=T)
##' 
##' eval$strat <- NA
##' eval$strat[grepl("13",eval$Model.name)] <- "Full"
##' eval$strat[!(grepl("11",eval$Model.name)|
##'              grepl("12",eval$Model.name)|
##'              grepl("13",eval$Model.name))] <- "Random"
##' eval$strat[grepl("11",eval$Model.name)|grepl("12",eval$Model.name)] <- "Strat"
##' 
##' boxplot(eval$Testing.data~ eval$strat, ylab="ROC AUC")
##' }


BIOMOD_cv <-
  function (data, k = 5, repetition = 5, do.full.models = TRUE, 
            stratified.cv = FALSE, stratify = "both", balance = "pres") 
  {
    DataSplitTable.y <- DataSplitTable.x <- DataSplitTable <- NULL
    if (stratified.cv) {
      repetition <- 1
      if (balance == "absences") {
        balance <- data@data.species == 1 | data@data.species == 
          0
      }
      else {
        balance <- data@data.species == 1
      }
      if (stratify == "x" | stratify == "both") {
        DataSplitTable.x <- matrix(NA, nrow(data@coord), 
                                   k)
        bands <- quantile(data@coord[balance, 1], probs = seq(0, 
                                                              100, 100/k)/100)
        bands[1] <- bands[1] - 1
        bands[k + 1] <- bands[k + 1] + 1
        for (i in 1:k) {
          DataSplitTable.x[, i] <- data@coord[, 1] >= bands[i] & 
            data@coord[, 1] < bands[i + 1]
        }
        if (stratify == "x") {
          DataSplitTable <- DataSplitTable.x
        }
      }
      if (stratify == "y" | stratify == "both") {
        DataSplitTable.y <- matrix(NA, nrow(data@coord), 
                                   k)
        bands <- quantile(data@coord[balance, 2], probs = seq(0, 
                                                              100, 100/k)/100)
        bands[1] <- bands[1] - 1
        bands[k + 1] <- bands[k + 1] + 1
        for (i in 1:k) {
          DataSplitTable.y[, i] <- data@coord[, 2] >= bands[i] & 
            data@coord[, 2] < bands[i + 1]
        }
        if (stratify == "y") {
          DataSplitTable <- DataSplitTable.y
        }
      }
      if (stratify == "both") {
        DataSplitTable <- cbind(DataSplitTable.x, DataSplitTable.y)
      }
    if (stratify == "block") {
      DataSplitTable <- as.data.frame(matrix(NA, nrow(data@coord), 
                                             4))

      blocks<-ENMeval::get.block(data@coord[data@data.species==1,],
                                 data@coord[data@data.species==0,])
      
      for(i in 1:4){
        DataSplitTable[data@data.species == 1,i] <-  blocks[[1]]!=i     
        DataSplitTable[data@data.species == 0,i] <-  blocks[[2]]!=i     
      }
    }      
      if (stratify != "block" & stratify != "x" & stratify != 
          "y" & stratify != "both") {
        DataSplitTable2 <- as.data.frame(matrix(NA, nrow(data@coord), k))
        bands <- quantile(data@data.env.var[balance, stratify], 
                          probs = seq(0, 100, 100/k)/100)
        bands[1] <- bands[1] - 1
        bands[k + 1] <- bands[k + 1] + 1
        for (i in 1:k) {
          DataSplitTable2[, i] <- data@data.env.var[balance, 
                                                    stratify] <= bands[i] | data@data.env.var[balance, 
                                                                                              stratify] > bands[i + 1]
        }
      }
    }
    else {
      for (rep in 1:repetition) {
        fold <- dismo::kfold(data@data.species, by = data@data.species, 
                             k = k)
        for (i in 1:k) {
          DataSplitTable <- cbind(DataSplitTable, fold != 
                                    i)
        }
      }
    }
    if(stratify != "block"){
      colnames(DataSplitTable) <- paste("RUN", 1:(k * repetition), 
                                        sep = "")
      if (do.full.models == TRUE) {
        DataSplitTable <- cbind(DataSplitTable, T)
        colnames(DataSplitTable)[k * repetition + 1] <- "Full"
      }
    }else{
      colnames(DataSplitTable) <- paste("RUN", 1:4, 
                                        sep = "")    
      if (do.full.models == TRUE) {
        DataSplitTable <- cbind(DataSplitTable, T)
        colnames(DataSplitTable)[5] <- "Full"
      }
    }
    
    return(DataSplitTable)
  }

# BIOMOD_RangeSize documentation ----------------------------------------------
##' @name BIOMOD_RangeSize
##' @author Wilfried Thuiller, Damien Georges, Bruno Lafourcade
##' 
##' @title Analyze the range size differences between projections of species distribution models
##' 
##' @description This function allows to calculate the absolute number of locations (pixels) lost, 
##' stable and gained, as well as the corresponding relative proportions, between two (or more) 
##' binary projections of (ensemble) species distribution models (\emph{which can represent new 
##' time scales or environmental scenarios for example}).
##' 
##' 
##' @param proj.current a \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
##' or \code{\link[terra:rast]{SpatRaster}} object containing the initial binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param proj.future a \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
##' or \code{\link[terra:rast]{SpatRaster}} object containing the final binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param thresholds
##' 
##' 
##' @return
##' 
##' A \code{list} containing two objects :
##' \describe{
##'   \item{Compt.By.Species}{a \code{data.frame} containing the summary of range change for each 
##'   comparison
##'   \itemize{
##'     \item \code{Loss} : number of pixels predicted to be lost
##'     \item \code{Stable0} : number of pixels not currently occupied and not predicted to be
##'     \item \code{Stable1} : number of pixels currently occupied and predicted to remain 
##'     occupied
##'     \item \code{Gain} : number of pixels predicted to be gained
##'     \item \code{PercLoss} : percentage of pixels currently occupied and predicted to be lost 
##'     (\code{Loss / (Loss + Stable1)})
##'     \item \code{PercGain} : percentage of pixels predicted to be gained compare to the 
##'     number of pixels currently occupied (\code{Gain / (Loss + Stable1)})
##'     \item \code{SpeciesRangeChange} : percentage of pixels predicted to change (loss or gain) 
##'     compare to the number of pixels currently occupied (\code{PercGain - PercLoss})
##'     \item \code{CurrentRangeSize} : number of pixels currently occupied
##'     \item \code{FutureRangeSize0Disp} : number of pixels predicted to be occupied, assuming 
##'     no migration
##'     \item \code{FutureRangeSize1Disp} : number of pixels predicted to be occupied, assuming 
##'     migration
##'   }
##'   }
##'   \item{Diff.By.Pixel}{an object in the same form than the input data (\code{proj.current} and 
##'   \code{proj.future}) and containing a value for each point/pixel of each comparison among :
##'   \itemize{
##'     \item \code{-2} : predicted to be lost
##'     \item \code{-1} : predicted to remain occupied
##'     \item \code{0} : predicted to remain unoccupied
##'     \item \code{1} : predicted to be gained
##'   }
##'   }
##' }
##' 
##' 
##' @details 
##' 
##' Note that \bold{this function is only relevant to compare binary projections, made on the 
##' same area with the same resolution}. \cr
##' 
##' Comparison between \code{proj.current} and \code{proj.future} depends 
##' on the number of projection in both objects:
##'| \code{proj.current}   | \code{proj.future} | \bold{Comparison}  |
##'| ------------------------- | ---------------------- | --------------------  |
##'| \bold{1 projection} (\emph{e.g. data.frame with 1 column, SpatRaster with 1 layer}) | \bold{1 projection}  (\emph{e.g. data.frame with 1 column, SpatRaster with 1 layer})  | comparison of both projection  (\emph{e.g. current vs future conditions for the same model ; current vs current condition for two different models}) |
##'| \bold{\code{n} projections}  (\emph{e.g. data.frame with n column, SpatRaster with n layer}) |  \bold{\code{n} projections}  (\emph{e.g. data.frame with n column, SpatRaster with n layer}) |  comparing projection \code{i} in \code{proj.current} to projection \code{i} in \code{proj.future}  (\emph{e.g. comparing current vs future condition for n models}) |
##'| \bold{\code{1} projection}   (\emph{e.g. data.frame with 1 column, SpatRaster with 1 layer}) |  \bold{\code{n} projections}  (\emph{e.g. data.frame with n column, SpatRaster with n layer}) |  comparing projection in \code{proj.current} to each projection in \code{proj.future}  (\emph{e.g. comparing current vs n different future condition (e.g. climate change scenario) for 1 model}) |
##' 
##' \code{Diff.By.Pixel} object is obtained by applying the simple following formula :
##' \deqn{proj.future - 2 * proj.current}
##' 
##' @md
##' 
##' @keywords "species range change" projections gain loss
##' 
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}, 
##' \code{\link{bm_PlotRangeSize}}
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
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                        expl.var = myExpl,
##'                                        resp.xy = myRespXY,
##'                                        resp.name = myRespName)
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' }
##' 
##' models.proj <- get_built_models(myBiomodModelOut, algo = "RF")
##'   # Project single models
##'   myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                     proj.name = 'CurrentRangeSize',
##'                                     new.env = myExpl,
##'                                     models.chosen = models.proj,
##'                                     metric.binary = 'all',
##'                                     build.clamping.mask = TRUE)
##' 
##' 
##' # --------------------------------------------------------------- #
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' data(bioclim_future)
##' myExplFuture <- terra::rast(bioclim_future)
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExplFuture <- terra::crop(myExplFuture, myExtent)
##' }
##' # Project onto future conditions
##' myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                               proj.name = 'FutureRangeSize',
##'                                               new.env = myExplFuture,
##'                                               models.chosen = models.proj,
##'                                               metric.binary = 'TSS')
##' 
##' # Load current and future binary projections
##' CurrentProj <- get_predictions(myBiomodProj,
##'                                metric.binary = "TSS",
##'                                model.as.col = TRUE)
##' FutureProj <- get_predictions(myBiomodProjectionFuture,
##'                                metric.binary = "TSS",
##'                                model.as.col = TRUE)
##' # Compute differences
##' myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj)
##' 
##' myBiomodRangeSize$Compt.By.Models
##' plot(myBiomodRangeSize$Diff.By.Pixel)
##' 
##' # Represent main results 
##' bm_PlotRangeSize(bm.range = myBiomodRangeSize)
##' 
##' 
##' @importFrom foreach foreach %do%
##' @importFrom terra rast nlyr `add<-` values
##' 
##' @export
##' 
##' 
## BIOMOD_RangeSize generic method ---------------------------------------------


setGeneric("BIOMOD_RangeSize",
           def = function(proj.current, proj.future, thresholds) {
             if (inherits(proj.current, "Raster") && inherits(proj.future, "Raster")) {
               return(
                 BIOMOD_RangeSize(rast(proj.current), rast(proj.future), thresholds)
               )
             } else {
               stop("'proj.current' and 'proj.future' must have the same class among 'data.frame', 'SpatRaster' and 'array'" )
             }
           })


## BIOMOD_RangeSize data.frame-data.frame Method ----------------------

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'data.frame', proj.future = 'data.frame'),
          function(proj.current, proj.future, thresholds)
          {
            .bm_cat("Do Range Size Computation")
            args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future, thresholds)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            if (ncol(proj.future) == ncol(proj.current)) {
              if (nonbinary) {
                Diff.By.Pixel <- as.data.frame((proj.future - proj.current)/proj.current) *100 
              } else {
                Diff.By.Pixel <- as.data.frame(proj.future - 2 * proj.current)
              }
              
              this_rownames <- colnames(proj.current)
            } else {
              Diff.By.Pixel <- foreach(thiscol = seq_len(ncol(proj.future)), .combine = 'cbind') %do% {
                if (nonbinary){
                  tmp <- as.data.frame((proj.future[,thiscol] - proj.current[,1])/proj.current[,1]) *100 
                } else {
                  tmp <- as.data.frame(proj.future[,thiscol] - 2 * proj.current[,1])
                }
                colnames(tmp) <- colnames(proj.future)[thiscol]
                tmp
              }
              this_rownames <- colnames(proj.future)
            }
            
            if (nonbinary){
              nbpixels <- sum(!is.na(Diff.By.Pixel[]))
              
              CBS[i, "Stable"] <- length(which(Diff.By.Pixel[] > -min(thresholds) & Diff.By.Pixel[] <= min(thresholds)))
              CBS[i, "Stableperc"] <- round(CBS[i, "Stable"] / nbpixels * 100, digits = 3)
              
              for (thre in thresholds){
                if (thre == max(thresholds)){
                  thre_max <- Inf
                } else {
                  thre_max <- thresholds[which(thresholds == thre) + 1]
                }
                
                CBS[i, paste0("loss", -thre)] <- length(which(Diff.By.Pixel[] > -thre_max & Diff.By.Pixel[] <= -thre))
                CBS[i, paste0("gain", thre)] <- length(which(Diff.By.Pixel[] > thre & Diff.By.Pixel[] <= thre_max))
                
                CBS[i, paste0("loss", -thre, "perc")] <- round(CBS[i, paste0("loss", -thre)] / nbpixels * 100, digits = 3)
                CBS[i, paste0("gain", thre, "perc")] <- round(CBS[i, paste0("gain", thre)] / nbpixels * 100, digits = 3)
              }
              
              CBS[i, "SpeciesRangeChange"] <- 100 - CBS[i, "Stableperc"]
              CBS[i, "CurrentRangeSize"] <- length(which(proj.current[] > 0))
              CBS[i, "FutureRangeSize.NoDisp"] <- length(which(proj.future[] > 0) & which(proj.current[] > 0))
              CBS[i, "FutureRangeSize.FullDisp"] <- length(which(proj.future[] > 0))
              
              
              names.res <- c(paste0("loss", sort(-thresholds)), "Stable", paste0("gain", thresholds))
              names.res <- c(names.res, paste0(names.res, "perc"))
              names.res <- c(names.res, c("SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
              
              Compt.By.Models <- CBS
              dimnames(Compt.By.Models) <- names.res
              
            } else {
              Compt.By.Models <- as.data.frame(.CompteurSp(Diff.By.Pixel, c(-2, 0, -1, 1)))
              Compt.By.Models[, seq(5,10)] <- NA
              Compt.By.Models[, 8] <- Compt.By.Models[, 1] + Compt.By.Models[, 3]
              Compt.By.Models[, 9] <- Compt.By.Models[, 3]
              Compt.By.Models[, 10] <- Compt.By.Models[, 4] + Compt.By.Models[, 3]
              
              Compt.By.Models[, 5] <- (100 * Compt.By.Models[, 1]) / Compt.By.Models[, 8]
              Compt.By.Models[, 6] <- (100 * Compt.By.Models[, 4]) / Compt.By.Models[, 8]
              Compt.By.Models[, 7] <- Compt.By.Models[, 6] - Compt.By.Models[, 5]
              
              dimnames(Compt.By.Models) <- list(this_rownames, c("Loss", "Stable0", "Stable1", "Gain"
                                                                 , "PercLoss", "PercGain", "SpeciesRangeChange"
                                                                 , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            }
            
            Output <- list(Compt.By.Models = Compt.By.Models, Diff.By.Pixel = Diff.By.Pixel)
            .bm_cat("Done")
            invisible(Output)
          })


## BIOMOD_RangeSize SpatRaster-SpatRaster Method ----------------------

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'SpatRaster', proj.future = 'SpatRaster'),
          function(proj.current, proj.future, thresholds)
          {
            .bm_cat("Do Range Size Computation")
            args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future, thresholds)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            if (nonbinary == FALSE){
              names.res = c("Loss", "Stable0", "Stable1", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            } else {
              names.res <- c(paste0("loss", sort(-thresholds)), "Stable", paste0("gain", thresholds))
              names.res <- c(names.res, paste0(names.res, "perc"))
              names.res <- c(names.res, c("SpeciesRangeChange", "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            }
            
            if(nlyr(proj.current) > 1){
              CBS <- matrix(ncol = length(names.res), nrow = nlyr(proj.current),
                          dimnames = list(names(proj.current), names.res))
            } else {
              CBS <- matrix(ncol = length(names.res), nrow = nlyr(proj.future),
                            dimnames = list(names(proj.future), names.res))
            }

            sp.rast <- rast()
            Cur <- proj.current[[1]]
            for (i in seq_len(nlyr(proj.future))) {
              if(nlyr(proj.current) > 1){
                Cur <- proj.current[[i]]
              }
              Fut <- proj.future[[i]]
              
              if (nonbinary){
                ## DiffByPixel
                Ras <- (Fut - Cur)/ Cur *100
                add(sp.rast) <- Ras
                
                ## ComptBySpecies
                nbpixels <- sum(!is.na(Ras[]))
                
                CBS[i, "Stable"] <- length(which(Ras[] > -min(thresholds) & Ras[] <= min(thresholds)))
                CBS[i, "Stableperc"] <- round(CBS[i, "Stable"] / nbpixels * 100, digits = 3)
                
                for (thre in thresholds){
                  if (thre == max(thresholds)){
                    thre_max <- Inf
                  } else {
                    thre_max <- thresholds[which(thresholds == thre) + 1]
                  }
                  
                  CBS[i, paste0("loss", -thre)] <- length(which(Ras[] > -thre_max & Ras[] <= -thre))
                  CBS[i, paste0("gain", thre)] <- length(which(Ras[] > thre & Ras[] <= thre_max))
                  
                  CBS[i, paste0("loss", -thre, "perc")] <- round(CBS[i, paste0("loss", -thre)] / nbpixels * 100, digits = 3)
                  CBS[i, paste0("gain", thre, "perc")] <- round(CBS[i, paste0("gain", thre)] / nbpixels * 100, digits = 3)
                }
                
                CBS[i, "SpeciesRangeChange"] <- 100 - CBS[i, "Stableperc"]
                CBS[i, "CurrentRangeSize"] <- length(which(Cur[] > 0))
                CBS[i, "FutureRangeSize.NoDisp"] <- length(which(Fut[] > 0) & which(Cur[] > 0))
                CBS[i, "FutureRangeSize.FullDisp"] <- length(which(Fut[] > 0))

              } else {
                ## DiffByPixel
                Ras <- Fut - (Cur + Cur)
                add(sp.rast) <- Ras
                
                ## ComptBySpecies
                CBS[i, 1] <- length(which(Ras[] == -2))
                CBS[i, 2] <- length(which(Ras[] == 0))
                CBS[i, 3] <- length(which(Ras[] == -1))
                CBS[i, 4] <- length(which(Ras[] == 1))
                
                CBS[i, 8] <- CBS[i, 1] + CBS[i, 3]
                CBS[i, 9] <- CBS[i, 3]
                CBS[i, 10] <- CBS[i, 3] + CBS[i, 4]
                
                CBS[i, 5] <- round(CBS[i, 1] / CBS[i, 8] * 100, digits = 3)
                CBS[i, 6] <- round(CBS[i, 4] / CBS[i, 8] * 100, digits = 3)
                CBS[i, 7] <- round(CBS[i, 10] / CBS[i, 8] * 100 - 100, digits = 3)
              }
              
            }
            names(sp.rast) <- rownames(CBS)
            
            .bm_cat("Done")
            return(list(Compt.By.Models = CBS, Diff.By.Pixel = sp.rast))
          })


# Argument Check ---------------------------------------------------------------

.BIOMOD_RangeSize.check.args <- function(proj.current, proj.future, thresholds) {
  
  ## dimensions checking ------------------------
  if (inherits(proj.current, "data.frame")) {
    dim_current <- nrow(proj.current)
    dim_future <- nrow(proj.future)
  } else { # SpatRaster case
    dim_current <- dim(proj.current)[c(1,2)]
    dim_future <- dim(proj.future)[c(1,2)]
  }
  
  if (any(dim_current != dim_future)) {
    stop(paste0("'proj.current' and 'proj.future' do not have the same dimensions ",
                "(data.frame must have either the same number of lines ; ",
                "raster must have the same spatial extent and resolution)."))
  }
  
  
  ## checking number of models to be compared ----------------------------
  if (inherits(proj.current, "data.frame")) {
    n_current <- ncol(proj.current)
    n_future <- ncol(proj.future)
  } else { # SpatRaster case
    n_current <- nlyr(proj.current)
    n_future <- nlyr(proj.future)
  }
  
  if (n_current == 1) {
    if (n_future == 1) {
      cat("\n Comparing 'proj.current' and 'proj.future'. ")
    } else if (n_future > 1) {
      cat("\n Comparing 'proj.current' with the ", n_future, " projections in 'proj.future'. ")
    } else {
      stop("'proj.future' require at least one layer or one column.")
    }
  } else if (n_current > 1) {
    if (n_future == n_current) {
      cat("\n Each projection in 'proj.current' will be compared once with its corresponding projection in 'proj.future'. ")
    } else {
      stop("When proj.current' have more than 1 projection, 'proj.future' must have the same number of projection.")
    }
  } else {
    stop("'proj.current' require at least one layer or one column.")
  }
  
  
  ## checking 0/1 ----------------------------
  if (inherits(proj.current, "data.frame")) {
    test_binary_current <- all(na.omit(unlist(proj.current)) %in% c(0,1))
    test_binary_future <- all(na.omit(unlist(proj.future)) %in% c(0,1))
  } else { # SpatRaster case
    test_binary_current <- all(na.omit(values(proj.current)) %in% c(0,1))
    test_binary_future <- all(na.omit(values(proj.future)) %in% c(0,1))
  }
  
  nonbinary <- ifelse((!test_binary_current | !test_binary_future), TRUE, FALSE)
  # if (!test_binary_current | !test_binary_future) {
  #   stop("'proj.current' and 'proj.future' must have only values among 0, 1 or NA.")
  # }
  
  # Thresholds
  .fun_testIfInherits(TRUE, "thresholds", thresholds, "numeric")
  thresholds <- sort(thresholds)
  
  
  return(list("proj.current" = proj.current,
              "proj.future" = proj.future,
              nonbinary = nonbinary,
              thresholds = thresholds))
}


# Additionnal tools -------------------------------------------------------

.CompteurSp <- function(Data, Value) {
  if (is.data.frame(Data)) {
    N <- dim(Data)[2]
    Compt <- as.data.frame(matrix(0, ncol = 4, nrow = dim(Data)[2]))
    i <- 1
    while(i <= N) {
      Compt[i, 1] <- length(Data[Data[, i] == Value[1], i])
      Compt[i, 2] <- length(Data[Data[, i] == Value[2], i])
      Compt[i, 3] <- length(Data[Data[, i] == Value[3], i])
      Compt[i, 4] <- length(Data[Data[, i] == Value[4], i])
      i <- i + 1
    }
    return(Compt)
  }
}

.CompteurThresh <- function(Data, Value) {
  if (is.data.frame(Data)) {
    N <- dim(Data)[2]
    V <- (length(Value)*2) +1
    Compt <- as.data.frame(matrix(FALSE, ncol = V, nrow = N))
    m <- (V+1)/2
    i <- 1
    bornes <- unique(sort(c(Value, -Value)))
    nb <- length(bornes) +1
    while(i <= N) {
      Compt[i,m] <- length(Data[Data[, i] == 0, i])
      for (v in 1:(m-1)){
        Compt[i, v] <- length(Data[(Data[, i] <= bornes[v] & Data[, i] > bornes[v-1]) , i])
        Compt[i, ncol(Compt) +1 - v] <- length(Data[(Data[, i] >= bornes[nb -v] & Data[, i] < bornes[nb-v+1]) , i])
      }
      i <- i + 1
    }
    colnames(Compt) <- c("verslinfini",bornes,"etaudela")
    return(Compt)
  }
}


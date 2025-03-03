###################################################################################################
##' @name bm_RangeSize
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
##' @param ordinal a \code{logical} indicating whether or not the projections should be considered as ordinal data
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
##'     \item \code{Stable_Abs} : number of pixels not currently occupied and not predicted to be
##'     \item \code{Stable_Pres} : number of pixels currently occupied and predicted to remain 
##'     occupied
##'     \item \code{Gain} : number of pixels predicted to be gained
##'     \item \code{PercLoss} : percentage of pixels currently occupied and predicted to be lost 
##'     (\code{Loss / (Loss + Stable_Pres)})
##'     \item \code{PercGain} : percentage of pixels predicted to be gained compare to the 
##'     number of pixels currently occupied (\code{Gain / (Loss + Stable_Pres)})
##'     \item \code{SpeciesRangeChange} : percentage of pixels predicted to change (loss or gain) 
##'     compare to the number of pixels currently occupied (\code{PercGain - PercLoss})
##'     \item \code{CurrentRangeSize} : number of pixels currently occupied
##'     \item \code{FutureRangeSize0Disp} : number of pixels predicted to be occupied, assuming 
##'     no migration
##'     \item \code{FutureRangeSize1Disp} : number of pixels predicted to be occupied, assuming 
##'     migration
##'   }
##'   }
##'   \item{loss.gain}{an object in the same form than the input data (\code{proj.current} and 
##'   \code{proj.future}) and containing a value for each point/pixel of each comparison among :
##'   \itemize{
##'     \item \code{-2} : predicted to be lost
##'     \item \code{-1} : predicted to remain occupied
##'     \item \code{0} : predicted to remain unoccupied
##'     \item \code{1} : predicted to be gained
##'   }
##'   }
##'   \item{Diff.By.Pixel}{an object in the same form than the input data (\code{proj.current} and 
##'   \code{proj.future}) and containing a value for each point/pixel of each comparison obtain with :
##'   \itemize{
##'   \item \code{Future - 2* Current} for binary data
##'   \item \code{Future - Current} for nonbinary
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
##' For binary data,\code{Diff.By.Pixel} object is obtained by applying the simple following formula :
##' \deqn{proj.future - 2 * proj.current}  
##' 
##' For non binary data, \code{Diff.By.Pixel} object is obtained by rescaling the two projection to a 0 to 1 scale
##' and applying the simple following formula :
##' \deqn{proj.future - proj.current}  
##' 
##' 
##' @keywords "species range change" projections gain loss
##' 
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}, 
##' \code{\link{bm_PlotRangeSize}}
##' @family Secondary functions
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
##'   myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
##'                                        resp.var = myResp,
##'                                        resp.xy = myRespXY,
##'                                        expl.var = myExpl)
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS', 'ROC'),
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
##' myBiomodRangeSize <- bm_RangeSize(proj.current = CurrentProj, proj.future = FutureProj)
##' 
##' myBiomodRangeSize$Compt.By.Models
##' plot(myBiomodRangeSize$Diff.By.Pixel)
##' 
##' 
##' 
##' @importFrom foreach foreach %do%
##' @importFrom terra rast nlyr `add<-` values `values<-`
##' @importFrom dplyr count
##' @importFrom scales rescale
##' 
##' @export
##' 
##' 
###################################################################################################


setGeneric("bm_RangeSize",
           def = function(proj.current, proj.future, ordinal = FALSE) {
             if (inherits(proj.current, "Raster") && inherits(proj.future, "Raster")) {
               return(
                 bm_RangeSize(rast(proj.current), rast(proj.future), ordinal)
               )
             } else {
               stop("'proj.current' and 'proj.future' must have the same class among 'data.frame', 'SpatRaster' and 'array'" )
             }
           })


## bm_RangeSize data.frame-data.frame Method ----------------------

##'
##' @rdname bm_RangeSize
##' @export
##'

setMethod('bm_RangeSize', signature(proj.current = 'data.frame', proj.future = 'data.frame'),
          function(proj.current, proj.future, ordinal)
          {
            .bm_cat("Do Range Size Computation")
            args <- .bm_RangeSize.check.args(proj.current, proj.future)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            
            if (ncol(proj.future) != ncol(proj.current)) {
              proj_current <- matrix(rep(proj.current[,1], ncol(proj.future)),
                                    ncol = ncol(proj.future),
                                    dimnames = )
              proj_current <- as.data.frame(proj_current)
              colnames(proj_current) <- colnames(proj.future)
            }  
            
            if (nonbinary && !ordinal) {
              ##rescale
              rescale_df <- rbind(proj.current, proj.future)
              rescale_df <- apply(rescale_df, 2, FUN = rescale, to = c(0,1))
              
              mid <- nrow(proj.current)
              end <- nrow(rescale_df)
              
              Diff.By.Pixel <- as.data.frame(rescale_df[(mid+1):end,] - rescale_df[1:mid,])
              loss.gain <- .loss.gain.nonbinary(current = proj.current, future = proj.future)
            } else if (!ordinal){
              Diff.By.Pixel <- as.data.frame(proj.future - 2 * proj.current)
            } else {
              proj.future <- apply(proj.future, 2, as.numeric)
              proj.current <- apply(proj.current, 2, as.numeric)
              Diff.By.Pixel <- as.data.frame(proj.future - proj.current)
              loss.gain <- .loss.gain.nonbinary(current = (proj.current - 1), future = (proj.future -1))
            }
            this_rownames <- colnames(proj.future)
            # } else {
            #   Diff.By.Pixel <- foreach(thiscol = seq_len(ncol(proj.future)), .combine = 'cbind') %do% {
            #     if (nonbinary && !ordinal){
            #       tmp <- as.data.frame(proj.future[,thiscol] - proj.current[,1])
            #       loss.gain <- .loss.gain.nonbinary(current = proj.current[,1], future = proj.future[,thiscol])
            #     } else if (!ordinal) {
            #       tmp <- as.data.frame(proj.future[,thiscol] - 2 * proj.current[,1])
            #     } else {
            #       tmp <- as.data.frame(proj.future[,thiscol] - proj.current[,1])
            #       loss.gain <- .loss.gain.nonbinary(current = (proj.current[,1] - 1), future = (proj.future[,thiscol] -1))
            #     }
            #     colnames(tmp) <- colnames(proj.future)[thiscol]
            #     tmp
            #   }
            #   this_rownames <- colnames(proj.future)
            # }
            # 
            
            if(nonbinary){
              comp <- as.data.frame(loss.gain)
            } else {
              comp <- Diff.By.Pixel
            }
            
          
            Compt.By.Models <- as.data.frame(.CompteurSp(comp, c(-2, 0, -1, 1)))
            Compt.By.Models[, seq(5,10)] <- NA
            Compt.By.Models[, 8] <- Compt.By.Models[, 1] + Compt.By.Models[, 3]
            Compt.By.Models[, 9] <- Compt.By.Models[, 3]
            Compt.By.Models[, 10] <- Compt.By.Models[, 4] + Compt.By.Models[, 3]
            
            Compt.By.Models[, 5] <- (100 * Compt.By.Models[, 1]) / Compt.By.Models[, 8]
            Compt.By.Models[, 6] <- (100 * Compt.By.Models[, 4]) / Compt.By.Models[, 8]
            Compt.By.Models[, 7] <- Compt.By.Models[, 6] - Compt.By.Models[, 5]
            
            dimnames(Compt.By.Models) <- list(this_rownames, c("Loss", "Stable_Abs", "Stable_Pres", "Gain"
                                                               , "PercLoss", "PercGain", "SpeciesRangeChange"
                                                               , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            
            
            Output <- list(Compt.By.Models = Compt.By.Models,
                           Diff.By.Pixel = Diff.By.Pixel,
                           loss.gain = comp,
                           data.type = ifelse(nonbinary, ifelse(ordinal, "ordinal", "nonbinary"), "binary")
            )
            .bm_cat("Done")
            invisible(Output)
          })


## bm_RangeSize SpatRaster-SpatRaster Method ----------------------

##'
##' @rdname bm_RangeSize
##' @export
##'

setMethod('bm_RangeSize', signature(proj.current = 'SpatRaster', proj.future = 'SpatRaster'),
          function(proj.current, proj.future, ordinal)
          {
            .bm_cat("Do Range Size Computation")
            args <- .bm_RangeSize.check.args(proj.current, proj.future)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            names.res = c("Loss", "Stable_Abs", "Stable_Pres", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            
            
            if(nlyr(proj.current) > 1){
              CBS <- matrix(ncol = length(names.res), nrow = nlyr(proj.current),
                          dimnames = list(names(proj.current), names.res))
            } else {
              CBS <- matrix(ncol = length(names.res), nrow = nlyr(proj.future),
                            dimnames = list(names(proj.future), names.res))
            }

            sp.rast <- rast()
            loss.gain_rast <- rast()
            Cur <- proj.current[[1]]
            for (i in seq_len(nlyr(proj.future))) {
              if(nlyr(proj.current) > 1){
                Cur <- proj.current[[i]]
              }
              Fut <- proj.future[[i]]
              
              if (nonbinary && !ordinal){
                ## DiffByPixel
                
                rescale_map <- c(Cur, Fut)
                values(rescale_map) <- rescale(values(rescale_map), to = c(0,1))
                Ras <- rescale_map[[2]] - rescale_map[[1]]
                add(sp.rast) <- Ras
                
                loss.gain <- .loss.gain.nonbinary(current = Cur - 1, future = Fut - 1)
                add(loss.gain_rast) <- loss.gain 
              } else if (!ordinal) {
                ## DiffByPixel
                Ras <- Fut - (Cur + Cur)
                add(sp.rast) <- Ras
                loss.gain <- Ras
                add(loss.gain_rast) <- loss.gain 
              } else {
                Ras <- Fut - Cur
                add(sp.rast) <- Ras
                loss.gain <- .loss.gain.nonbinary(current = proj.current, future = proj.future)
                add(loss.gain_rast) <- loss.gain 
              }
              
              ## ComptBySpecies
              CBS[i, 1] <- length(which(loss.gain[] == -2))
              CBS[i, 2] <- length(which(loss.gain[] == 0))
              CBS[i, 3] <- length(which(loss.gain[] == -1))
              CBS[i, 4] <- length(which(loss.gain[] == 1))
              
              CBS[i, 8] <- CBS[i, 1] + CBS[i, 3]
              CBS[i, 9] <- CBS[i, 3]
              CBS[i, 10] <- CBS[i, 3] + CBS[i, 4]
              
              CBS[i, 5] <- round(CBS[i, 1] / CBS[i, 8] * 100, digits = 3)
              CBS[i, 6] <- round(CBS[i, 4] / CBS[i, 8] * 100, digits = 3)
              CBS[i, 7] <- round(CBS[i, 10] / CBS[i, 8] * 100 - 100, digits = 3)
            } 

            names(sp.rast) <- rownames(CBS)
            Output <- list(Compt.By.Models = as.data.frame(CBS),
                          Diff.By.Pixel = sp.rast,
                          loss.gain = loss.gain_rast,
                          data.type = ifelse(nonbinary, ifelse(ordinal, "ordinal", "nonbinary"), "binary")
            )
            .bm_cat("Done")
            return(Output)
          })


###################################################################################################

.bm_RangeSize.check.args <- function(proj.current, proj.future)
{
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
  
  
  return(list("proj.current" = proj.current,
              "proj.future" = proj.future,
              nonbinary = nonbinary))
}


###################################################################################################

.CompteurSp <- function(Data, Value)
{
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


.loss.gain.nonbinary <- function(current, future){
  category <- current > 0
  category[] <- ifelse(current[] == future[], ifelse(category[], -1, 0),
                       ifelse(current[] < future[], 1, -2)) #gain, loss
  
  return(category)
}




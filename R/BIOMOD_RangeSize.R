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
##' @param proj.current an \code{array}, \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
##' or \code{\link[terra:rast]{SpatRaster}} object containing the initial binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param proj.future an \code{array}, \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
##' or \code{\link[terra:rast]{SpatRaster}} object containing the final binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' 
##' 
##' @return
##' 
##' A \code{list} containing two objects :
##' \describe{
##'   \item{Compt.By.Species}{a \code{data.frame} containing the summary of range change for each 
##'   comparison
##'   \itemize{
##'     \item{\code{Loss} : }{number of pixels predicted to be lost}
##'     \item{\code{Stable0} : }{number of pixels not currently occupied and not predicted to be}
##'     \item{\code{Stable1} : }{number of pixels currently occupied and predicted to remain 
##'     occupied}
##'     \item{\code{Gain} : }{number of pixels predicted to be gained}
##'     \item{\code{PercLoss} : }{percentage of pixels currently occupied and predicted to be lost 
##'     (\code{Loss / (Loss + Stable1)})}
##'     \item{\code{PercGain} : }{percentage of pixels predicted to be gained compare to the 
##'     number of pixels currently occupied (\code{Gain / (Loss + Stable1)})}
##'     \item{\code{SpeciesRangeChange} : }{percentage of pixels predicted to change (loss or gain) 
##'     compare to the number of pixels currently occupied (\code{PercGain - PercLoss})}
##'     \item{\code{CurrentRangeSize} : }{number of pixels currently occupied}
##'     \item{\code{FutureRangeSize0Disp} : }{number of pixels predicted to be occupied, assuming 
##'     no migration}
##'     \item{\code{FutureRangeSize1Disp} : }{number of pixels predicted to be occupied, assuming 
##'     migration}
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
##' same area with the same resolution}.
##' 
##' \cr Comparison between \code{proj.current} and \code{proj.future} depends 
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
##'   # Create default modeling options
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       bm.options = myBiomodOptions,
##'                                       nb.rep = 2,
##'                                       data.split.perc = 80,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       do.full.models = FALSE,
##'                                       seed.val = 42)
##' }
##' 
##' file.proj <- paste0(myRespName, "/proj_Current/", myRespName, ".Current.projection.out")
##' if (file.exists(file.proj)) {
##'   myBiomodProj <- get(load(file.proj))
##' } else {
##' 
##'   # Project single models
##'   myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                     proj.name = 'Current',
##'                                     new.env = myExpl,
##'                                     models.chosen = 'all',
##'                                     metric.binary = 'all',
##'                                     metric.filter = 'all',
##'                                     build.clamping.mask = TRUE)
##' }
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
##'                                               proj.name = 'Future',
##'                                               new.env = myExplFuture,
##'                                               models.chosen = 'all',
##'                                               metric.binary = 'TSS',
##'                                               build.clamping.mask = TRUE)
##' 
##' # Load current and future binary projections
##' proj.current <- terra::rast("GuloGulo/proj_Current/proj_Current_GuloGulo_TSSbin.grd")
##' proj.future <- terra::rast("GuloGulo/proj_Future/proj_Future_GuloGulo_TSSbin.grd")
##' 
##' # Compute differences
##' myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = proj.current, proj.future = proj.future)
##' 
##' myBiomodRangeSize$Compt.By.Models
##' plot(myBiomodRangeSize$Diff.By.Pixel)
##' 
##' # Represent main results 
##' bm_PlotRangeSize(bm.range = myBiomodRangeSize)
##' 
##' 
##' @importFrom terra rast nlyr `add<-`
##' 
##' @export
##' 
##' 
## BIOMOD_RangeSize generic method ---------------------------------------------


setGeneric("BIOMOD_RangeSize",
           def = function(proj.current, proj.future) {
             standardGeneric("BIOMOD_RangeSize")
           })


## BIOMOD_RangeSize data.frame-data.frame Method ----------------------

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'data.frame', proj.future = 'data.frame'),
          function(proj.current, proj.future)
          {
            .bm_cat("Do Range Size Computation")
            args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future)

            if (ncol(proj.future) == ncol(proj.current)) {
              Diff.By.Pixel <- as.data.frame(proj.future - 2 * proj.current)
              this_rownames <- colnames(proj.current)
            } else {
              Diff.By.Pixel <- foreach(thiscol = seq_len(ncol(proj.future)), .combine = 'cbind') %do% {
                tmp <- as.data.frame(proj.future[,thiscol] - 2 * proj.current[,1])
                colnames(tmp) <- colnames(proj.future)[thiscol]
                tmp
              }
              this_rownames <- colnames(proj.future)
            }
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
            
            Output <- list(Compt.By.Models = Compt.By.Models, Diff.By.Pixel = Diff.By.Pixel)
            .bm_cat("Done")
            invisible(Output)
          })

## BIOMOD_RangeSize array-array Method ----------------------
##'
##' @rdname BIOMOD_RangeSize
##' @export
##'


setMethod('BIOMOD_RangeSize', signature(proj.current = 'array', proj.future = 'array'),
          function(proj.current, proj.future)
          {
            
            ## Transform Current array into data.frame
            proj.current <- as.data.frame(proj.current)
            names(proj.current) <- unlist(lapply(strsplit(names(proj.current), ".", fixed = TRUE), function(x)
            { return(paste(x[3], x[2], x[1], sep = "_")) }))
            names(proj.current) <- unlist(lapply(strsplit(names(proj.current), ".", fixed = TRUE), function(x)
            {
              x.rev <- rev(x) ## reverse the order of splitted vector to have algo at the end
              data.set.id <- x.rev[1]
              cross.valid.id <- x.rev[2]
              algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
              model.id <- paste(data.set.id, cross.valid.id, algo.id, sep = "_")
              return(model.id)
            }))

            ## Transform Future array into data.frame
            proj.future <- as.data.frame(proj.future)
            names(proj.future) <- unlist(lapply(strsplit(names(proj.future), ".", fixed = TRUE), function(x)
            {
              x.rev <- rev(x) ## reverse the order of splitted vector to have algo at the end
              data.set.id <- x.rev[1]
              cross.valid.id <- x.rev[2]
              algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
              model.id <- paste(data.set.id, cross.valid.id, algo.id, sep = "_")
              return(model.id)
            }))
            
            return(BIOMOD_RangeSize(proj.current, proj.future))
          })


## BIOMOD_RangeSize SpatRaster-SpatRaster Method ----------------------

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'SpatRaster', proj.future = 'SpatRaster'),
          function(proj.current, proj.future)
          {
            .bm_cat("Do Range Size Computation")
            args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future)
            
            names.res = c("Loss", "Stable0", "Stable1", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            if(nlyr(proj.current) > 1){
              CBS <- matrix(ncol = 10, nrow = nlyr(proj.current),
                          dimnames = list(names(proj.current), names.res))
            } else {
              CBS <- matrix(ncol = 10, nrow = nlyr(proj.future),
                            dimnames = list(names(proj.future), names.res))
            }

            sp.rast <- rast()
            Cur <- proj.current[[1]]
            for (i in seq_len(nlyr(proj.future))) {
              ## DiffByPixel
              if(nlyr(proj.current) > 1){
                Cur <- proj.current[[i]]
              }
              Fut <- proj.future[[i]]
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
            names(sp.rast) <- rownames(CBS)
            
            .bm_cat("Done")
            return(list(Compt.By.Models = CBS, Diff.By.Pixel = sp.rast))
          })


## BIOMOD_RangeSize Raster-Raster Method ----------------------
##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize',
          signature(proj.current = 'Raster', proj.future = 'Raster'),
          function(proj.current, proj.future) {
            BIOMOD_RangeSize(proj.current = rast(proj.current), proj.future = rast(proj.future))
          })

# Argument Check ---------------------------------------------------------------

.BIOMOD_RangeSize.check.args <- function(proj.current, proj.future) {
  
  ## dimensions checking ------------------------
  if(inherits(proj.current, "data.frame")){
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
  if(inherits(proj.current, "data.frame")){
    n_current <- ncol(proj.current)
    n_future <- ncol(proj.future)
  } else { # SpatRaster case
    n_current <- nlyr(proj.current)
    n_future <- nlyr(proj.future)
  }
  
  if (n_current == 1) {
    if(n_future == 1) {
      cat("\n Comparing 'proj.current' and 'proj.future'. ")
    } else if (n_future > 1) {
      cat("\n Comparing 'proj.current' with the ", n_future, " projections in 'proj.future'. ")
    } else {
      stop("'proj.future' require at least one layer or one column.")
    }
  } else if (n_current > 1) {
    if(n_future == n_current) {
      cat("\n Each projection in 'proj.current' will be compared once with its corresponding projection in 'proj.future'. ")
    } else {
      stop("When proj.current' have more than 1 projection, 'proj.future' must have the same number of projection.")
    }
  } else {
    stop("'proj.current' require at least one layer or one column.")
  }
  
  
  ## checking 0/1 ----------------------------
  if(inherits(proj.current, "data.frame")){
    test_binary_current <- all(na.omit(unlist(proj.current)) %in% c(0,1))
    test_binary_future <- all(na.omit(unlist(proj.future)) %in% c(0,1))
  } else { # SpatRaster case
    test_binary_current <- all(na.omit(values(proj.current)) %in% c(0,1))
    test_binary_future <- all(na.omit(values(proj.future)) %in% c(0,1))
  }
  
  if(!test_binary_current | !test_binary_future){
    stop("'proj.current' and 'proj.future' must have only values among 0, 1 or NA.")
  }
  
}


# Additionnal tools -------------------------------------------------------

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

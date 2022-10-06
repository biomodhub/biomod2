###################################################################################################
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
##' or \code{\link[raster:stack]{RasterStack}} object containing the initial binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param proj.future an \code{array}, \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the final binary projection(s) 
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
##' \code{Diff.By.Pixel} object is obtained by applying the simple following formula :
##' \deqn{proj.future - 2 * proj.current}
##' 
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
##' myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl <- raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExpl <- raster::stack(raster::crop(myExpl, myExtent))
##' }
##' 
##' # ---------------------------------------------------------------
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
##' # ---------------------------------------------------------------
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0('external/bioclim/future/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExplFuture = raster::stack(system.file(myFiles, package = 'biomod2'))
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExplFuture <- raster::stack(raster::crop(myExplFuture, myExtent))
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
##' proj.current <- raster::stack("GuloGulo/proj_Current/proj_Current_GuloGulo_TSSbin.grd")
##' proj.future <- raster::stack("GuloGulo/proj_Future/proj_Future_GuloGulo_TSSbin.grd")
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
##' @importFrom raster stack addLayer
##' 
##' @export
##' 
##' 
###################################################################################################


setGeneric("BIOMOD_RangeSize",
           def = function(proj.current, proj.future) {
             standardGeneric("BIOMOD_RangeSize")
           })


###################################################################################################

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'data.frame', proj.future = 'data.frame'),
          function(proj.current, proj.future)
          {
            .bm_cat("Do Range Size Computation")
            args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future)

            CompteurSp <- function(Data, Value)
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
            
            Diff.By.Pixel <- as.data.frame(proj.future - 2 * proj.current)
            Compt.By.Models <- as.data.frame(CompteurSp(Diff.By.Pixel, c(-2, 0, -1, 1)))
            Compt.By.Models[, 8] <- Compt.By.Models[, 1] + Compt.By.Models[, 3]
            Compt.By.Models[, 9] <- Compt.By.Models[, 3]
            Compt.By.Models[, 10] <- Compt.By.Models[, 4] + Compt.By.Models[, 3]
            
            Compt.By.Models[, 5] <- (100 * Compt.By.Models[, 1]) / Compt.By.Models[, 8]
            Compt.By.Models[, 6] <- (100 * Compt.By.Models[, 4]) / Compt.By.Models[, 8]
            Compt.By.Models[, 7] <- Compt.By.Models[, 6] - Compt.By.Models[, 5]

            dimnames(Compt.By.Models) <- list(colnames(proj.current), c("Loss", "Stable0", "Stable1", "Gain"
                                                                       , "PercLoss", "PercGain", "SpeciesRangeChange"
                                                                       , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            
            Output <- list(Compt.By.Models = Compt.By.Models, Diff.By.Pixel = Diff.By.Pixel)
            .bm_cat("Done")
            invisible(Output)
          })

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


###################################################################################################

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'RasterStack', proj.future = 'RasterStack'),
          function(proj.current, proj.future)
          {
            .bm_cat("Do Range Size Computation")
            args <- .BIOMOD_RangeSize.check.args(proj.current, proj.future)
            
            names.res = c("Loss", "Stable0", "Stable1", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            CBS <- matrix(ncol = 10, nrow = length(proj.current@layers),
                          dimnames = list(names(proj.current), names.res))

            sp.stack <- stack()
            for (i in 1:length(proj.current@layers)) {
              
              ## DiffByPixel
              Cur <- proj.current@layers[[i]]
              Fut <- proj.future@layers[[i]]
              Ras <- Fut - (Cur + Cur)
              sp.stack <- addLayer(sp.stack, Ras)
              
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
            names(sp.stack) <- rownames(CBS)
            
            .bm_cat("Done")
            return(list(Compt.By.Models = CBS, Diff.By.Pixel = sp.stack))
          })

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize', signature(proj.current = 'RasterLayer', proj.future = 'RasterStack'),
          function(proj.current, proj.future)
          {
            .bm_cat("Do Range Size Computation")
            
            names.res = c("Loss", "Stable0", "Stable1", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            CBS <- matrix(ncol = 10, nrow = length(proj.future@layers),
                          dimnames = list(names(proj.future), names.res))
            
            sp.stack <- stack()
            for (i in 1:length(proj.future@layers)) {
              
              ## DiffByPixel
              Cur <- proj.current
              Fut <- proj.future@layers[[i]]
              Ras <- Fut - (Cur + Cur)
              sp.stack <- addLayer(sp.stack, Ras)
              
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
            names(sp.stack) <- rownames(CBS)
            
            .bm_cat("Done")
            return(list(Compt.By.Models = CBS, Diff.By.Pixel = sp.stack))
          })

##'
##' @rdname BIOMOD_RangeSize
##' @export
##'

setMethod('BIOMOD_RangeSize',
          signature(proj.current = 'RasterLayer', proj.future = 'RasterLayer'),
          function(proj.current, proj.future) {
            BIOMOD_RangeSize(proj.current = proj.current, proj.future = stack(proj.future))
          })



###################################################################################################

.BIOMOD_RangeSize.check.args <- function(proj.current, proj.future)
{
  # dimensions checking
  if (sum(!(dim(proj.current) == dim(proj.future)) > 0)) {
    stop(paste0("'proj.current' and 'proj.future' do not have the same dimensions "
                , "('proj.current' must have either 1 projection or the same number of projections as 'proj.future')"))
  }
}

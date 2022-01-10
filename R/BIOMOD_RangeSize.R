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
##' @param CurrentProj an \code{array}, \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the initial binary projection(s) 
##' of the (ensemble) species distribution model(s)
##' @param FutureProj an \code{array}, \code{data.frame}, \code{\link[raster:stack]{RasterLayer}} 
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
##'   \item{Diff.By.Pixel}{an object in the same form than the input data (\code{CurrentProj} and 
##'   \code{FutureProj}) and containing a value for each point/pixel of each comparison among :
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
##' \deqn{FutureProj - 2 * CurrentProj}
##' 
##' 
##' @keywords species range change, projections, gain, loss
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_ModelingOptions}}, 
##' \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_EnsembleForecasting}}, \code{\link{level.plot}}
##' 
##'   
##' @examples
##' 
##' # species occurrences
##' myFile <- system.file("external/species/mammals_table.csv", package="biomod2")
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[, c("X_WGS84", "Y_WGS84")]
##' 
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0("external/bioclim/current/bio", c(3, 4, 7, 11, 12), ".grd")
##' myExpl = raster::stack(system.file(myFiles[1], package = "biomod2"),
##'                        system.file(myFiles[2], package = "biomod2"),
##'                        system.file(myFiles[3], package = "biomod2"),
##'                        system.file(myFiles[4], package = "biomod2"),
##'                        system.file(myFiles[5], package = "biomod2"))
##' 
##' # 1. Formating Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' myBiomodOption <- BIOMOD_ModelingOptions()
##' 
##' # 3. Doing Modelisation
##' myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
##'                                     models = c('SRE','RF'),
##'                                     models.options = myBiomodOption,
##'                                     NbRunEval = 2,
##'                                     DataSplit = 80,
##'                                     VarImport = 0,
##'                                     models.eval.meth = c('TSS','ROC'),
##'                                     do.full.models = FALSE,
##'                                     modeling.id = 'test')
##' 
##' # 4.1 Projecting on current environmental conditions
##' myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
##'                                         new.env = myExpl,
##'                                         proj.name = 'current',
##'                                         chosen.models = 'all',
##'                                         binary.meth = 'TSS',
##'                                         compress = FALSE,
##'                                         build.clamping.mask = FALSE)
##' 
##' # 4.2 Projecting on future environmental conditions
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0("external/bioclim/future/bio", c(3, 4, 7, 11, 12), ".grd")
##' myExplFuture = raster::stack(system.file(myFiles[1], package = "biomod2"),
##'                              system.file(myFiles[2], package = "biomod2"),
##'                              system.file(myFiles[3], package = "biomod2"),
##'                              system.file(myFiles[4], package = "biomod2"),
##'                              system.file(myFiles[5], package = "biomod2"))
##' 
##' myBiomodProjectionFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
##'                                               new.env = myExplFuture,
##'                                               proj.name = 'future',
##'                                               chosen.models = 'all',
##'                                               binary.meth = 'TSS',
##'                                               compress = FALSE,
##'                                               build.clamping.mask = TRUE)
##' 
##' # 5. Detecting predicted species range change
##' 
##' # load binary projections
##' CurrentProj <- raster::stack("GuloGulo/proj_current/proj_current_GuloGulo_TSSbin.grd")
##' futurePred <- raster::stack("GuloGulo/proj_future/proj_future_GuloGulo_TSSbin.grd")
##' 
##' myBiomodRangeSize <- BIOMOD_RangeSize(CurrentProj = CurrentProj, FutureProj = futurePred)
##' 
##' # print summary and visualize changes
##' myBiomodRangeSize$Compt.By.Models
##' plot(myBiomodRangeSize$Diff.By.Pixel)
##' 
##' 
##' @importFrom raster stack addLayer
##' 
##' @export
##' 
##' 
###################################################################################################


setGeneric("BIOMOD_RangeSize",
           def = function(CurrentProj, FutureProj, ...) {
             standardGeneric("BIOMOD_RangeSize")
           })


###################################################################################################

setMethod('BIOMOD_RangeSize', signature(CurrentProj = 'data.frame', FutureProj = 'data.frame'),
          function(CurrentProj, FutureProj)
          {
            args <- .BIOMOD_RangeSize.check.args(CurrentProj, FutureProj)

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
            
            Diff.By.Pixel <- as.data.frame(FutureProj - 2 * CurrentProj)
            Compt.By.Models <- as.data.frame(CompteurSp(Diff.By.Pixel, c(-2, 0, -1, 1)))
            Compt.By.Models[, 8] <- Compt.By.Models[, 1] + Compt.By.Models[, 3]
            Compt.By.Models[, 9] <- Compt.By.Models[, 3]
            Compt.By.Models[, 10] <- Compt.By.Models[, 4] + Compt.By.Models[, 3]
            
            Compt.By.Models[, 5] <- (100 * Compt.By.Models[, 1]) / Compt.By.Models[, 8]
            Compt.By.Models[, 6] <- (100 * Compt.By.Models[, 4]) / Compt.By.Models[, 8]
            Compt.By.Models[, 7] <- Compt.By.Models[, 6] - Compt.By.Models[, 5]

            dimnames(Compt.By.Models) <- list(colnames(CurrentProj), c("Loss", "Stable0", "Stable1", "Gain"
                                                                       , "PercLoss", "PercGain", "SpeciesRangeChange"
                                                                       , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp"))
            
            Output <- list(Compt.By.Models = Compt.By.Models, Diff.By.Pixel = Diff.By.Pixel)
            invisible(Output)
          })


setMethod('BIOMOD_RangeSize', signature(CurrentProj = 'array', FutureProj = 'array'),
          function(CurrentProj, FutureProj)
          {
            
            ## Transform Current array into data.frame
            CurrentProj <- as.data.frame(CurrentProj)
            names(CurrentProj) <- unlist(lapply(strsplit(names(CurrentProj), ".", fixed = TRUE), function(x)
            { return(paste(x[3], x[2], x[1], sep = "_")) }))
            names(CurrentProj) <- unlist(lapply(strsplit(names(CurrentProj), ".", fixed = TRUE), function(x)
            {
              x.rev <- rev(x) ## reverse the order of splitted vector to have algo at the end
              data.set.id <- x.rev[1]
              cross.valid.id <- x.rev[2]
              algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
              model.id <- paste(data.set.id, cross.valid.id, algo.id, sep = "_")
              return(model.id)
            }))

            ## Transform Future array into data.frame
            FutureProj <- as.data.frame(FutureProj)
            names(FutureProj) <- unlist(lapply(strsplit(names(FutureProj), ".", fixed = TRUE), function(x)
            {
              x.rev <- rev(x) ## reverse the order of splitted vector to have algo at the end
              data.set.id <- x.rev[1]
              cross.valid.id <- x.rev[2]
              algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
              model.id <- paste(data.set.id, cross.valid.id, algo.id, sep = "_")
              return(model.id)
            }))
            
            return(BIOMOD_RangeSize(CurrentProj, FutureProj))
          })


###################################################################################################

setMethod('BIOMOD_RangeSize', signature(CurrentProj = 'RasterStack', FutureProj = 'RasterStack'),
          function(CurrentProj, FutureProj)
          {
            args <- .BIOMOD_RangeSize.check.args(CurrentProj, FutureProj)
            
            names.res = c("Loss", "Stable0", "Stable1", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            CBS <- matrix(ncol = 10, nrow = length(CurrentProj@layers),
                          dimnames = list(names(CurrentProj), names.res))

            sp.stack <- stack()
            for (i in 1:length(CurrentProj@layers)) {
              
              ## DiffByPixel
              Cur <- CurrentProj@layers[[i]]
              Fut <- FutureProj@layers[[i]]
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
            return(list(Compt.By.Models = CBS, Diff.By.Pixel = sp.stack))
          })


setMethod('BIOMOD_RangeSize', signature(CurrentProj = 'RasterLayer', FutureProj = 'RasterStack'),
          function(CurrentProj, FutureProj)
          {
            names.res = c("Loss", "Stable0", "Stable1", "Gain"
                          , "PercLoss", "PercGain", "SpeciesRangeChange"
                          , "CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
            CBS <- matrix(ncol = 10, nrow = length(FutureProj@layers),
                          dimnames = list(names(FutureProj), names.res))
            
            sp.stack <- stack()
            for (i in 1:length(FutureProj@layers)) {
              
              ## DiffByPixel
              Cur <- CurrentProj
              Fut <- FutureProj@layers[[i]]
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
            return(list(Compt.By.Models = CBS, Diff.By.Pixel = sp.stack))
          })


setMethod('BIOMOD_RangeSize',
          signature(CurrentProj = 'RasterLayer', FutureProj = 'RasterLayer'),
          function(CurrentProj, FutureProj) {
            BIOMOD_RangeSize(CurrentProj = CurrentProj, FutureProj = stack(FutureProj))
          })



###################################################################################################

.BIOMOD_RangeSize.check.args <- function(CurrentProj, FutureProj)
{
  # dimensions checking
  if (sum(!(dim(CurrentProj) == dim(FutureProj)) > 0)) {
    stop("CurrentProj & FutureProj dimensions mismatched!")
  }
}

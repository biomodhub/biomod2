\name{BIOMOD_RangeSize}
\alias{BIOMOD_RangeSize}
\alias{BIOMOD_RangeSize-methods}
\alias{BIOMOD_RangeSize,data.frame,data.frame-method}
\alias{BIOMOD_RangeSize,array,array-method}
\alias{BIOMOD_RangeSize,RasterStack,RasterStack-method}
\alias{BIOMOD_RangeSize,RasterLayer,RasterLayer-method}
\alias{BIOMOD_RangeSize,RasterLayer,RasterStack-method}



\title{ Analysis of the range size changes }
\description{
  This function allows to estimate the proportion and relative number of pixels (or habitat) lost, gained
and stable for the time slice considered in species-climate modeling under future scenarios.
}
\usage{

\S4method{BIOMOD_RangeSize}{data.frame,data.frame}(CurrentPred, FutureProj,  SpChange.Save=NULL)

\S4method{BIOMOD_RangeSize}{array,array}(CurrentPred, FutureProj,  SpChange.Save=NULL)

\S4method{BIOMOD_RangeSize}{RasterStack,RasterStack}(CurrentPred, FutureProj,  SpChange.Save=NULL)

\S4method{BIOMOD_RangeSize}{RasterStack,RasterStack}(CurrentPred, FutureProj,  SpChange.Save=NULL)

\S4method{BIOMOD_RangeSize}{RasterLayer,RasterLayer}(CurrentPred, FutureProj,  SpChange.Save=NULL)

}

\arguments{
  \item{CurrentPred}{ a data.frame of n (number of models) columns OR an array formatted as an output of \code{\link[biomod2]{BIOMOD_Projection}} function OR a \code{RasterStack} (a layer by model) giving the current state of the species in binary  }
  \item{FutureProj}{ a data.frame of n (number of models) columns OR an array formatted as an output of \code{\link[biomod2]{BIOMOD_Projection}} function OR a \code{RasterStack} (a layer by model) giving the future projections of the species in binary according to a future scenario}
  \item{SpChange.Save}{ the name given to the new object storing the results }
}
\details{
  Note that this function is only relevant if you make projections on the same area with the same resolution and
  for a different time slice (past or future) than for the current data.
}
\value{
  A list of two items is created: Compt.By.Species and Diff.By.Pixel ; Compt.By.Species is the summary of this function; Diff.By.Pixel is in the same form than your input data.

  \item{Compt.By.Species}{ stores the summary of range change for each species (sorted by rows).
  The first four columns are absolute values whereas the next 3 ones are relative values:}
  \item{    Disa}{ represents the number of pixels predicted to be lost by the given species.}
  \item{    Stable0}{ is the number of pixels which are not currently occupied by the given species and not predicted to be.}
  \item{    Stable1}{ represents the number of pixels currently occupied by the given species, and predicted to remain occupied into the future.}
  \item{    Gain}{ represents the number of pixels which are currently not occupied by the given species but predicted to be into the future.}
  \item{    PercLoss}{corresponds to the percentage of currently occupied sites to be lost (Disa/(Disa+Stable1)}
  \item{    PercGain}{corresponds to the percentage of new sites considering the species' current distribution size (Gain/(Disa+Stable1).
  For example, if the there are 30 sites currently occupied and 15 new sites are projected to be occupied in future, it makes PercGain=+50(\%).}
  \item{    SpeciesRangeChange}{ it is the overall projection outcome, equal to PercGain-PercLoss. It does not assess for any migration shifts as it
  strictly compares the range sizes between current and future states.}
  \item{    CurrentRangeSize}{ represents the modeled current range size (number of pixels occupied) of the given species.}
  \item{    FutureRangeSize0Disp}{ represents the future modeled range size assuming no migration of the given species.}
  \item{    FutureRangeSize1Disp}{ represents the future modeled range size assuming migration of the given species (depending on the datasets given in input, if Migration has been used or not).}
  
  \item{Diff.By.Pixel}{ the summary of range change for each species (sorted by columns and with the pixel
  in rows). For each species, a pixel could have four different values :
  -2 if the given pixel is predicted to be lost by the species. 
  -1 if the given pixel is predicted to be stable for the species.
   0 is the given pixel was not occupied, and will not be in the future.
   1 if the given pixel was not occupied, and is predicted to be into the future.
  This table could be easily plotted into GIS software in order to represent the pattern of change for the selected species (or even with the level.plot() function).
  }

}


\author{ Wilfried Thuiller, Damien Georges, Bruno Lafourcade }

\examples{
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)


# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('CTA','RF'), 
                                       models.options = myBiomodOption, 
                                       models.eval.meth ='TSS',
                                       rescal.all.models=FALSE)


# 4.1 Projection on current environemental conditions

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = 'all',
                                          binary.meth = 'TSS',
                                          compress = FALSE,
                                          build.clamping.mask = FALSE)

# 4.2 Projection on future environemental conditions

myExplFuture = raster::stack(system.file("external/bioclim/future/bio3.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio4.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio7.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio11.grd",package="biomod2"),
                     system.file("external/bioclim/future/bio12.grd",package="biomod2"))

myBiomodProjectionFuture <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                              new.env = myExplFuture,
                                              proj.name = 'future',
                                              selected.models = 'all',
                                              binary.meth = 'TSS',
                                              compress = FALSE,
                                              build.clamping.mask = TRUE)

# 5. Detect where our species occurances state is forecasted to change

# load binary projections
# here is rasters objects ('.grd')
currentPred <- raster::stack("GuloGulo/proj_current/proj_current_GuloGulo_TSSbin.grd")
futurePred <- raster::stack("GuloGulo/proj_future/proj_future_GuloGulo_TSSbin.grd")


# call the Range size function
myBiomodRangeSize <- BIOMOD_RangeSize(
  CurrentPred=currentPred,
  FutureProj=futurePred)

# see the results
myBiomodRangeSize$Compt.By.Models
plot(myBiomodRangeSize$Diff.By.Pixel)

}

\keyword{ IO }

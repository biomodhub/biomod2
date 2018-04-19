\name{BIOMOD_Projection}
\alias{BIOMOD_Projection}

\title{ Project the calibrated models within \pkg{biomod2} into new space or time}
\description{
  For all the models currently implemented, \pkg{biomod2} is able to project potential distributions of species in other areas, other resolutions or other time scales.
}
\usage{
BIOMOD_Projection(modeling.output,
                  new.env,
                  proj.name,
                  xy.new.env = NULL,
                  selected.models = 'all',
                  binary.meth = NULL,
                  filtered.meth = NULL,
                  compress = TRUE,
                  build.clamping.mask = TRUE,
                  ...)
}

\arguments{

\item{modeling.output}{  \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} object produced by a \code{\link[biomod2]{BIOMOD_Modeling}} run}

\item{new.env}{ A set of explanatory variables onto which models will be projected . It could be a \code{data.frame}, a \code{matrix}, or a \code{rasterStack} object. Make sure the column names (\code{data.frame} or \code{matrix}) or layer Names (\code{rasterStack}) perfectly match with the names of variables used to build the models in the previous steps.}

\item{proj.name}{a character defining the projection name (a new folder will be created with this name)}

\item{xy.new.env}{optional coordinates of new.env data. Ignored if new.env is a \code{rasterStack}}

\item{selected.models}{'all' when all models have to be used to render projections or a subset vector of modelling.output models computed (accessing with the slot @models.computed of your \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} object)}

\item{binary.meth}{a vector of a subset of models evaluation method computed before (see \code{\link[biomod2]{BIOMOD_Modeling}}). If NULL then no binary transformation computed, else the given binary techniques will be used to transform the projection into 0/1 data. }

\item{filtered.meth}{a vector of a subset of models evaluation method computed before (see \code{\link[biomod2]{BIOMOD_Modeling}}). if NULL then no filtering transformation computed, else the given binary techniques will be used to transform the projection by stetting to 0 the probability values below a specific threshold.}

\item{compress}{Boolean or character, the compression format of objects stored on your hard drive. May be one of \sQuote{TRUE}, \sQuote{FALSE}, \sQuote{xz} or \sQuote{gzip} (see \code{\link{save}})}

\item{build.clamping.mask}{ if TRUE, a clamping mask will be saved on hard drive different (see details)}

\item{...}{ Additional arguments (see details section) }
}


\details{
Projections are done for all selected models, that means (by default) for all evaluation run, and pseudo absences selections if applicable. This projections may be used later to compute ensemble forecasting.

If \code{build.clamping.mask} is set to \code{TRUE} a file (same type than \code{new.env} arg) will be saved in your projection folder. This mask will identifies locations where predictions are uncertain because the values of the variables are outside the range used for calibrating the models. The \sQuote{build.clamping.mask} values correspond to the number of variables that are out of their calibrating/training range. (see vignette for more details)

\ldots may be :
\itemize{
  \item{\code{silent}:}{logical, if TRUE, console outputs are turned off}
  \item{\code{do.stack}: }{logical, if TRUE, attempt to save all projections in a unique object {i.e \code{RasterStack}}. If FALSE or if objects are too heavy to be load all together in memory, projections will be stored into separated files.}
  \item{\code{keep.in.memory}:}{logical, if FALSE only the link pointing to a hard drive copy of projections are stored in output object. That can be useful to prevent memory issues.}
  \item{\code{output.format}:}{whether \sQuote{.RData}, \sQuote{.grd} or \sQuote{.img} defining projections saving format (on hard drive). If \code{new.env} argument is under table format (\code{data.frame} or \code{matrix}), the only choice you have is \sQuote{.RData} }
  \item{\code{omit.na}:}{logical, if TRUE (default), all not fully referenced environmental points will get a NA as prediction. If FALSE, models that can produce predictions with incomplete data will return a prediction value for this points. }
  \item{\code{on_0_1000}:}{logical, if TRUE (default), 0 - 1 probabilities are converted into a 0 - 1000 integer scale. This implies a lot of memory saving. User that want to comeback on a 0 - 1 scale latter will just have to divide all projections by 1000  }
}

}

\value{
Returns the projections for all selected model (\code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} object), and stored in the hard drive on the specific directory names by the name of the projection. The data is a 4-dimensions array (see ...) if new.env is a \code{matrix} or a \code{data.frame}. It is a rasterStack if new.env is a \code{rasterStack} and or several rasterLayers if the \code{rasterStack} is too large.

A new folder is also created on your hard drive. This folder contains the created projection object (basic one and binary and filtered ones if selected). The object are loaded with the \code{\link{load}} function. The loaded object can be then plotted and analyzed.

}

\author{
Wilfried Thuiller, Damien Georges
}

\seealso{
\code{\link{BIOMOD_Modeling}, \link{BIOMOD_FormatingData}},  \code{\link{BIOMOD_ModelingOptions}}
}

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
                                       models = c('SRE','RF'),
                                       models.options = myBiomodOption,
                                       NbRunEval=1,
                                       DataSplit=70,
                                       models.eval.meth = c('TSS'),
                                       do.full.models = FALSE)


# 4.1 Projection on current environemental conditions

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = myExpl,
                                          proj.name = 'current',
                                          selected.models = 'all',
                                          binary.meth = 'TSS',
                                          compress = FALSE,
                                          build.clamping.mask = FALSE)

\dontrun{
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

# print summary and plot projections
myBiomodProjectionFuture
plot(myBiomodProjectionFuture)
}

}
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }


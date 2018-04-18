\name{BIOMOD_LoadModels}
\alias{BIOMOD_LoadModels}

\title{ Load models built within BIOMOD_Modeling function}

\description{
  This function was implemented to help \pkg{biomod2} users to load individual models built during the \code{\link[biomod2]{BIOMOD_Modeling}} step.
}
\usage{
  BIOMOD_LoadModels(bm.out, ... )
}

\arguments{
  \item{bm.out}{ a \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} returned by \code{\link[biomod2]{BIOMOD_Modeling}} }
  \item{...}{ additional arguments (see details)  }
}

\details{
  This function is particularly useful when you plane to make some response plot analyses. It will induce models, built at \code{\link[biomod2]{BIOMOD_Modeling}} step, loading in your working space.

  If you run this function referencing only \code{bm.out} argument, all models built will be loaded. However, you can make a models subselection using the additional arguments (see below).

  \bold{Additional arguments (\ldots{}) : }
  All the following arguments are optional.

  \itemize{

    \item{\bold{models:}}{ a character vector defining the names of models (e.g c(\sQuote{GLM}, \sQuote{GAM}, \sQuote{RF}) ) you want to load (models subselection) }
    \item{\bold{run.eval:}}{ a character vector defining the names of evaluation run (e.g c(\sQuote{RUN1}, \sQuote{Full}) ) you want to load (repetition subselection) }
    \item{\bold{data.set:}}{ a character vector defining the names of data.set (e.g c(\sQuote{PA1}, \sQuote{PA2})) you want to load (pseudo absences subselection) }
    \item{\bold{path:}}{ the path to file where the species folder is. To be filled if species folder is different from your working directory) }

  }

}

\note{
  \code{SRE} models are not supported yet. They will be automatically excluded even if they are selected.
}

\value{
  A character vector filled with the loaded models names.
}


\author{ Damien Georges }

\seealso{ \code{\link{BIOMOD_Modeling}}}
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
                                       models = c('RF'),
                                       models.options = myBiomodOption,
                                       NbRunEval=2,
                                       DataSplit=70,
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       do.full.models = FALSE)

# 4. Loading some models built

myLoadedModels <- BIOMOD_LoadModels(myBiomodModelOut, models='RF')

myLoadedModels

}

\keyword{ models }
\keyword{ datasets }

\name{BIOMOD_EnsembleForecasting}
\alias{BIOMOD_EnsembleForecasting}

\title{ Ensemble projections of species over space and time }
\description{
This function use projections of \sQuote{individual models} and ensemble models from \code{\link[biomod2]{BIOMOD_EnsembleModeling}} to build an ensemble of species' projections over space and time.
}
\usage{
BIOMOD_EnsembleForecasting( EM.output,
                            projection.output = NULL,
                            new.env = NULL,
                            xy.new.env = NULL,
                            selected.models = 'all',
                            proj.name = NULL,
                            binary.meth = NULL,
                            filtered.meth = NULL,
                            compress = TRUE,
                            ...)
}


\arguments{
  \item{EM.output}{ a \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} returned by \code{\link[biomod2]{BIOMOD_EnsembleModeling}}}
  \item{projection.output}{a \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} returned by \code{\link[biomod2]{BIOMOD_Projection}} }
  \item{new.env}{ a \code{RasterStack} or a \code{data.frame} object containing explanatory data for the studied area. Needed only if \code{projection.output} is \code{NULL}. Prefer to use \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} if you have ever done the calculations.}
  \item{xy.new.env}{the matching coordinates of \code{new.env} if \code{new.env} is defined and if it is a \code{data.frame}}
  \item{selected.models}{ if not 'all', a character vector containing a subset of ensemble-models you want make projection}
  \item{proj.name}{the projection name (results will be saved within proj_proj.name directory). Only needed if \code{projection.output} is \code{NULL} }
  \item{binary.meth}{vector specifying the names of evaluation metrics and associated thresholds to transform the probabilities of presence into presence and absence (binary transformation).  }
  \item{filtered.meth}{vector specifying the names of evaluation metrics and associated thresholds to transform into 0 the probabilities of presence lower than the thresholds.}
  \item{compress}{boolean or character, the compression format of objects stored on your hard drive. May be one of \sQuote{TRUE}, \sQuote{FALSE}, \sQuote{xz} or \sQuote{gzip} (see \code{\link{save}})}
  \item{\ldots}{ further arguments (see details)}
}

\details{
This function requires to have successfully run \pkg{biomod2} modeling, ensemble-modeling and projection steps. Ensemble projections will be created in respect to \code{projection.output} projections, which are combined following \code{EM.output} ensemble-modeling rules.

The \sQuote{total.consensus} projection is basically the mean of all projections (for having only one output).  

\ldots may be :
\itemize{
  \item{\code{on_0_1000}:}{logical, if TRUE (default), 0 - 1 probabilities are converted into a 0 - 1000 integer scale. This implies a lot of memory saving. User that want to comeback on a 0 - 1 scale latter will just have to divide all projections by 1000  }
}

} % end of detail 

\value{
Nothing returned but specific \sQuote{projection files} are saved on the hard drive projection folder. This files are either an \code{array} or a \code{RasterStack} depending the original projections data type.
Load these created files to plot and work with them.

}% end of value


\author{ 
Wilfried Thuiller, Damien Georges, Robin Engler
}

\seealso{
\code{\link[biomod2]{BIOMOD_EnsembleModeling}}, \code{\link[biomod2]{BIOMOD_Projection}}
}

\examples{
# 0. Load data & Selecting Data
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

# 3. Running the models 
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                       models = c('RF'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=2, 
                                       DataSplit=60, 
                                       Yweights=NULL, 
                                       VarImport=0, 
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)
                                       
# 4. Creating the ensemble models 
myBiomodEM <- BIOMOD_EnsembleModeling( 
                 modeling.output = myBiomodModelOut,
                 chosen.models = grep('_RF', get_built_models(myBiomodModelOut), 
                                  value=TRUE),
                 em.by = 'algo',
                 eval.metric = c('TSS'),
                 eval.metric.quality.threshold = c(0.7),
                 prob.mean = TRUE,
                 prob.cv = FALSE,
                 prob.ci = FALSE,
                 prob.ci.alpha = 0.05,
                 prob.median = FALSE,
                 committee.averaging = FALSE,
                 prob.mean.weight = FALSE,
                 prob.mean.weight.decay = 'proportional' )
                                       
# 5. Individual models projections on current environmental conditions
myBiomodProjection <- BIOMOD_Projection(
                        modeling.output = myBiomodModelOut,
                        new.env = myExpl,
                        proj.name = 'current',
                        selected.models = grep('_RF', get_built_models(
                                                myBiomodModelOut), value=TRUE),
                        compress = FALSE,
                        build.clamping.mask = FALSE)
                          

# 4. Creating the ensemble projections
BIOMOD_EnsembleForecasting( projection.output = myBiomodProjection,
                            EM.output = myBiomodEM)
}

\keyword{ models }

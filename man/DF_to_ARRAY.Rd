\name{DF_to_ARRAY}
\alias{DF_to_ARRAY}
\alias{LIST_to_ARRAY}

\title{Convert a biomod2 data.frame (or list) into array }

\description{ DF_to_ARRAY  and LIST_to_ARRAY are \pkg{biomod2} internal functions that can be useful to help users to transform a standard \pkg{biomod2} \code{data.frame} or \code{list} output into the standard \code{array} one.
}

\usage{
DF_to_ARRAY(df)
LIST_to_ARRAY(ll)
}



\arguments{
  \item{df}{  a standard \pkg{biomod2} \code{data.frame} output }
  \item{ll}{  a standard \pkg{biomod2} \code{list} output (e.g \code{evaluate} function output )}
}


\details{
This functions can be useful when you want to re-convert \pkg{biomod2} \code{data.frame} or \code{list} (e.g. projections, evaluations) into their initial \code{array} format (see \code{\link[biomod2]{BIOMOD_Projection}} for further details) 
}


\value{
A standard \pkg{biomod2} \code{array} output.
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_Projection}} }
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
                             
# Keep only points where we have info                       
myExpl <- raster::extract(myExpl, myRespXY)

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
                                       Yweights=NULL, 
                                       VarImport=0, 
                                       models.eval.meth = c('ROC'),
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)
                                       
                                       
# 4 Projection on current environemental conditions

myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                          new.env = data.frame(myExpl),
                                          proj.name = 'current',
                                          selected.models = 'all')
                                          


# 5. Get projection under data.frame format
myProjDF <- get_predictions(myBiomodProjection, as.data.frame=TRUE)
class(myProjDF)
dim(myProjDF)
dimnames(myProjDF)

# 6. Transform data.frame into array
myProjArray <- DF_to_ARRAY(myProjDF)
class(myProjArray)
dim(myProjArray)
dimnames(myProjArray)

}

\keyword{models}
\keyword{formula}
\keyword{options}

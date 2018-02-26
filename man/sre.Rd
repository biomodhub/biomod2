\name{sre}
\alias{sre}

\title{ Surface Range Envelope }
\description{
  Run a rectilinear surface range envelop (equivalent to BIOCLIM) using the extreme percentiles as recommended by Nix or Busby.
  The SRE performs a simple analysis of within which range of each variable the data is recorded and renders predictions.
}
\usage{
sre(Response = NULL, 
    Explanatory = NULL, 
    NewData = NULL, 
    Quant=0.025, 
    return_extremcond = FALSE)
}

\arguments{
  \item{Response}{ a vector (resp. a matrix/data.frame, a SpatialPointsDataFrame or a RasterLayer ) giving your species' presences and absences data }
  \item{Explanatory}{ a matrix (resp. a data.frame, a SpatialPointsDataFrame or a RasterStack ) containing the environmental variables for the sites given in Response. It must have as many rows as there are elements in Response. }
  \item{NewData}{ The data for which you want to render predictions with the sre. It must be a matrix (resp. a data.frame, a SpatialPointsDataFrame or a RasterStack ) of the same type as the one given in Explanatory and with precisely the same variable names. }
  \item{Quant}{ the value defines the most extreme values for each variable not to be taken into account
   for determining the tolerance boundaries for the considered species.  }
   \item{return_extremcond}{ boolean, if TRUE a matrix containing extrem conditions supported is returned}
}
\details{

  The more variables you put in, the more restrictive your model will be (if non-colinear variables).

  This method is very much influenced by the data input, and more specifically by the extremes.
  Where a linear model can discriminate the extreme values from the main tendency, the SRE considers it
  equal as any other data point which leads to notable changes in predictions.
  Note that, as a consequence of its functioning, the predictions are directly given in binary, a site being
  either potentially suitable for all the variables, either out of bounds for at least one variable and therefore
  considered unsuitable. 
  
  The quants argument determines the threshold at which the data will be taken into account for calibration : the default
  of 0.05 induces that the 5\% most extreme values will be avoided for each variable on each side of its distribution along
  the gradient. So it in fact takes 5\% away at each end of the variables distribution, giving a total of 10\% of data not considered.
}
\value{
  A vector (resp. a RasterLayer ) of the same length as there are rows in NewData giving the prediction in binary (1=presence, 0=absence)
}


\author{ Wilfried Thuiller, Bruno Lafourcade, Damien Georges }
\seealso{ \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{BIOMOD_ModelingOptions}}, \code{\link[biomod2]{BIOMOD_Projection}} }
\examples{
require(raster)

# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'#'PteropusGiganteus'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[which(myResp==1),c("X_WGS84","Y_WGS84")]

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

# we build a raster layer based on environmental rasters for our response variable
myResp <- reclassify(subset(myExpl,1,drop=TRUE), c(-Inf,Inf,0))
myResp[cellFromXY(myResp,myRespXY)] <- 1


# Compute some SRE for several quantile values
g <- sre(Response = myResp, Explanatory = myExpl, NewData=myExpl, Quant=0)
gg <- sre(Response = myResp, Explanatory = myExpl, NewData=myExpl, Quant=0.025)
ggg <- sre(Response = myResp, Explanatory = myExpl, NewData=myExpl, Quant=0.05)

# visualise results
par(mfrow=c(2,2),mar=c(6, 5, 5, 3))
plot(myResp, main=paste(myRespName,"original distrib."))
plot(g, main="full data calibration")
plot(gg, main="Perc025")
plot(ggg, main="Perc05")


}

\keyword{ models }
\keyword{ multivariate }

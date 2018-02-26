\name{multiple.plot}
\alias{multiple.plot}

\title{ Plot and compare prediction maps within BIOMOD }
\description{
  This function allows a direct comparison of the predictions across models. The
  maps produced are the same than with the \code{\link{level.plot}} function 
}
\usage{
multiple.plot(Data,
              coor, 
              color.gradient='red', 
              plots.per.window=9, 
              cex=1, 
              save.file="no", 
              name="multiple plot", 
              ImageSize = "small", 
              AddPresAbs = NULL, 
              PresAbsSymbol = c(cex * 0.8, 16, 4))
}

\arguments{
  \item{Data}{ the data to visualise }
  \item{coor}{ a 2 columns matrix of the same length as data.in giving the coordinates of the data points }  
  \item{color.gradient}{ available : 'red', 'grey' and 'blue' }
  \item{plots.per.window}{ the number of plots to visualise per graphic window }
  \item{cex}{ to change the point size : >1 will increase size, <1 will decrease it }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
  \item{ImageSize}{ The image size for JPEG and TIFF files saved to disk. Available : 'small', 'standard' and 'large' }
  \item{AddPresAbs}{ Optional: adds the presences and absences used for calibration to the plot. Data must be entered as a matrix/dataframe with 3 columns (in this order): X-coordinate, Y-coordinate, Presence(1) or Absence(0). X and Y coordinates must be in the same system as the plot. }
  \item{PresAbsSymbol}{ Optional: a 3 element vector giving the symbols to be used by the AddPresAbs argument for plotting. The elements of the vector must be in this order: size of presence/absence symbols given as a multiplication factor of the 'cex' value entered in the function (e.g. a value of 0.5 means that the symbols will be drawn at a size = 0.5*cex value entered in the function), symbol (in PCH code) to be used for presences, symbol (in PCH code) to be used for absences. An example of input vector for this parameter is 'c(0.4,16,4)' }
}


\author{ Wilfried Thuiller, Bruno Lafourcade }

\seealso{ \code{\link{level.plot}} }
\examples{
\dontrun{
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)

# the name of studied species
myRespName <- c("ConnochaetesGnou", "GuloGulo", "PantheraOnca", 
                "PteropusGiganteus", "TenrecEcaudatus", "VulpesVulpes")

# the presence/absences data for our species 
myResp <- DataSpecies[,myRespName]

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]

multiple.plot(Data = myResp,
              coor = myRespCoord )
}
}

\keyword{ plot }

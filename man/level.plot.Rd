\name{level.plot}
\alias{level.plot}

\title{ Plot 2-dimensional data for visualizing distribution of species or environment}
\description{
  Enables to plot data with 2-dimensional geographical coordinates
}
\usage{
level.plot( data.in, 
            XY,
            color.gradient = "red",
            cex = 1,
            level.range = c(min(data.in), max(data.in)),
            show.scale = TRUE,
            title = "level plot",
            SRC=FALSE,
            save.file="no",
            ImageSize="small",
            AddPresAbs=NULL,
            PresAbsSymbol=c(cex*0.8,16,4),
            ...)
}

\arguments{
  \item{data.in}{ the data you want to visualize }
  \item{XY}{ a 2 columns matrix of the same length as data.in giving the coordinates of the data points}
  \item{color.gradient}{ available : red, grey and blue }
  \item{cex}{ to change the point size : >1 will increase size, <1 will decrease it }
  \item{level.range}{ the range of values for which you want the color gradient to be used, useful to increase the resolution
  of the graph especially if you have extreme values (see examples section) }
  \item{show.scale}{ a feature for just keeping the graph without the scale (if set to False) }
  \item{title}{ the title wanted for the graph }
  \item{SRC}{ if TRUE, the function recognizes the values as being 'Species Range Change' (see XXX) values and associates the appropriate colors to them }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{ImageSize}{ The image size for JPEG and TIFF files saved to disk. Available : 'small', 'standard' and 'large' }
  \item{AddPresAbs}{ Optional: adds the presences and absences used for calibration to the plot. Data must be entered as a matrix/dataframe with 3 columns (in this order): X-coordinate, Y-coordinate, Presence(1) or Absence(0). X and Y coordinates must be in the same system as the plot. }
  \item{PresAbsSymbol}{ Optional: a 3 element vector giving the symbols to be used by the AddPresAbs argument for plotting. The elements of the vector must be in this order: size of presence/absence symbols given as a multiplication factor of the 'cex' value entered in the function (e.g. a value of 0.5 means that the symbols will be drawn at a size = 0.5*cex value entered in the function), symbol (in PCH code) to be used for presences, symbol (in PCH code) to be used for absences. An example of input vector for this parameter is 'c(0.4,16,4)' }
  \item{...}{extra args}
}

\author{ Bruno Lafourcade }

\seealso{ \code{\link{multiple.plot}} }
\examples{
\dontrun{
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species 
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


level.plot(data.in=myResp, XY=myRespXY)
}
}

\keyword{ plot }

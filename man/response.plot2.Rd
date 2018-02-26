\name{response.plot2}
\alias{response.plot2}

\title{ Function for for plotting predicted responses from species distribution models in 2 or 3 dimensions}
\description{ Adaptation of the Evaluation Strip method proposed by Elith et al.(2005). This function enables to plot the response curves of a model independently of the algorithm used for
building the model. It therefore permits a direct comparison of predicted responses from the different statistical approaches on the same data. 
  
}
\usage{
response.plot2( models,
                Data,
                show.variables=seq(1:ncol(Data)),
                do.bivariate = FALSE,
                fixed.var.metric = 'mean',
                save.file="no",
                name="response_curve",
                ImageSize=480,
                plot=TRUE,
                ...)
}

\arguments{
  \item{models}{ a character vector specifying the models for which the response curves have to be plotted. Compatible with GAM, GBM, GLM, ANN, CTA, RF, FDA, MARS and MAXENT.}
  \item{Data}{ a dataframe containing the variables for which the response curves have to be plotted. They must have the same names as the ones used to calibrate the model. RasterStack are also supported. }
  \item{show.variables}{ the names or the column numbers of 'Data' for selecting the variables to be plotted. By default all columns are taken }
  \item{do.bivariate}{ 'logical', if FALSE (default), the predicted response curves are plotted for every singe variable independently (2-D dimension). If TRUE, the predicted response curves are represented in 3-D for all pairs of variables}
  \item{fixed.var.metric}{ either 'mean' (default), 'median', 'min' or 'max' specifying the statistic used to fix as constant the remaining variables when the predicted response is estimated for one of the variables}
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
  \item{ImageSize}{ the size of the image in pixels if save.file is different to "no". Affects "jpeg" and "tiff" outputs only. Default if 480 pixels which is the R default. }
  \item{plot}{ if TRUE (the default) then a plot is produced}
  \item{\ldots}{ further arguments :
    \itemize{
      \item{\code{data_species} : vector containing data species occurances. Have to match with \code{Data}. If given, the statistic used to fix variables value will be calculated only over presences points. (Considered only if \code{Data} is under table format)}
      \item{\code{col} : vector of colors to be used (only for univariate case)}
      \item{\code{lty} : vector of lines types to be used}
      \item{\code{main} : character, the title of the graph (default one based on first model class is automatically produced if not refered)}
      \item{\code{display_title} : logical, display or not graph title}
      \item{\code{legend} : logical, add legend to plot (only for univariate case)}
    }}
}
\details{
  For building the predicted response curves, n-1 variables are set constant to a fixed value (mean, median, min or max i.e \code{fixed.var.metric} arg) and only the remaining one (resp. 2 for 3D response plot) is varying across its whole range (given by \code{Data}). n the case of categorical variable, the most represented class is taken. The variations observed and the curve thus obtained shows the sensibility of the model to that specific variable. This method does not account for interactions between variables.
In the evaluation strip initially proposed by Elith et al. 2005 the remaining variables are set to the mean. 
}

\value{
    a 4D array is returned. It contains the necessary outputs to produce the plots. This is useful to make your own custom response plot graphics.  
    
  Array returned structure : 
  \itemize{
    \item{First dimension:}{ the dimension of the predictions}
    \item{Second dimension:}{ 2 or 3 columns: The first one (resp. the first two) is (are) the explanatory variable(s) to plot, the last one, the probability of occurrence}
    \item{Third dimension:}{ The set of environmental variables for which the response.plot was asked to run.}
    \item{Fourth dimension:}{ the selected models}
  
  }

  
  
}

\author{ Wilfried Thuiller, Damien Georgs}

\references{ 
Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting 
predicted responses from species distribution models. Ecological Modelling 186, 280-289.
}

\seealso{ \code{\link{BIOMOD_Modeling}}  }
\examples{
\dontrun{
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'VulpesVulpes'

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
                                       models = c('RF','GLM'), 
                                       models.options = myBiomodOption, 
                                       NbRunEval=2, 
                                       DataSplit=80, 
                                       Prevalence=0.5,
                                       VarImport=0, 
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       do.full.models=FALSE)


# 4. Plot response curves

# 4.1 Load the models for which we want to extract the predicted response curves
myGLMs <- BIOMOD_LoadModels(myBiomodModelOut, models='GLM')

# 4.2 plot 2D response plots
myRespPlot2D <- response.plot2(models  = myGLMs,
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'median',
                               col = c("blue", "red"),
                               legend = TRUE,
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'))

# 4.2 plot 3D response plots
## here only for a lone model (i.e "VulpesVulpes_PA1_AllData_GLM")
myRespPlot3D <- response.plot2(models  = myGLMs[1],
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = TRUE,
                               fixed.var.metric = 'median',
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                               display_title=FALSE)

### all the values used to produce this plot are stored into the returned object
### you can redo plots by yourself and customised them
dim(myRespPlot2D)
dimnames(myRespPlot2D)

dim(myRespPlot3D)
dimnames(myRespPlot3D)
}
}

\keyword{ dplot }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }

\name{response.plot}
\alias{response.plot}

\title{ Analysis of the response curves of a model within Biomod}
\description{ Depreciated function, please use \code{\link[biomod2]{response.plot2}} instead
}

\usage{
response.plot(model, 
              Data, 
              show.variables=seq(1:ncol(Data)), 
              save.file="no", 
              name="response_curve", 
              ImageSize=480, 
              plot=TRUE)
}

\arguments{
  \item{model}{ the model for which you want the response curves to be plotted. Compatible with GAM, GBM, GLM, ANN, CTA, RF, FDA and MARS.}
  \item{Data}{ the variables for which you want the response curves to be plotted. A data frame is wanted with one column per vriable. They have to have the same names 
   as the ones used to calibrate the model.}
  \item{show.variables}{ give in the column numbers of 'Data' for selecting the variables that are wanted for plotting }
  \item{save.file}{ can be set to "pdf", "jpeg" or "tiff" to save the plot. Pdf options can be changed by setting the default values of pdf.options().}
  \item{name}{ the name of the file produced if save.file is different to "no" (extensions are already included) }
  \item{ImageSize}{ the size of the image in pixels if save.file is different to "no". Affects "jpeg" and "tiff" outputs only. Default if 480 pixels which is the R default. }
  \item{plot}{ if TRUE (the default) then a plot is produced. If not, an array containing predictions is returned (see details)}
}

\details{

  Depreciated function, please use \code{\link[biomod2]{response.plot2}} instead.
  
}

\author{ Wilfried Thuiller}

\references{ 
Elith, J., Ferrier, S., Huettmann, FALSE. & Leathwick, J. R. 2005 The evaluation strip: A new and robust method for plotting 
predicted responses from species distribution models. Ecological Modelling 186, 280-289.
}

\keyword{ dplot }
\keyword{ models }
\keyword{ regression }
\keyword{ nonlinear }
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{ tree }

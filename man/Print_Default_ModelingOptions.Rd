\name{Print_Default_ModelingOptions}
\alias{Print_Default_ModelingOptions}

\title{ Get default values of BIOMOD inner models' options }
\description{
  This function print the default parameters used to build each model.
}

\usage{
Print_Default_ModelingOptions()
}

\details{
This function is useful to parameterize the selected models. It gives a formatted list of all parameters the user can modify. You can copy this function output, modify all parameters you want in a text editor and paste the modified string as argument to (\link{BIOMOD_ModelingOptions}) function.
}

\value{ 
 Nothing returned.
}

\author{ Wilfried Thuiller, Damien Georges }

\examples{
# print default models options
Print_Default_ModelingOptions()
}

\keyword{ models }
\keyword{ options }
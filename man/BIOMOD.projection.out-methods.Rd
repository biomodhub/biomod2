\name{BIOMOD.projection.out-method}

%\alias{get_predictions}
%\alias{get_predictions,BIOMOD.projection.out-method}
\alias{get_projected_models}
\alias{get_projected_models,BIOMOD.projection.out-method}
\alias{free}
\alias{free,BIOMOD.projection.out-method}

\title{ BIOMOD.projection.out getters }

\description{
Functions to get attributes of \code{\link[biomod2]{BIOMOD_Projection}} outputs
}

\usage{
  %get_predictions(obj, ...)
  get_projected_models(obj, ...)
  free(obj, ...)
}


\arguments{
  \item{obj}{ \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} object }
  \item{...}{ extra arguments (see details) }
}

\note{
\code{\link{get_predictions}} is also available for \code{"\link[=BIOMOD.projection.out-class]{BIOMOD.projection.out}"} object.
}


\details{

  \ldots{ } available values :

  \itemize{
    \item{get_predictions}
    \itemize{
      \item{\code{as.data.frame}:}{"logical", return projections into a \code{data.frame} (TRUE) or an \code{array} (FALSE, default)}
      \item{\code{full.name}:}{NULL or "character", full names of models you want to extract projections }
      \item{\code{model}:}{NULL or "character", \emph{type name} of models you want to extract projections }
      \item{\code{run.eval}:}{NULL or "character", \emph{run name} of models you want to extract projections }
      \item{\code{data.set}:}{NULL or "character", \emph{dataset name} of models you want to extract projections }

    }
  }
}


\value{

  \enumerate{
    \item{\bold{get_predictions: }}{ an \code{array} or a \code{data.frame} containing models projections}
    \item{\bold{get_projected_models: }}{ \code{character} containing names of models that have been projected}
    \item{\bold{free: }}{remove projection value from R memory and only keep the path to corresponding file on hard drive.}
  }
}

\author{ Wilfried Thuiller, Damien Georges}

\seealso{ \code{\link{help}}}

\keyword{models}
\keyword{projection}
\keyword{project}
\keyword{forecast}

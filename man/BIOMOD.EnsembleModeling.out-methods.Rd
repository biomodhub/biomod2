\name{BIOMOD.EnsembleModeling.out-method}

%%\alias{get_built_models}
%%\alias{get_built_models,BIOMOD.EnsembleModeling.out-method}
%%\alias{get_evaluations}
%%\alias{get_evaluations,BIOMOD.EnsembleModeling.out-method}
\alias{get_kept_models}
\alias{get_kept_models,BIOMOD.EnsembleModeling.out-method}
\alias{get_needed_models}
\alias{get_needed_models,BIOMOD.EnsembleModeling.out-method}

\title{ BIOMOD.EnsembleModeling.out getters }

\description{
Functions to get attributes of \code{\link[biomod2]{BIOMOD_EnsembleModeling}} outputs
}

\usage{
  %get_built_models(obj, ...)
  %get_evaluations(obj, ...)
  get_kept_models(obj, ...)
  get_needed_models(obj, ...)
}



\arguments{
  \item{obj}{ \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} object }
  \item{...}{ extra arguments (see details) }
}

\note{
\code{\link{get_built_models}} & \code{\link{get_evaluations}} are also available for \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} object.
}


\details{

  \ldots{ } available values :

  \itemize{

    \item{\bold{get_evaluations}}{
      \itemize{
        \item{\code{as.data.frame}:}{"logical", ( FALSE by default ) if TRUE, a standardized \code{data.frame} will be produced else a list is returned }
      }
    }

    \item{\bold{get_kept_models}}{
      \itemize{
        \item{model}{"character" or "numeric" referring to model names (get_built_models()) or model id }
      }
    }


  }
}


\value{

  \enumerate{
    \item{\bold{get_built_models: }}{ a \code{character} vector indicating set of ensemble-modeling algorithms ran}

    \item{\bold{get_built_models: }}{ a \code{character} vector indicating the names of ensemble models computed }

    \item{\bold{get_evaluations: }}{ an \code{array}, a \code{data.frame} or a \code{list} containing ensemble models evaluation scores}

    \item{\bold{get_kept_models: }}{ an \code{character} vector indicating names of selected models for ensemble-models building}

    \item{\bold{get_needed_models: }}{ an \code{character} vector indicating names of all needed models required to build all ensemble-models}

  }
}

\author{ Wilfried Thuiller, Damien Georges}

\seealso{ \code{\link{help}}}

\keyword{ models }
\keyword{ ensemble }
\keyword{ get }

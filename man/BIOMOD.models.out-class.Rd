\name{BIOMOD.models.out-class}

\docType{class}

\alias{BIOMOD.models.out-class}

% BIOMOD.models.out-method
\alias{show,BIOMOD.models.out-method}


\title{ BIOMOD_modelling() outputs objects class}

\description{
The BIOMOD.models.out objects are created, used and returned by BIOMOD functions.

\itemize{
  \item{output of: }{\code{\link[biomod2]{BIOMOD_Modeling}}}
  \item{input of: }{\code{\link[biomod2]{BIOMOD_EnsembleModeling}}, \code{\link[biomod2]{BIOMOD_Projection}}}
}

}


\section{Slots}{

BIOMOD.models.out objects
  \describe{
    \item{\code{modeling.id}:}{ "character", id of modeling procedure}
    \item{\code{sp.name}:}{ "character", species name }
    \item{\code{expl.var.names}:}{ "character", explanatory variables names }
    \item{\code{has.evaluation.data}:}{"logical", is some data are reserved for evaluating models?}
    \item{\code{models.computed}:}{ "character", names of computed models }
    \item{\code{models.failed}:}{ "character", names of failed models }
    \item{\code{models.evaluation}:}{ "BIOMOD.stored.array", evaluations of each model computed according to selected evaluation methods }
    \item{\code{variables.importances}:}{ "BIOMOD.stored.array", models variable importances }
    \item{\code{models.prediction}:}{ "BIOMOD.stored.array", predictions of each models on their own calibration + validation dataset }
    \item{\code{models.prediction.eval}:}{ "BIOMOD.stored.array", predictions of each models on evaluation dataset (if defined) }
    \item{\code{formated.input.data}:}{ "BIOMOD.stored.formated.data", input data }
    \item{\code{calib.lines}:}{ "BIOMOD.stored.array", calibrations and evaluation lines selected for each run }
    \item{\code{models.options}:}{ "BIOMOD.stored.models.options", options used to build each model }
    \item{\code{rescal.all.models}:}{"logical", is all models scaling done?}
    \item{\code{link}:}{'character', the path to object hard drive copy}
  }

}

\author{ Damien Georges }


\examples{
showClass("BIOMOD.models.out")
}

\keyword{models}
\keyword{option}


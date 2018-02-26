\name{BIOMOD.EnsembleModeling.out-class}

\docType{class}

\alias{BIOMOD.EnsembleModeling.out-class}
\alias{BIOMOD.EnsembleModeling.out}

% BIOMOD.models.out-method
\alias{show,BIOMOD.EnsembleModeling.out-method}


\title{ BIOMOD_EnsembleModeling() outputs objects class}

\description{
EnsembleModeling objects are created, used and returned by BIOMOD functions. It's contains information relative to an \pkg{biomod2} ensemble modeling procedure. 

\itemize{
  \item{output of: }{\code{\link[biomod2]{BIOMOD_EnsembleModeling}}}
  \item{input of: }{\code{\link[biomod2]{BIOMOD_EnsembleForecasting}}}
}

}

\section{Slots}{

  \describe{
    \item{\code{sp.name}:}{ "character", species name }
    \item{\code{expl.var.names}:}{ "character", explanatory variables names }
    \item{\code{models.out.obj}:}{"BIOMOD.stored.models.out", object which contains information on individuals models that have been combined}
    \item{\code{eval.metric}:}{ "character", evaluation metrics choosed for models selection }
    \item{\code{eval.metric.quality.threshold}:}{ "numeric", thresholds defined for models selection }
    \item{\code{em.computed}:}{ "character", ensemble models built names}
    \item{\code{em.by}:}{ "character", way models are combined}
    \item{\code{em.models}:}{ "ANY", list of built biomod2.ensemble.models objects }
    \item{\code{modeling.id}:}{ "character", the id of the whole modelling process}
    \item{\code{link}:}{ "character", the path to correspunding hard drive saved object}
  }

}

\seealso{
  \code{\link[biomod2]{BIOMOD_Projection}}, \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{BIOMOD_EnsembleModeling}}, \code{\link[biomod2]{BIOMOD_EnsembleForecasting}}
}

\author{ Damien Georges }

\examples{
showClass("BIOMOD.EnsembleModeling.out")
}

\keyword{models, ensemble}

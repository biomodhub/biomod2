\name{BIOMOD.projection.out-class}

\docType{class}

\alias{BIOMOD.projection.out-class}
\alias{BIOMOD.projection.out}

% BIOMOD.projection.out-method
\alias{show,BIOMOD.projection.out-method}
\alias{plot,BIOMOD.projection.out,missing-method}

\title{ BIOMOD_Projection() outputs objects class}

\description{
BIOMOD.projection.out object is created, used and returned by \pkg{biomod2} functions. It contains information relative to projections of \code{\link[biomod2]{BIOMOD_Modeling}} outputs over a defined area. This object may be reused latter by \code{\link[biomod2]{BIOMOD_EnsembleForecasting}} function.

\itemize{
  \item{output of: }{\code{\link[biomod2]{BIOMOD_Projection}}}
  \item{input of: }{\code{\link[biomod2]{BIOMOD_EnsembleForecasting}}}
}

}


\section{Slots}{

  \describe{
    \item{\code{proj.names}:}{ "character", projection name (name of folder where projections are stored) }
    \item{\code{sp.name}:}{ "character", species name }
    \item{\code{expl.var.names}:}{ "character", explanatory variables names }
    \item{\code{models.projected}:}{"character", models projected names }
    \item{\code{scaled.models}:}{ "logical", was projected models scaled? }
    \item{\code{modeling.object}:}{"BIOMOD.stored.data", associated modeling objects (link to)} 
    \item{\code{modeling.object.id}:}{"character",id of associated modeling objects}
    \item{\code{type}:}{ "character", type of projection done (\code{array}, \code{data.frame}, \code{RasterStack}, \code{RasterLayer} \ldots) }
    \item{\code{proj}:}{ "BIOMOD.stored.data", object containing projected values}
    \item{\code{xy.coord}:}{ "matrix", coordinates of projected points}
  }

}

                        
\author{ Damien Georges }

\seealso{
  \code{\link[biomod2]{BIOMOD_Projection}}, \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{BIOMOD_EnsembleModeling}}, \code{\link[biomod2]{BIOMOD_EnsembleForecasting}}

}

\examples{
showClass("BIOMOD.projection.out")
}

\keyword{models, projection, ensemble, forecast}

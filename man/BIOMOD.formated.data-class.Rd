\name{BIOMOD.formated.data-class}

\docType{class}

\alias{BIOMOD.formated.data-class}
\alias{BIOMOD.formated.data}
\alias{BIOMOD.formated.data.PA-class}
\alias{BIOMOD.formated.data.PA}

\alias{BIOMOD.formated.data,data.frame,ANY-method}
\alias{BIOMOD.formated.data,numeric,RasterStack-method}
\alias{BIOMOD.formated.data,numeric,data.frame-method}
\alias{BIOMOD.formated.data,numeric,matrix-method}

% BIOMOD.models.out-method
\alias{show,BIOMOD.formated.data-method}
\alias{plot,BIOMOD.formated.data,missing-method}
\alias{show,BIOMOD.formated.data.PA-method}
\alias{plot,BIOMOD.formated.data.PA,missing-method}

\title{ BIOMOD_FormatingData() outputs objects class}

\description{
BIOMOD.formated.data objects are created, used and returned by BIOMOD functions. It's contains the minimal set of data \pkg{biomod2} needs to work. Input data given to \code{\link[biomod2]{BIOMOD_FormatingData}} are rearranged to fit with \code{\link[biomod2]{BIOMOD_Modeling}} input format. All data are stored into matrix (even environmental raster) explaining why some objects appears to be quite heavy.

If you ask for pseudo absences selection in \code{\link[biomod2]{BIOMOD_FormatingData}}, you will get a \code{BIOMOD.formated.data.PA}, else you will get a \code{BIOMOD.formated.data.PA} object.

\itemize{
  \item{output of: }{\code{\link[biomod2]{BIOMOD_FormatingData}}}
  \item{input of: }{\code{\link[biomod2]{BIOMOD_Modeling}}}
}

}


\section{Slots}{

BIOMOD.formated.data

  \describe{
    \item{\code{sp.name}:}{ "character", species name }
    \item{\code{coord}:}{ "data.frame", species points XY coordinates }
    \item{\code{data.species}:}{"numeric", species presences, absences and no information kept points}
    \item{\code{data.env.var}:}{ "data.frame", explanatory variables associated to species points }
    \item{\code{has.data.eval}:}{ "logical", was specific models evaluation dataset given ? }
    \item{\code{eval.coord}:}{ "data.frame", species models evaluation points XY coordinates}
    \item{\code{eval.data.species}:}{ "numeric", species presences, absences and no information models evaluation kept points}
    \item{\code{eval.data.env.var}:}{ "data.frame", explanatory variables associated to models evaluation species points}
  }

BIOMOD.formated.data

  \describe{
    \item{\code{sp.name}:}{ "character", species name }
    \item{\code{coord}:}{ "data.frame", species points XY coordinates }
    \item{\code{data.species}:}{"numeric", species presences, absences and no information kept points}
    \item{\code{data.env.var}:}{ "data.frame", explanatory variables associated to species points }
    \item{\code{has.data.eval}:}{ "logical", was specific models evaluation dataset given ? }
    \item{\code{eval.coord}:}{ "data.frame", species models evaluation points XY coordinates}
    \item{\code{eval.data.species}:}{ "numeric", species presences, absences and no information models evaluation kept points}
    \item{\code{eval.data.env.var}:}{ "data.frame", explanatory variables associated to models evaluation species points}
    \item{\code{PA}:}{"data.frame", each column contains 1 or 0 indicating for each species points if it was select or not in associated PA dataset }
  }

}


\author{ Damien Georges }

\examples{
showClass("BIOMOD.formated.data")
}

\keyword{models}
\keyword{data}
\keyword{formating}

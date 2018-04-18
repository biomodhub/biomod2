\name{BIOMOD.models.out-methods}

\alias{get_predictions}
\alias{get_predictions,BIOMOD.models.out-method}
\alias{get_predictions,BIOMOD.EnsembleModeling.out-method}
\alias{get_predictions,BIOMOD.projection.out-method}
\alias{get_evaluations}
\alias{get_calib_lines,BIOMOD.models.out-method}
\alias{get_calib_lines}
\alias{get_evaluations,BIOMOD.models.out-method}
\alias{get_evaluations,BIOMOD.EnsembleModeling.out-method}
\alias{get_variables_importance}
\alias{get_variables_importance,BIOMOD.models.out-method}
\alias{get_variables_importance,BIOMOD.EnsembleModeling.out-method}
\alias{get_options}
\alias{get_options,BIOMOD.models.out-method}
\alias{get_predictions}
\alias{get_predictions,BIOMOD.models.out-method}
\alias{get_built_models}
\alias{get_built_models,BIOMOD.models.out-method}
\alias{get_built_models,BIOMOD.EnsembleModeling.out-method}
\alias{get_formal_data}
\alias{get_formal_data,BIOMOD.models.out-method}

\title{ BIOMOD.models.out getters }

\description{
Functions to get attributes of \code{\link[biomod2]{BIOMOD_Modeling}} outputs
}

\usage{
  get_predictions(obj, ...)
  get_formal_data(obj, ...)
  get_evaluations(obj, ...)
  get_calib_lines(obj, ...)
  get_variables_importance(obj, ...)
  get_options(obj, ...)
  get_built_models(obj, ...)
}



\arguments{
  \item{obj}{ \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} object }
  \item{...}{ extra arguments (see details) }
}

\details{

  \ldots{ } available values :

  \itemize{

    \item{\bold{get_evaluations}}{
      \itemize{
        \item{\code{as.data.frame}:}{"logical", ( FALSE by default ) if TRUE, a standardized \code{data.frame} will be produced else a 4-dimension \code{array} is returned }
      }
    }

    \item{\bold{get_calib_lines}}{
    }

    \item{\bold{get_predictions (for BIOMOD_Modeling() outputs only)}}{
      \itemize{
        \item{\code{as.data.frame}:}{logical(default FALSE). If TRUE, models predictions will be returned as \code{data.frame} rather than \code{array}}
        \item{\code{evaluation}:}{logical (default FALSE). If TRUE, model prediction over evaluation data will be returned }}
    }

    \item{\bold{get_formal_data}}{
      \itemize{
        \item{\code{subinfo}:}{character. Flag defining a specific information to extract from \code{"\link[=BIOMOD.formated.data-class]{BIOMOD.formated.data}"} object. Supported values are: }
          \itemize{
            \item{\code{NULL}:}{ (default) the whole \code{"\link[=BIOMOD.formated.data-class]{BIOMOD.formated.data}"} object is returned }
            \item{\code{'MinMax'}:}{ All explanatory variables ranges returned}
            \item{\code{'resp.var'}:}{ Response variables \code{vector} returned}
            \item{\code{'eval.resp.var'}:}{ Evaluation response variables \code{vector} returned}
            \item{\code{'expl.var'}:}{ Explanatory variables \code{data.frame} returned}
            \item{\code{'eval.expl.var'}:}{ Evaluation explanatory variables \code{data.frame} returned}
            \item{\code{'expl.var.names'}:}{ Explanatory variables names returned}


          }
      }
    }
  }
}


\value{

  \enumerate{
    \item{\bold{get_predictions: }}{ an \code{array} (or a \code{data.frame}) containing models predictions over calibrating and testing data (those used for evaluate models)}

    \item{\bold{get_calib_lines: }}{ an \code{array} (or a \code{data.frame}) having the same dimension than the output of \bold{get_predictions()} of logical values. All lines containing TRUE have been used to calibrate the model}

    \item{\bold{get_evaluations: }}{ an \code{array}, a \code{data.frame} or a \code{list} containing models evaluation scores}


    \item{\bold{get_variables_importance: }}{ an \code{array} containing models variables importance}

    \item{\bold{get_options: }}{ a \code{"\link[=BIOMOD.Model.Options-class]{BIOMOD.Model.Options}"} reporting options used to build individual models}

    \item{\bold{get_formal_data: }}{a \code{"\link[=BIOMOD.formated.data-class]{BIOMOD.formated.data}"} object containing data used for models building and evaluation, or a part of this object}

    \item{\bold{get_built_models: }}{a character vector giving the names of models successfully computed}

  }
}

\author{ Wilfried Thuiller, Damien Georges}

\seealso{ \code{\link{help}}}

\keyword{ models }

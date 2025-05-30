% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomod2_classes_5.R
\docType{class}
\name{biomod2_ensemble_model}
\alias{biomod2_ensemble_model}
\alias{biomod2_ensemble_model-class}
\alias{EMmean_biomod2_model-class}
\alias{EMmedian_biomod2_model-class}
\alias{EMcv_biomod2_model-class}
\alias{EMci_biomod2_model-class}
\alias{EMca_biomod2_model-class}
\alias{EMwmean_biomod2_model-class}
\alias{EMmode_biomod2_model-class}
\alias{EMfreq_biomod2_model-class}
\alias{show,biomod2_ensemble_model-method}
\title{Ensemble model output object class (when running \code{BIOMOD_EnsembleModeling()})}
\usage{
\S4method{show}{biomod2_ensemble_model}(object)
}
\arguments{
\item{object}{a \code{\link{biomod2_ensemble_model}} object}
}
\description{
Class created by \code{\link{BIOMOD_EnsembleModeling}}
}
\details{
\code{biomod2_model} is the basic object for \pkg{biomod2} ensemble species distribution models. 
\cr All listed classes below are derived from \code{biomod2_model}, and have a 
\code{model_class} slot specific value :

\itemize{
  \item \code{biomod2_ensemble_model} : \code{model_class} is \code{EM}
  \item \code{EMmean_biomod2_model} : \code{model_class} is \code{EMmean}
  \item \code{EMmedian_biomod2_model} : \code{model_class} is \code{EMmedian}
  \item \code{EMcv_biomod2_model} : \code{model_class} is \code{EMcv}
  \item \code{EMci_biomod2_model} : \code{model_class} is \code{EMci}
  \item \code{EMca_biomod2_model} : \code{model_class} is \code{EMca}
  \item \code{EMwmean_biomod2_model} : \code{model_class} is \code{EMwmean}
  \item \code{EMmode_biomod2_model} : \code{model_class} is \code{EMmode}
  \item \code{EMfreq_biomod2_model} : \code{model_class} is \code{EMfreq}
}
}
\section{Slots}{

\describe{
\item{\code{modeling.id}}{a \code{character} corresponding to the name (ID) of the simulation set}

\item{\code{model_name}}{a \code{character} corresponding to the model name}

\item{\code{model_class}}{a \code{character} corresponding to the model class}

\item{\code{model_options}}{a \code{list} containing the model options}

\item{\code{model}}{the corresponding model object}

\item{\code{scaling_model}}{the corresponding scaled model object}

\item{\code{dir_name}}{a \code{character} corresponding to the modeling folder}

\item{\code{resp_name}}{a \code{character} corresponding to the species name}

\item{\code{expl_var_names}}{a \code{vector} containing names of explanatory variables}

\item{\code{expl_var_type}}{a \code{vector} containing classes of explanatory variables}

\item{\code{expl_var_range}}{a \code{list} containing ranges of explanatory variables}

\item{\code{model_evaluation}}{a \code{data.frame} containing the model evaluations}

\item{\code{model_variables_importance}}{a \code{data.frame} containing the model variables importance}
}}

\examples{

showClass("biomod2_ensemble_model")
showClass("EMmean_biomod2_model")
showClass("EMmedian_biomod2_model")
showClass("EMcv_biomod2_model")
showClass("EMci_biomod2_model")
showClass("EMca_biomod2_model")
showClass("EMwmean_biomod2_model")
showClass("EMmode_biomod2_model")
showClass("EMfreq_biomod2_model")


}
\seealso{
\code{\link{biomod2_model}}, \code{\link{BIOMOD_EnsembleModeling}}

Other Toolbox objects: 
\code{\link{BIOMOD.ensemble.models.out}},
\code{\link{BIOMOD.formated.data}},
\code{\link{BIOMOD.formated.data.PA}},
\code{\link{BIOMOD.models.options}},
\code{\link{BIOMOD.models.out}},
\code{\link{BIOMOD.options.dataset}},
\code{\link{BIOMOD.options.default}},
\code{\link{BIOMOD.projection.out}},
\code{\link{BIOMOD.rangesize.out}},
\code{\link{BIOMOD.stored.data}},
\code{\link{biomod2_model}}
}
\author{
Damien Georges
}
\concept{Toolbox objects}

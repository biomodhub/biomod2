\name{biomod2_model-class}

\docType{class}

\alias{biomod2_model-class}
\alias{biomod2_model}
\alias{biomod2_ensemble_model-class}
\alias{biomod2_ensemble_model}

% biomod2_model-method
\alias{show,biomod2_model-method}
\alias{get_formal_model}
\alias{get_formal_model,biomod2_model-method}
\alias{get_scaling_model}
\alias{get_scaling_model,biomod2_model-method}
\alias{check_data_range}
\alias{get_var_range}
\alias{get_var_type}


\alias{ANN_biomod2_model-class}
\alias{ANN_biomod2_model}
\alias{predict,ANN_biomod2_model-method}

\alias{CTA_biomod2_model-class}
\alias{CTA_biomod2_model}
\alias{predict,CTA_biomod2_model-method}

\alias{FDA_biomod2_model-class}
\alias{FDA_biomod2_model}
\alias{predict,FDA_biomod2_model-method}

\alias{GAM_biomod2_model-class}
\alias{GAM_biomod2_model}
\alias{predict,GAM_biomod2_model-method}

\alias{GLM_biomod2_model-class}
\alias{GLM_biomod2_model}
\alias{predict,GLM_biomod2_model-method}

\alias{GBM_biomod2_model-class}
\alias{GBM_biomod2_model}
\alias{predict,GBM_biomod2_model-method}

\alias{MARS_biomod2_model-class}
\alias{MARS_biomod2_model}
\alias{predict,MARS_biomod2_model-method}

\alias{MAXENT.Phillips_biomod2_model-class}
\alias{MAXENT.Phillips_biomod2_model}
\alias{predict,MAXENT.Phillips_biomod2_model-method}

\alias{MAXENT.Tsuruoka_biomod2_model-class}
\alias{MAXENT.Tsuruoka_biomod2_model}
\alias{predict,MAXENT.Tsuruoka_biomod2_model-method}

\alias{RF_biomod2_model-class}
\alias{RF_biomod2_model}
\alias{predict,RF_biomod2_model-method}

\alias{SRE_biomod2_model-class}
\alias{SRE_biomod2_model}
\alias{predict,SRE_biomod2_model-method}

\alias{EMca_biomod2_model-class}
\alias{EMca_biomod2_model}
\alias{predict,EMca_biomod2_model-method}

\alias{EMci_biomod2_model-class}
\alias{EMci_biomod2_model}
\alias{predict,EMci_biomod2_model-method}

\alias{EMcv_biomod2_model-class}
\alias{EMcv_biomod2_model}
\alias{predict,EMcv_biomod2_model-method}

\alias{EMmean_biomod2_model-class}
\alias{EMmean_biomod2_model}
\alias{predict,EMmean_biomod2_model-method}

\alias{EMmedian_biomod2_model-class}
\alias{EMmedian_biomod2_model}
\alias{predict,EMmedian_biomod2_model-method}

\alias{EMwmean_biomod2_model-class}
\alias{EMwmean_biomod2_model}
\alias{predict,EMwmean_biomod2_model-method}


\title{ biomod2 models objects class and functions}

\description{
  This objects represent all biomod2 models as it was a 'classical' \R model. It is produce during \code{\link[biomod2]{BIOMOD_Modeling}} step. It can be printed. \code{\link{predict}} function is also implemented for this object and supports \code{\link{data.frame}} and \code{{RasterStack}} inputs. You cen rise formal models with \code{get_formal_model} and the associated scaling GLM (if exists) with \code{get_scaling_model}. All this object have also their own \code{predict} function.

}


\section{Slots}{

  \describe{
    \item{\code{model_name}:}{ "character", model name }
    \item{\code{model_class}:}{"character", type of model (e.g. 'GLM', 'RF', 'MAXENT.Phillips','EMmean')}
    \item{\code{model_options}:}{ "list", list of options used for model building }
    \item{\code{model}:}{"ANY", the formal R model (if exists)}
    \item{\code{scaling_model}:}{"ANY", the associated glm to scale prediction on 0-1 (if exists)}
    \item{\code{resp_name}:}{"character", response variable name}
    \item{\code{expl_var_names}:}{"character", explanatory variables names}
    \item{\code{model_evaluation}:}{"matrix", model evaluation scores (if exist)}
    \item{\code{model_variables_importance}:}{"matrix", model importance of variables (if exists)}
    
    \item{\code{model_output_dir}:}{"character", path to model output directory}

  }


}


\author{ Damien Georges }

\examples{
showClass("ANN_biomod2_model")
}

\keyword{models}
\keyword{predict}

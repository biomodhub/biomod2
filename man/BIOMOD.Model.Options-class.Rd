\name{BIOMOD.Model.Options-class}

\docType{class}

\alias{BIOMOD.Model.Options-class}

% BIOMOD.Model.Options-method
\alias{show,BIOMOD.Model.Options-method}

\title{ BIOMOD_ModelingOptions() outputs objects class}

\description{
BIOMOD.Model.Options objects are created, used and returned by BIOMOD functions. These objects will contains for each model support within \pkg{biomod2}, a set of options that users can change. Please refer to \code{\link[biomod2]{BIOMOD_ModelingOptions}} for futher details. 

\itemize{
  \item{output of: }{\code{\link[biomod2]{BIOMOD_ModelingOptions}}}
  \item{input of: }{\code{\link[biomod2]{BIOMOD_Modeling}}}
}

}


\section{Slots}{

Please refer to \code{\link[biomod2]{BIOMOD_ModelingOptions}} for each model arguments supported.
  \describe{
    \item{\code{GLM}:}{ "list", list of GLM supported options }
    \item{\code{GBM}:}{ "list", list of GBM supported options }
    \item{\code{GAM}:}{ "list", list of GAM supported options }
    \item{\code{CTA}:}{ "list", list of CTA supported options }
    \item{\code{ANN}:}{ "list", list of ANN supported options }
    \item{\code{SRE}:}{ "list", list of SRE supported options }
    \item{\code{FDA}:}{ "list", list of FDA supported options }
    \item{\code{MARS}:}{ "list", list of MARS supported options }
    \item{\code{RF}:}{ "list", list of RF supported options }
    \item{\code{MAXENT.Phillips}:}{ "list", list of MAXENT.Phillips supported options }
    \item{\code{MAXENT.Tsuruoka}:}{ "list", list of MAXENT.Tsuruoka supported options }
  }

}

\author{ Damien Georges }

\seealso{
  \code{\link[biomod2]{BIOMOD_ModelingOptions}}.
}

\examples{
showClass("BIOMOD.Model.Options")
}

\keyword{models, options}

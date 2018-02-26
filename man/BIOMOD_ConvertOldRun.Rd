\name{BIOMOD_ConvertOldRun}
\alias{BIOMOD_ConvertOldRun}

\title{ Convert objects and outputs from BIOMOD.xx into biomod2.xx objects and outputs}

\description{
  This function converts workspace, modelling outputs, results and objects created with version xx of \pkg{BIOMOD} into \pkg{biomod2} objects and re-organized the directories to be used with \pkg{biomod2}
}
\usage{
BIOMOD_ConvertOldRun(savedObj,
                     path = NULL)
}

\arguments{
  \item{savedObj}{ a BIOMOD.1.xx workspace image. It's a .Rdata file named 'Biomod_run.RData' for plurispecific run or Your_Species_Name.RData if you have done monospecific modelling }
  \item{path}{ Optional path to your 'savedObj' if you don't have give the full path to the object }
}

\details{
  This function is useful to convert former \pkg{BIOMOD} runs into the new \pkg{biomod2} object structure. This is mostly interesting in the case users want to relaunched some projections or analyses within the \pkg{biomod2} new structure. Returned 'BIOMOD.models.out' objects can be then used as classic object for making projections for instance (\code{\link{BIOMOD_Projection}}).
  Be aware that because \pkg{biomod2} has strongly changed between the first and second version, some new additional functions and information could not be used with converted objects (i.e. Calibration Lines access, Maxent run, SRE projections...).   

}

\value{
  A list of 'BIOMOD.models.out' (one per species modeled) containing information of your old run. Specific directories are also created on your hard drive (see \code{\link{BIOMOD_Modeling}}) 
}


\author{ Damien Georges }

\seealso{ \code{\link{BIOMOD_Modeling}}}
\examples{
## Not Done Yet ##
##
##
}
                 
\keyword{ models }
\keyword{ datasets }

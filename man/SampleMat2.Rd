\name{SampleMat2}
\alias{SampleMat2}

\title{ Sample binary vector  }

\description{ \code{SampleMat2} is an internal \pkg{biomod2} function that can help user to sample a binary vector keeping the same proportion of 0s and 1s than in the initial vector.
}

\usage{
  SampleMat2(ref, ratio, as.logi=FALSE)
}

\arguments{
  \item{ref}{ a binary vector }
  \item{ratio}{ the proportion of \code{ref} to sample }
  \item{as.logi}{ logical, if FALSE (default) id of cell will be return; if TRUE, logical vector of same length than ref will be return }
}

\details{
  This function can be useful to help users to select a part of initial dataset that will be only kept for all validation procedures. 

}
\value{
  A list of 2 elements is returned :
  \itemize{
    \item{calibration}{ Ids of cells selected for calibration (1st sample)}
    \item{evaluation}{ Ids of cells selected for evaluation (1st sample complementary)}
  }
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_FormatingData}}}

\examples{
  a <- sample(c(0,1),100, replace=TRUE)
  SampleMat2(ref=a, ratio=0.7)
}

\keyword{ models, formula, options }

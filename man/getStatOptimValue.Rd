\name{getStatOptimValue}
\alias{getStatOptimValue}

\title{get the optimal score of evaluation statistical metrics }

\description{ \code{getStatOptimValue} is an internal \pkg{biomod2} function to get the best value that some referenced evaluation statistical metrics.
}

\usage{
  getStatOptimValue(stat)
}

\arguments{
  \item{stat}{  either  'TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD', 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR' or 'ORSS'}
}

\details{
  Please refer to \code{\link[biomod2]{BIOMOD_Modeling}} to get more information about this metrics.

}
\value{
  The best value that \code{stat} could rise.
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{calculate.stat}}, \code{\link[biomod2]{Find.Optim.Stat}} }

\examples{
  getStatOptimValue('TSS')
  getStatOptimValue('KAPPA')
  getStatOptimValue('POFD')
}

\keyword{models}
\keyword{formula}
\keyword{options}

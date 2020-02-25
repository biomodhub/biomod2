\name{calculate.stat}
\alias{calculate.stat}

\title{Calculate evaluation metrics based on a misclassification table }

\description{ \code{calculate.stat} is an internal \pkg{biomod2} function to get scores, based on a misclassification table, of some referenced evaluation metrics.
}

\usage{
  calculate.stat( Misc,
                  stat='TSS')
}

\arguments{
  \item{Misc}{ a misclassification table }
  \item{stat}{  either 'TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD', 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR' or 'ORSS'}
}

\details{
  Please refer to \code{\link[biomod2]{BIOMOD_Modeling}} to get more information about this metrics.

}
\value{
  The \code{stat} score for the \code{Misc} table.
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link[biomod2]{getStatOptimValue}}, \code{\link[biomod2]{Find.Optim.Stat}} }

\examples{
  a <- sample(c(0,1),100, replace=TRUE)
  b <- sample(c(0,1),100, replace=TRUE)
  
  miscTab_aa <- table(a,a)
  miscTab_ab <- table(a,b)
  
  # perfect score
  calculate.stat( miscTab_aa, stat='TSS')
  # random score
  calculate.stat( miscTab_ab, stat='TSS')
  
}

\keyword{models}
\keyword{formula}
\keyword{options}


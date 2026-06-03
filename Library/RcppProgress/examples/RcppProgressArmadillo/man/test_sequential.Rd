\name{test_multithreaded}
\alias{test_multithreaded}
\title{
  Test/example of using RcppProgress in multithreaded code
}
\description{
  Test function for the package rngOpenMP in a multithreaded mode.
}
\usage{
  test_multithreaded(max=100, nb=1000, threads=0, display_progress=TRUE)
}
\arguments{
  \item{max}{ number of long computation jobs to perform }
  \item{nb}{ number of interruptible steps a long computaton job is made of }
  \item{threads}{ number of threads to use, or 0 to use the default OpenMP number, that is usually the number of cores }
  \item{display_progress}{ whether to display the progress bar or not }
}
\details{
 Will perform max long jobs made of nb interruptible steps.
 The jobs are parallelized. The computation can be interrupted by the user.
}
\value{
  None
}
\author{
  Karl Forner
}



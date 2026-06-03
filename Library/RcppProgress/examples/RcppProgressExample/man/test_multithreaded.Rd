\name{test_sequential}
\alias{test_sequential}
\title{
  Test/example of using RcppProgress in sequential code, i.e not in multithreaded code
}
\description{
  Test function for the package rngOpenMP in a sequential mode.
}
\usage{
  test_sequential(max=100, nb=1000, display_progress=TRUE)
}
\arguments{
  \item{max}{ number of long computation jobs to perform }
  \item{nb}{ number of interruptible steps a long computaton job is made of }
  \item{display_progress}{ whether to display the progress bar or not }
}
\details{
 Will perform max long jobs made of nb interruptible steps.
 The jobs are performed sequentially. The computation can be interrupted by the user.
}
\value{
  None
}
\author{
  Karl Forner
}



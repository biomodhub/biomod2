\name{randomise_data}
\alias{randomise_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 data set shuffling tool
%%  ~~function to do ... ~~
}
\description{
  This function is wrapper to shuffle elegantly data-set. This function shuffle selected columns of a given data.set according to a defined shuffling algorithm.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
randomise_data(data,variable,method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ a \code{data.frame} you want to shuffle}
  \item{variable}{ the variables names or columns id you want to shuffle}
  \item{method}{ the randomisation method (only 'full_rand' available yet)}
}

\value{
  a \code{data.frame} with selected columns shuffled compared to the original table.
}

\author{
  Damien Georges
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{variables_importance}}, \code{\link{full_suffling}}
}

\examples{
xx <- data.frame(a=1:10,b=11:20,c=21:30)
randomise_data(data=xx, variable='b', method='full_rand')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ suffle }
\keyword{ random }
\keyword{ importance }% __ONLY ONE__ keyword per line

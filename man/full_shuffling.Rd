\name{full_suffling}
\alias{full_suffling}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 data set shuffling tool
%%  ~~function to do ... ~~
}
\description{
  This function is developer tool to shuffle elegantly data-set. This function shuffle one column of a given data.set
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
full_suffling(x,id=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a \code{data.frame} you want to shuffle}
  \item{id}{ the columns you want to shuffle}
}


\value{
  a \code{data.frame} with selected columns shuffled compared to the original table.
}

\author{
  Damien Georges
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{variables_importance}}
}

\examples{
xx <- matrix(rep(1:10,3),10,3)
full_suffling(xx,c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ shuffle }
\keyword{ random }
\keyword{ importance }% __ONLY ONE__ keyword per line

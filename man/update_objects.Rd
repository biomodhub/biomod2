\name{update_objects}
\alias{update_objects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 update biomod2 objects
%%  ~~function to do ... ~~
}
\description{
  This function is wrapper to update objects construct with a old version of biomod2 to a current one
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{
  update_objects(obj, recursive=TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ a \code{biomod2} object you want to update (e.g. 'BIOMOD.formated.data', 'biomod2_model' object)}
  \item{recursive}{ logical, if TRUE all objects on which updated object is based will be also updtaed }
}

\details{
  This function will add/change/delete all object slots that have evolved between 2 versions of the package. 
  If necessary, objects stored on hard drive will also be updated.
}


\value{
  the updated version of the  biomod2 object is return
}

\author{
  Damien Georges
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{variables_importance}}, \code{\link{full_suffling}}
}

\examples{
## not done yet ##
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


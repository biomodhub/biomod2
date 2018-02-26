\name{BIOMOD.stored.objects-class}

\docType{class}

\alias{BIOMOD.stored.data}
\alias{BIOMOD.stored.data-class}
\alias{BIOMOD.stored.files}
\alias{BIOMOD.stored.files-class}
\alias{BIOMOD.stored.data.frame}
\alias{BIOMOD.stored.data.frame-class}
\alias{BIOMOD.stored.array}
\alias{BIOMOD.stored.array-class}
\alias{BIOMOD.stored.formated.data}
\alias{BIOMOD.stored.formated.data-class}
\alias{BIOMOD.stored.models.options}
\alias{BIOMOD.stored.models.options-class}
\alias{BIOMOD.stored.models.out}
\alias{BIOMOD.stored.models.out-class}
\alias{BIOMOD.stored.raster.stack}
\alias{BIOMOD.stored.raster.stack-class}

\title{ BIOMOD.stored.xxx objects class}

\description{
All BIOMOD.stored.xxx objects are there to make some RAM savings by loading some objects only when they are needed. We just keep a link to the place were those objects are stored when we don't need them anymore. They are kind of pointers to objects created within \pkg{biomod2} and having a copy on hard drive. It can be consider just as a code trick.
}



\section{Slots}{

Each BIOMOD.stored.xxx objects has 3 slot:
  \describe{
    \item{\code{inMemory}:}{ "logical", is the object already loaded on RAM? }
    \item{\code{link}:}{ "character", path to the hard drive copy of the object }
    \item{\code{val}:}{ NULL or the object if loaded.. The object class is different for each BIOMOD.stored.xxx objects }
  }

}


\author{ Damien Georges }

\examples{
showClass("BIOMOD.stored.files")
}

\keyword{models}
\keyword{ensemble}
\keyword{object}
\keyword{storing}

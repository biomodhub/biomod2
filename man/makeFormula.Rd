\name{makeFormula}
\alias{makeFormula}

\title{Standardized formula maker }

\description{ makeFormula is an internal \pkg{biomod2} function that can be useful to help users to build easily some standardized formula used later by statistical models.
}

\usage{
makeFormula(respName, 
            explVar, 
            type = 'simple', 
            interaction.level = 0,
            ...)
}



\arguments{
  \item{respName}{  a \code{character} indicating the response variable name }
  \item{explVar}{ a \code{matrix} or a \code{data.frame}, the explanatory variables table that will be considered at modelling step}
  \item{type}{ either 'simple', 'quadratic', 'polynomial' or 's_smoother' defining the type of formula you want to build}
  \item{interaction.level}{ an \code{integer}, the interaction level depth between explanatory variables}
  \item{...}{some additional arguments (see details)}
}


\details{
It's advise to give only a subset of \code{explVar} table to avoid useless memory consuming. If some explanatory variables are factorial ones, you have to give a \code{data.frame} for \code{explVar} where associated columns are define as \code{factor}.

\code{...} argument availables values are :
\itemize{
  \item{k}{ the smooting parameter value (used only if \code{type = 's_smoother'}) corespunding to \code{k} parameter of \pkg{mgcv} \code{\link[mgcv]{s}}} or \code{df} \pkg{gam} \code{\link[gam]{s}} arguments.

}
}


\value{
A \code{link[stats]{formula}} class object that can be directly given to most of \R statistical models.   
}


\author{ Damien Georges }

\seealso{ \code{\link[biomod2]{BIOMOD_ModelingOptions}}, \code{link[stats]{formula}} }
\examples{

# create simulated data
myResp <- sample(c(0,1),20, replace=TRUE)
myExpl <- matrix(runif(60), ncol=3, dimnames=list(NULL,c('var1','var2','var3')))

# create a formula
myFormula <- makeFormula( respName = 'myResp', 
                          explVar = head(myExpl), 
                          type = 'quadratic', 
                          interaction.level = 0)
# show formula created
myFormula

}

\keyword{ models, formula, options }

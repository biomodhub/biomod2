\name{variables_importance}
\alias{variables_importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Variables importance calculation
%%  ~~function to do ... ~~
}
\description{
 This function will return a variable importance value for each variable involved within your model.
}
\usage{
variables_importance(model, data, method="full_rand", nb_rand=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ the model you want to study variables importance (one of the models supported within biomod2, ensemble models are also supported)}
  \item{data}{ the \code{data.set} on which you want to perform analyses}
  \item{method}{ the randomisation method (only 'full_rand' available yet)}
  \item{nb_rand}{ the number of permutation done for each variable}
  \item{...}{ additional args (not implemented yet)}
}

\details{
It's more or less base on the same principle than \code{\link[randomForest]{randomForest}} variables importance algorithm. The principle is to shuffle a single variable of the given data. Make model prediction with this 'shuffled' data.set. Then we compute a simple correlation (Pearson's by default) between references predictions and the 'shuffled' one. The return score is 1-cor(pred_ref,pred_shuffled). The highest the value, the more influence the variable has on the model. A value of this 0 assumes no influence of that variable on the model. Note that this technique does not account for interactions between the variables.
}

\value{
  a \code{list} of class "BIOMOD_variables_importances" which contains:
  \itemize{
    \item{mat}{a \code{data.frame} containing variables importance scores for each permutation run. }
  }
  
}

\author{
  Damien Georges
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link[biomod2]{randomise_data}}, \code{\link[biomod2]{full_suffling}}
}

\examples{
xx <- data.frame( a=sample(c(0,1),100,replace=TRUE),
                  b=rnorm(100),
                  c=1:100)

mod <- glm(a~b+c, data=xx)

variables_importance(model=mod, data=xx[,c('b','c')], method="full_rand", nb_rand=3)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ suffle }
\keyword{ random }
\keyword{ importance }% __ONLY ONE__ keyword per line

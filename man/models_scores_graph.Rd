\name{models_scores_graph}
\alias{models_scores_graph}

\title{ Produce models evaluation bi-dimensional graph }

\description{
  This function is a graphic tool to represent evaluation scores of models produced with  biomod2 according to 2 different evaluation methods. Models can be grouped in several ways (by algo, by CV run, ...) to highlight potential differences in models quality due to chosen models, cross validation sampling bias,...
  Each point represents the average evaluation score across each group. Lines represents standard deviation of evaluation scores of the group.
}

\usage{
  models_scores_graph(obj, 
                      metrics = NULL, 
                      by = 'models', 
                      plot = TRUE, 
                      ... )
}

\arguments{
  \item{obj}{ a \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} ( returned by \code{\link[biomod2]{BIOMOD_Modeling}} ) or a \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} (returned by \code{\link[biomod2]{BIOMOD_EnsembleModeling}}) }
  \item{metrics}{ character vector of 2 chosen metrics (e.g c("ROC", "TSS")); if not filled the two first evaluation methods computed at modeling stage will be selected.}
  \item{by}{ character ('models'), the way evaluation scores are grouped. Should be one of 'models', 'algos', 'CV_run' or 'data_set' (see detail section)  }
  \item{plot}{ logical (TRUE), does plot should be produced }
  \item{...}{ additional graphical arguments (see details)  }
}

\details{
  \bold{ \code{by} argument description : }
  
  \code{by} arg refers to the way models scores will be combined to compute mean and sd. It should take the following values:

  \itemize{
    \item{\code{models}} : group evaluation scores according to top level models. Top level models should be for instance GLM, GAM, RF, SRE... in \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} input case whereas it should be EMcaByTSS (committee averaging using TSS score), EMwmeanByROC (weighted mean using ROC scores),... or whatever in \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} input case.
    \item{\code{algos}} : If you work with \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} then \code{algos} is equivalent to \code{models}. If you work with  \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"} then it refer to formal models i.e. GLM, GAM, RF, SRE... (should also be mergedAlgo in this case).
    \item{\code{cv_run}} : the cross validation run e.g. run1, run2,... Should be mergedRun in EnsembleModels input case.
    \item{\code{data_set}} : the data set (resp. pseudo absences data set) used to group scores e.g PA1, PA2,... if pseudo absences sampling have been computed or AllData inf not. Should also be mergedData in EnsembleModels case.
  
  }
  

  \bold{Additional arguments (\ldots{}) : }
  
  Additional graphical parameters should be.
  
  \itemize{
    \item{\code{xlim}}{ the graphical range represent for the first evaluation metric }
    \item{\code{ylim}}{ the graphical range represent for the second evaluation metric }
    \item{\code{main}}{ main plot title }
  }
  
}

\note{
  This function have been instigate by \cite{Elith*, J., H. Graham*, C., P. Anderson, R., Dudik, M., Ferrier, S., Guisan, A., J. Hijmans, R., Huettmann, F., R. Leathwick, J., Lehmann, A., Li, J., G. Lohmann, L., A. Loiselle, B., Manion, G., Moritz, C., Nakamura, M., Nakazawa, Y., McC. M. Overton, J., Townsend Peterson, A., J. Phillips, S., Richardson, K., Scachetti-Pereira, R., E. Schapire, R., Soberon, J., Williams, S., S. Wisz, M. and E. Zimmermann, N. (2006), Novel methods improve prediction of species distributions from occurrence data. Ecography, 29: 129-151. doi: 10.1111/j.2006.0906-7590.04596.x} (fig 3)
}

\value{
  A ggplot2 plotting object is return. It means that user should then easily customize this plot (see example) 
}


\author{ Damien Georges }

\seealso{ \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}}

\examples{

## this example is based on BIOMOD_Modeling function example
example(BIOMOD_Modeling)

## we will need ggplot2 package to produce our custom version of the graphs
require(ggplot2)

## plot evaluation models score graph

### by models
gg1 <- models_scores_graph( myBiomodModelOut,
                            by = 'models',
                            metrics = c('ROC','TSS') )
## we see a influence of model selected on models capabilities
## e.g. RF are much better than SRE

### by cross validation run
gg2 <- models_scores_graph( myBiomodModelOut,
                            by = 'cv_run',
                            metrics = c('ROC','TSS') )
## there is no difference in models quality if we focus on 
## cross validation sampling

### some graphical customisations
gg1_custom <- 
  gg1 + 
  ggtitle("Diff between RF and SRE evaluation scores") + ## add title
  scale_colour_manual(values=c("green", "blue")) ## change colors

gg1_custom

}
                 
\keyword{ evaluation }
\keyword{ scores }
\keyword{ graph }

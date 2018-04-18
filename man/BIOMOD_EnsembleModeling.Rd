\name{BIOMOD_EnsembleModeling}
\alias{BIOMOD_EnsembleModeling}

\title{ Create and evaluate an ensemble set of models and predictions}
\description{
\code{BIOMOD_EnsembleModeling} combines models and make ensemble predictions built with \code{\link[biomod2]{BIOMOD_Modeling}}. The ensemble predictions can also be evaluated against the original data given to  \code{\link[biomod2]{BIOMOD_Modeling}}. Biomod2 proposes a range of options to build ensemble models and predictions and to assess the modeling uncertainty. The created ensemble models can then be used to project distributions over space and time as classical \pkg{biomod2} models.
}

\usage{
BIOMOD_EnsembleModeling( modeling.output,
                         chosen.models = 'all',
                         em.by = 'PA_dataset+repet',
                         eval.metric = 'all',
                         eval.metric.quality.threshold = NULL,
                         models.eval.meth = c('KAPPA','TSS','ROC'),
                         prob.mean = TRUE,
                         prob.cv = FALSE,
                         prob.ci = FALSE,
                         prob.ci.alpha = 0.05,
                         prob.median = FALSE,
                         committee.averaging = FALSE,
                         prob.mean.weight = FALSE,
                         prob.mean.weight.decay = 'proportional',
                         VarImport = 0)
}

\arguments{
  \item{modeling.output}{a \code{"\link[=BIOMOD.models.out-class]{BIOMOD.models.out}"} returned by \code{\link[biomod2]{BIOMOD_Modeling}}}
  \item{chosen.models}{a character vector (either 'all' or a sub-selection of model names) that defines the models kept for building the ensemble models (might be useful for removing some non-preferred models)}
  \item{em.by}{Character. Flag defining the way the models will be combined to build the ensemble models. Available values are \code{'PA_dataset+repet'} (default), \code{'PA_dataset+algo'}, \code{'PA_dataset'}, \code{'algo'} and  \code{'all'}}
  \item{eval.metric}{vector of names of evaluation metric used to build ensemble models. It is involved for formal models exclusion if \code{eval.metric.quality.threshold} is defined and/or for building ensemble models that are dependent of formal models evaluation scores (e.g. weighted mean and committee averaging). If 'all', the same evaluation metrics than those of \code{modeling.output} will be automatically selected}
  \item{eval.metric.quality.threshold}{If not \code{NULL}, the minimum scores below which models will be excluded of the ensemble-models building.}
  \item{models.eval.meth}{ the evaluation methods used to evaluate ensemble models ( see \code{"\link[=BIOMOD_Modeling]{BIOMOD_Modeling}"}  \bold{models.eval.meth} section for more detailed informations )}
  \item{prob.mean}{Logical. Estimate the mean probabilities across predictions }
  \item{prob.cv}{Logical. Estimate the coefficient of variation across predictions}
  \item{prob.ci}{Logical . Estimate the confidence interval around the \code{prob.mean}}
  \item{prob.ci.alpha}{Numeric. Significance level for estimating the confidence interval. Default = 0.05}
  \item{prob.median}{Logical. Estimate the median of probabilities}
  \item{committee.averaging}{Logical. Estimate the committee averaging across predictions}
  \item{prob.mean.weight}{Logical. Estimate the weighted sum of probabilities}
  \item{prob.mean.weight.decay}{Define the relative importance of the weights. A high value will strongly discriminate the 'good' models from the 'bad' ones (see the details section). If the value of this parameter is set to 'proportional' (default), then the attributed weights are proportional to the evaluation scores given by 'weight.method'(\code{eval.metric})}
  \item{VarImport}{Number of permutation to estimate variable importance}
}

\details{
  \enumerate{
    \item{\bold{Models sub-selection (\code{chosen.models})}}{

      Useful to exclude some models that have been selected in the previous steps (\code{modeling.output}). This vector of model names can be access applying \code{get_built_models} to your \code{modeling.output} data. It makes easier the selection of models. The default value (i.e. \sQuote{all}) will kept all available models.}

    \item{\bold{Models assembly rules (\code{em.by})}}{

    Please refer to \href{../doc/EnsembleModelingAssembly.pdf}{EnsembleModelingAssembly} vignette that is dedicated to this parameter.

      5 different ways to combine models can be considered. You can make ensemble models considering :
      \itemize{
        \item{ Dataset used for models building (Pseudo Absences dataset and repetitions done)}{: \code{'PA_dataset+repet'}}
        \item{ Dataset used and statistical models }{: \code{'PA_dataset+algo'}}
        \item{ Pseudo-absences selection dataset }{: \code{'PA_dataset'}}
        \item{ Statistical models }{: \code{'algo'}}
        \item{ A total consensus model  }{: \code{'all'}}
      }

      The value chosen for this parameter will control the number of ensemble models built.
      If no evaluation data was given the at \code{\link[biomod2]{BIOMOD_FormatingData}} step, some ensemble models evaluation may be a bit unfair because the data that will be used for evaluating ensemble models could differ from those used for evaluate \code{\link[biomod2]{BIOMOD_Modeling}} models (in particular, some data used for 'basal models' calibration can be re-used for ensemble models evaluation). You have to keep it in mind ! (\href{../doc/EnsembleModelingAssembly.pdf}{EnsembleModelingAssembly} vignette for extra details)
    }

    \item{\bold{Evaluation metrics}}{
      \itemize{
        \item{\bold{\code{eval.metric}}}{

          The selected metrics here are necessary the ones chosen at the \code{\link[biomod2]{BIOMOD_Modeling}} step.
          If you select several, ensembles will be built according to each of them.
          The chosen metrics will be used at different stages in this function :
          \enumerate{
            \item{to remove \sQuote{bad models} (having a score lower than \code{eval.metric.quality.threshold} (see bellow))}
            \item{to make the binary transformation needed for committee averaging computation}
            \item{to weight the models in the probability weighted mean model}
            \item{to test (and/or evaluate) your ensemble-models forecasting ability (at this step, each ensemble-model (ensemble will be evaluated according to each evaluation metric)}
          } % end of enumerate
        }


        \item{\bold{\code{eval.metric.quality.threshold}}}{

          You have to give as many threshold as \code{eval.metric} you have selected.  If you have selected  several evaluation metrics , take care to ensure that the order of thresholds matches the order of \code{eval.metric}.
          All models having a score lower than these quality thresholds will not be kept for building ensemble-models.
        }
      }
    }

    \item{\bold{Ensemble-models algorithms}}{
      \enumerate{
        \item{\bold{Mean of probabilities (\code{prob.mean}) }}{

          This ensemble-model corresponds to the mean probabilities over the selected models.
        }

        \item{\bold{Coefficient of variation of Probabilities (\code{prob.cv}) }}{

          This ensemble-model corresponds to the coefficient of variation (i.e. sd / mean) of the probabilities over the selected models. This model is not scaled. It will be evaluated like all other ensemble-models although this interpretation is obviously different. CV is a measure of uncertainty rather a measure of probability of occurrence. If the CV gets a high evaluation score it means that the uncertainty is high where the species is observed (which might not be a good feature of the models). The lower is the score, the better are the models. CV is a nice complement to the mean probability.
        }

        \item{\bold{Confidence interval (\code{prob.ci} & \code{prob.ci.alpha}) }}{
         This is the confidence interval around the mean probability (see above). This is also a nice complement to the mean probability.
          Two ensemble-models will be built if \code{prob.ci} is \code{TRUE} :
          \itemize{
            \item{The upper one (there is less than a 100*\code{prob.ci.alpha}/2 \% of chance to get probabilities upper than the given ones)}
                        \item{The lower one (there is less than a 100*\code{prob.ci.alpha}/2 \% of chance to get probabilities lower the than given ones)}
          }
          These intervals are calculated with the following function :

          \deqn{
            I_c = [ \bar{x} -  \frac{t_\alpha sd }{ \sqrt{n} }; \bar{x} +  \frac{t_\alpha sd }{ \sqrt{n} }]
          }
        }

        \item{\bold{ Median of probabilities (\code{prob.median})}}{

          This ensemble-model corresponds to the median probability over the selected models. The median is less sensitive to outliers than the mean. In practical terms, calculating the median requires more time and memory than the mean (or even weighting mean) as it asks to load all predictions to then extract the median. It may need to be considered in case of large dataset.
        }

        \item{\bold{ Models committee averaging (\code{committee.averaging})}}{

          To do this model, the probabilities from the selected models are first transformed into binary data according to the thresholds defined at \code{\link[biomod2]{BIOMOD_Modeling}} step (maximizing evaluation metric score over \sQuote{testing dataset} ). The committee averaging score is then the average of binary predictions. It is built on the analogy of a simple vote. Each model vote for the species being ether present or absent. For each site, the sum of 1 is then divided by the number of models. The interesting feature of this measure is that it gives both a prediction and a measure of uncertainty. When the prediction is close to 0 or 1, it means that all models agree to predict 0 and 1 respectively. When the prediction is around 0.5, it means that half the models predict 1 and the other half 0.
        }

        \item{\bold{ Weighted mean of probabilities (\code{prob.mean.weight} & \code{prob.mean.weight.decay})}}{

          This algorithm return the mean weighted (or more precisely this is the weighted sum) by the selected evaluation method scores (better a model is, more importance it has in the ensemble). The scores come from \code{\link[biomod2]{BIOMOD_Modeling}} step.

          The \code{prob.mean.weight.decay} is the ratio between a weight and the following or prior one. The formula is : W = W(-1) * \code{prob.mean.weight.decay}. For example, with the value of 1.6 and 4 weights wanted, the relative importance of the weights will be 1 /1.6/2.56(=1.6*1.6)/4.096(=2.56*1.6) from the weakest to the strongest, and gives 0.11/0.17/0.275/0.445 considering that the sum of the weights is equal to one. The lower the \code{prob.mean.weight.decay}, the smoother the differences between the weights enhancing a weak discrimination between models.

          The value 'proportional' (default) is also possible for the \code{prob.mean.weight.decay}: the weights are awarded for each method proportionally to their evaluation scores. The advantage is that the discrimination is more fair than with the \code{prob.mean.weight.decay}. In the latter case, close scores can strongly diverge in the weights they are awarded, when the proportional method will consider them as being fairly similar in prediction quality and award them a similar weight.
          It is also possible to define a function as \code{prob.mean.weight.decay} argument. In this case the given function will be applied to models scores to transform them into weights that will be used for weighted mean ensemble model building. For instance if you specified \code{ function(x){x^2} } as \code{prob.mean.weight.decay}, the squared of evaluation score of each model will be used to weight formal models predictions.
        }

      }
    }

  } % end of itemize

} % end of detail

\value{
A \code{"\link[=BIOMOD.EnsembleModeling.out-class]{BIOMOD.EnsembleModeling.out}"}. This object will be later given to \code{\link[biomod2]{BIOMOD_EnsembleForecasting}} if you want to make some projections of this ensemble-models.

You can access to evaluation scores with the \code{get_evaluations} function and to the built models names with the \code{get_built_models} function (see example).
}

\note{
Models are now combined by repetition, other way to combine them (e.g. by Models, all together...) will be available soon
}


\author{
Damien Georges & Wilfried Thuiller with participation of Robin Engler
}

\seealso{
\code{\link[biomod2]{BIOMOD_Modeling}}, \code{\link{BIOMOD_Projection}}, \code{\link[biomod2]{BIOMOD_EnsembleForecasting}}
}

\examples{
# species occurrences
DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
                                    package="biomod2"), row.names = 1)
head(DataSpecies)

# the name of studied species
myRespName <- 'GuloGulo'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]


# Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
myExpl = raster::stack( system.file( "external/bioclim/current/bio3.grd",
                     package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio7.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio11.grd",
                             package="biomod2"),
                system.file( "external/bioclim/current/bio12.grd",
                             package="biomod2"))

# 1. Formatting Data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# 2. Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()

# 3. Doing Modelisation

myBiomodModelOut <- BIOMOD_Modeling( myBiomodData,
                                       models = c('SRE','CTA','RF'),
                                       models.options = myBiomodOption,
                                       NbRunEval=1,
                                       DataSplit=80,
                                       Yweights=NULL,
                                       VarImport=3,
                                       models.eval.meth = c('TSS'),
                                       SaveObj = TRUE,
                                       rescal.all.models = FALSE,
                                       do.full.models = FALSE)

# 4. Doing Ensemble Modelling
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                           chosen.models = 'all',
                           em.by = 'all',
                           eval.metric = c('TSS'),
                           eval.metric.quality.threshold = c(0.7),
                           models.eval.meth = c('TSS','ROC'),
                           prob.mean = TRUE,
                           prob.cv = FALSE,
                           prob.ci = FALSE,
                           prob.ci.alpha = 0.05,
                           prob.median = FALSE,
                           committee.averaging = FALSE,
                           prob.mean.weight = TRUE,
                           prob.mean.weight.decay = 'proportional' )

# print summary
myBiomodEM

# get evaluation scores
get_evaluations(myBiomodEM)


}

\keyword{ models }

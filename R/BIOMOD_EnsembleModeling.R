## BIOMOD_EnsembleModeling documentation ---------------------------------------
##' @name BIOMOD_EnsembleModeling
##' @author Wilfried Thuiller, Damien Georges, Robin Engler
##' 
##' @title Create and evaluate an ensemble set of models and predictions
##' 
##' @description This function allows to combine a range of models built with the 
##' \code{\link{BIOMOD_Modeling}} function in one (or several) ensemble model. Modeling 
##' uncertainty can be assessed as well as variables importance, ensemble predictions can be 
##' evaluated against original data, and created ensemble models can be projected over new 
##' conditions (see Details).
##' 
##' 
##' @param bm.mod a \code{\link{BIOMOD.models.out}} object returned by the 
##' \code{\link{BIOMOD_Modeling}} function
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' @param em.by a \code{character} corresponding to the way kept models will be combined to build 
##' the ensemble models, must be among \code{PA_dataset+repet}, \code{PA_dataset+algo}, 
##' @param em.algo (\emph{optional, default} \code{c('prob.mean')}) a \code{vector}
##'  corresponding to the ensemble models that should  be computed. 
##'  May contain up to six values among \code{'prob.mean'}, \code{'prob.median'}, 
##'  \code{'prob.cv'}, \code{'prob.ci'}, \code{'committee.averaging'}, 
##'  \code{'prob.mean.weight'}.
##' @param metric.select a \code{vector} containing evaluation metric names to be used together with 
##' \code{metric.select.thresh} to exclude single models based on their evaluation scores 
##' (for ensemble methods like probability weighted mean or committee averaging). Must be among  
##' \code{all} (same evaluation metrics than those of \code{bm.mod}), \code{user.defined} 
##' (and defined through \code{metric.select.table}) or \code{ROC}, \code{TSS}, \code{KAPPA}, 
##' \code{ACCURACY}, \code{BIAS}, \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, 
##' \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, \code{ORSS}
##' @param metric.select.thresh (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} of \code{numeric} values corresponding to the minimum scores (one for each 
##' \code{metric.select}) below which single models will be excluded from the ensemble model 
##' building
##' @param metric.select.table (\emph{optional, default} \code{NULL}) \cr 
##' If \code{metric.select = 'user.defined'}, a \code{data.frame} containing evaluation scores 
##' calculated for each single models and that will be compared to \code{metric.select.thresh} 
##' values to exclude some of them from the ensemble model building, with \code{metric.select} 
##' rownames, and \code{models.chosen} colnames
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param prob.ci.alpha (\emph{optional, default} \code{0.05}) \cr 
##' A \code{numeric} value corresponding to the significance level to estimate confidence interval
##' @param prob.mean.weight.decay (\emph{optional, default} \code{proportional}) \cr 
##' A value defining the relative importance of the weights (if \code{prob.mean.weight = TRUE}). 
##' A high value will strongly discriminate \emph{good} models from the \emph{bad} ones (see Details), while \code{proportional} will 
##' attribute weights proportionally to the models evaluation scores
##' @param prob.mean (\emph{optional, default} \code{TRUE}) \cr A \code{logical}
##'   value defining whether to compute the mean probabilities across
##'   predictions or not
##' @param prob.median (\emph{obsolete}) \cr A \code{logical} value defining
##'   whether to compute the median probabilities across predictions or not
##' @param prob.cv (\emph{obsolete}) \cr A \code{logical} value defining whether
##'   to compute the coefficient of variation across predictions or not
##' @param prob.ci (\emph{obsolete}) \cr A \code{logical} value defining whether
##'   to compute te confidence interval around the \code{prob.mean} ensemble
##'   model or not
##' @param committee.averaging (\emph{obsolete}) \cr A \code{logical} value
##'   defining whether to compute the committee averaging across predictions or
##'   not
##' @param prob.mean.weight (\emph{obsolete}) \cr A \code{logical} value
##'   defining whether to compute the weighted sum of probabilities across
##'   predictions or not
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' @param do.progress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the progress bar is to be rendered or not
##' 
##' 
##' @return
##' 
##' A \code{BIOMOD.ensemble.models.out} object containing models outputs, or links to saved 
##' outputs. \cr Models outputs are stored out of \R (for memory storage reasons) in 2 different 
##' folders created in the current working directory :
##' \enumerate{
##'   \item a \emph{models} folder, named after the \code{resp.name} argument of 
##'   \code{\link{BIOMOD_FormatingData}}, and containing all ensemble models
##'   \item a \emph{hidden} folder, named \code{.BIOMOD_DATA}, and containing outputs related 
##'   files (original dataset, calibration lines, pseudo-absences selected, predictions, 
##'   variables importance, evaluation values...), that can be retrieved with 
##'   \href{https://biomodhub.github.io/biomod2/reference/getters.out.html}{\code{get_[...]}} or \code{\link{load}} functions, and used by other 
##'   \pkg{biomod2} functions, like \code{\link{BIOMOD_EnsembleForecasting}}
##' }
##' 
##' 
##' @details 
##' 
##' \describe{
##'   \item{Models sub-selection (\code{models.chosen})}{Applying \code{\link{get_built_models}} 
##'   function to the \code{bm.mod} object gives the names of the single models created 
##'   with the \code{\link{BIOMOD_Modeling}} function. The \code{models.chosen} argument can take 
##'   either a sub-selection of these single model names, or the \code{all} default value, to 
##'   decide which single models will be used for the ensemble model building.}
##' 
##'   \item{Models assembly rules (\code{em.by})}{Single models built with the 
##'   \code{\link{BIOMOD_Modeling}} function can be combined in 5 different ways to obtain 
##'   ensemble models :
##'   \itemize{
##'     \item{\code{PA_dataset+repet} : }{each combination of pseudo-absence and repetition 
##'     datasets is done, \emph{merging} algorithms together}
##'     \item{\code{PA_dataset+algo} : }{each combination of pseudo-absence and algorithm datasets 
##'     is done, \emph{merging} repetitions together}
##'     \item{\code{PA_dataset} : }{pseudo-absence datasets are considered individually, 
##'     \emph{merging} algorithms and repetitions together}
##'     \item{\code{algo} : }{algorithm datasets are considered individually, \emph{merging} 
##'     pseudo-absence and repetitions together}
##'     \item{\code{all} : }{all models are combined into one}
##'   }
##'   Hence, depending on the chosen method, the number of ensemble models built will vary. \cr
##'   \emph{Be aware that if no evaluation data was given to the 
##'   \code{\link{BIOMOD_FormatingData}} function, some ensemble model evaluations may be biased 
##'   due to difference in data used for single model evaluations.}}
##' 
##'   \item{Evaluation metrics}{
##'   \itemize{
##'     \item{\bold{\code{metric.select}} : }{the selected metrics must be chosen among the ones used 
##'     within the \code{\link{BIOMOD_Modeling}} function to build the \code{model.output} object, 
##'     unless \code{metric.select = 'user.defined'} and therefore values will be provided through 
##'     the \code{metric.select.table} parameter. \cr In the case of the selection of several 
##'     metrics, they will be used at different steps of the ensemble modeling function : 
##'     \enumerate{
##'       \item remove \emph{low quality} single models, having a score lower than 
##'       \code{metric.select.thresh}
##'       \item perform the binary transformation needed if \code{committee.averaging = TRUE}
##'       \item weight models if \code{prob.mean.weight = TRUE}
##'       \item test and/or evaluate the ensemble models built
##'     }
##'     }
##'     \item{\bold{\code{metric.select.thresh}} : }{as many values as evaluation metrics 
##'     selected with the \code{metric.select} parameter, and defining the corresponding quality 
##'     thresholds below which the single models will be excluded from the ensemble model 
##'     building.}
##'     \item{\bold{\code{metric.select.table}} : }{a \code{data.frame} must be given if 
##'     \code{metric.select = 'user.defined'} to allow the use of evaluation metrics other than 
##'     those calculated within \pkg{biomod2}. The \code{data.frame} must contain as many columns 
##'     as \code{models.chosen} with matching names, and as many rows as evaluation metrics to be 
##'     used. The number of rows must match the length of the \code{metric.select.thresh} 
##'     parameter. The values contained in the \code{data.frame} will be compared to those defined 
##'     in \code{metric.select.thresh} to remove \emph{low quality} single models from 
##'     the ensemble model building.}
##'   }
##'   }
##' 
##'   \item{Ensemble-models algorithms}{The set of models to be calibrated on the data. \cr 
##'   10 modeling techniques are currently available :
##'   \itemize{
##'     \item{\bold{\code{prob.mean}} : }{Mean of probabilities over the selected models}
##'     
##'     \item{\bold{\code{prob.median}} : }{Median of probabilities over the selected models \cr 
##'     The median is less sensitive to outliers than the mean, however it requires more 
##'     computation time and memory as it loads all predictions (on the contrary to the mean or 
##'     the weighted mean).}
##'     
##'     \item{\bold{\code{prob.cv}} : }{Coefficient of variation (sd / mean) of probabilities 
##'     over the selected models \cr 
##'     This model is not scaled. It will be evaluated like all other ensemble models although its 
##'     interpretation will be obviously different. CV is a measure of uncertainty rather a 
##'     measure of probability of occurrence. If the CV gets a high evaluation score, it means 
##'     that the uncertainty is high where the species is observed (which might not be a good 
##'     feature of the model). \emph{The lower is the score, the better are the models.} 
##'     CV is a nice complement to the mean probability.}
##'     
##'     \item{\bold{\code{prob.ci}} & \bold{\code{prob.ci.alpha}} : }{Confidence interval around 
##'     the mean of probabilities of the selected models \cr 
##'     It is also a nice complement to the mean probability. It creates 2 ensemble models : 
##'     \itemize{
##'       \item \emph{LOWER} : there is less than \code{100 * prob.ci.alpha / 2} \% of chance to 
##'       get probabilities lower than the given ones
##'       \item \emph{UPPER} : there is less than \code{100 * prob.ci.alpha / 2} \% of chance to 
##'       get probabilities upper than the given ones
##'     }
##'     These intervals are calculated with the following function :
##'     \deqn{I_c = [ \bar{x} -  \frac{t_\alpha sd }{ \sqrt{n} }; 
##'     \bar{x} +  \frac{t_\alpha sd }{ \sqrt{n} }]}
##'     }
##'     
##'     \item{\bold{\code{committee.averaging}} : }{Probabilities from the selected models are 
##'     first transformed into binary data according to the thresholds defined when building the 
##'     \code{model.output} object with the \code{BIOMOD_Modeling} function, maximizing the 
##'     evaluation metric score over the testing dataset. The committee averaging score is 
##'     obtained by taking the average of these binary predictions. It is built on the analogy 
##'     of a simple vote :
##'     \itemize{
##'       \item each single model votes for the species being either present (\code{1}) or absent 
##'       (\code{0})
##'       \item the sum of \code{1} is then divided by the number of single models \emph{voting}
##'     }
##'     The interesting feature of this measure is that it gives both a prediction and a measure 
##'     of uncertainty. When the prediction is close to \code{0} or \code{1}, it means that all 
##'     models agree to predict \code{0} or \code{1} respectively. When the prediction is around 
##'     \code{0.5}, it means that half the models predict \code{1} and the other half \code{0}. 
##'     \cr}
##'     
##'     \item{\bold{\code{prob.mean.weight}} & \bold{\code{prob.mean.weight.decay}} : }{
##'     Probabilities from the selected models are weighted according to their evaluation scores 
##'     obtained when building the \code{model.output} object with the \code{BIOMOD_Modeling} 
##'     function (\emph{better a model is, more importance it has in the ensemble}) and summed.}
##'   }
##'   
##' The \code{prob.mean.weight.decay} is the ratio between a weight and the next or previous one. 
##' The formula is : \code{W = W(-1) * prob.mean.weight.decay}. \emph{For example, with the value 
##' of \code{1.6} and \code{4} weights wanted, the relative importance of the weights will be 
##' \code{1/1.6/2.56(=1.6*1.6)/4.096(=2.56*1.6)} from the weakest to the strongest, and gives 
##' \code{0.11/0.17/0.275/0.445} considering that the sum of the weights is equal to one. The 
##' lower the \code{prob.mean.weight.decay}, the smoother the differences between the weights 
##' enhancing a weak discrimination between models.}
##' 
##' If \code{prob.mean.weight.decay = 'proportional'}, the weights are assigned to each model 
##' proportionally to their evaluation scores. The discrimination is fairer than using the 
##' \emph{decay} method where close scores can have strongly diverging weights, while the 
##' proportional method would assign them similar weights.
##' 
##' It is also possible to define the \code{prob.mean.weight.decay} parameter as a function that 
##' will be applied to single models scores and transform them into weights. \emph{For example, 
##' if \code{prob.mean.weight.decay = function(x) {x^2}}, the squared of evaluation score of each 
##' model will be used to weight the models predictions.}}
##' }
##' 
##' 
##' @keywords models ensemble weights
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_ModelingOptions}}, 
##' \code{\link{BIOMOD_CrossValidation}}, \code{\link{bm_VariablesImportance}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleForecasting}},
##' \code{\link{bm_PlotEvalMean}}, \code{\link{bm_PlotEvalBoxplot}}, 
##' \code{\link{bm_PlotVarImpBoxplot}}, \code{\link{bm_PlotResponseCurves}}
##' @family Main functions
##' 
##'   
##' @examples
##' 
##' library(terra)
##' # Load species occurrences (6 species available)
##' data(DataSpecies)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- 'GuloGulo'
##' 
##' # Get corresponding presence/absence data
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # Get corresponding XY coordinates
##' myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
##' 
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' data(bioclim_current)
##' myExpl <- terra::rast(bioclim_current)
##' 
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExpl <- terra::crop(myExpl, myExtent)
##' }
##' 
##' ## ----------------------------------------------------------------------- #
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                        expl.var = myExpl,
##'                                        resp.xy = myRespXY,
##'                                        resp.name = myRespName)
##' 
##'   # Create default modeling options
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       bm.options = myBiomodOptions,
##'                                       nb.rep = 2,
##'                                       data.split.perc = 80,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       do.full.models = FALSE,
##'                                       seed.val = 42)
##' }
##' 
##' ## ----------------------------------------------------------------------- #
##' # Model ensemble models
##' myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
##'                                       models.chosen = 'all',
##'                                       em.by = 'all',
##'                                       metric.select = c('TSS'),
##'                                       metric.select.thresh = c(0.7),
##'                                       metric.eval = c('TSS', 'ROC'),
##'                                       var.import = 3,
##'                                       prob.mean = TRUE,
##'                                       prob.median = FALSE,
##'                                       prob.cv = FALSE,
##'                                       prob.ci = FALSE,
##'                                       prob.ci.alpha = 0.05,
##'                                       committee.averaging = TRUE,
##'                                       prob.mean.weight = FALSE,
##'                                       prob.mean.weight.decay = 'proportional',
##'                                       seed.val = 42)
##' myBiomodEM
##' 
##' # Get evaluation scores & variables importance
##' get_evaluations(myBiomodEM, as.data.frame = TRUE)
##' get_variables_importance(myBiomodEM, as.data.frame = TRUE)
##' 
##' # Represent evaluation scores
##' bm_PlotEvalMean(bm.out = myBiomodEM, group.by = 'model')
##' bm_PlotEvalBoxplot(bm.out = myBiomodEM, group.by = c('model', 'model'))
##' 
##' # Represent variables importance
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'model', 'model'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('expl.var', 'model', 'dataset'))
##' # bm_PlotVarImpBoxplot(bm.out = myBiomodEM, group.by = c('model', 'expl.var', 'dataset'))
##' 
##' # Represent response curves
##' # bm_PlotResponseCurves(bm.out = myBiomodEM, 
##' #                       models.chosen = get_built_models(myBiomodEM)[c(1, 2)],
##' #                       fixed.var = 'median')
##' # bm_PlotResponseCurves(bm.out = myBiomodEM, 
##' #                       models.chosen = get_built_models(myBiomodEM)[c(1, 2)],
##' #                       fixed.var = 'min')
##' # bm_PlotResponseCurves(bm.out = myBiomodEM, 
##' #                       models.chosen = get_built_models(myBiomodEM)[2],
##' #                       fixed.var = 'median',
##' #                       do.bivariate = TRUE)
##' 
##' 
##' @export
##' @importFrom terra rast  
##' 
## BIOMOD_EnsembleModeling function ------------------------------------------- 

BIOMOD_EnsembleModeling <- function(bm.mod,
                                    models.chosen = 'all',
                                    em.by = 'PA_dataset+repet',
                                    em.algo,
                                    metric.select = 'all',
                                    metric.select.thresh = NULL,
                                    metric.select.table = NULL,
                                    metric.eval = c('KAPPA', 'TSS', 'ROC'),
                                    var.import = 0,
                                    prob.ci.alpha = 0.05,
                                    prob.mean.weight.decay = 'proportional',
                                    nb.cpu = 1,
                                    seed.val = NULL,
                                    do.progress = TRUE,
                                    prob.mean, # deprecated
                                    prob.median, # deprecated
                                    prob.cv, # deprecated
                                    prob.ci, # deprecated  
                                    committee.averaging, # deprecated
                                    prob.mean.weight) { # deprecated
  .bm_cat("Build Ensemble Models")
  
  ## 0. Check arguments --------------------------------------------------------
  args <- .BIOMOD_EnsembleModeling.check.args(bm.mod = bm.mod,
                                              models.chosen = models.chosen,
                                              em.by = em.by,
                                              em.algo = em.algo,
                                              metric.select = metric.select,
                                              metric.select.thresh = metric.select.thresh,
                                              metric.select.table = metric.select.table,
                                              metric.eval = metric.eval,
                                              prob.ci.alpha = prob.ci.alpha,
                                              prob.mean.weight.decay = prob.mean.weight.decay,
                                              prob.mean = prob.mean,
                                              prob.cv = prob.cv,
                                              prob.ci = prob.ci,
                                              prob.median = prob.median,
                                              committee.averaging = committee.averaging,
                                              prob.mean.weight = prob.mean.weight)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## Get selected options 
  em.options <- list(em.by = em.by)
  expl_var_type = get_var_type(get_formal_data(bm.mod, 'expl.var'))
  expl_var_range = get_var_range(get_formal_data(bm.mod, 'expl.var'))
  
  
  ## 1. Create output object ---------------------------------------------------
  EM <- new('BIOMOD.ensemble.models.out',
            modeling.id = bm.mod@modeling.id,
            dir.name = bm.mod@dir.name,
            sp.name = bm.mod@sp.name,
            expl.var.names = bm.mod@expl.var.names)
  EM@models.out@link <- bm.mod@link
  EM@em.by = em.by
  
  ## Various objects will be stored (models, predictions, evaluation)
  name.BIOMOD_DATA = file.path(EM@dir.name, EM@sp.name, ".BIOMOD_DATA", EM@modeling.id, "ensemble.models")
  
  
  ## 2. Do Ensemble modeling ---------------------------------------------------
  em.out <- foreach(assemb = names(em.mod.assemb)) %do%
  {
    cat("\n\n  >", assemb, "ensemble modeling")
    models.kept <- em.mod.assemb[[assemb]]
    
    #### defined data that will be used for models performances calculation 
    if (bm.mod@has.evaluation.data) {
      eval.obs <- get_formal_data(bm.mod, 'eval.resp.var')
      eval.expl <- get_formal_data(bm.mod, 'eval.expl.var')
    }
    
    ### subselection of observations according to dataset used to produce ensemble models ---------
    obs <-  get_formal_data(bm.mod, 'resp.var')
    expl <- get_formal_data(bm.mod, 'expl.var')
    if (em.by %in% c("PA_dataset", 'PA_dataset+algo', 'PA_dataset+repet') &&
        unlist(strsplit(assemb, "_"))[3] != 'AllData') {
      if (inherits(get_formal_data(bm.mod), "BIOMOD.formated.data.PA")) {
        kept_cells <- get_formal_data(bm.mod)@PA.table[, unlist(strsplit(assemb, "_"))[3]]
      } else {
        kept_cells <- rep(TRUE, length(obs))
      }
    } else if (em.by %in% c("algo","all")) { # no other option should be possible
      if (inherits(get_formal_data(bm.mod), "BIOMOD.formated.data.PA")) {
        # get the union of pseudo absences
        kept_cells <- apply(get_formal_data(bm.mod)@PA.table, 1, any) 
      } else {
        kept_cells <- rep(TRUE, length(obs))
      }
    } else { # in case 'AllData'
      kept_cells <- rep(TRUE, length(obs))
    }
    
    obs <- obs[kept_cells]
    expl <- expl[kept_cells, , drop = FALSE]
    obs[is.na(obs)] <- 0
    
    ## get needed models predictions ----------------------------------------
    needed_predictions <- .get_needed_predictions(bm.mod, em.by, models.kept
                                                  , metric.select, metric.select.thresh
                                                  , metric.select.user, metric.select.table
                                                  , nb.cpu)
    
    ## if no prediction selected => switch to next model
    if (length(needed_predictions) == 0) next
    
    ## LOOP over evaluation metrics ------------------------------------------
    em.out.eval <- foreach (eval.m = metric.select) %do%
    {
      models.kept <- needed_predictions$models.kept[[eval.m]]
      models.kept.scores <- needed_predictions$models.kept.scores[[eval.m]]
      
      ## LOOP over em.algo ---------------------------------------------------
      em.out.algo <- foreach (algo = em.algo) %do%
      {
        ListOut <- list(model = NULL,
                        calib.failure = NULL,
                        pred = NULL,
                        pred.eval = NULL,
                        evaluation = NULL,
                        var.import = NULL)
        
        algo.1 <- algo.2 <- algo.3 <- NULL
        models.kept.tmp = models.kept
        if (algo == 'prob.mean') {
          algo.1 <- "Mean of probabilities"
          algo.2 <- algo.3 <- "EMmean"
        } else if (algo == 'prob.cv') {
          algo.1 <- "Coef of variation of probabilities"
          algo.2 <- algo.3 <- "EMcv"
        } else if (algo == 'prob.median') {
          algo.1 <- "Median of probabilities"
          algo.2 <- algo.3 <- "EMmedian"
        } else if (algo == 'prob.ci.inf') {
          algo.1 <- "Confidence Interval"
          algo.2 <- "EMciInf"
          algo.3 <- "EMci"
        } else if (algo == 'prob.ci.sup') {
          algo.1 <- "Confidence Interval"
          algo.2 <- "EMciSup"
          algo.3 <- "EMci"
        } else if (algo == 'committee.averaging') {
          algo.1 <- "Committee averaging"
          algo.2 <- algo.3 <- "EMca"
          
          ## remove models if some thresholds are undefined
          models.kept.thresh <- unlist(lapply(models.kept.tmp, function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(get_evaluations(bm.mod, data.set = dat, run.eval = run, algo = mod, Metric.eval = eval.m)[, "Cutoff"])
          }))
          names(models.kept.thresh) <- models.kept.tmp
          models.kept.tmp = models.kept.tmp[is.finite(models.kept.thresh)]
          models.kept.thresh.tmp = models.kept.thresh[is.finite(models.kept.thresh)]
        } else if (algo == 'prob.mean.weight') {
          algo.1 = "Probabilities weighting mean"
          algo.2 = algo.3 = "EMwmean"
          
          # remove SRE models if ROC
          models.kept.scores.tmp <- models.kept.scores
          if (eval.m == 'ROC') {
            sre.id <- grep("_SRE", models.kept.tmp)
            if (length(sre.id) > 0) {
              cat("\n\n     !! SRE modeling cannot be used with prob.mean.weight by ROC and will be switched off for", assemb, " selected by ROC.")
              models.kept.tmp <- models.kept.tmp[-sre.id]
              models.kept.scores.tmp <- models.kept.scores[-sre.id]
              
              if(length(models.kept.tmp) == 1){
                cat("\n     !! due to SRE switched off, ensemble models for prob.mean.weight in ", assemb, "will be based on only one single model.")
                cat("\n     !! Please make sure this is intended or review your selection metrics and threshold.")
              }
            }
          }
          
          ## remove models if score is not defined
          models.kept.tmp <- models.kept.tmp[is.finite(models.kept.scores.tmp)]
          models.kept.scores.tmp <- models.kept.scores.tmp[is.finite(models.kept.scores.tmp)]
          
          # weights are "decay" times decreased for each subsequent model in model quality order.
          models.kept.scores.tmp <- round(models.kept.scores.tmp, 3) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.
          
          # dealing with numerical decay
          cat("\n\t\t", " original models scores = ", models.kept.scores.tmp)
          if (is.numeric(prob.mean.weight.decay)) {
            DecayCount <- sum(models.kept.scores.tmp > 0)
            WOrder <- order(models.kept.scores.tmp, decreasing = TRUE)
            Dweights <- models.kept.scores.tmp
            for (J in 1:DecayCount) {
              Dweights[WOrder[J]] <- I(prob.mean.weight.decay ^ (DecayCount - J + 1))
            }
            # If 2 or more scores are identical -> make a mean weight between the ones concerned
            for (J in 1:length(models.kept.scores.tmp)) {
              comp = models.kept.scores.tmp[J] == models.kept.scores.tmp
              if (sum(comp) > 1) {
                Dweights[which(comp == TRUE)] <- mean(Dweights[which(comp == TRUE)])
              }
            }
            models.kept.scores.tmp <- round(Dweights, digits = 3)
            rm(list = c('Dweights', 'DecayCount', 'WOrder'))
          } else if (is.function(prob.mean.weight.decay)) { # dealing with function decay
            models.kept.scores.tmp <- sapply(models.kept.scores.tmp, prob.mean.weight.decay)
          }
          
          ### Standardise model weights
          models.kept.scores.tmp <- round(models.kept.scores.tmp / sum(models.kept.scores.tmp, na.rm = TRUE)
                                          , digits = 3)
          cat("\n\t\t", " final models weights = ", models.kept.scores.tmp)
        }
        
        
        # Models building --------------------------------------------------
        cat("\n   >", algo.1, "...")
        model_name <- paste0(bm.mod@sp.name, "_", algo.2, "By", eval.m, "_", assemb)
        model.bm <- new(paste0(algo.3, "_biomod2_model"),
                        model = models.kept.tmp,
                        model_name = model_name,
                        model_class = algo.3,
                        model_options = em.options,
                        dir_name = bm.mod@dir.name,
                        resp_name = bm.mod@sp.name,
                        expl_var_names = bm.mod@expl.var.names,
                        expl_var_type = expl_var_type,
                        expl_var_range = expl_var_range,
                        modeling.id = bm.mod@modeling.id)
        
        if (algo == 'prob.ci.inf') {
          model.bm@alpha <- prob.ci.alpha
          model.bm@side <- 'inferior'
        } else if (algo == 'prob.ci.sup') {
          model.bm@alpha <- prob.ci.alpha
          model.bm@side <- 'superior'
        } else if (algo == 'committee.averaging') {
          model.bm@thresholds <- models.kept.thresh.tmp
        } else if (algo == 'prob.mean.weight') {
          model.bm@penalization_scores <- models.kept.scores.tmp
        }
        
        # Models Predictions --------------------------------------------------
        
        ## create the suitable directory architecture
        pred.bm.name <- paste0(model_name, ".predictions")
        pred.bm.outfile <- file.path(name.BIOMOD_DATA, "ensemble.models.predictions", pred.bm.name)
        dir.create(dirname(pred.bm.outfile), showWarnings = FALSE, recursive = TRUE)
        pred.newdata <- needed_predictions$predictions[which(needed_predictions$predictions$full.name %in% model.bm@model), c("full.name", "Points", "pred")]
        pred.newdata <- tapply(X = pred.newdata$pred, INDEX = list(pred.newdata$Points, pred.newdata$full.name), FUN = mean)
        pred.newdata <- as.data.frame(pred.newdata)
        
        ## store models prediction on the hard drive
        pred.bm <- try(predict(model.bm
                               , newdata = pred.newdata
                               , data_as_formal_predictions = TRUE
                               , on_0_1000 = TRUE
                               , seedval = seed.val))
        
        if (inherits(pred.bm, "try-error")) {
          # keep the name of uncompleted modelisations
          cat("\n   ! Note : ", model_name, "failed!\n")
          ListOut$calib.failure = model_name
          return(ListOut) ## end of function.
        } else {
          ListOut$model <- model_name
          ListOut$pred <- pred.bm
          assign(pred.bm.name, pred.bm)
          save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
          rm(list = pred.bm.name)
          
          ## do the same for evaluation data
          if (exists('eval.obs') && exists('eval.expl') && !inherits(pred.bm, "try-error")){
            
            pred.bm.eval.outfile <- paste0(pred.bm.outfile,"Eval")
            pred.bm.name <- paste0(model_name, ".predictionsEval")
            eval_pred.bm <- predict(model.bm, newdata = eval.expl, seedval = seed.val)
            ListOut$pred.eval <- eval_pred.bm
            assign(pred.bm.name, eval_pred.bm)
            save(list = pred.bm.name, file = pred.bm.eval.outfile, compress = TRUE)
            rm(list = pred.bm.name)
          }
          
          # Models Evaluation ----------------------------------------------------
          
          if (length(metric.eval) > 0 ) {
            cat("\n\t\t\tEvaluating Model stuff...")
            if (!(algo %in% c('prob.cv', 'prob.ci.inf','prob.ci.sup'))) {
              if (em.by == "PA_dataset+repet") {
                ## select the same evaluation data than formal models
                ## get formal models calib/eval lines
                calib_lines <- get_calib_lines(bm.mod)
                ## get info on wich dataset and which repet this ensemble model is based on
                pa_dataset_id <- paste0("_", unlist(strsplit(assemb, "_"))[3])
                repet_id <- paste0("_", unlist(strsplit(assemb, "_"))[2])
                ## define and extract the subset of points model will be evaluated on
                if (repet_id == "_Full") {
                  eval_lines <- rep(TRUE, length(pred.bm))
                } else {
                  ## trick to detect when it is a full model but with a non common name
                  ## i.e. all lines used for calib => full model
                  eval_lines <- !na.omit(calib_lines[, repet_id, pa_dataset_id])
                  if (all(!eval_lines)) { eval_lines <- !eval_lines }
                }
              } else {
                eval_lines <- rep(TRUE, length(pred.bm))
              }
              
              cross.validation <- foreach(xx = metric.eval, .combine = "rbind") %do% {
                bm_FindOptimStat(metric.eval = xx,
                                 obs = obs[eval_lines],
                                 fit = pred.bm[eval_lines])
              }
              colnames(cross.validation)[which(colnames(cross.validation) == "Best.stat")] <- "Testing.data"
              
              if (exists('eval_pred.bm')) {
                true.evaluation <- foreach(xx = metric.eval, .combine = "rbind") %do% {
                  bm_FindOptimStat(metric.eval = xx,
                                   obs = eval.obs,
                                   fit = eval_pred.bm * 1000,
                                   threshold = cross.validation["Cutoff", xx])
                }
                cross.validation$Evaluating.data <- true.evaluation$Best.stat
              } else {
                cross.validation$Evaluating.data <- NA
              }
              
              ## store results
              for (col.i in 2:ncol(cross.validation)) {
                cross.validation[, col.i] <- round(cross.validation[, col.i], digits = 3)
              }
              ListOut$evaluation <- cross.validation
              model.bm@model_evaluation <- cross.validation
            }
          }
          
          # Models Variable Importance -------------------------------------------
          if (var.import > 0 ) {
            cat("\n\t\t\tEvaluating Predictor Contributions...", "\n")
            variables.importance <- 
              bm_VariablesImportance(bm.model = model.bm
                                     , expl.var = expl
                                     , nb.rep = var.import
                                     , seed.val = seed.val
                                     , do.progress = do.progress)
            ListOut$var.import <- variables.importance
            model.bm@model_variables_importance <- variables.importance
            
          }
          
          # Models saving --------------------------------------------------------
          assign(model_name, model.bm)
          save(list = model_name, file = file.path(bm.mod@dir.name, bm.mod@sp.name, "models",
                                                   bm.mod@modeling.id, model_name))
          return(ListOut)
        }
      }
      names(em.out.algo) <- em.algo
      return(em.out.algo)
    }
    names(em.out.eval) <- metric.select
    return(em.out.eval)
  }
  names(em.out) <- names(em.mod.assemb)
  
  
  ### check at least one model was computed -----------------------------------
  EM@em.computed <- .transform_outputs_list("em", em.out, out = "model")
  EM@em.failed <- .transform_outputs_list("em", em.out, out = "calib.failure")
  
  if(length(EM@em.computed) == 1 && EM@em.computed == "none"){
    cat("\n! All models failed")
    return(EM)
  }
  
  ### SAVE EM outputs ---------------------------------------------------------
  models.evaluation <- .transform_outputs_list("em", em.out, out = "evaluation")
  EM = .fill_BIOMOD.models.out("models.evaluation", models.evaluation, EM
                               , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
  if (var.import > 0) {
    variables.importance <- .transform_outputs_list("em", em.out, out = "var.import")
    EM = .fill_BIOMOD.models.out("variables.importance", variables.importance, EM
                                 , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
  }
  models.prediction <- .transform_outputs_list("em", em.out, out = "pred")
  EM = .fill_BIOMOD.models.out("models.prediction", models.prediction, EM
                               , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
  if (bm.mod@has.evaluation.data) {
    models.prediction.eval <- .transform_outputs_list("em", em.out, out = "pred.eval")
    EM = .fill_BIOMOD.models.out("models.prediction.eval", models.prediction.eval, EM
                                 , inMemory = TRUE, nameFolder = name.BIOMOD_DATA)
  }
  
  #### fix models names ---------------------------------------------------------
  name.OUT <- paste0(EM@sp.name, '.', EM@modeling.id, '.ensemble.models.out')
  EM@link <- file.path(EM@dir.name, EM@sp.name, name.OUT)
  assign(x = name.OUT, value = EM)
  save(list = name.OUT, file = EM@link)
  
  .bm_cat("Done")
  return(EM)
}


# Argument check function  -----------------------------------------------------

.BIOMOD_EnsembleModeling.check.args <- function(bm.mod,
                                                models.chosen,
                                                em.by,
                                                em.algo,
                                                metric.select,
                                                metric.select.thresh,
                                                metric.select.table,
                                                metric.eval,
                                                prob.ci.alpha,
                                                prob.mean.weight.decay,
                                                prob.mean, # deprecated
                                                prob.cv, # deprecated
                                                prob.ci, # deprecated
                                                prob.median, # deprecated
                                                committee.averaging, # deprecated
                                                prob.mean.weight) { # deprecated
  ## 1. Check bm.mod ----------------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
  
  ## 2. Check models.chosen ---------------------------------------------------
  if ( is.null(models.chosen) | (length(models.chosen) == 1 && models.chosen[1] == 'all')) {
    cat("\n   ! all models available will be included in ensemble.modeling")
    models.chosen <- bm.mod@models.computed
  } else {
    .fun_testIfIn(TRUE, "models.chosen", models.chosen, bm.mod@models.computed)
  }
  
  # 3. check argument em.algo ----------------------------------------------
  em.avail.check <- c('prob.mean', 'prob.cv', 'prob.ci',
                      'prob.median', 'committee.averaging', 'prob.mean.weight')
  em.avail <- c('prob.mean', 'prob.cv', 'prob.ci.inf', 'prob.ci.sup',
                'prob.median', 'committee.averaging', 'prob.mean.weight')
  
  if (missing(em.algo)) {
    if (all(c(missing(prob.mean), missing(prob.cv),
              missing(prob.ci), missing(prob.median),
              missing(committee.averaging), missing(prob.mean.weight)))) {
      em.algo <- 'prob.mean'
      cat("\n! setting em.algo to its default value c('prob.mean')")
    } else {
      if (missing(prob.mean)) prob.mean <- FALSE
      if (missing(prob.cv)) prob.cv <- FALSE
      if (missing(prob.ci)) prob.ci <- FALSE
      if (missing(prob.median)) prob.median <- FALSE
      if (missing(committee.averaging)) committee.averaging <- FALSE
      if (missing(prob.mean.weight)) prob.mean.weight <- FALSE
      
      if (!is.logical(prob.mean) | !is.logical(prob.median) |
          !is.logical(prob.cv) | !is.logical(prob.ci) | 
          !is.logical(committee.averaging) | !is.logical(prob.mean.weight)) {
        stop("prob.mean, prob.cv, prob.ci, prob.median, committee.averaging and prob.mean.weight arguments must be logical")
      }
      
      em.algo <- em.avail[c(prob.mean, prob.cv,  prob.ci,  prob.ci,
                            prob.median, committee.averaging, prob.mean.weight)]
      
      em.algo.user <- em.avail.check[c(prob.mean, prob.cv, 
                                       prob.ci, prob.median,
                                       committee.averaging, prob.mean.weight)]
      # warning for obsolete arguments
      cat(paste0("\n   ! arguments ", paste0(em.avail.check[-6], collapse = ", "),
                 " and ",em.avail.check[6], " are obsolete. Please use `em.algo = c('",paste0(em.algo.user, collapse = "', '"),"')` instead."))
    }
  } else {
    .fun_testIfIn(TRUE, "em.algo", em.algo, em.avail.check)
    em.algo <- unique(em.algo)
    testCI <- grepl(pattern = "prob.ci", x = em.algo)
    if(any(testCI)){
      em.algo <- em.algo[-which(testCI)]
      em.algo <- c(em.algo, 'prob.ci.inf', 'prob.ci.sup')
    }
  }
  
  
  ## 4. Check metric.select ---------------------------------------------------
  if (!is.null(metric.select)) {
    if (!is.character(metric.select)) {
      stop("metric.select must be a character vector or NULL")
    }
    metric.select.user = FALSE
    if ('user.defined' %in% metric.select) {
      metric.select.user = TRUE
      if (!is.null(metric.select.table)) {
        .fun_testIfIn(TRUE, "models.chosen", models.chosen, colnames(metric.select.table))
        metric.select.table <- metric.select.table[, models.chosen, drop = FALSE]
        metric.select <- rownames(metric.select.table)
      } else {
        stop("metric.select.table must be a data.frame or NULL")
      }
    } else {
      if ('all' %in% metric.select) {
        metric.select <- unique(get_evaluations(bm.mod)$Metric.eval)
      }
      .fun_testIfIn(TRUE, "metric.select", metric.select, unique(get_evaluations(bm.mod)$Metric.eval))
    }
  }
  
  ## 5. Check metric.select.thresh --------------------------------------------
  if (!is.null(metric.select)) {
    if (!is.null(metric.select.thresh)) {
      if (!is.numeric(metric.select.thresh)) {
        stop("metric.select.thresh must be NULL or a numeric vector")
      }
      if (length(metric.select) != length(metric.select.thresh)) {
        stop("you must specify as many metric.select.thresh as metric.select (if you specify some)")
      }
      cat("\n   > Evaluation & Weighting methods summary :\n")
      cat(paste(metric.select, metric.select.thresh, sep = " over ", collapse = "\n      ")
          , fill = TRUE, labels = "     ")
    } else {
      cat("\n   ! No metric.select.thresh -> All models will be kept for Ensemble Modeling")
      metric.select.thresh <- rep(0, length(metric.select))
    }
  } else {
    metric.select <- 'none'
  }
  
  ## 6. Check metric.eval -----------------------------------------------------
  metric.eval <- unique(metric.eval)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC')
  .fun_testIfIn(TRUE, "metric.eval", metric.eval, avail.eval.meth.list)
  
  ## 7. Check selected EM algo ------------------------------------------------
  
  if (is.null(metric.select) && 
      any(c("committee.averaging", "prob.mean.weight") %in% em.algo)) {
    stop("You must choose metric.select if you want to compute Committee Averaging or Probability Weighted Mean algorithms")
  }
  
  ## 7.1 Check alpha for Confident interval
  if ("prob.ci" %in% em.algo) {
    .fun_testIfPosNum(TRUE, "prob.ci.alpha", prob.ci.alpha)
    if (prob.ci.alpha <= 0 | prob.ci.alpha >= 0.5) {
      stop("prob.ci.alpha must be a numeric between 0 and 0.5")
    }
  }
  
  ## 7.2 Check decay for wmean
  if ("prob.mean.weight" %in% em.algo) {
    if ((!is.numeric(prob.mean.weight.decay) &&
         !is.character(prob.mean.weight.decay) &&
         !is.function(prob.mean.weight.decay)) ||
        (is.numeric(prob.mean.weight.decay) && prob.mean.weight.decay < 0) ||
        (is.character(prob.mean.weight.decay) && prob.mean.weight.decay != 'proportional')) {
      stop("'prob.mean.weight.decay' should be either 'proportional', a numeric value > 0 or a function")
    }
  }
  
  ## 8. Check em.by -----------------------------------------------------------
  if(length(em.by) != 1){
    stop("\nem.by should be of length 1")
  }
  .fun_testIfIn(TRUE, "em.by", em.by, c('PA_dataset', 'algo', 'all', 'PA_dataset+repet', 'PA_dataset+algo'))
  
  # check that repetition are note merged with full models
  if(any(grepl(pattern = "RUN",  x = models.chosen)) &&
     any(grepl(pattern = "Full", x = models.chosen)) &&
     em.by != 'PA_dataset+repet') {
    cat("\n!!! Removed models using the Full dataset as ensemble models cannot merge repetition dataset (RUN1, RUN2, ...) with Full dataset unless em.by = 'PA_dataset+repet'.")
    models.chosen <- models.chosen[!grepl(pattern = "Full", x = models.chosen)]
  }
  
  ## 9. Check that ensemble model have > 1 model to run -------------
  ## make a list of models names that will be combined together according to em.by argument
  em.mod.assemb <- .get_models_assembling(models.chosen, em.by)
  ### Check that all EM have > 1 model selected ----------------------------
  out.check <- foreach(assemb = names(em.mod.assemb), .combine = 'rbind') %do% {
    out <- .get_kept_models(
      bm.mod, em.mod.assemb[[assemb]], 
      metric.select, metric.select.thresh,
      metric.select.user, metric.select.table
    )$models.kept
    data.frame(models.kept = sapply(out, length),
               metric.select = names(out),
               assemb = assemb,
               row.names = NULL)
  }
  
  for (thismetric in metric.select)  {
    out.check.sub <-
      out.check[which(out.check$metric.select == thismetric),]
    assemb.1 <- 
      out.check.sub[which(out.check.sub$models.kept == 1), "assemb"]
    assemb.0 <- 
      out.check.sub[which(out.check.sub$models.kept == 0), "assemb"]
    
    if(length(assemb.0) > 0 || length(assemb.1) > 0){
      cat("\n")
      if(length(assemb.0) > 0){
        cat("\n     !! Ensemble Model", assemb.0, "selected by", thismetric, "have no model selected and will fail.")
      }   
      if(length(assemb.1) > 0){
        cat("\n     !! Ensemble Model", assemb.1, "selected by", thismetric, "have only one single model selected.")
      }
      cat("\n     !! Please make sure this is intended or review your selection metrics and threshold.")
    }
  }
  
  
  
  
  ## Return Args ------------------------------------------------
  return(list(bm.mod = bm.mod,
              models.chosen = models.chosen,
              em.algo = em.algo,
              metric.select = metric.select,
              metric.select.thresh = metric.select.thresh,
              metric.select.user = metric.select.user,
              metric.select.table = metric.select.table,
              metric.eval = metric.eval,
              prob.ci.alpha = prob.ci.alpha,
              prob.mean.weight.decay = prob.mean.weight.decay,
              em.by = em.by,
              em.mod.assemb = em.mod.assemb))
  
}


# .get_models_assembling --------------------------------------------------

.get_models_assembling <- function(models.chosen, em.by)
{
  assembl.list = list()
  if (em.by == 'PA_dataset') {
    for (dat in .extract_modelNamesInfo(models.chosen, info = 'data.set')) {
      assembl.list[[paste0("mergedAlgo_mergedRun_", dat)]] <- models.chosen[grep(paste0("_", dat, "_"), models.chosen)]
    }
  } else if (em.by == 'algo') {
    for (algo in .extract_modelNamesInfo(models.chosen, info = 'models')) {
      # grep use a regexp to match end of model name with algo name
      assembl.list[[paste0(algo, "_mergedRun_mergedData")]] <- 
        models.chosen[grep(paste0("*\\_", algo,"$"), models.chosen)]
    }
  } else if (em.by == 'all') {
    assembl.list[["mergedAlgo_mergedRun_mergedData"]] <- models.chosen
  } else if (em.by == 'PA_dataset+repet') {
    for (dat in .extract_modelNamesInfo(models.chosen, info = 'data.set')) {
      for (repet in .extract_modelNamesInfo(models.chosen, info = 'run.eval')) {
        mod.tmp <- intersect(x = grep(paste0("_", dat, "_"), models.chosen)
                             , y = grep(paste0("_", repet, "_"), models.chosen))
        if (length(mod.tmp) > 0) {
          assembl.list[[paste0("mergedAlgo_", repet, "_", dat)]] <- models.chosen[mod.tmp]
        }
      }
    }
  } else if (em.by == 'PA_dataset+algo') {
    for (dat in .extract_modelNamesInfo(models.chosen, info = 'data.set')) {
      for (algo in .extract_modelNamesInfo(models.chosen, info = 'models')) {
        # grep use a regexp to match end of model name with algo name
        mod.tmp <- intersect(x = grep(paste0("_", dat, "_"), models.chosen)
                             , y = grep(paste0("*\\_", algo,"$"), models.chosen))
        if (length(mod.tmp) > 0) {
          assembl.list[[paste0(algo, "_mergedRun_", dat)]] <- models.chosen[mod.tmp]
        }
      }
    }
  }
  return(assembl.list)
}


# .get_needed_predictions -------------------------------------------------

.get_needed_predictions <- function(bm.mod, em.by,  models.kept, metric.select
                                    , metric.select.thresh, metric.select.user
                                    , metric.select.table, nb.cpu) {
  
  out <- .get_kept_models(bm.mod, models.kept, 
                          metric.select, metric.select.thresh,
                          metric.select.user, metric.select.table)
  
  models.kept.union <- unique(unlist(out$models.kept))
  
  if (length(models.kept.union) > 0) {
    ## load prediction on each PA_dataset
    if (em.by %in% c("PA_dataset", 'PA_dataset+algo', 'PA_dataset+repet') || 
        !inherits(get_formal_data(bm.mod), "BIOMOD.formated.data.PA") ||
        ncol(get_formal_data(bm.mod)@PA.table) == 1) {
      out$predictions <- get_predictions(bm.mod, full.name = models.kept.union)
    } else {
      ## only for models with pseudo absence
      ## some prediction should be added for the union of pseudo-absences
      
      cat("\n   ! Additional projection required for ensemble models merging several pseudo-absence dataset...")
      
      # temp folders for additionnal predictions
      temp_name <- paste0('tmp_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE))
      PA.table <- get_formal_data(bm.mod)@PA.table
      # data that are kept at least in one PA dataset
      kept_data <- apply(PA.table, 1, any)
      # tells which PA dataset is used by which model
      models.kept.PA <-   sapply(models.kept.union, function(x){
        .extract_modelNamesInfo(x, info = "data.set")
      })
      
      out$predictions <-
        foreach(thisPA = unique(models.kept.PA), .combine = "rbind") %do% {
          ## model kept for this PA dataset
          thismodels <- names(models.kept.PA)[which(models.kept.PA == thisPA)]
          ## index of data to predict and data already predicted
          index_to_predict <- which(!PA.table[, thisPA] & kept_data)
          index_current <- which(PA.table[, thisPA])
          
          ## retrieve predictions for this PA dataset
          current_prediction <- get_predictions(bm.mod, full.name = thismodels)
          current_prediction$Points <- index_current
          
          # subsetting environment and coord
          env_to_predict <- get_formal_data(bm.mod)@data.env.var[index_to_predict,]
          coord_to_predict <- get_formal_data(bm.mod)@coord[index_to_predict,]
          
          # prediction on the other PA datasets
          new_prediction <-
            get_predictions(
              BIOMOD_Projection(
                bm.mod = bm.mod,
                new.env = env_to_predict,
                proj.name = temp_name,
                xy.new.env = coord_to_predict,
                models.chosen = thismodels,
                compress = TRUE,
                build.clamping.mask = FALSE,
                do.stack = TRUE,
                nb.cpu = nb.cpu
              ))
          new_prediction$Points <- index_to_predict
          
          ## combining old and new predictions
          index_full <- c(index_current, index_to_predict)
          
          res <- rbind(current_prediction, new_prediction)
          res <- res[, c("full.name", "data.set", "run.eval", "algo", "Points", "pred")]
          res <- res[order(res$full.name, res$Points), ]
          return(res)
        }
      
      # delete temporary directory
      unlink(file.path(bm.mod@dir.name, bm.mod@sp.name, paste0("proj_", temp_name))
             , recursive = TRUE, force = TRUE)
      cat("\n")
    }
    return(out)
  } else {
    cat("\n   ! No models kept due to threshold filtering... Ensemble Modeling was skipped!")
    return(NULL)
  }
}



.get_kept_models <- function(bm.mod, models.kept, 
                             metric.select, metric.select.thresh,
                             metric.select.user, metric.select.table){
  out <- list(predictions = NULL, models.kept = NULL, models.kept.scores = NULL)
  for (eval.m in metric.select) {
    if (eval.m != 'none') {
      if (metric.select.user) {
        models.kept.scores <- metric.select.table[eval.m, models.kept]
      } else {
        models.kept.scores <- unlist(lapply(models.kept, function(x) {
          mod <- tail(unlist(strsplit(x, "_")), 3)[3]
          run <- tail(unlist(strsplit(x, "_")), 3)[2]
          dat <- tail(unlist(strsplit(x, "_")), 3)[1]
          # select evaluations scores obtained for Evaluation Data if exists or CV if not
          out <- get_evaluations(bm.mod, data.set = dat, run.eval = run, algo = mod, Metric.eval = eval.m)
          if (bm.mod@has.evaluation.data) {
            return(out[, "Evaluating.data"])
          } else {
            return(out[, "Testing.data"])
          }
        }))
      }
      ## set NA to -1
      if (!is.null(models.kept.scores)) {
        models.kept.scores[is.na(models.kept.scores)] <- -1
      }
      thresh = metric.select.thresh[which(metric.select == eval.m)]
      out$models.kept[[eval.m]] <- models.kept[models.kept.scores > thresh]
      out$models.kept.scores[[eval.m]] <- models.kept.scores[models.kept.scores > thresh]
    } else {
      out$models.kept[[eval.m]] <- models.kept
    }
  }
  out
}

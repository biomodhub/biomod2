# BIOMOD_ModelingOptions documentation -----------------------------------------
##' @name BIOMOD_ModelingOptions
##' @aliases BIOMOD_ModelingOptions
##' @aliases bm_DefaultModelingOptions
##' @author Damien Georges, Wilfried Thuiller
##' 
##' @title Configure the modeling options for each selected model
##'
##' @description Parametrize and/or tune \pkg{biomod2}'s single models options.
##'
##'
##' @param GLM (\emph{optional, default} \code{NULL}) \cr A \code{list} containing GLM options
##' @param GBM (\emph{optional, default} \code{NULL}) \cr A \code{list} containing GBM options
##' @param GAM (\emph{optional, default} \code{NULL}) \cr A \code{list} containing GAM options
##' @param CTA (\emph{optional, default} \code{NULL}) \cr A \code{list} containing CTA options
##' @param ANN (\emph{optional, default} \code{NULL}) \cr A \code{list} containing ANN options
##' @param SRE (\emph{optional, default} \code{NULL}) \cr A \code{list} containing SRE options
##' @param FDA (\emph{optional, default} \code{NULL}) \cr A \code{list} containing FDA options
##' @param MARS (\emph{optional, default} \code{NULL}) \cr A \code{list} containing MARS options
##' @param RF (\emph{optional, default} \code{NULL}) \cr A \code{list} containing RF options
##' @param MAXENT (\emph{optional, default} \code{NULL}) \cr A \code{list} 
##' containing MAXENT options
##'
##'
##' @return 
##' 
##' A \code{\link{BIOMOD.models.options}} object that can be used to build species distribution 
##' model(s) with the \code{\link{BIOMOD_Modeling}} function.
##' 
##' 
##' @details
##' 
##' This function allows advanced user to change some default parameters of \pkg{biomod2} inner 
##' models. \cr 10 single models are available within the package, and their options can be set 
##' with this function through \code{list} objects.
##' 
##' The \code{\link{bm_DefaultModelingOptions}} function prints all default parameter values for 
##' all available models. \cr This output can be copied and pasted to be used as is (with wanted 
##' changes) as function arguments (see \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}).
##' 
##' Below is the detailed list of all modifiable parameters for each available model.
##'
##' @section GLM : (\code{\link[stats]{glm}})
##' \itemize{
##'   \item{\code{myFormula}}{ : a typical \code{formula} object (see 
##'   \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}). \cr If not \code{NULL}, \code{type} 
##'   and \code{interaction.level} parameters are switched off. \cr 
##'   You can choose to either :
##'   \itemize{
##'     \item{generate automatically the GLM formula with the following parameters :
##'     \itemize{
##'       \item{\code{type = 'quadratic'}}{ : formula given to the model, must 
##'       be \code{simple}, \code{quadratic} or \code{polynomial}}
##'       \item{\code{interaction.level = 0}}{ : an \code{integer} corresponding 
##'       to the interaction level between considered variables considered (\emph{be aware that 
##'       interactions quickly enlarge the number of effective variables used into the GLM !})}
##'     }
##'     }
##'     \item{or construct specific formula}
##'   }
##'   }
##'   \item{\code{test = 'AIC'}}{ : information criteria for the stepwise 
##'   selection procedure, must be \code{AIC} (\emph{Akaike Information Criteria}, \code{BIC} 
##'   (\emph{Bayesian Information Criteria}) or \code{none} (\emph{consider only the full model, 
##'   no stepwise selection, but this can lead to convergence issue and strange results !})}
##'   \item{\code{family = binomial(link = 'logit')}}{ : a \code{character} 
##'   defining the error distribution and link function to be used in the model, mus be a family 
##'   name, a family function or the result of a call to a family function (see \link{family}) 
##'   (\emph{so far, \pkg{biomod2} only runs on presence-absence data, so binomial family is the 
##'   default !})}
##'   \item{\code{control}}{ : a \code{list} of parameters to control the fitting process (passed to 
##'   \code{\link{glm.control}})}
##' }
##'
##' @section GBM : (default \code{\link[gbm]{gbm}})
##' 
##' \emph{Please refer to \code{\link[gbm]{gbm}} help file for more details.}
##' \itemize{
##'   \item{\code{distribution = 'bernoulli'}}
##'   \item{\code{n.trees = 2500}}
##'   \item{\code{interaction.depth = 7}}
##'   \item{\code{n.minobsinnode = 5}}
##'   \item{\code{shrinkage = 0.001}}
##'   \item{\code{bag.fraction = 0.5}}
##'   \item{\code{train.fraction = 1}}
##'   \item{\code{cv.folds = 3}}
##'   \item{\code{keep.data = FALSE}}
##'   \item{\code{verbose = FALSE}}
##'   \item{\code{perf.method = 'cv'}}
##'   \item{\code{n.cores = 1}}
##' }
##'
##' @section GAM : (\code{\link[gam]{gam}} or \code{\link[mgcv]{gam}})
##' \itemize{
##'   \item{\code{algo = 'GAM_gam'}}{ : a \code{character} defining the chosen GAM function, must 
##'   be \code{GAM_gam} (see \code{\link[gam]{gam}}), \code{GAM_mgcv} (see \code{\link[mgcv]{gam}}) 
##'   or \code{BAM_mgcv} (see \code{\link[mgcv]{bam}})}
##'   \item{\code{myFormula}}{ : a typical \code{formula} object (see 
##'   \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}). \cr If not \code{NULL}, \code{type} 
##'   and \code{interaction.level} parameters are switched off. \cr 
##'   You can choose to either :
##'   \itemize{
##'     \item{generate automatically the GAM formula with the following parameters :
##'     \itemize{
##'       \item{\code{type = 's_smoother'}}{ : the smoother used to generate the formula}
##'       \item{\code{interaction.level = 0}}{ : an \code{integer} corresponding 
##'       to the interaction level between considered variables considered (\emph{be aware that 
##'       interactions quickly enlarge the number of effective variables used into the GLM !})}
##'     }
##'     }
##'     \item{or construct specific formula}
##'   }
##'   }
##'   \item{\code{k = -1}}{a smooth term in a formula argument to gam, must be \code{-1} or 
##'   \code{4} (see \pkg{gam} \code{\link[gam]{s}} or \pkg{mgcv} \code{\link[mgcv]{s}})}
##'   \item{\code{family = binomial(link = 'logit')}}{ : a \code{character} defining 
##'   the error distribution and link function to be used in the model, mus be a family name, a 
##'   family function or the result of a call to a family function (see \link{family}) 
##'   (\emph{so far, \pkg{biomod2} only runs on presence-absence data, so binomial family is the 
##'   default !})}
##'   \item{\code{control}}{ : a \code{list} of parameters to control the fitting process (passed to 
##'   \code{\link[mgcv]{gam.control}} or \code{\link[gam]{gam.control}})}
##'   \item{some options specific to \code{GAM_mgcv} (\emph{ignored if \code{algo = 'GAM_gam'}})}{
##'   \itemize{
##'     \item{\code{method = 'GCV.Cp'})}
##'     \item{\code{optimizer = c('outer','newton')}}
##'     \item{\code{select = FALSE}}
##'     \item{\code{knots = NULL}}
##'     \item{\code{paramPen = NULL}}
##'   }
##'   }
##' }
##'
##' @section CTA : (\code{\link[rpart]{rpart}})
##' 
##' \emph{Please refer to \code{\link[rpart]{rpart}} help file for more details.}
##' \itemize{
##'   \item{\code{method = 'class'}}
##'   \item{\code{parms = 'default'}}{ : if \code{'default'}, default \pkg{rpart} 
##'   \code{parms} value are kept}
##'   \item{\code{cost = NULL}}
##'   \item{\code{control}}{ : see \code{\link[rpart]{rpart.control}}}
##' }
##'
##' @section ANN : (\code{\link[nnet]{nnet}})
##' \itemize{
##'   \item{\code{NbCV = 5}}{ : an \code{integer} corresponding to the number of cross-validation 
##'   repetitions to find best size and decay parameters}
##'   \item{\code{size = NULL}}{ : an \code{integer} corresponding to the number of units in the 
##'   hidden layer. If \code{NULL} then size parameter will be optimized by cross-validation based 
##'   on model AUC (\code{NbCv} cross-validations ; tested size will be the following : 
##'   \code{c(2, 4, 6, 8)}). It is also possible to give a \code{vector} of size values to be tested, 
##'   and the one giving the best model AUC will be kept.}
##'   \item{\code{decay = NULL}}{ : a \code{numeric} corresponding to weight decay. If \code{NULL} 
##'   then decay parameter will be optimized by cross-validation based on model AUC (\code{NbCv} 
##'   cross-validations ; tested size will be the following : \code{c(0.001, 0.01, 0.05, 0.1)}). 
##'   It is also possible to give a \code{vector} of decay values to be tested, and the one giving 
##'   the best model AUC will be kept.}
##'   \item{\code{rang = 0.1}}{ : a \code{numeric} corresponding to the initial random weights on 
##'   \code{[-rang, rang]}}
##'   \item{\code{maxit = 200}}{ : an \code{integer} corresponding to the maximum number of 
##'   iterations}
##' }
##'
##' @section SRE : (\code{\link{bm_SRE}})
##' \itemize{
##'   \item{\code{quant = 0.025}}{ : a \code{numeric} corresponding to the quantile of 
##'   '\emph{extreme environmental variable}' removed to select species envelops}
##' }
##'
##' @section FDA : (\code{\link[mda]{fda}})
##' 
##' \emph{Please refer to \code{\link[mda]{fda}} help file for more details.}
##' \itemize{
##'   \item{\code{method = 'mars'}}
##'   \item{\code{add_args = NULL}}{ : a \code{list} of additional parameters to \code{method} and 
##'   given to the \code{...} options of \code{\link[mda]{fda}} function}
##' }
##'
##' @section MARS : (\code{\link[earth]{earth}})
##' 
##' \emph{Please refer to \code{\link[earth]{earth}} help file for more details.}
##' \itemize{
##'   \item{\code{myFormula}}{ : a typical \code{formula} object (see 
##'   \href{https://biomodhub.github.io/biomod2/reference/BIOMOD_ModelingOptions.html#examples}{Examples}). \cr If not \code{NULL}, \code{type} 
##'   and \code{interaction.level} parameters are switched off. \cr 
##'   You can choose to either :
##'   \itemize{
##'     \item{generate automatically the MARS formula with the following parameters :
##'     \itemize{
##'       \item{\code{type = 'simple'}}{ : formula given to the model, must 
##'       be \code{simple}, \code{quadratic} or \code{polynomial}}
##'       \item{\code{interaction.level = 0}}{ : an \code{integer} corresponding 
##'       to the interaction level between considered variables considered (\emph{be aware that 
##'       interactions quickly enlarge the number of effective variables used into the MARS !})}
##'     }
##'     }
##'     \item{or construct specific formula}
##'   }
##'   }
##'   \item{\code{nk = NULL}}{ : an \code{integer} corresponding to the maximum number of model 
##'   terms. \cr 
##'   If \code{NULL} default MARS function value is used : \code{max(21, 2 * nb_expl_var + 1)}}
##'   \item{\code{penalty = 2}}
##'   \item{\code{thresh = 0.001}}
##'   \item{\code{nprune = NULL}}
##'   \item{\code{pmethod = 'backward'}}
##' }
##'
##' @section RF : (\code{\link[randomForest]{randomForest}})
##' \itemize{
##'   \item{\code{do.classif = TRUE}}{ : if \code{TRUE} \emph{random.forest classification} will 
##'   be computed, otherwise \emph{random.forest regression} will be done}
##'   \item{\code{ntree = 500}}
##'   \item{\code{mtry = 'default'}}
##'   \item{\code{sampsize = NULL}}
##'   \item{\code{nodesize = 5}}
##'   \item{\code{maxnodes = NULL}}
##' }
##'
##' @section MAXENT : (\url{https://biodiversityinformatics.amnh.org/open_source/maxent/})
##' \itemize{
##'   \item{\code{path_to_maxent.jar = getwd()}}{ : a \code{character}
##'   corresponding to \pkg{maxent.jar} file link} 
##'   
##'   \item{\code{memory_allocated = 512}}{ : an \code{integer} corresponding to
##'   the amount of memory (in Mo) reserved for \code{java} to run
##'   \code{MAXENT}, must be \code{64}, \code{128}, \code{256},
##'   \code{512}, \code{1024}... or \code{NULL} to use default \code{java}
##'   memory limitation parameter}
##'   
##'   \item{\code{initial_heap_size = NULL}}{ : a \code{character} initial heap
##'   space (shared memory space) allocated to java. Argument transmitted to
##'   \code{-Xms} when calling java. Used in \code{\link{BIOMOD_Projection}} but
##'   not in \code{\link{BIOMOD_Modeling}}. Values can be \code{1024K},
##'   \code{4096M}, \code{10G} ... or \code{NULL} to use default \code{java}
##'   parameter}
##'   
##'   \item{\code{max_heap_size = NULL}}{ : a \code{character} initial heap
##'   space (shared memory space) allocated to java. Argument transmitted to
##'   \code{-Xmx} when calling java. Used in \code{\link{BIOMOD_Projection}} but
##'   not in \code{\link{BIOMOD_Modeling}}. Must be larger than
##'   \code{initial_heap_size}. Values can be \code{1024K}, \code{4096M},
##'   \code{10G} ... or \code{NULL} to use default \code{java} parameter}
##'   
##'   \item{\code{background_data_dir}}{ : a \code{character} corresponding to
##'   directory path where explanatory variables are stored as \code{ASCII}
##'   files (raster format). If specified, \code{MAXENT} will generate
##'   its own background data from explanatory variables rasters (as usually
##'   done in \code{MAXENT} studies). Otherwise \pkg{biomod2} pseudo-absences
##'   will be used (see \code{\link{BIOMOD_FormatingData}})}
##'   
##'   \item{\code{maximumbackground}}{ : an \code{integer} corresponding to the
##'   maximum number of background data to sample if the
##'   \code{background_data_dir} parameter has been set}
##'   
##'   \item{\code{maximumiterations = 200}}{ : an \code{integer} corresponding
##'   to the maximum number of iterations to do} 
##'   
##'   \item{\code{visible = FALSE}}{ : a \code{logical} to make the
##'   \code{MAXENT} user interface available}
##'   
##'   \item{\code{linear = TRUE}}{ : a \code{logical} to allow linear features
##'   to be used} \item{\code{quadratic = TRUE}}{ : a \code{logical} to allow
##'   quadratic features to be used} 
##'   
##'   \item{\code{product = TRUE}}{ : a \code{logical} to allow product features
##'   to be used}
##'   
##'   \item{\code{threshold = TRUE}}{ : a \code{logical} to allow threshold
##'   features to be used}
##'   
##'   \item{\code{hinge = TRUE}}{ : a \code{logical} to allow hinge features to
##'   be used} \item{\code{lq2lqptthreshold = 80}}{ : an \code{integer}
##'   corresponding to the number of samples at which product and threshold
##'   features start being used} 
##'   
##'   \item{\code{l2lqthreshold = 10}}{ : an
##'   \code{integer} corresponding to the number of samples at which quadratic
##'   features start being used} 
##'   
##'   \item{\code{hingethreshold = 15}}{ : an
##'   \code{integer} corresponding to the number of samples at which hinge
##'   features start being used}
##'   
##'   \item{\code{beta_threshold = -1.0}}{ : a
##'   \code{numeric} corresponding to the regularization parameter to be applied
##'   to all threshold features (\emph{negative value enables automatic
##'   setting})} 
##'   
##'   \item{\code{beta_categorical = -1.0}}{ : a \code{numeric}
##'   corresponding to the regularization parameter to be applied to all
##'   categorical features (\emph{negative value enables automatic setting})}
##'   
##'   \item{\code{beta_lqp = -1.0}}{ : a \code{numeric} corresponding to the
##'   regularization parameter to be applied to all linear, quadratic and
##'   product features (\emph{negative value enables automatic setting})}
##'   
##'   \item{\code{beta_hinge = -1.0}}{ : a \code{numeric} corresponding to the
##'   regularization parameter to be applied to all hinge features
##'   (\emph{negative value enables automatic setting})}
##'   
##'   \item{\code{betamultiplier = 1}}{ : a \code{numeric} to multiply all
##'   automatic regularization parameters \cr (\emph{higher number gives a more
##'   spread-out distribution})} 
##'   
##'   \item{\code{defaultprevalence = 0.5}}{ : a
##'   \code{numeric} corresponding to the default prevalence of the species \cr
##'   (\emph{probability of presence at ordinary occurrence points})}
##' }
##'
##' % @section \bold{MAXENT.Tsuruoka (\code{\link[maxent]{maxent}})} :
##' % \itemize{
##' %   \item{\code{l1_regularizer = 0.0}}{ : a \code{numeric} turning on L1 regularization and setting 
##' %   the regularization parameter (\emph{a value of \code{0} will disable L1 regularization})}
##' %   \item{\code{l2_regularizer = 0.0}}{ : a \code{numeric} turning on L2 regularization and setting 
##' %   the regularization parameter (\emph{a value of \code{0} will disable L2 regularization})}
##' %   \item{\code{use_sgd = FALSE}}{ : a \code{logical} to use SGD parameter estimation}
##' %   \item{\code{set_heldout = 0}}{ : an \code{integer} corresponding to the number of documents to 
##' %   hold out (\emph{to test against and prevent overfitting, use carefully in case of dataset with 
##' %   low number of occurrences})}
##' %   \item{\code{verbose = FALSE}}{ : a \code{logical} specifying whether to provide descriptive 
##' %   output about the training process}
##' % }
##'
##'
##' @keywords models options
##' 
##' 
##' @seealso \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_Modeling}}
##' @family Main functions
##' 
##' 
##' @examples
##' library(terra)
##' 
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
##' # ---------------------------------------------------------------#
##' # Print default modeling options
##' bm_DefaultModelingOptions()
##' 
##' # Create default modeling options
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##' myBiomodOptions
##' 
##' # # Part (or totality) of the print can be copied and customized
##' # # Below is an example to compute quadratic GLM and select best model with 'BIC' criterium
##' # myBiomodOptions <- BIOMOD_ModelingOptions(
##' #   GLM = list(type = 'quadratic',
##' #              interaction.level = 0,
##' #              myFormula = NULL,
##' #              test = 'BIC',
##' #              family = 'binomial',
##' #              control = glm.control(epsilon = 1e-08,
##' #                                    maxit = 1000,
##' #                                    trace = FALSE)))
##' # myBiomodOptions
##' # 
##' # # It is also possible to give a specific GLM formula
##' # myForm <- 'Sp277 ~ bio3 + log(bio10) + poly(bio16, 2) + bio19 + bio3:bio19'
##' # myBiomodOptions <- BIOMOD_ModelingOptions(GLM = list(myFormula = formula(myForm)))
##' # myBiomodOptions
##'
##'
## @importFrom gam gam.control
## @importFrom mgcv gam.control
##' @importFrom methods as new validObject
##' 
##' 
##' @export
##'
##'
## -------------------------------------------------------------------------- ##

BIOMOD_ModelingOptions <- function(GLM = NULL,
                                   GBM = NULL,
                                   GAM = NULL,
                                   CTA = NULL,
                                   ANN = NULL,
                                   SRE = NULL,
                                   FDA = NULL,
                                   MARS = NULL,
                                   RF = NULL,
                                   MAXENT = NULL)
{
  # .bm_cat("Build Modeling Options")
  
  ## 1. create a defaut BIOMOD.models.options object --------------------------
  opt <- new('BIOMOD.models.options')
  
  ## 2. modify it if necessary ------------------------------------------------
  fam_GLM = fam_GAM = c('binomial', 'gaussian', 'Gamma', 'inverse.gaussian'
                        , 'poisson', 'quasi', 'quasibinomial', 'quasipoisson')
  
  ## 2.1 GLM ------------------------------------------------------------------
  if (!is.null(GLM)) {
    if (!is.null(GLM$type)) { opt@GLM$type <- GLM$type }
    if (!is.null(GLM$interaction.level)) { opt@GLM$interaction.level <- GLM$interaction.level }
    if (!is.null(GLM$myFormula)) { opt@GLM$myFormula <- GLM$myFormula }
    if (!is.null(GLM$test)) { opt@GLM$test <- GLM$test }
    if (!is.null(GLM$family)) {
      fam.test <- TRUE
      if (inherits(GLM$family, 'family')) {
        opt@GLM$family <- GLM$family
      } else if (is.character(GLM$family)) {
        if (! unlist(strsplit(GLM$family, "[/(]"))[1] %in% fam_GLM) { fam.test <- FALSE }
        if (grepl(')', GLM$family)) { # check string formalisation to add () if necessary
          opt@GLM$family <- eval(parse(text = GLM$family))
        } else {
          opt@GLM$family <- eval(parse(text = paste0(GLM$family, "()")))
        }
      } else{ fam.test <- FALSE }
      if (!fam.test) {
        cat("\n!!! invalid GLM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GLM$family <- binomial(link = 'logit')
      }
    }
    if (!is.null(GLM$mustart)) { opt@GLM$mustart <- GLM$mustart }
    if (!is.null(GLM$control)) { opt@GLM$control <- GLM$control }
  }
  
  ## 2.2 GBM ------------------------------------------------------------------
  if (!is.null(GBM)){
    if (!is.null(GBM$distribution)) { opt@GBM$distribution <- GBM$distribution }
    if (!is.null(GBM$n.trees)) { opt@GBM$n.trees <- GBM$n.trees }
    if (!is.null(GBM$interaction.depth)) { opt@GBM$interaction.depth <- GBM$interaction.depth }
    if (!is.null(GBM$n.minobsinnode)) { opt@GBM$n.minobsinnode <- GBM$n.minobsinnode }
    if (!is.null(GBM$shrinkage)) { opt@GBM$shrinkage <- GBM$shrinkage }
    if (!is.null(GBM$bag.fraction)) { opt@GBM$bag.fraction <- GBM$bag.fraction }
    if (!is.null(GBM$train.fraction)) { opt@GBM$train.fraction <- GBM$train.fraction }
    if (!is.null(GBM$cv.folds)) { opt@GBM$cv.folds <- GBM$cv.folds }
    if (!is.null(GBM$keep.data)) { opt@GBM$keep.data <- GBM$keep.data }
    if (!is.null(GBM$verbose)) { opt@GBM$verbose <- GBM$verbose }
    if (!is.null(GBM$perf.method)) { opt@GBM$perf.method <- GBM$perf.method }
    if (!is.null(GBM$n.cores)) { opt@GBM$n.cores <- GBM$n.cores } else { opt@GBM$n.cores <- NULL }
  }
  
  ## 2.3 GAM ------------------------------------------------------------------
  if (!is.null(GAM)) {
    if (!is.null(GAM$algo)) { opt@GAM$algo <- GAM$algo }
    if (!is.null(GAM$type)) { opt@GAM$type <- GAM$type }
    if (!is.null(GAM$k)) {
      opt@GAM$k <- GAM$k
    } else if (opt@GAM$algo == 'GAM_gam') {
      opt@GAM$k <- 4
    } else {
      opt@GAM$k <- -1
    }
    if (!is.null(GAM$interaction.level)) { opt@GAM$interaction.level <- GAM$interaction.level }
    if (!is.null(GAM$myFormula)) { opt@GAM$myFormula <- GAM$myFormula }
    if (!is.null(GAM$family)) {
      fam.test <- TRUE
      if (inherits(GAM$family, 'family')) {
        opt@GAM$family <- GAM$family
      } else if ( is.character(GAM$family)){
        if (! unlist(strsplit(GAM$family,"[/(]"))[1] %in% fam_GAM){ fam.test <- FALSE}
        if (grepl(')', GAM$family)) { # check string formalisation to add () if necessary
          opt@GAM$family <- eval(parse(text = GAM$family))
        } else {
          opt@GAM$family <- eval(parse(text = paste(GAM$family, "()")))
        }
      } else { fam.test <- FALSE }
      if (!fam.test) {
        cat("\n!!! invalid GAM$family given -> binomial(link = 'logit') was automatically set up !!!")
        opt@GAM$family <- binomial(link = 'logit')
      }
    }
    
    if (is.null(GAM$control)) {
      if (opt@GAM$algo == 'GAM_gam') {
        if(!requireNamespace('gam', quietly = TRUE)) stop("Package 'gam' not found")
        opt@GAM$control <- gam::gam.control()
      } else {
        if(!requireNamespace('mgcv', quietly = TRUE)) stop("Package 'mgcv' not found")
        opt@GAM$control <- mgcv::gam.control()
      }
    } else {
      user.control.list <- GAM$control
      if (opt@GAM$algo == 'GAM_gam') {
        if(!requireNamespace('gam', quietly = TRUE)) stop("Package 'gam' not found")
        default.control.list <- gam::gam.control()
      } else {
        if(!requireNamespace('mgcv', quietly = TRUE)) stop("Package 'mgcv' not found")
        default.control.list <- mgcv::gam.control()
      }
      control.list <- lapply(names(default.control.list), function(x) {
        if (x %in% names(user.control.list)) {
          return(user.control.list[[x]])
        } else {
          return(default.control.list[[x]])
        }
      })
      names(control.list) <- names(default.control.list)
      opt@GAM$control <- control.list
    }
    
    if (!is.null(GAM$method)) { opt@GAM$method <- GAM$method }
    if (!is.null(GAM$optimizer)) { opt@GAM$optimizer <- GAM$optimizer }
    if (!is.null(GAM$select)) { opt@GAM$select <- GAM$select }
    if (!is.null(GAM$knots)) { opt@GAM$knots <- GAM$knots }
    if (!is.null(GAM$paraPen)) { opt@GAM$paraPen <- GAM$paraPen }
  } else {
    if (opt@GAM$algo == 'GAM_gam') {
      opt@GAM$control <- gam::gam.control()
      opt@GAM$k <- 4
    } else {
      opt@GAM$control <- mgcv::gam.control()
      opt@GAM$k <- -1
    }
  }
  
  ## 2.4 CTA ------------------------------------------------------------------
  if (!is.null(CTA)) {
    if (!is.null(CTA$method)) { opt@CTA$method <- CTA$method }
    if (!is.null(CTA$parms)) { opt@CTA$parms <- CTA$parms }
    if (!is.null(CTA$control)) { opt@CTA$control <- CTA$control }
    if (!is.null(CTA$cost)) { opt@CTA$cost <- CTA$cost }
  }
  
  ## 2.5 ANN ------------------------------------------------------------------
  if (!is.null(ANN)) {
    if (!is.null(ANN$NbCV)) { opt@ANN$NbCV <- ANN$NbCV }
    if (!is.null(ANN$size)) { opt@ANN$size <- ANN$size }
    if (!is.null(ANN$decay)) { opt@ANN$decay <- ANN$decay }
    if (!is.null(ANN$rang)) { opt@ANN$rang <- ANN$rang }
    if (!is.null(ANN$maxit)) { opt@ANN$maxit <- ANN$maxit }
  }
  
  ## 2.6 SRE ------------------------------------------------------------------
  if (!is.null(SRE)) {
    if (!is.null(SRE$quant)) { opt@SRE$quant <- SRE$quant }
  }
  
  ## 2.7 FDA ------------------------------------------------------------------
  if (!is.null(FDA)) {
    if (!is.null(FDA$method)) { opt@FDA$method <- FDA$method }
    if (!is.null(FDA$add_args)) { opt@FDA$add_args <- FDA$add_args } ## additional args such as degree, nk
  }
  
  ## 2.8 MARS -----------------------------------------------------------------
  if (!is.null(MARS)) {
    if (!is.null(MARS$type)) { opt@MARS$type <- MARS$type }
    if (!is.null(MARS$interaction.level)) { opt@MARS$interaction.level <- MARS$interaction.level }
    if (!is.null(MARS$myFormula)) { opt@MARS$myFormula <- MARS$myFormula }
    if (!is.null(MARS$nk)) { opt@MARS$nk <- MARS$nk }
    if (!is.null(MARS$penalty)) { opt@MARS$penalty <- MARS$penalty }
    if (!is.null(MARS$thresh)) { opt@MARS$thresh <- MARS$thresh }
    if (!is.null(MARS$nprune)) { opt@MARS$nprune <- MARS$nprune }
    if (!is.null(MARS$pmethod)) { opt@MARS$pmethod <- MARS$pmethod }
  }
  
  ## 2.9 RF -------------------------------------------------------------------
  if (!is.null(RF)) {
    if (!is.null(RF$type)) { opt@RF$type <- RF$type }
    if (!is.null(RF$ntree)) { opt@RF$ntree <- RF$ntree }
    if (!is.null(RF$mtry)) { opt@RF$mtry <- RF$mtry }
    if (!is.null(RF$sampsize)) { opt@RF$sampsize <- RF$sampsize }
    if (!is.null(RF$nodesize)) { opt@RF$nodesize <- RF$nodesize }
    if (!is.null(RF$maxnodes)) { opt@RF$maxnodes <- RF$maxnodes }
  }
  
  ## 2.10 MAXENT -----------------------------------------------------
  if (!is.null(MAXENT)) {
    if (!is.null(MAXENT$path_to_maxent.jar)) {
      opt@MAXENT$path_to_maxent.jar <- normalizePath(sub("maxent.jar", "", MAXENT$path_to_maxent.jar)) # ensure path format validity
    } else {
      opt@MAXENT$path_to_maxent.jar <- getwd()
    }
    if (!is.null(MAXENT$memory_allocated)) {
      opt@MAXENT$memory_allocated <- MAXENT$memory_allocated
    }
    if (!is.null(MAXENT$initial_heap_size)) {
      opt@MAXENT$initial_heap_size <- MAXENT$initial_heap_size
    }
    if (!is.null(MAXENT$max_heap_size)) {
      opt@MAXENT$max_heap_size <- MAXENT$max_heap_size
    }
    if (!is.null(MAXENT$background_data_dir)) {
      opt@MAXENT$background_data_dir <- MAXENT$background_data_dir
    }
    if (!is.null(MAXENT$maximumbackground)) {
      opt@MAXENT$maximumbackground <- MAXENT$maximumbackground
    }
    if (!is.null(MAXENT$maximumiterations)) {
      opt@MAXENT$maximumiterations <- MAXENT$maximumiterations
    }
    if (!is.null(MAXENT$visible)) {
      opt@MAXENT$visible <- MAXENT$visible
    }
    if (!is.null(MAXENT$linear)) {
      opt@MAXENT$linear <- MAXENT$linear
    }
    if (!is.null(MAXENT$quadratic)) {
      opt@MAXENT$quadratic <- MAXENT$quadratic
    }
    if (!is.null(MAXENT$product)) {
      opt@MAXENT$product <- MAXENT$product
    }
    if (!is.null(MAXENT$threshold)) {
      opt@MAXENT$threshold <- MAXENT$threshold
    }
    if (!is.null(MAXENT$hinge)) {
      opt@MAXENT$hinge <- MAXENT$hinge
    }
    if (!is.null(MAXENT$lq2lqptthreshold)) {
      opt@MAXENT$lq2lqptthreshold <- MAXENT$lq2lqptthreshold
    }
    if (!is.null(MAXENT$l2lqthreshold)) {
      opt@MAXENT$l2lqthreshold <- MAXENT$l2lqthreshold
    }
    if (!is.null(MAXENT$hingethreshold)) {
      opt@MAXENT$hingethreshold <- MAXENT$hingethreshold
    }
    if (!is.null(MAXENT$beta_threshold)) {
      opt@MAXENT$beta_threshold <- MAXENT$beta_threshold
    }
    if (!is.null(MAXENT$beta_categorical)) {
      opt@MAXENT$beta_categorical <- MAXENT$beta_categorical
    }
    if (!is.null(MAXENT$beta_lqp)) {
      opt@MAXENT$beta_lqp <- MAXENT$beta_lqp
    }
    if (!is.null(MAXENT$beta_hinge)) {
      opt@MAXENT$beta_hinge <- MAXENT$beta_hinge
    }
    if (!is.null(MAXENT$betamultiplier)) {
      opt@MAXENT$betamultiplier <- MAXENT$betamultiplier
    }
    if (!is.null(MAXENT$defaultprevalence)) {
      opt@MAXENT$defaultprevalence <- MAXENT$defaultprevalence
    }
  } else {
    opt@MAXENT$path_to_maxent.jar <- getwd()
  }
  
  # if (!is.null(MAXENT.Tsuruoka)) {
  #   if (!is.null(MAXENT.Tsuruoka$l1_regularizer)) { opt@MAXENT.Tsuruoka$l1_regularizer <- MAXENT.Tsuruoka$l1_regularizer }
  #   if (!is.null(MAXENT.Tsuruoka$l2_regularizer)) { opt@MAXENT.Tsuruoka$l2_regularizer <- MAXENT.Tsuruoka$l2_regularizer }
  #   if (!is.null(MAXENT.Tsuruoka$use_sgd)) { opt@MAXENT.Tsuruoka$use_sgd <- MAXENT.Tsuruoka$use_sgd }
  #   if (!is.null(MAXENT.Tsuruoka$set_heldout)) { opt@MAXENT.Tsuruoka$set_heldout <- MAXENT.Tsuruoka$set_heldout }
  #   if (!is.null(MAXENT.Tsuruoka$verbose)) { opt@MAXENT.Tsuruoka$verbose <- MAXENT.Tsuruoka$verbose }
  # }
  
  ## 3. test validity ---------------------------------------------------------
  test <- as.logical(validObject(object = opt, test = TRUE, complete = FALSE))
  
  if (!test) {
    cat("\n\n!!! NULL object returned because of invalid parameters given !!!")
    return(NULL)
  }
  # .bm_cat("Done")
  return(opt)
}


##'
##' @rdname BIOMOD_ModelingOptions
##' @export
##'

bm_DefaultModelingOptions <- function()
{
  cat('\n Defaut modeling options. Copy, change what you want, and paste it as arg to BIOMOD_ModelingOptions().\n\n')
  opt_tmp <- BIOMOD_ModelingOptions()
  print(opt_tmp)
}

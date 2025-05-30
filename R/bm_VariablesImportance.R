###################################################################################################
##' @name bm_VariablesImportance
##' @author Damien Georges
##' 
##' @title Variables' importance calculation
##' 
##' @description This internal \pkg{biomod2} function allows the user to compute a variable 
##' importance value for each variable involved in the given model.
##' 
##' 
##' @param bm.model a \code{biomod2_model} object (or \code{nnet}, \code{rpart}, \code{fda}, 
##' \code{gam}, \code{glm}, \code{lm}, \code{gbm}, \code{mars}, \code{randomForest}, 
##' \code{xgb.Booster}) that can be obtained with the \code{\link{get_formal_model}} function 
##' @param expl.var a \code{data.frame} containing the explanatory variables that will be used to 
##' compute the variables importance
##' @param variables (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the names of the explanatory variables that will be considered
##' @param method a \code{character} corresponding to the randomization method to be used, must be 
##' \code{full_rand} (\emph{only method available so far})
##' @param nb.rep an \code{integer} corresponding to the number of permutations to be done for 
##' each variable
## @param nb.cpu (\emph{optional, default} \code{1}) \cr 
## An \code{integer} value corresponding to the number of computing resources to be used to 
## parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' @param do.progress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the progress bar is to be rendered or not
##' @param temp.workdir (\emph{optional, default} \code{NULL}) \cr 
##' A \code{character} value corresponding to the folder name containing temporal prediction files 
##' when using \code{MAXENT}
##' 
##' 
##' @return  
##' 
##' A \code{3} columns \code{data.frame} containing variable's importance scores for each 
##' permutation run :
##' \itemize{
##'   \item \code{expl.var} : the considered explanatory variable (the one permuted)
##'   \item \code{rand} : the ID of the permutation run
##'   \item \code{var.imp} : the variable's importance score
##' }
##' 
##' @details
##' 
##' For each variable to be evaluated :
##' \enumerate{
##'   \item shuffle the original variable
##'   \item compute model prediction with shuffled variable
##'   \item calculate Pearson's correlation between reference and shuffled predictions
##'   \item return score as \code{1 - cor}
##' }
##' The highest the value, the less reference and shuffled predictions are correlated, and the 
##' more influence the variable has on the model. A value of \code{0} assumes no influence of 
##' the variable on the model.
##' 
##' \emph{Note that this calculation does not account for variables' interactions.}
##' 
##' The same principle is used in \code{\link[randomForest]{randomForest}}.
##' 
##' 
##' @keywords shuffle random importance "Pearson correlation"
##' 
##' 
##' @seealso \code{\link[randomForest]{randomForest}}, 
##' \code{\link{bm_RunModelsLoop}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{BIOMOD_EnsembleModeling}}, \code{\link{bm_PlotVarImpBoxplot}}, 
##' \code{\link{get_variables_importance}}
##' @family Secondary functions
##' 
##' 
##' @examples
##' ## Create simple simulated data
##' myResp.s <- sample(c(0, 1), 20, replace = TRUE)
##' myExpl.s <- data.frame(var1 = sample(c(0, 1), 100, replace = TRUE),
##'                        var2 = rnorm(100),
##'                        var3 = 1:100)
##' 
##' ## Compute variables importance
##' mod <- glm(var1 ~ var2 + var3, data = myExpl.s)
##' bm_VariablesImportance(bm.model = mod, 
##'                        expl.var = myExpl.s[, c('var2', 'var3')],
##'                        method = "full_rand",
##'                        nb.rep = 3)
##' 
##' 
##' @importFrom foreach foreach %do% %:% 
##' @importFrom stats cor
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


bm_VariablesImportance <- function(bm.model, 
                                   expl.var,
                                   variables = NULL,
                                   method = "full_rand", 
                                   nb.rep = 1,
                                   seed.val = NULL,
                                   do.progress = TRUE,
                                   temp.workdir = NULL)
{
  args <- .bm_VariablesImportance.check.args(bm.model, expl.var, variables, method, temp.workdir)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## Test if prediction is computable
  ref <- try(
    predict(bm.model,
            newdata = expl.var,
            temp_workdir = temp.workdir,
            seedval = seed.val)
  )
  
  if (inherits(ref, "try-error")) { stop("Unable to make model prediction") }
  
  if (model_type == "ordinal") {
    if (!(bm.model@model_name %in% c("GLM", "GAM", "XGBOOST"))) { ## keep numeric values for these 3 models
      ref <- as.numeric(factor(ref, ordered = TRUE))
    }
  }
  
  if (model_type == "multiclass") {
    if (!(bm.model@model_name == "XGBOOST")) { 
      ref <- as.numeric(factor(ref))
    }
  }
  
  ## Make randomisation
  cat('\n')
  if (do.progress) {
    PROGRESS = txtProgressBar(min = 0, max = nb.rep * length(variables), style = 3)
    i.iter = 0
  }
  out = foreach (r = 1:nb.rep, .combine = "rbind") %:%
    foreach (v = variables, .combine = "rbind") %do%
    {
      data_rand <- .randomise_data(expl.var, v, method) #, seedval = seed.val)
      shuffled.pred <- predict(bm.model, data_rand, temp_workdir = temp.workdir, seedval = seed.val)
      
      if (model_type == "ordinal") {
        if (!(bm.model@model_name %in% c("GLM", "GAM", "XGBOOST"))) { ## keep numeric values for these 3 models
          shuffled.pred <- as.numeric(factor(shuffled.pred, ordered = TRUE))
        }
      }
      if (model_type == "multiclass") {
        if (!(bm.model@model_name == "XGBOOST")) { 
          shuffled.pred <- as.numeric(factor(shuffled.pred))
        }
      }
      method_cor <- ifelse(model_type %in% c("ordinal", "multiclass"), "pearson", "spearman")
      out_vr <- 1 - max(round(
        cor(x = ref, y = shuffled.pred, use = "pairwise.complete.obs", method = method_cor)
        , digits = 6), 0, na.rm = TRUE)
      if (do.progress) {
        i.iter = i.iter + 1
        setTxtProgressBar(pb = PROGRESS, value = i.iter)
      }
      return(data.frame(expl.var = v, rand = r, var.imp = out_vr))
    }
  if (do.progress) { close(PROGRESS) }
  
  return(out)
}


###################################################################################################

.bm_VariablesImportance.check.args <- function(bm.model, expl.var, variables, method, temp.workdir)
{
  # test that input data is supported
  .fun_testIfInherits(TRUE, "bm.model", bm.model, c("biomod2_model", "nnet", "rpart", "fda", "gam"
                                                    , "glm", "lm", "gbm", "mars", "randomForest"))
  
  model_type <- ifelse(inherits(bm.model, "biomod2_model"), bm.model@model_type, "binary")
  
  # get variables names
  if (is.null(variables)) { variables <- colnames(expl.var) }
  
  # test method is supported
  .fun_testIfIn(TRUE, "method", method, c('full_rand'))
  
  if(missing(temp.workdir)){ temp.workdir <- NULL }
  
  return(list(bm.model = bm.model
              , expl.var = expl.var
              , method = method
              , variables = variables
              , temp.workdir = temp.workdir
              , model_type = model_type))
}


###################################################################################################

.randomise_data <- function(expl.var, variable, method, seedval = NULL)
{
  if (method == 'full_rand') {
    return(.full_shuffling(expl.var, variable, seedval))
  }
}

.full_shuffling <- function(x, id = NULL, seedval = NULL)
{
  if (!(is.vector(x) | is.matrix(x) | is.data.frame(x))) {
    stop("x must be a 1 or 2 dimension odject")
  }
  
  ## Set a new random seed to ensure that sampling is random
  ## (issue when CTA is involved and seed needs to be set to a fix number)
  if (is.null(seedval)) {
    set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6")) * 1000000)
  }
  
  out <- NULL
  if (is.null(id)) {
    out <- x[sample.int(length(x))]
  } else {
    out <- x
    for (idd in id) { out[, idd] <- out[sample.int(nrow(x)), idd]  }
  }
  
  return(out)
}


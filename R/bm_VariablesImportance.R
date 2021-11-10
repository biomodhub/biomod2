##' ###############################################################################################
##' @name bm_VariablesImportance
##' @aliases bm_VariablesImportance
##' @author Damien Georges
##' 
##' @title Variables' importance calculation
##' 
##' @description
##' 
##' This internal \pkg{biomod2} function allows the user to compute a variable importance value 
##' for each variable involved in the given model.
##' 
##' @param model a \code{BIOMOD.models.out} object (coming either from 
##' \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} function) and for 
##' which variables importance is to be computed
##' @param data a \code{data.frame} containing the explanatory variables to use to perform 
##' variables importance calculation
##' @param method a \code{character} corresponding to the randomisation method, must be 
##' \code{full_rand} (\emph{only method available so far})
##' @param nb_rand an \code{integer} corresponding to the number of permutations to be done for 
##' each variable
##' 
##' 
##' @return  
##' 
##' A \code{BIOMOD.variables.importance} object containing :
##' \describe{
##'   \item{mat}{a \code{data.frame} containing variables importance scores for each 
##'   permutation run}
##' }
##' 
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
##' @keywords shuffle, random, importance, Pearson correlation
##' 
##' 
##' @examples
##' 
##' xx <- data.frame(a = sample(c(0, 1), 100, replace = TRUE),
##'                  b = rnorm(100),
##'                  c = 1:100)
##'   
##' mod <- glm(a ~ b + c, data = xx)
##' 
##' bm_VariablesImportance(model = mod, 
##'                        data = xx[, c('b', 'c')],
##'                        method = "full_rand",
##'                        nb_rand = 3)
##' 
##' 
##' ###############################################################################################

bm_VariablesImportance <- function(model, 
                                   data, 
                                   method = "full_rand", 
                                   nb_rand = 1)
{
  args <- .bm_VariablesImportance.check.args(model = model, data = data, method = method)
  
  ## Test if prediction is computable
  ref <- try(predict(args$model, args$data))
  if (inherits(ref, "try-error")) { stop("Unable to make model prediction") }
  
  ## Prepare output matrix
  out <- list()
  out$mat <- matrix(0, nrow = length(args$variables), ncol = nb_rand
                    , dimnames = list(args$variables, paste0('rand', 1:nb_rand)))
  
  ## Make randomisation
  for (r in 1:nb_rand) {
    for (v in args$variables) {
      data_rand <- .randomise_data(args$data, v, method)
      shuffled.pred <- predict(args$model, data_rand)
      out$mat[v, r] <- 1 - max(round(
        cor(x = ref, y = shuffled.pred, use = "pairwise.complete.obs", method = "pearson")
        , digits = 6), 0, na.rm = TRUE)
    }
  }
  
  class(out) <- "BIOMOD.variables.importance"
  return(out)
}


###################################################################################################


.bm_VariablesImportance.check.args <- function(...)
{
  args <- list(...)
  
  # test that input data is supported
  supported_models <- c("biomod2_model", "nnet", "rpart", "fda", "gam", "glm", "lm", "gbm", "mars", "randomForest")
  if (!inherits(args$model, supported_models)) { stop("Model class unsupported") }
  
  # test method is supported
  supported_methods <- c("full_rand")
  if (!(args$method %in% supported_methods)) { stop("Unknown method") }
  
  # get variables names
  if (is.null(args$variables)) { args$variables <- colnames(args$data) }
  
  return(args)
}


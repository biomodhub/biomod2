###################################################################################################
##' @name bm_MakeFormula
##' @author Damien Georges
##' 
##' @title Standardized formula maker
##' 
##' @description
##' 
##' This internal \pkg{biomod2} function allows the user to create easily a standardized formula 
##' that can be used later by statistical models.
##' 
##' @param respName a \code{character} defining the response variable name
##' @param explVar a \code{matrix} or \code{data.frame} corresponding to the explanatory 
##' variables table that will be considered at modelling step
##' @param type a \code{character} defining the wanted type of formula, must be \code{simple}, 
##' \code{quadratic}, \code{polynomial} or \code{s_smoother}
##' @param interaction.level an \code{integer} corresponding to the interaction level depth 
##' between explanatory variables
##' @param \ldots some additional arguments (see \href{bm_MakeFormula.html#details}{Details})
##' 
##' 
##' @return  
##' 
##' A \code{\link[stats]{formula}} class object that can be directly given to most of \R 
##' statistical models.
##' 
##' 
##' @details
##' 
##' It is advised to give only a subset of \code{explVar} table to avoid useless memory consuming. 
##' \cr If some explanatory variables are factorial, \code{explVar} must be a \code{data.frame} 
##' whose corresponding columns are defined as \code{factor}. \cr \cr
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item{\code{k}}{ : an \code{integer} corresponding to the smoothing parameter value of 
##'   \code{\link[mgcv]{s}} or \code{\link[gam]{s}} arguments (\emph{used only if 
##'   \code{type = 's_smoother'}})}
##' }
##' 
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}, 
##' \code{\link[stats]{formula}}
##' 
##' 
##' @keywords models, formula, options
##' 
##' 
##' @examples
##' 
##' ## Create simulated data
##' myResp <- sample(c(0, 1), 20, replace = TRUE)
##' myExpl <- matrix(runif(60), 
##'                  ncol = 3,
##'                  dimnames = list(NULL, c('var1', 'var2', 'var3')))
##' 
##' ## Create a formula
##' myFormula <- bm_MakeFormula(respName = 'myResp',
##'                             explVar = head(myExpl),
##'                             type = 'quadratic',
##'                             interaction.level = 0)
##'   
##' ## Show formula created
##' myFormula
##' 
##' 
##' @export
##' 
##'
###################################################################################################


bm_MakeFormula <- function(respName, 
                           explVar, 
                           type = 'simple', 
                           interaction.level = 0, 
                           ...)
{
  ## 1. Check parameters --------------------------------------------------------------------------
  sup_args <- list(...)
  
  if (!is.character(respName) || length(respName)!=1) {
    stop("respName must be a unique response variable name")
  }
  if (!is.data.frame(explVar) &&  !is.matrix(explVar)) {
    stop("explVar must be a data.frame or matrix")
  }
  .fun_testIfIn(TRUE, "type", type, c("simple", "quadratic", "polynomial", "s_smoother"))
  
  explVarNames <- colnames(explVar)
  if (respName %in% explVarNames) { # remove the response variable if given in explVar
    explVar <- explVar[, -which(explVarNames == respName), drop = FALSE]
    explVarNames <- colnames(explVar)
  }
  
  interaction.level <- min(interaction.level, ncol(explVar))
  
  ## 2. Create the formula ------------------------------------------------------------------------
  junk <- c(1)
  switch(EXPR = type
         , "simple" = {
           junk <- paste(junk, paste(explVarNames, collapse = " + "), sep = " + ")
         }
         , "quadratic" = {
           for (v in 1:ncol(explVar)) {
             if (is.numeric(explVar[, v])) {
               junk <- paste(junk, paste0(explVarNames[v], "+I(", explVarNames[v], "^2)"), sep = " + ")
             } else { junk <- paste(junk, explVarNames[v], sep = " + ") }
           }
         }
         , "polynomial" = {
           for (v in 1:ncol(explVar)) {
             if (is.numeric(explVar[, v])) {
               junk <- paste(junk, paste0(explVarNames[v], "+I(", explVarNames[v], "^2)+I(", explVarNames[v], "^3)"), sep = " + ")
             } else { junk <- paste(junk, explVarNames[v], sep = " + ") }
           }
         }
         , "s_smoother" = {
           for (v in 1:ncol(explVar)) {
             if (is.numeric(explVar[, v])) {
               if (is.null(sup_args$k)) {
                 junk <- paste(junk, paste0("s(", explVarNames[v], ")"), sep = " + ")
               } else {
                 junk <- paste(junk, paste0("s(", explVarNames[v], ",k=", sup_args$k, ")"), sep = " + ")
               }
             } else { junk <- paste(junk, explVarNames[v], sep = " + ") }
           }
         })
                
  
  ## 3. Add interactions if requested -------------------------------------------------------------
  junk.inter <- NULL
  if (interaction.level > 0) {
    for (i.l in 1:interaction.level) {
      inter.tab <- combn(explVarNames, i.l + 1)
      junk.inter <- paste(junk.inter, paste(apply(inter.tab, 2, paste, collapse = ":"), collapse = " + "), sep = " + ")
    }
  }
  if (length(junk.inter)) { junk <- paste0(junk, junk.inter) }
  
  ## 4. Return the formula ------------------------------------------------------------------------
  return(as.formula(paste0(respName," ~ ", junk)))
}

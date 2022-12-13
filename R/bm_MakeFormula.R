###################################################################################################
##' @name bm_MakeFormula
##' @author Damien Georges
##' 
##' @title Standardized formula maker
##' 
##' @description This internal \pkg{biomod2} function allows the user to create easily a 
##' standardized formula that can be used later by statistical models.
##' 
##' 
##' @param resp.name a \code{character} corresponding to the response variable name
##' @param expl.var a \code{matrix} or \code{data.frame} containing the explanatory variables that 
##' will be used at the modeling step
##' @param type a \code{character} corresponding to the wanted type of formula, must be 
##' \code{simple}, \code{quadratic}, \code{polynomial} or \code{s_smoother}
##' @param interaction.level an \code{integer} corresponding to the interaction level depth 
##' between explanatory variables
##' @param k (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the smoothing parameter value of \code{\link[mgcv]{s}} 
##' or \code{\link[gam]{s}} arguments (\emph{used only if \code{type = 's_smoother'}})
##' 
##' @return  
##' 
##' A \code{\link[stats]{formula}} class object that can be directly given to most of 
##' \R statistical models.
##' 
##' 
##' @details
##' 
##' It is advised to give only a subset of \code{expl.var} table to avoid useless memory consuming. 
##' \cr If some explanatory variables are factorial, \code{expl.var} must be a \code{data.frame} 
##' whose corresponding columns are defined as \code{factor}.
##' 
##' 
##' @keywords models formula options
##' 
##' 
##' @seealso \code{\link[stats]{formula}}, \code{\link[mgcv]{s}}, \code{\link[gam]{s}}, 
##' \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Tuning}}, 
##' \code{\link{bm_RunModelsLoop}}
##' @family Secundary functions
##' 
##' 
##' @examples
##' 
##' ## Create simple simulated data
##' myResp.s <- sample(c(0, 1), 20, replace = TRUE)
##' myExpl.s <- data.frame(var1 = sample(c(0, 1), 100, replace = TRUE),
##'                        var2 = rnorm(100),
##'                        var3 = 1:100)
##' 
##' ## Generate automatic formula
##' bm_MakeFormula(resp.name = 'myResp.s',
##'                expl.var = head(myExpl.s),
##'                type = 'quadratic',
##'                interaction.level = 0)
##' 
##' 
##' @importFrom utils combn head read.csv setTxtProgressBar tail txtProgressBar write.table
##' 
##' 
##' @export
##' 
##'
###################################################################################################

bm_MakeFormula <- function(resp.name, 
                           expl.var, 
                           type = 'simple', 
                           interaction.level = 0, 
                           k = NULL)
{
  ## 1. Check parameters --------------------------------------------------------------------------
  if (!is.character(resp.name) || length(resp.name)!=1) {
    stop("resp.name must be a unique response variable name")
  }
  if (!is.data.frame(expl.var) &&  !is.matrix(expl.var)) {
    stop("expl.var must be a data.frame or matrix")
  }
  .fun_testIfIn(TRUE, "type", type, c("simple", "quadratic", "polynomial", "s_smoother"))
  
  explVarNames <- colnames(expl.var)
  if (resp.name %in% explVarNames) { # remove the response variable if given in expl.var
    expl.var <- expl.var[, -which(explVarNames == resp.name), drop = FALSE]
    explVarNames <- colnames(expl.var)
  }
  
  interaction.level <- min(interaction.level, ncol(expl.var))
  
  ## 2. Create the formula ------------------------------------------------------------------------
  junk <- c(1)
  switch(EXPR = type
         , "simple" = {
           junk <- paste(junk, paste(explVarNames, collapse = " + "), sep = " + ")
         }
         , "quadratic" = {
           for (v in 1:ncol(expl.var)) {
             if (is.numeric(expl.var[, v])) {
               junk <- paste(junk, paste0(explVarNames[v], "+I(", explVarNames[v], "^2)"), sep = " + ")
             } else { junk <- paste(junk, explVarNames[v], sep = " + ") }
           }
         }
         , "polynomial" = {
           for (v in 1:ncol(expl.var)) {
             if (is.numeric(expl.var[, v])) {
               junk <- paste(junk, paste0(explVarNames[v], "+I(", explVarNames[v], "^2)+I(", explVarNames[v], "^3)"), sep = " + ")
             } else { junk <- paste(junk, explVarNames[v], sep = " + ") }
           }
         }
         , "s_smoother" = {
           for (v in 1:ncol(expl.var)) {
             if (is.numeric(expl.var[, v])) {
               if (is.null(k)) {
                 junk <- paste(junk, paste0("gam::s(", explVarNames[v], ")"), sep = " + ")
               } else {
                 junk <- paste(junk, paste0("gam::s(", explVarNames[v], ",k=", k, ")"), sep = " + ")
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
  if (length(junk.inter) > 0) { junk <- paste0(junk, junk.inter) }
  
  ## 4. Return the formula ------------------------------------------------------------------------
  return(as.formula(paste0(resp.name," ~ ", junk)))
}

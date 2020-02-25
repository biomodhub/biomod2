##' @name makeFormula
##' @title Standardized formula maker
##' 
##' @description
##' makeFormula is an internal \pkg{biomod2} function that can be useful
##' to help users to build easily some standardized formula used later by
##' statistical models.
##' 
##' @param respName a \code{character} indicating the response variable
##' name
##' @param explVar a \code{matrix} or a \code{data.frame}, the 
##' explanatory variables table that will be considered at modelling step
##' @param type either 'simple', 'quadratic', 'polynomial' or 's_smoother'
##' defining the type of formula you want to build
##' @param interaction.level an \code{integer}, the interaction level
##' depth between explanatory variables
##' @param \ldots some additional arguments (see details)
##' 
##' @details
##' It is advised to give only a subset of \code{explVar} table to avoid
##' useless memory consuming. If some explanatory variables are factorial
##' ones, you have to give a \code{data.frame} for \code{explVar} where
##' associated columns are define as \code{factor}.
##' 
##' \code{...} argument available values are :
##' 
##' - `k` the smoothing parameter value (used only if 
##' \code{type = 's_smoother'}) corresponding to \code{k} parameter 
##' of \pkg{mgcv} \code{\link[mgcv]{s}}  or \code{df} \pkg{gam} 
##' \code{\link[gam]{s}} arguments.
##'
##' @return a \code{link[stats]{formula}} class object that can be
##' directly given to most of \R statistical models.
##' 
##' @author Damien Georges
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}, 
##' \code{link[stats]{formula}}
##' 
##' @keywords models
##' @keywords formula
##' @keywords options
##' 
##' @examples
##' ##' create simulated data
##' myResp <- sample(c(0, 1), 20, replace = TRUE)
##' myExpl <- 
##'   matrix(
##'     runif(60), 
##'     ncol = 3, 
##'     dimnames=list(NULL, c('var1', 'var2', 'var3'))
##'   )
##' 
##' ##' create a formula
##' myFormula <- 
##'   makeFormula( 
##'     respName = 'myResp',
##'     explVar = head(myExpl),
##'     type = 'quadratic',
##'     interaction.level = 0
##'   )
##'   
##' ##' show formula created
##' myFormula
##'
makeFormula <- function(
  respName, 
  explVar, 
  type = 'simple', 
  interaction.level = 0, 
  ...
){
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  # This function return a string in a well formated way. May be give as formula argument to a "basic"
  # statistical model.
  # Several types of models are available
  #
  # D.GEORGES 12/2011
  # =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  sup_args <- list(...)

  # 0. Supported Types =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  availableTypes = c("simple", "quadratic", "polynomial", "s_smoother")

  # 1. Check Given Args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(!is.character(respName) || length(respName)!=1){
    stop("You must give a unique response variable name")
  }

  if(!is.data.frame(explVar) &&  !is.matrix(explVar)){
    stop("You must give explanatory variable table")
  }

  if(!(type %in% availableTypes)){
    stop(paste("Formula type must be one of : ", toString(availableTypes), sep=""))
  }

  explVarNames <- colnames(explVar)
  if(respName %in% explVarNames){ # remove the response variable data if it's given
    explVar <- explVar[, - which(explVarNames == respName), drop=FALSE]
    explVarNames <- colnames(explVar)
  }

  interaction.level <- min(interaction.level, ncol(explVar))

  # 2. Create the formula =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  junk <- c(1)

  switch(EXPR=type,
         "simple" = {junk <- paste(junk,paste(explVarNames, collapse=" + "), sep=" + ") },

         "quadratic" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
#                junk <- c(junk, paste(explVarNames[v], "+I(", explVarNames[v],
#                                       "^2)+I(",explVarNames[v],"^3)", sep="") )
                junk <- paste(junk, paste(explVarNames[v], "+I(", explVarNames[v],
                                      "^2)", sep=""), sep=" + " )
             } else { junk <- paste(junk, explVarNames[v], sep=" + ") }
           } },

         "polynomial" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
               junk <- paste(junk, paste(explVarNames[v],
                                      "+I(", explVarNames[v],
                                      "^2)+I(",explVarNames[v],
                                      "^3)",sep=""), sep=" + " )
#                   junk <- c(junk, paste(explVarNames[v],
#                                       "+poly(",explVarNames[v],
#                                       ",3)",sep="") )
             } else { junk <- paste(junk, explVarNames[v], sep=" + ") }
           } },

         "s_smoother" = {
           for (v in 1:ncol(explVar) ){
             if(is.numeric(explVar[,v])){
               if(is.null(sup_args$k)){
                 junk <- paste(junk, paste("s(",explVarNames[v],
                                       ")",sep=""), sep=" + " )
               } else{
                 junk <- paste(junk, paste("s(",explVarNames[v],
                                       ",k=",sup_args$k,")",sep=""), sep=" + " )
               }

             } else { junk <- paste(junk, explVarNames[v], sep=" + ") }
           } })


  # interactions
  junk.inter <- NULL
  if(interaction.level > 0){
    for(i.l in 1:interaction.level){
      inter.tab <- combn(explVarNames,i.l+1)
      junk.inter <- paste(junk.inter, paste(apply(inter.tab,2,paste,collapse=":"),collapse=" + "), sep= " + ")
    }
  }

  if(length(junk.inter)){
    junk <- paste(junk,junk.inter,sep = "")
  }

  # 2. Return the formula =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  return(as.formula(paste(respName," ~ ", junk, sep="")))

}

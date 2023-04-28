
## --------------------------------------------------------------------------- #
## 1. BIOMOD.options.default -------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.options.default
##' @aliases BIOMOD.options.default-class
##' @author Maya Gueguen
##' 
## @title \code{BIOMOD_FormatingData()} output object class
##' 
## @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
## \code{\link{BIOMOD_Tuning}}, \code{\link{bm_CrossValidation}} and 
## \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @slot model a \code{character} corresponding to the model
##' @slot type a \code{character} corresponding to the data type 
##' (\code{binary}, \code{binary.PA}, \code{abundance}, \code{compositional})
##' @slot package a \code{character} corresponding to the package containing 
##' the model function to be called
##' @slot func a \code{character} corresponding to the model function name 
##' to be called
##' @slot args.names a \code{vector} containing \code{character} corresponding 
##' to the model function arguments
##' @slot args.default a \code{list} containing for each dataset the default 
##' values for all arguments listed in \code{args.names}
##' 
##' 
## @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_Tuning}}, 
## \code{\link{bm_CrossValidation}}, \code{\link{BIOMOD_Modeling}}, 
## \code{\link{bm_RunModelsLoop}}
## @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.options.default")
##' 
##' ## ----------------------------------------------------------------------- #
##' 
##' 
NULL

##' @name BIOMOD.options.default-class
##' @rdname BIOMOD.options.default
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.options.default",
         representation(model = 'character',
                        type = 'character',
                        package = "character",
                        func = "character",
                        args.names = "character",
                        args.default = "list"),
         validity = function(object){ return(TRUE) })


# 1.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.options.default", def = function(mod, typ, pkg, fun) { standardGeneric("BIOMOD.options.default") })

.BIOMOD.options.default.check.args <- function(mod, typ, pkg, fun)
{
  ## check if model is supported
  avail.models.list <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS'
                         , 'RF', 'MAXENT', 'MAXNET', 'XGBOOST')
  .fun_testIfIn(TRUE, "mod", mod, avail.models.list)
  
  ## check if type is supported
  avail.types.list <- c('binary', 'binary.PA', 'abundance', 'compositional')
  .fun_testIfIn(TRUE, "typ", typ, avail.types.list)
  
  ## check package exists
  if (!isNamespaceLoaded(pkg)) { requireNamespace(pkg, quietly = TRUE) } ## ATTENTION n'a pas l'air de suffire
  
  ## check function exists
  avail.functions.list <- lsf.str(paste0("package:", pkg)) ## Du coup ne marche pas si le package n'est pas explicitement loadé
  .fun_testIfIn(TRUE, "fun", fun, avail.functions.list)
}


## BIOMOD.options.default -------------------------------------------------------
##' 
##' @rdname BIOMOD.options.default
##' @export
##' 

setMethod('BIOMOD.options.default', signature(mod = 'character', typ = 'character'),
          function(mod, typ, pkg, fun) 
          {
            .BIOMOD.options.default.check.args(mod, typ, pkg, fun)
            BOM <- new(
              'BIOMOD.options.default',
              model = mod,
              type = typ,
              package = pkg,
              func = fun,
              args.names = formalArgs(fun),
              args.default = as.list(formals(fun))
            )
            return(BOM)
          }
)


# 1.3 Other Functions -----------------------------------------------------------------------------



## --------------------------------------------------------------------------- #
## 2. BIOMOD.options.dataset -------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.options.dataset
##' @aliases BIOMOD.options.dataset-class
##' @author Maya Gueguen
##' 
## @title \code{BIOMOD_FormatingData()} output object class
##' 
## @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
## \code{\link{BIOMOD_Tuning}}, \code{\link{bm_CrossValidation}} and 
## \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @slot ... slots from \code{\link{BIOMOD.options.default}} object
##' @slot args.values a \code{list} containing for each dataset the to-be-used  
##' values for all arguments listed in \code{args.names}
##' 
##' 
## @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_Tuning}}, 
## \code{\link{bm_CrossValidation}}, \code{\link{BIOMOD_Modeling}}, 
## \code{\link{bm_RunModelsLoop}}
## @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.options.dataset")
##' 
##' ## ----------------------------------------------------------------------- #
##' 
##' 
NULL

##' @name BIOMOD.options.dataset-class
##' @rdname BIOMOD.options.dataset
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.options.dataset",
         contains = "BIOMOD.options.default",
         representation(args.values = "list"),
         validity = function(object){ return(TRUE) })


# 1.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.options.dataset", def = function(strategy, val = NULL, bm.format = NULL, calib.lines = NULL, ...) {
  standardGeneric("BIOMOD.options.dataset") })

.BIOMOD.options.dataset.check.args <- function(strategy, val = NULL, bm.format = NULL, calib.lines = NULL)
{
  ## check if strategy is supported
  avail.strategy.list <- c('default', 'bigboss', 'user.defined', 'tuned')
  .fun_testIfIn(TRUE, "strategy", strategy, avail.strategy.list)
  
  ## USER DEFINED parameterisation --------------
  if (strategy == "user.defined") {
    .fun_testIfInherits(TRUE, "val", val, c("list"))
  }
  
  ## TUNING with bm_Tuning parameterisation -----
  if (strategy == "tuned") {
    .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
    
    if (!is.null(calib.lines)) {
      .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
      
      expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun")
      if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
        expected_CVnames <- c(expected_CVnames
                              , sapply(1:ncol(bm.format@PA.table)
                                       , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
                                                             , paste0("_PA", this_PA, "_allRun"))))
      } 
      .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
    }
  }
}


## BIOMOD.options.dataset -------------------------------------------------------
##' 
##' @rdname BIOMOD.options.dataset
##' @export
##' 

setMethod('BIOMOD.options.dataset', signature(strategy = 'character'),
          function(mod, typ, pkg, fun, strategy, val = NULL, bm.format = NULL, calib.lines = NULL)
          {
            .BIOMOD.options.dataset.check.args(strategy, val, bm.format, calib.lines)
            
            BOM <- BIOMOD.options.default(mod, typ, pkg, fun)
            
            ## GET parameter values according to strategy -------------------------------
            if (is.null(calib.lines)) {
              argsval <- list("AllData_AllRun" = BOM@args.default)
            } else {
              argsval <- lapply(1:ncol(calib.lines), function(xx) { BOM@args.default })
              names(argsval) <- colnames(calib.lines)
            }
            if (strategy == "bigboss") {
              ## create and load specific dataset
              ## TODO
            } else if (strategy == "user.defined") {
              if (!("..." %in% BOM@args.names)) {
                .fun_testIfIn(TRUE, "names(val)", names(val), BOM@args.names)
              } else {
                ## ???
              }
              argsval <- val
            } else if (strategy == "tuned") {
              ## Call to bm_Tuning for one model at a time
              argsval <- bm_Tuning(obj = BOM, bm.format = bm.format, calib.lines = calib.lines)
            }
            
            ## SPECIFIC case of formula -------------------------------------------------
            if ("formula" %in% BOM@args.names && !is.null(bm.format)) {
              for (ii in 1:length(argsval)) {
                if (is.null(argsval[[ii]]$formula) || nchar(argsval[[ii]]$formula) == 0) {
                  argsval[[ii]]$formula <- bm_MakeFormula(resp.name = bm.format@sp.name
                                                          , expl.var = head(bm.format@data.env.var)
                                                          , type = 'simple'
                                                          , interaction.level = 0)
                }
              }
            }
            ## ATTENTION : si on ne donne pas bm.format, on n'a pas de formula du coup
            
            BOD <- new('BIOMOD.options.dataset', BOM, args.values = argsval)
            return(BOD)
          }
)


# 1.3 Other Functions -----------------------------------------------------------------------------


###################################################################################################

### BIOMOD.stored.options ------------------------------------------------
##' @name BIOMOD.stored.options-class
##' @rdname BIOMOD.stored.data
##' 
setClass("BIOMOD.stored.options",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.options.dataset'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) })


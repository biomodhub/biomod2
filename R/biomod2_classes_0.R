
## --------------------------------------------------------------------------- #
## 1. BIOMOD.options.default -------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.options.default
##' @aliases BIOMOD.options.default-class
##' @author Maya Gueguen
##' 
##' @title \code{\link{bm_ModelingOptions}} output object class
##' 
##' @description Class returned by \code{\link{bm_ModelingOptions}} (a 
##' \code{list} of \code{BIOMOD.options.dataset} more exactly), and used by 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @param mod a \code{character} corresponding to the model name to be computed, must be either 
##' \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, 
##' \code{MARS}, \code{RF}, \code{MAXENT}, \code{MAXNET}, \code{XGBOOST}
##' @param typ a \code{character} corresponding to the data type to be used, must be either 
##' \code{binary}, \code{binary.PA}, \code{abundance}, \code{compositional}
##' @param pkg a \code{character} corresponding to the package containing 
##' the model function to be called
##' @param fun a \code{character} corresponding to the model function name 
##' to be called
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
##' @seealso \code{\link{BIOMOD.options.dataset}}, \code{\link{bm_ModelingOptions}}, 
##' \code{\link{bm_Tuning}}, \code{\link{BIOMOD_Modeling}}, \code{\link{bm_RunModelsLoop}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.options.default")
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


# 1.2 Constructors --------------------------------------------------------------------------------
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
  
  if (mod != 'MAXENT') {
    ## check package exists
    # if (!isNamespaceLoaded(pkg)) { requireNamespace(pkg, quietly = TRUE) } ## ATTENTION n'a pas l'air de suffire
    eval(parse(text = paste0("require(", pkg, ")")))
    
    ## check function exists
    avail.functions.list <- lsf.str(pos = paste0("package:", pkg)) ## Du coup ne marche pas si le package n'est pas explicitement loadé
    .fun_testIfIn(TRUE, "fun", fun, avail.functions.list)
  }
}


## BIOMOD.options.default -----------------------------------------------------
##' 
##' @rdname BIOMOD.options.default
##' @export
##' 

setMethod('BIOMOD.options.default', signature(mod = 'character', typ = 'character'),
          function(mod, typ, pkg, fun) 
          {
            .BIOMOD.options.default.check.args(mod, typ, pkg, fun)
            if (mod != 'MAXENT') {
              BOM <- new(
                'BIOMOD.options.default',
                model = mod,
                type = typ,
                package = pkg,
                func = fun,
                args.names = formalArgs(fun),
                args.default = as.list(formals(fun))
              )
            } else {
              params.MAXENT = list(path_to_maxent.jar = getwd(),
                                   memory_allocated = 512,
                                   initial_heap_size = NULL,
                                   max_heap_size = NULL,
                                   background_data_dir = 'default',
                                   maximumbackground = 'default',
                                   maximumiterations = 200,
                                   visible = FALSE,
                                   linear = TRUE,
                                   quadratic = TRUE,
                                   product = TRUE,
                                   threshold = TRUE,
                                   hinge = TRUE,
                                   lq2lqptthreshold = 80,
                                   l2lqthreshold = 10,
                                   hingethreshold = 15,
                                   beta_threshold = -1.0,
                                   beta_categorical = -1.0,
                                   beta_lqp = -1.0,
                                   beta_hinge = -1.0,
                                   betamultiplier = 1,
                                   defaultprevalence = 0.5)
              
              BOM <- new(
                'BIOMOD.options.default',
                model = mod,
                type = typ,
                package = pkg,
                func = fun,
                args.names = names(params.MAXENT),
                args.default = params.MAXENT
              )
            }
            return(BOM)
          }
)



## --------------------------------------------------------------------------- #
## 2. BIOMOD.options.dataset -------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.options.dataset
##' @aliases BIOMOD.options.dataset-class
##' @author Maya Gueguen
##' 
##' @title \code{\link{bm_ModelingOptions}} output object class
##' 
##' @description Class returned by \code{\link{bm_ModelingOptions}} (a 
##' \code{list} of \code{BIOMOD.options.dataset} more exactly), and used by 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @inheritParams BIOMOD.options.default
##' @param strategy a \code{character} corresponding to the method to 
##' select models' parameters values, must be either \code{default}, 
##' \code{bigboss}, \code{user.defined}, \code{tuned}
##' @param val (\emph{optional, default} \code{NULL}) \cr
##' A \code{list} containing parameters values
##' @param bm.format (\emph{optional, default} \code{NULL}) \cr
##' A \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' A \code{data.frame} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{bm_CrossValidation}} functions, to explore the distribution of calibration 
##' and validation datasets
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
##' @slot args.values a \code{list} containing for each dataset the to-be-used  
##' values for all arguments listed in \code{args.names}
##' 
##' 
##' @seealso \code{\link{BIOMOD.options.default}}, \code{\link{bm_ModelingOptions}}, 
##' \code{\link{bm_Tuning}}, \code{\link{BIOMOD_Modeling}}, \code{\link{bm_RunModelsLoop}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.options.dataset")
##' 
##' 
NULL

##' @name BIOMOD.options.dataset-class
##' @rdname BIOMOD.options.dataset
##' @export
##' 

# 2.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.options.dataset",
         contains = "BIOMOD.options.default",
         representation(args.values = "list"),
         validity = function(object){ return(TRUE) })


# 2.2 Constructors --------------------------------------------------------------------------------
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
  
  # if (!is.null(MAXENT$path_to_maxent.jar)) {
  #   opt@MAXENT$path_to_maxent.jar <- normalizePath(sub("maxent.jar", "", MAXENT$path_to_maxent.jar)) # ensure path format validity
  # } else {
  #   opt@MAXENT$path_to_maxent.jar <- getwd()
  # }
  
  if (!is.null(bm.format)) {
    .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  }
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
  #   allPA <- sapply(colnames(calib.lines), function(xx) strsplit(xx, "_")[[1]][2])
  #   allrun <- sapply(colnames(calib.lines), function(xx) strsplit(xx, "_")[[1]][3])
  # } else if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
  #   allPA <- colnames(bm.format@PA.table)
  #   allrun <- rep(NA, length(allPA))
  # }
  # 
  # return(list(allPA = allPA, allrun = allrun))
}


## BIOMOD.options.dataset -----------------------------------------------------
##' 
##' @rdname BIOMOD.options.dataset
##' @export
##' 

setMethod('BIOMOD.options.dataset', signature(strategy = 'character'),
          function(mod, typ, pkg, fun, strategy, val = NULL, bm.format = NULL, calib.lines = NULL)
          {
            cat('\n\t> ', mod, 'options (datatype:', typ, ', package:', pkg, ', function:', fun, ')...')
            
            .BIOMOD.options.dataset.check.args(strategy, val, bm.format, calib.lines)
            
            BOM <- BIOMOD.options.default(mod, typ, pkg, fun)
            
            ## GET parameter values according to strategy -------------------------------
            if (strategy %in% c("default", "bigboss")) {
              argstmp <- BOM@args.default
              if (mod == "GLM") { ## SHOULD BE MOVED to place when testing values !!
                argstmp$control = list()
              }
              
              if (strategy == "bigboss") {
                if (mod == "ANN") {
                  argstmp$size = NULL
                  argstmp$decay = 5
                  argstmp$trace = FALSE
                  argstmp$rang = 0.1
                  argstmp$maxit = 200
                  argstmp$nbCV = 5
                } else if (mod == "CTA") {
                  argstmp$method = "class"
                  argstmp$control = list(xval = 5, 
                                         minbucket = 5, 
                                         minsplit = 5,
                                         cp = 0.001, 
                                         maxdepth = 25)
                  argstmp$cost = NULL
                } else if (mod == "FDA") {
                  argstmp$method = "mars"
                } else if (mod == "GAM" && pkg == "mgcv") {
                  if (!is.null(bm.format)) {
                    argstmp$formula = bm_MakeFormula(resp.name = bm.format@sp.name
                                                     , expl.var = head(bm.format@data.env.var)
                                                     , type = 's_smoother'
                                                     , interaction.level = 0
                                                     , k = NULL)
                  }
                  argstmp$family = binomial(link = 'logit')
                  argstmp$method = "GCV.Cp"
                  argstmp$control = list(epsilon = 1e-06, trace = FALSE, maxit = 100)
                } else if (mod == "GBM") {
                  argstmp$n.trees = 2500
                  argstmp$interaction.depth = 7
                  argstmp$n.minobsinnode = 5
                  argstmp$shrinkage = 0.001
                  argstmp$cv.folds = 3
                  argstmp$keep.data = FALSE
                  argstmp$n.cores = 1
                } else if (mod == "GLM") {
                  if (!is.null(bm.format)) {
                    argstmp$formula = bm_MakeFormula(resp.name = bm.format@sp.name
                                                     , expl.var = head(bm.format@data.env.var)
                                                     , type = 'quadratic'
                                                     , interaction.level = 0)
                  }
                  argstmp$family = binomial(link = 'logit')
                  argstmp$mustart = 0.5
                  argstmp$control = glm.control(maxit = 50)
                } else if (mod == "MARS") {
                  argstmp$glm = list(family = binomial)
                  argstmp$ncross = 0
                  argstmp$nk = NULL
                  argstmp$penalty = 2
                  argstmp$thresh = 0.001
                  argstmp$nprune = NULL
                  argstmp$pmethod = 'backward'
                } else if (mod == "RF") {
                  argstmp$do.classif = TRUE
                  argstmp$ntree = 500
                  argstmp$mtry = NULL
                  argstmp$strata = factor(c(0, 1))
                  argstmp$sampsize = NULL
                  argstmp$nodesize = 5
                  argstmp$maxnodes = NULL
                } else if (mod == "SRE") {
                  argstmp$do.extrem = TRUE
                } else if (mod == "XGBOOST") {
                  argstmp$max.depth = 2
                  argstmp$eta = 1
                  argstmp$nthread = 2
                  argstmp$nrounds = 4
                  argstmp$objective = "binary:logistic"
                }
              }
                
              if (is.null(calib.lines)) {
                argsval <- list("_allData_allRun" = argstmp)
                if (inherits(bm.format, "BIOMOD.formated.data.PA")) {
                  for (PA in colnames(bm.format@PA.table)) {
                    argsval[[paste0("_", PA, "_allRun")]] <- argsval[["_allData_allRun"]]
                  }
                }
              } else {
                argsval <- lapply(1:ncol(calib.lines), function(xx) { argstmp })
                names(argsval) <- colnames(calib.lines)
              }
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
                if (is.null(argsval[[ii]]$formula) || (length(argsval[[ii]]$formula) == 1 && nchar(argsval[[ii]]$formula) == 0)) {
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


# 2.3 Other Functions -----------------------------------------------------------------------------

### show BIOMOD.options.dataset -----------------------------------------------
##'
##' @rdname BIOMOD.options.dataset
##' @importMethodsFrom methods show
##' @export
##'


## Attention ! print only values for _allData_allRun for now
setMethod('show', signature('BIOMOD.options.dataset'),
          function(object)
          {
            cat('\n\t> ', object@model, 'options (datatype:', object@type, ', package:', object@package, ', function:', object@func, ') :')
            # for (arg in object@args.names) { ## NOT working for bigboss for example, if new parameters
            for (arg in names(object@args.values[["_allData_allRun"]])) {
              val.def = capture.output(object@args.default[[arg]])
              val.used = capture.output(object@args.values[["_allData_allRun"]][[arg]])
              
              cat('\n\t\t- ', arg, "=", sub("\\[1\\] ", "", val.used))
              if (!is.null(val.used) && !is.null(val.def) && 
                  (length(val.used) != length(val.def) || val.used != val.def)) {
                cat('   (default:', sub("\\[1\\] ", "", val.def), ')')
              }
            }
            cat("\n")
          }
)


## --------------------------------------------------------------------------- #
## 3. BIOMOD.models.options --------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.models.options
##' @aliases BIOMOD.models.options-class
##' @author Maya Gueguen
##' 
##' @title \code{\link{bm_ModelingOptions}} output object class
##' 
##' @description Class returned by \code{\link{bm_ModelingOptions}} and used by 
##' \code{\link{BIOMOD_Modeling}
##' 
##' 
##' @slot models a \code{vector} containing model names for which options have 
##' been retrieved and defined, must be \code{algo.datatype.package.function}
##' @slot options a \code{list} containing \code{\link{BIOMOD.options.dataset}} 
##' object for each model
##' 
##' @param object a \code{\link{BIOMOD.models.options}} object
##' 
##' 
##' @seealso \code{\link{BIOMOD.options.default}}, 
##' \code{\link{BIOMOD.options.dataset}}, 
##' \code{\link{bm_ModelingOptions}}, \code{\link{bm_Tuning}}, 
##' \code{\link{BIOMOD_Modeling}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.models.options")
##' 
##' 
NULL

##' @name BIOMOD.models.options-class
##' @rdname BIOMOD.models.options
##' @export
##' 

# 3.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.models.options",
         representation(models = "character",
                        options = "list"),
         validity = function(object){ return(TRUE) })

# 3.3 Other Functions -----------------------------------------------------------------------------

### show BIOMOD.models.options ------------------------------------------------
##'
##' @rdname BIOMOD.models.options
##' @importMethodsFrom methods show
##' @export
##'

setMethod('show', signature('BIOMOD.models.options'),
          function(object)
          {
            .bm_cat("BIOMOD.models.options")
            for (mod in object@models) {
              show(object@options[[mod]])
            }
            .bm_cat()
          }
)

# test <- .fun_testIfIn(test, "GLM$type", object@GLM$type, c("simple", "quadratic", "polynomial", "user.defined")) ## MOVE in formula ?
# test <- .fun_testIfIn(test, "GLM$test", object@GLM$test, c("AIC", "BIC", "none"))
# test <- .fun_testIfIn(test, "GBM$distribution", object@GBM$distribution, c("bernoulli", "huberized", "multinomial", "adaboost"))
# test <- .fun_testIfIn(test, "GBM$perf.method", object@GBM$perf.method, c('OOB', 'test', 'cv'))
# test <- .fun_testIfIn(test, "GAM$type", object@GAM$type, c('s_smoother', 's', 'lo', 'te')) ## MOVE in formula ?
# test <- .fun_testIfIn(test, "GAM$method", object@GAM$method, c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML"))
# if (any(!object@GAM$optimizer %in%  c("perf", "outer", "newton", "bfgs", "optim", "nlm", "nlm.fd"))) {
#   cat("\nGAM$optimizer bad definition (see ?mgcv::gam)")
# }
# test <- .fun_testIfIn(test, "CTA$method", object@CTA$method, c("anova", "poisson", "class", "exp"))
# test <- .fun_testIfIn(test, "FDA$method", object@FDA$method, c('polyreg', 'mars', 'bruto'))
# test <- .fun_testIfIn(test, "MARS$type", object@MARS$type, c("simple", "quadratic", "polynomial", "user.defined")) ## MOVE in formula ?
# supported.pmethod <- c('backward', 'none', 'exhaustive', 'forward', 'seqrep', 'cv')
# if(!is.element(object@MARS$pmethod, supported.pmethod)){
#   cat("\nMARS$pmethod must be a one of", supported.pmethod);
# }
# 
# ## MAXENT ------------------------------------------------------------
# if (!is.character(object@MAXENT$path_to_maxent.jar)) {
#   cat("\nMAXENT$path_to_maxent.jar must be a character")
#   test <- FALSE
# }
# if (!is.null(object@MAXENT$memory_allocated)) {
#   if (!is.numeric(object@MAXENT$memory_allocated)) {
#     cat("\nMAXENT$memory_allocated must be a positive integer or NULL for unlimited memory allocation")
#     test <- FALSE
#   }
# }
# if (!is.character(object@MAXENT$background_data_dir)) {
#   cat("\nMAXENT$background_data_dir must be 'default' (=> use the same pseudo absences than other models as background) or a path to the directory where your environmental layer are stored")
#   test <- FALSE
# }
# tt <- is.character(object@MAXENT$maximumbackground) | is.numeric(object@MAXENT$maximumbackground)
# if (is.character(object@MAXENT$maximumbackground)) if (object@MAXENT$maximumbackground != "default") tt <- FALSE
# if (!tt) {
#   cat("\nMAXENT$maximumbackground must be 'default' or numeric")
# }
# test <- .fun_testIfPosInt(test, "MAXENT$maximumiterations", object@MAXENT$maximumiterations)
# if (!is.logical(object@MAXENT$visible)) {
#   cat("\nMAXENT$visible must be a logical")
# }
# if (!is.logical(object@MAXENT$linear)) {
#   cat("\nMAXENT$linear must be a logical")
# }
# if (!is.logical(object@MAXENT$quadratic)) {
#   cat("\nMAXENT$quadratic must be a logical")
# }
# if (!is.logical(object@MAXENT$product)) {
#   cat("\nMAXENT$product must be a logical")
# }
# if (!is.logical(object@MAXENT$threshold)) {
#   cat("\nMAXENT$threshold must be a logical")
# }
# if (!is.logical(object@MAXENT$hinge)) {
#   cat("\nMAXENT$hinge must be a logical")
# }
# if (!is.numeric(object@MAXENT$lq2lqptthreshold)) {
#   cat("\nMAXENT$lq2lqptthreshold must be a numeric")
# }
# if (!is.numeric(object@MAXENT$l2lqthreshold)) {
#   cat("\nMAXENT$l2lqthreshold must be a numeric")
# }
# if (!is.numeric(object@MAXENT$lq2lqptthreshold)) {
#   cat("\nMAXENT$lq2lqptthreshold must be a numeric")
# }
# if (!is.numeric(object@MAXENT$hingethreshold)) {
#   cat("\nMAXENT$hingethreshold must be a numeric")
# }
# if (!is.numeric(object@MAXENT$beta_threshold)) {
#   cat("\nMAXENT$beta_threshold must be a numeric")
# }
# if (!is.numeric(object@MAXENT$beta_categorical)) {
#   cat("\nMAXENT$beta_categorical must be a numeric")
# }
# if (!is.numeric(object@MAXENT$beta_lqp)) {
#   cat("\nMAXENT$beta_lqp must be a numeric")
# }
# if (!is.numeric(object@MAXENT$beta_hinge)) {
#   cat("\nMAXENT$beta_hinge must be a numeric")
# }
# if (!is.numeric(object@MAXENT$betamultiplier)) {
#   cat("\nMAXENT$betamultiplier must be a numeric")
# }
# if (!is.numeric(object@MAXENT$defaultprevalence)) {
#   cat("\nMAXENT$defaultprevalence must be a numeric")
# }
# 
# if(!is.null(object@MAXENT$initial_heap_size)){
#   test <- .check_bytes_format(test,
#                               object@MAXENT$initial_heap_size,
#                               "initial_heap_size")
# }
# if(!is.null(object@MAXENT$max_heap_size)){
#   test <- .check_bytes_format(test,
#                               object@MAXENT$max_heap_size,
#                               "max_heap_size")
# }


# ### show.BIOMOD.models.options -------------------------------------------------
# 
# ## GLM options : control = glm.control(", .print_control(object@GLM$control), ") ),", sep = "", fill = .Options$width)
# ## FDA options : paste(.print_control(object@FDA$add_args), collapse = "")
# cat("\n MAXNET = list( myFormula = ", .print_formula(object@MAXNET$myFormula), ",", sep = "")




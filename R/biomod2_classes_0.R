
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
##' \code{ANN}, \code{CTA}, \code{FDA}, \code{GAM}, \code{GBM}, \code{GLM}, \code{MARS}, 
##' \code{MAXENT}, \code{MAXNET}, \code{RF}, \code{SRE}, \code{XGBOOST}
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
##' @importFrom utils lsf.str
##' @importFrom methods formalArgs
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
  avail.models.list <- c('ANN', 'CTA', 'FDA', 'GAM', 'GBM', 'GLM', 'MARS', 'MAXENT', 'MAXNET', 'RF', 'SRE', 'XGBOOST')
  .fun_testIfIn(TRUE, "mod", mod, avail.models.list)
  
  ## check if type is supported
  # avail.types.list <- c('binary', 'binary.PA', 'abundance', 'compositional')
  avail.types.list <- c('binary')
  .fun_testIfIn(TRUE, "typ", typ, avail.types.list)
  
  if (mod != 'MAXENT') {
    ## check package exists
    # lsf.str can be used only with attached package so we are
    # forced to load the packages used in the modeling with attachNamespace
    # attachNamespace(pkg) ==> ERROR : l'espace de noms est déjà attaché
    eval(parse(text = paste0("require(", pkg, ")")))
    
    ## check function exists
    avail.functions.list <- lsf.str(pos = paste0("package:", pkg))
    .fun_testIfIn(TRUE, "fun", fun, avail.functions.list)
  } else {
    if (!file.exists(file.path(getwd(), "maxent.jar"))) {
      warning(paste0("'maxent.jar' file is missing in current working directory ("
                     , getwd(), ").\n"
                     , "It must be downloaded (https://biodiversityinformatics.amnh.org/open_source/maxent/) "
                     , "and put in the working directory."))
    }
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
                args.names = eval(parse(text = paste0("formalArgs(", pkg, "::", fun, ")"))),
                args.default = eval(parse(text = paste0("as.list(formals(", pkg, "::", fun, "))")))
              )
            } else {
              params.MAXENT = list(path_to_maxent.jar = getwd(),
                                   memory_allocated = 512,
                                   initial_heap_size = NULL,
                                   max_heap_size = NULL,
                                   background_data_dir = 'default',
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
##' @param user.val (\emph{optional, default} \code{NULL}) \cr
##' A \code{list} containing parameters values
##' @param user.base (\emph{optional, default} \code{NULL}) \cr A character, 
##' \code{default} or \code{bigboss} used when \code{strategy = 'user.defined'}. 
##' It sets the bases of parameters to be modified by user defined values.
##' @param tuning.fun (\emph{optional, default} \code{NULL}) \cr
##' A \code{character} corresponding to the model function name 
##' to be called through \code{\link[caret]{train}} function for tuning parameters
##' @param bm.format (\emph{optional, default} \code{NULL}) \cr
##' A \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' A \code{data.frame} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{bm_CrossValidation}} functions, to explore the distribution of calibration 
##' and validation datasets
##' 
##' @param object a \code{\link{BIOMOD.options.dataset}} object
##' @param x a \code{\link{BIOMOD.options.dataset}} object
##' @param dataset a \code{character} corresponding to the name of a dataset contained 
##' in the \code{arg.values} slot
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
##' @importFrom stats rnorm runif na.omit na.exclude
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
setGeneric("BIOMOD.options.dataset",
           def = function(strategy, user.val = NULL, user.base = NULL, 
                          tuning.fun = NULL, bm.format = NULL, calib.lines = NULL, ...) {
             standardGeneric("BIOMOD.options.dataset")
           })

.BIOMOD.options.dataset.check.args <- function(strategy, user.val = NULL, user.base = NULL, tuning.fun = NULL, bm.format = NULL, calib.lines = NULL)
{
  ## check if strategy is supported
  avail.strategy.list <- c('default', 'bigboss', 'user.defined', 'tuned')
  .fun_testIfIn(TRUE, "strategy", strategy, avail.strategy.list)
  
  ## USER DEFINED parameterisation --------------
  if (strategy == "user.defined") {
    avail.user.base <- c('default', 'bigboss')
    .fun_testIfIn(TRUE, "user.base", user.base, avail.user.base)
    
    if (!is.null(user.val)) {
      .fun_testIfInherits(TRUE, "user.val", user.val, c("list"))
      
      
    } else if (user.base == "bigboss") {
      strategy <- "bigboss" # revert to bigboss options 
    } else {
      strategy <- "default" # revert to default options
    }
  }
  
  ## TUNING parameterisation --------------------
  if (strategy == "tuned") {
    all.fun <- c('avNNet', 'rpart', 'rpart2', 'fda', 'gamSpline', 'bam', 'gam', 'gbm', 'glm', 'earth', 'rf', 'xgbTree')
    .fun_testIfIn(TRUE, "tuning.fun", tuning.fun, c(all.fun, "bm_SRE", "ENMevaluate", "maxnet"))
    .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  }
  
  ## check bm.format, bm.format@PA.table and calib.lines
  if (!is.null(bm.format)) {
    .fun_testIfInherits(TRUE, "bm.format", bm.format, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  }
  expected_CVnames <- "_allData_allRun"
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), expected_CVnames)
    
    if (!is.null(bm.format) && inherits(bm.format, "BIOMOD.formated.data.PA")) {
      expected_CVnames <- c(expected_CVnames
                            , sapply(1:ncol(bm.format@PA.table)
                                     , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
                                                           , paste0("_PA", this_PA, "_allRun"))))
    } 
    .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
    expected_CVnames <- colnames(calib.lines)
  } else {
    if (!is.null(bm.format) && inherits(bm.format, "BIOMOD.formated.data.PA")) {
      expected_CVnames <- c(expected_CVnames
                            , sapply(1:ncol(bm.format@PA.table)
                                     , function(this_PA) paste0("_PA", this_PA, "_allRun")))
    }
  }
  if (strategy == "user.defined" && !is.null(user.val)) {
    .fun_testIfIn(TRUE, "names(user.val)", names(user.val), expected_CVnames)
    if (length(names(user.val)) != length(expected_CVnames)) {
      warning(paste0("Options will be changed only for a subset of datasets ("
                     , paste0(names(user.val), collapse = ", ")
                     , ") and not the others ("
                     , paste0(setdiff(expected_CVnames, names(user.val)), collapse = ", ")
                     , "). \nPlease update 'user.val' argument if this is not wanted."))
    }
  }
  
  return(list(strategy = strategy,
              expected_CVnames = expected_CVnames))
}

# .BIOMOD.options.dataset.test <- function(bm.opt)
# {
#   ## create false dataset
#   mySp = "Hariba"
#   myResp = c(rep(1, 100), rep(0, 100))
#   myRespXY = data.frame(x = sample(1:100, 200, replace = TRUE)
#                         , y = sample(1:100, 200, replace = TRUE))
#   myExpl = data.frame(var1 = rnorm(200, mean = 3, sd = 0.4)
#                       , var2 = runif(200)
#                       , var3 = sample(1:5, 200, replace = TRUE))
#   mySpExpl = cbind(myResp, myExpl)
#   colnames(mySpExpl)[1] = mySp
#   
#   ## get options for specific model
#   if (bm.opt@model == "GAM") {
#     subclass_name <- paste0(bm.opt@model, "_", bm.opt@type, "_", bm.opt@package)
#     .load_gam_namespace(model_subclass = subclass_name)
#   }
#   
#   ## 
#   for (xx in 1:length(bm.opt@args.values)) { ## SHOULD BE MOVED to place when testing values ??
#     if ('...' %in% names(bm.opt@args.values[[xx]])) {
#       bm.opt@args.values[[xx]][['...']] <- NULL
#     }
#   }
#   print(bm.opt)
#   
#   ## run test for each dataset
#   for (dataset_name in names(bm.opt@args.values)) {
#     bm.opt.val <- bm.opt@args.values[[dataset_name]]
#     
#     ## 2. CREATE MODELS -----------------------------------------------------------------------------
#     set.seed(42)
#     
#     if (bm.opt@model != "MAXENT") { ## ANY MODEL BUT MAXENT ------------------------------------------------
#       ## PRELIMINAR ---------------------------------------------------
#       bm.opt.val$formula <- bm_MakeFormula(resp.name = mySp
#                                            , expl.var = head(myExpl)
#                                            , type = 'simple'
#                                            , interaction.level = 0)
#       
#       if (bm.opt@model == "RF" && !is.null(bm.opt.val$type) && bm.opt.val$type == 'classification') {
#         # defining occurences as factor for doing classification and not regression in RF
#         mySpExpl <- mySpExpl %>% mutate_at(mySp, factor)
#       }
#       
#       ## FILL data parameter ------------------------------------------
#       if (bm.opt@model %in% c("ANN", "CTA", "FDA", "GAM", "GBM", "MARS", "RF", "GLM")) {
#         bm.opt.val$data <- mySpExpl
#       } else if (bm.opt@model == "MAXNET") {
#         bm.opt.val$p <- myResp
#         bm.opt.val$data <- myExpl
#       } else if (bm.opt@model == "SRE") {
#         bm.opt.val$resp.var <- myResp
#         bm.opt.val$expl.var <- myExpl
#       } else if (bm.opt@model == "XGBOOST") {
#         bm.opt.val$label <- myResp
#         bm.opt.val$data <- as.matrix(myExpl)
#       }
#       
#       ## RUN model ----------------------------------------------------
#       print("you")
#       model.sp <- do.call(bm.opt@func, bm.opt.val)
#       print("pi")
#     } else {
#       ## STUFF MAXENT
#     }
#   }
# }


## BIOMOD.options.dataset -----------------------------------------------------
##' 
##' @rdname BIOMOD.options.dataset
##' @export
##' 

setMethod('BIOMOD.options.dataset', signature(strategy = 'character'),
          function(mod, typ, pkg, fun, strategy
                   , user.val = NULL, user.base = NULL
                   , tuning.fun = NULL
                   , bm.format = NULL, calib.lines = NULL)
          {
            cat('\n\t> ', mod, 'options (datatype:', typ, ', package:', pkg, ', function:', fun, ')...')
            
            args <- .BIOMOD.options.dataset.check.args(strategy = strategy, user.val = user.val, user.base = user.base, tuning.fun = tuning.fun
                                                       , bm.format = bm.format, calib.lines = calib.lines)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            BOM <- BIOMOD.options.default(mod, typ, pkg, fun)
            
            argstmp <- BOM@args.default
            ## NEEDED TO WORK !!!! ------------------------------------------------------
            ## SHOULD BE MOVED to place when testing values !! ??
            if (mod == "ANN") { 
              argstmp[["x"]] = NULL
              argstmp$size = 2
            }
            
            if (mod == "CTA") { argstmp$method <- "class" }
            
            if (mod == "FDA") {
              argstmp$dimension = NULL
              argstmp$keep.fitted = NULL
            }
            if (mod == "GAM") {
              argstmp[["x"]] = NULL
              argstmp[["y"]] = NULL
              argstmp$family = binomial(link = 'logit')
              if (pkg == "gam") { argstmp$control = gam::gam.control() }
              if (pkg == "mgcv") {
                argstmp$method = "GCV.Cp"
                argstmp$control = mgcv::gam.control()
              }
            }
            if (mod == "GLM") {
              argstmp$family = binomial(link = 'logit')
              argstmp$control = list()
            }
            if (mod == "MAXNET") { argstmp[["f"]] = NULL }
            if (mod == "RF") {
              argstmp[["x"]] = NULL
              argstmp$mtry = 1
              argstmp$type <- "classification"
            }
            if (mod == "XGBOOST") { argstmp$nrounds = 4 }
            
            argstmp[["..."]] = NULL
            BOM@args.default <- argstmp
            ## SHOULD BE MOVED to place when testing values !! ??
            
            ## SPECIFIC case of formula -------------------------------------------------
            if ("formula" %in% BOM@args.names) {
              if (!is.null(bm.format)) {
                if (is.null(argstmp$formula) || 
                    (length(argstmp$formula) == 1 && nchar(argstmp$formula) == 0) ||
                    argstmp$formula == "formula(data)") {
                  argstmp$formula <- bm_MakeFormula(resp.name = bm.format@sp.name
                                                    , expl.var = head(bm.format@data.env.var)
                                                    , type = ifelse(mod == "GAM" && pkg == "mgcv", 's_smoother'
                                                                    , ifelse(mod == "GLM", 'quadratic'
                                                                             , 'simple'))
                                                    , interaction.level = 0)
                }
              } else {
                warning("No bm.format provided. No definition of formula through bm_MakeFormula.")
              }
            }
            ## ATTENTION : si on ne donne pas bm.format, on n'a pas de formula du coup
            
            ## GET parameter values according to strategy -------------------------------
            if (strategy %in% c("default", "bigboss") || (strategy == "user.defined" && user.base == "bigboss")) {
              if (strategy == "bigboss" || (strategy == "user.defined" && user.base == "bigboss")) {
                # data(OptionsBigboss) # internal data is readily available
                
                val <- OptionsBigboss@options[[paste0(c(mod, typ, pkg, fun), collapse = ".")]]@args.values[['_allData_allRun']]
                for (ii in names(val)) {
                  if (ii != "formula") { argstmp[[ii]] <- val[[ii]] }
                }
              }
              
              argsval <- lapply(expected_CVnames, function(xx) { argstmp })
              names(argsval) <- expected_CVnames
            } 
            
            if (strategy == "user.defined") {
              
              if (!("..." %in% BOM@args.names)) {
                for (CVname in names(user.val)) {
                  .fun_testIfIn(TRUE, paste0("names(user.val[['", CVname, "']])"), names(user.val[[CVname]]), BOM@args.names)
                }
              } else {
                ## ???
              }
              
              argsval <- lapply(expected_CVnames, function(xx) { argstmp })
              names(argsval) <- expected_CVnames
              
              for (CVname in names(user.val)) {
                val <- user.val[[CVname]]
                for (ii in names(val)) { 
                  argsval[[CVname]][[ii]] <- val[[ii]] 
                }
              }
              
            } else if (strategy == "tuned") {
              argsval <- bm_Tuning(model = mod, tuning.fun = tuning.fun, do.formula = TRUE, do.stepAIC = TRUE
                                   , bm.options = BOM, bm.format = bm.format, calib.lines = calib.lines
                                   , metric.eval = ifelse(mod == "MAXENT", "auc.val.avg", "TSS"))
            }
            BOD <- new('BIOMOD.options.dataset', BOM, args.values = argsval)
            
            ## TEST that all given options do not produce errors ------------------------
            # .BIOMOD.options.dataset.test(bm.opt = BOD)
            
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


## Attention ! print only values for _allData_allRun
setMethod('show', signature('BIOMOD.options.dataset'),
          function(object)
          {
            cat('\n\t> ', object@model, 'options (datatype:', object@type, ', package:', object@package, ', function:', object@func, ') :')
            # for (arg in object@args.names) { ## NOT working for bigboss for example, if new parameters
            dataset <- ifelse("_allData_allRun" %in% names(object@args.values)
                              , "_allData_allRun", names(object@args.values)[1])
            cat('\n\t   ( dataset', dataset, ')')

            for (arg in names(object@args.values[[dataset]])) {
              val.def = capture.output(object@args.default[[arg]])
              val.used = capture.output(object@args.values[[dataset]][[arg]])

              cat('\n\t\t- ', arg, "=", sub("\\[1\\] ", "", val.used))
              if (!is.null(val.used) && !is.null(val.def) &&
                  (length(val.used) != length(val.def) || any(val.used != val.def))) {
                cat('   (default:', sub("\\[1\\] ", "", val.def), ')')
              }
            }
            cat("\n")
          }
)

##'
##' @rdname BIOMOD.options.dataset
##' @export
##'

setMethod('print', signature('BIOMOD.options.dataset'),
          function(x, dataset = '_allData_allRun')
          {
            object = x
            cat('\n\t> ', object@model, 'options (datatype:', object@type, ', package:', object@package, ', function:', object@func, ') :')
            dataset <- ifelse(dataset %in% names(object@args.values)
                              , dataset, ifelse("_allData_allRun" %in% names(object@args.values)
                                                , "_allData_allRun", names(object@args.values)[1]))
            cat('\n\t   ( dataset', dataset, ')')
            
            for (arg in names(object@args.values[[dataset]])) {
              val.def = capture.output(object@args.default[[arg]])
              val.used = capture.output(object@args.values[[dataset]][[arg]])
              
              cat('\n\t\t- ', arg, "=", sub("\\[1\\] ", "", val.used))
              if (!is.null(val.used) && !is.null(val.def) && 
                  (length(val.used) != length(val.def) || any(val.used != val.def))) {
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
##' \code{\link{BIOMOD_Modeling}}
##' 
##' @param object a \code{\link{BIOMOD.models.options}} object
##' @param x a \code{\link{BIOMOD.models.options}} object
##' @param dataset a \code{character} corresponding to the name of a dataset contained 
##' in the \code{arg.values} slot of the \code{\link{BIOMOD.options.dataset}} object 
##' for each model
##' 
##' @slot models a \code{vector} containing model names for which options have 
##' been retrieved and defined, must be \code{algo.datatype.package.function}
##' @slot options a \code{list} containing \code{\link{BIOMOD.options.dataset}} 
##' object for each model
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
         representation(models = "character", options = "list"),
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

##'
##' @rdname BIOMOD.models.options
##' @export
##'

setMethod('print', signature('BIOMOD.models.options'),
          function(x, dataset = '_allData_allRun')
          {
            object = x
            .bm_cat("BIOMOD.models.options")
            for (mod in object@models) {
              print(object@options[[mod]], dataset = dataset)
            }
            .bm_cat()
          }
)

# test <- .fun_testIfIn(test, "GLM$test", object@GLM$test, c("AIC", "BIC", "none"))
# test <- .fun_testIfIn(test, "GBM$distribution", object@GBM$distribution, c("bernoulli", "huberized", "multinomial", "adaboost"))
# test <- .fun_testIfIn(test, "CTA$method", object@CTA$method, c("anova", "poisson", "class", "exp"))
# test <- .fun_testIfIn(test, "FDA$method", object@FDA$method, c('polyreg', 'mars', 'bruto'))


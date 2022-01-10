
##' 
##' @importFrom raster subset
##' 

###################################################################################################
## 0. Generic Functions definition 
## Used for different classes 
##    A = BIOMOD.models.out, B = BIOMOD.projection.out, C = BIOMOD.EnsembleModeling.out
###################################################################################################

setGeneric("get_options", function(obj, ...) { standardGeneric("get_options") }) ## A
setGeneric("get_calib_lines", function(obj, ...) { standardGeneric("get_calib_lines") }) ## A
setGeneric("get_formal_data", function(obj, ...) { standardGeneric("get_formal_data") }) ## A

setGeneric("get_projected_models", function(obj, ...) { standardGeneric("get_projected_models") }) ## B
setGeneric("free", function(obj, ...) { standardGeneric("free") }) ## B

setGeneric("get_predictions", function(obj, ...) { standardGeneric("get_predictions") }) ## ABC

setGeneric("get_needed_models", function(obj, ...) { standardGeneric("get_needed_models") }) ## C
setGeneric("get_kept_models", function(obj, ...) { standardGeneric("get_kept_models") }) ## C

setGeneric("get_built_models", function(obj, ...) { standardGeneric("get_built_models") }) ## AC
setGeneric("get_evaluations", function(obj, ...) { standardGeneric("get_evaluations") }) ## AC
setGeneric("get_variables_importance", function(obj, ...) { standardGeneric("get_variables_importance") }) ## AC


###################################################################################################
## 4. BIOMOD.models.out
###################################################################################################

# 4.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.models.out",
         representation(modeling.id = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.computed = 'character',
                        models.failed = 'character',
                        has.evaluation.data = 'logical',
                        rescal.all.models = 'logical',
                        models.evaluation = 'BIOMOD.stored.array',
                        variables.importance = 'BIOMOD.stored.array',
                        models.prediction = 'BIOMOD.stored.array',
                        models.prediction.eval = 'BIOMOD.stored.array',
                        formated.input.data = 'BIOMOD.stored.formated.data',
                        calib.lines = 'BIOMOD.stored.array',
                        models.options = 'BIOMOD.stored.models.options',
                        link = 'character'),
         prototype(modeling.id = as.character(format(Sys.time(), "%s")),
                   sp.name = '',
                   expl.var.names = '',
                   models.computed = '',
                   models.failed = '',
                   has.evaluation.data = FALSE,
                   rescal.all.models = TRUE,
                   models.evaluation = new('BIOMOD.stored.array'),
                   variables.importance = new('BIOMOD.stored.array'),
                   models.prediction = new('BIOMOD.stored.array'),
                   models.prediction.eval = new('BIOMOD.stored.array'),
                   formated.input.data = new('BIOMOD.stored.formated.data'),
                   calib.lines = new('BIOMOD.stored.array'),
                   models.options = new('BIOMOD.stored.models.options'),
                   link = ''),
         validity = function(object){ return(TRUE) } )

setClass("BIOMOD.stored.models.out",
         contains = "BIOMOD.stored.data",
         representation(val = 'BIOMOD.models.out'),
         prototype(val = NULL),
         validity = function(object){ return(TRUE) } )


# 4.3 Other functions -----------------------------------------------------------------------------
# #' @rdname BIOMOD.models.out-objects
# #' @docType method
# #' @aliases show, BIOMOD.models.out-method
setMethod('show', signature('BIOMOD.models.out'),
          function(object)
          {
            .bm_cat("BIOMOD.models.out")
            cat("\nModeling id :", object@modeling.id, fill = .Options$width)
            cat("\nSpecies modeled :", object@sp.name, fill = .Options$width)
            cat("\nConsidered variables :", object@expl.var.names, fill = .Options$width)
            cat("\n\nComputed Models : ", object@models.computed, fill = .Options$width)
            cat("\n\nFailed Models : ", object@models.failed, fill = .Options$width)
            .bm_cat()
          }
)

## ------------------------------------------------------------------------------------------------

setMethod("get_options", "BIOMOD.models.out",
          function(obj) {
            if (obj@models.options@inMemory) {
              return(obj@models.options@val)
            } else if (obj@models.options@link != '') {
              return(get(load(obj@models.options@link)))
            } else { return(NA) }
          }
)

setMethod("get_calib_lines", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, ...) {
            calib_lines <- load_stored_object(obj@calib.lines)
            return(calib_lines)
          }
)

setMethod("get_formal_data", "BIOMOD.models.out",
          function(obj, subinfo = NULL) {
            if (is.null(subinfo)) {
              if (obj@formated.input.data@inMemory) {
                return(obj@formated.input.data@val)
              } else if (obj@formated.input.data@link != '') {
                data <- get(load(obj@formated.input.data@link))
                return(data)
              } else { return(NA) }
            } else if (subinfo == 'MinMax') {
              MinMax = foreach::foreach(i = 1:ncol(env)) %do% {
                x = env[, i]
                if (is.numeric(x)) {
                  return(list(min = min(x, na.rm = T), max = max(x, na.rm = T)))
                } else if (is.factor(x)) {
                  return(list(levels = levels(x)))
                }
              }
              names(MinMax) = colnames(env)
              return(MinMax)
            } else if (subinfo == 'expl.var') {
              return(as.data.frame(get_formal_data(obj)@data.env.var))
            } else if (subinfo == 'expl.var.names') {
              return(obj@expl.var.names)
            } else if (subinfo == 'resp.var') {
              return(as.numeric(get_formal_data(obj)@data.species))
            } else if (subinfo == 'eval.resp.var') {
              return(as.numeric(get_formal_data(obj)@eval.data.species))
            } else if (subinfo == 'eval.expl.var') {
              return(as.data.frame(get_formal_data(obj)@eval.data.env.var))
            } else { stop("Unknown subinfo tag")}
          }
)

## ------------------------------------------------------------------------------------------------

setMethod("get_predictions", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE, evaluation = FALSE)
          {
            # check evaluation data availability
            if (evaluation & (!obj@has.evaluation.data)) {
              warning("calibration data returned because no evaluation data available")
              evaluation = FALSE
            }
            
            # select calibration or eval data
            if (evaluation) { pred <- obj@models.prediction.eval } else { pred <- obj@models.prediction }
            
            if (!as.data.frame) {
              if (pred@inMemory) {
                return(pred@val)
              } else if (pred@link != '') {
                return(get(load(pred@link)))
              } else {
                return(NULL)
              }
            } else {
              if (pred@inMemory) {
                mod.pred <- as.data.frame(pred@val)
              } else if (pred@link != '') {
                mod.pred <- as.data.frame(get(load(pred@link)))
              } else { return(NULL) }
              
              names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred), ".", fixed = TRUE), function(x)
              {
                x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo at the end
                data.set.id <- x.rev[1]
                cross.valid.id <- x.rev[2]
                algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
                model.id <- paste(obj@sp.name,
                                  data.set.id,
                                  cross.valid.id,
                                  algo.id, sep = "_")
                return(model.id)
              }))
              return(mod.pred)
            }
          }
)

setMethod("get_built_models", "BIOMOD.models.out", function(obj, ...) { return(obj@models.computed) })

setMethod("get_evaluations", "BIOMOD.models.out",
          function(obj, ...)
          {
            args <- list(...)
            
            ## fill some additional parameters
            as.data.frame <- ifelse(!is.null(args$as.data.frame), args$as.data.frame, FALSE)
            
            out <- NULL
            if (obj@models.evaluation@inMemory) {
              out <- obj@models.evaluation@val
            } else if(obj@models.evaluation@link != '') {
              out <- get(load(obj@models.evaluation@link))
            }
            
            ## transform into data.frame object if needed
            if(as.data.frame)
            {
              tmp <- reshape::melt.array(out, varnames = c("eval.metric", "test", "m", "r", "d"))
              model_names <- unique(apply(tmp[, c("m", "r", "d"), drop = F], 1, paste, collapse = "_"))
              out <- data.frame() #NULL
              for (mod in model_names)
              {
                m = unlist(strsplit(mod, "_"))[1]
                r = unlist(strsplit(mod, "_"))[2]
                d = unlist(strsplit(mod, "_"))[3]
                ind.mrd = which(tmp$m == m & tmp$r == r & tmp$d == d)
                eval.met = as.character(unique(tmp[ind.mrd, "eval.metric", drop = T]))
                for(em in eval.met)
                {
                  ind.em = intersect(ind.mrd, which(tmp$eval.metric == em))
                  out <- rbind(out, data.frame( Model.name = mod,
                                                Eval.metric = em,
                                                Testing.data = as.numeric( tmp[intersect(ind.em, which(tmp$test == "Testing.data")), "value", drop = T]), 
                                                Evaluating.data = ifelse("Evaluating.data" %in% tmp$test
                                                                         , as.numeric(tmp[intersect(ind.em, which(tmp$test == "Evaluating.data")), "value", drop = T]), NA ), 
                                                Cutoff = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Cutoff")), "value", drop = T]), 
                                                Sensitivity = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Sensitivity")), "value", drop = T]), 
                                                Specificity = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Specificity")), "value", drop = T]))
                  )
                } # end loop on eval metric
              } # end loop on models names
            }
            return(out)
          }
)

setMethod("get_variables_importance", "BIOMOD.models.out",
          function(obj, ...) {
            if (obj@variables.importance@inMemory) {
              return(obj@variables.importance@val)
            } else if (obj@variables.importance@link != '') {
              return(get(load(obj@variables.importance@link)))
            } else { return(NA) }
          }
)



###################################################################################################
## 5. BIOMOD.projection.out
###################################################################################################

# 5.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.projection.out",
         representation(proj.names = 'character',
                        sp.name = 'character',
                        expl.var.names = 'character',
                        models.projected = 'character',
                        scaled.models = 'logical',
                        modeling.object = 'BIOMOD.stored.data',
                        modeling.object.id = 'character',
                        type = 'character',
                        proj = 'BIOMOD.stored.data',
                        xy.coord = 'matrix'),
         prototype(proj.names = '',
                   sp.name='',
                   expl.var.names='',
                   models.projected='',
                   scaled.models=TRUE,
                   modeling.object.id='',
                   type='',
                   xy.coord = matrix()),
         validity = function(object){ return(TRUE) })


# 5.3 Other functions -----------------------------------------------------------------------------
setMethod('plot', signature(x = 'BIOMOD.projection.out', y = "missing"),
          function(x, col = NULL, str.grep = NULL)
          {
            models_selected <- x@models.projected
            if (length(str.grep)) { 
              models_selected <- grep(paste(str.grep, collapse = "|"), models_selected, value = T)
            }
            if (!length(models_selected)) { stop("invalid str.grep arg") }
            
            if(inherits(x@proj, "BIOMOD.stored.raster.stack")){
              requireNamespace("rasterVis")
              
              my.at <- seq(0, 1000, by = 100) ## breaks of color key
              my.labs.at <- seq(0, 1000, by = 250) ## labels placed vertically centered
              my.lab <- seq(0, 1000, by = 250) ## labels
              my.col <- colorRampPalette(c("grey90", "yellow4", "green4"))(100) ## colors
              
              ## try to use levelplot function
              try_plot <- try(levelplot(get_predictions(x, full.name = models_selected),
                                        at = my.at,
                                        margin = T,
                                        col.regions = my.col,
                                        main = paste(x@sp.name, x@proj.names, "projections"),
                                        colorkey = list(labels = list(labels = my.lab, at = my.labs.at))
              ))
              if (!inherits(try_plot, "try-error")) { ## produce plot
                print(try_plot)
              } else { ## try classical plot
                cat("\nrasterVis' levelplot() function failed. Try to call standard raster plotting function.",
                    "It can lead to unoptimal representations.",
                    "You should try to do it by yourself extracting predicions (see : get_predictions() function).",
                    fill = options()$width)
                try_plot <- try(plot(get_predictions(x, full.name = models_selected)))
                if (inherits(try_plot,"try-error")) {
                  cat("\n Plotting function failed.. You should try to do it by yourself!")
                }
              }
            } else if (inherits(x@proj, "BIOMOD.stored.array")) {
              if (ncol(x@xy.coord) != 2) {
                cat("\n ! Impossible to plot projections because xy coordinates are not available !")
              } else {
                multiple.plot(Data = get_predictions(x, full.name = models_selected, as.data.frame = T),
                              coor = x@xy.coord)
              }
            } else {
              cat("\n !  Biomod Projection plotting issue !", fill = .Options$width)
            }
          }
)

##' @rdname BIOMOD.projection.out-objects
##' @docType method
##' @aliases show, BIOMOD.projection.out-method
setMethod('show', signature('BIOMOD.projection.out'),
          function(object){
            .bm_cat("'BIOMOD.projection.out'")
            cat("\nProjection directory :", paste0(object@sp.name, "/", object@proj.names), fill = .Options$width)
            cat("\n")
            cat("\nsp.name :", object@sp.name, fill = .Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill = .Options$width)
            cat("\n")
            cat("\nmodeling id :", object@modeling.object.id , "(", object@modeling.object@link , ")", fill = .Options$width)
            cat("\nmodels projected :", toString(object@models.projected), fill = .Options$width)
            .bm_cat()
          })

## ------------------------------------------------------------------------------------------------

setMethod("get_projected_models", "BIOMOD.projection.out", function(obj){ return(obj@models.projected) })

setMethod('free', signature('BIOMOD.projection.out'),
          function(obj) {
            if (inherits(obj@proj, "BIOMOD.stored.array")) {
              obj@proj@val = array()
            } else if (inherits(obj@proj, "BIOMOD.stored.raster.stack")) {
              obj@proj@val = stack()
            } else {
              obj@proj@val = NULL
            }
            obj@proj@inMemory = FALSE
            return(obj)
          })

setMethod("get_predictions", "BIOMOD.projection.out",
          function(obj, as.data.frame = FALSE, full.name = NULL, model = NULL, run.eval = NULL, data.set = NULL)
          {
            models_selected <- get_projected_models(obj)
            if(length(full.name)){
              models_selected <- intersect(full.name, models_selected)
            } else if(length(model) | length(run.eval) | length(data.set)){
              # models subselection according to model, run.eval and sata.set parameters
              grep_model = grep_run.eval = grep_data.set = "*"
              if(length(model)) { grep_model <- paste0("(", paste(model, collapse = "|"), ")") }
              if (length(run.eval)) { grep_run.eval <- paste0("(", paste(run.eval, collapse = "|"), ")") }
              if (length(data.set)) { grep_data.set <- paste0("(", paste(data.set, collapse = "|"), ")") }
              grep_full <- paste0(grep_data.set, "_", grep_run.eval, "_", grep_model, "$")
              models_selected <- grep(pattern = grep_full, models_selected, value = T)
            }
            
            if (length(models_selected))
            {
              proj <- load_stored_object(obj@proj)
              names(proj) <- obj@models.projected
              if (inherits(proj, 'Raster')) {
                proj <- subset(proj, models_selected, drop = FALSE)
              } else if (length(dim(proj)) == 4) { ## 4D arrays
                proj <- proj[, .extractModelNamesInfo(model.names = models_selected, info = 'models'),
                             .extractModelNamesInfo(model.names = models_selected, info = 'run.eval'),
                             .extractModelNamesInfo(model.names = models_selected, info = 'data.set'), drop = FALSE]
              } else { ## matrix (e.g. from ensemble models projections)
                proj <- proj[, models_selected, drop = FALSE]
              }
              
              if(as.data.frame)
              {
                proj <- as.data.frame(proj)
                ## set correct names
                if (obj@type == 'array' & sum(!(names(proj) %in% models_selected)) > 0) { ## from array & not valid names
                  names(proj) <- unlist(lapply(strsplit(names(proj), ".", fixed = TRUE), function(x)
                  {
                    x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                    data.set.id <- x.rev[1]
                    cross.valid.id <- x.rev[2]
                    algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
                    model.id <- paste(obj@sp.name,
                                      data.set.id,
                                      cross.valid.id,
                                      algo.id, sep = "_")
                    return(model.id)
                  }))
                }
                proj <- proj[, models_selected] # reorder the data.frame
              }
            } else { proj <- NULL }
            return(proj)
          }
)


###################################################################################################
## 6. BIOMOD.EnsembleModeling.out
###################################################################################################

# 6.1 Class Definition ----------------------------------------------------------------------------
##' @name BIOMOD.EnsembleModeling.out-class
##' @rdname BIOMOD.EnsembleModeling.out-objects
##' @docType class
##' 
##' @aliases BIOMOD.EnsembleModeling.out-class
##' @aliases BIOMOD.EnsembleModeling.out
##' 
##' @title BIOMOD_EnsembleModeling() outputs objects class
##' 
##' @description
##' EnsembleModeling objects are created, used and returned by BIOMOD
##' functions. It's contains information relative to an \pkg{biomod2}
##' ensemble modeling procedure.
##' 
##' - output of: \code{\link{BIOMOD_EnsembleModeling}}
##' - input of: \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' @slot sp.name "character", species name
##' @slot expl.var.names "character", explanatory variables
##' names
##' @slot models.out.obj "BIOMOD.stored.models.out", object which
##' contains information on individuals models that have been combined
##' @slot eval.metric "character", evaluation metrics chose for
##' models selection
##' @slot eval.metric.quality.threshold "numeric", thresholds
##' defined for models selection 
##' @slot em.computed "character", ensemble models built names
##' @slot em.by "character", way models are combined
##' @slot em.models "ANY", list of built biomod2.ensemble.models
##' objects
##' @slot modeling.id "character", the id of the whole
##' modelling process
##' @slot link "character", the path to corresponding hard drive
##' saved object
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, 
##' \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{BIOMOD_EnsembleForecasting}}
##' 
##' @keywords models
##' @keywords ensemble
##' @author Damien Georges 
##' 
##' @examples
##' showClass("BIOMOD.EnsembleModeling.out")
##' 
setClass("BIOMOD.EnsembleModeling.out",
         representation(sp.name = 'character',
                        expl.var.names = 'character',
                        models.out.obj = 'BIOMOD.stored.models.out',
                        eval.metric = 'character',
                        eval.metric.quality.threshold = 'numeric',
                        em.computed = 'character',
                        em.by = 'character',
                        em.models = 'ANY',
                        modeling.id = 'character',
                        link = 'character'),
         #                         em.models.kept = 'list',
         #                         em.prediction = 'BIOMOD.stored.array',
         #                         em.evaluation = 'BIOMOD.stored.array',
         #                         em.res = 'list',
         #                         em.ci.alpha = 'numeric',
         #                         em.weight = 'list',
         #                         em.bin.tresh = 'list'),
         prototype(sp.name = '',
                   expl.var.names = '',
                   models.out.obj = new('BIOMOD.stored.models.out'),
                   eval.metric = '',
                   eval.metric.quality.threshold = 0,
                   em.models = NULL,
                   em.computed = character(),
                   modeling.id = '.',
                   link = ''),
         #                     em.models.kept = NULL,
         #                     em.prediction = NULL,
         #                     #                     em.evaluation = NULL,
         #                     em.res = list(),
         #                     em.ci.alpha = 0.05,
         #                     em.weight = list(),
         #                     em.bin.tresh = list()),
         validity = function(object){ return(TRUE) })


# 6.3 Other functions -----------------------------------------------------------------------------
##' @rdname BIOMOD.EnsembleModeling.out-objects
##' @docType method
##' @param object a BIOMOD.EnsembleModeling.out object
setMethod('show', signature('BIOMOD.EnsembleModeling.out'),
          function(object){
            .bm_cat("'BIOMOD.EnsembleModeling.out'")
            cat("\nsp.name :", object@sp.name, fill=.Options$width)
            cat("\nexpl.var.names :", object@expl.var.names, fill=.Options$width)
            cat("\n")
            cat("\nmodels computed:", toString(object@em.computed), fill=.Options$width)
            
            .bm_cat()
          })

## ------------------------------------------------------------------------------------------------

# setMethod("get_needed_models", "BIOMOD.EnsembleModeling.out",
#           function(obj, subset = 'all', ...) {
#             add.args <- list(...)
#             needed_models <- lapply(obj@em.models, function(x) { return(x@model) })
#             needed_models <- unique(unlist(needed_models))
#             return(needed_models)
#           })

setMethod("get_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, selected.models = 'all', ...) {
            add.args <- list(...)
            if (selected.models[[1]] == "all") {
              selected.index <- c(1:length(obj@em.models))
            } else {
              selected.index <- which(names(obj@em.models) %in% selected.models)
            }
            needed_models <- lapply(obj@em.models[selected.index], function(x) { return(x@model) })
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          })


setMethod("get_kept_models", "BIOMOD.EnsembleModeling.out",
          function(obj, model, ...) {
            if (is.character(model) | is.numeric(model)) {
              return(obj@em.models[[model]]@model)
            } else {
              kept_mod <- lapply(obj@em.models, function(x) { return(x@model) })
              names(kept_mod) <- names(obj@em.models)
              return(kept_mod)
            }
          })

## ------------------------------------------------------------------------------------------------

setMethod("get_predictions", "BIOMOD.EnsembleModeling.out",
          function(obj, ...)
          {
            ## note: ensemble models predicitons are stored within the directory
            ##  <sp.name>/.BIOMOD_DATA/<modelling.id>/ensemble.models/ensemble.models.projections/
            ##  This function is just a friendly way to load this data
            
            ## get the path to projections files we want to load
            files.to.load <- file.path(obj@sp.name, ".BIOMOD_DATA", obj@modeling.id, "ensemble.models",
                                       "ensemble.models.predictions", paste0(obj@em.computed, ".predictions"))
            ## load and merge projection files within a data.frame
            bm.pred <- do.call(cbind, lapply(files.to.load, function(ftl) get(load(ftl))))
            colnames(bm.pred) <- obj@em.computed
            return(bm.pred)
          })

setMethod("get_built_models", "BIOMOD.EnsembleModeling.out", function(obj, ...){ return(obj@em.computed) })

setMethod("get_evaluations", "BIOMOD.EnsembleModeling.out",
          function(obj, ...)
          {
            args <- list(...)
            
            ## fill some additional parameters
            as.data.frame <- ifelse(!is.null(args$as.data.frame), args$as.data.frame, FALSE)
            
            ## extract evaluation scores as a list
            out <- list()
            models <- obj@em.computed ## list of computed models
            for (mod in models) {
              out[[mod]] <- obj@em.models[[mod]]@model_evaluation[, , drop = F]
            }
            
            ## transform into data.frame object if needed
            if(as.data.frame)
            {
              tmp <- melt(out, varnames = c("eval.metric", "test"))
              tmp$model.name <- sapply(tmp$L1, function(x) { paste(unlist(strsplit(x, "_"))[-1], collapse = "_") })
              out <- data.frame() #NULL
              for (mod in unique(tmp$model.name))
              {
                eval.met = as.character(unique(tmp[which(tmp$model.name == mod), "eval.metric", drop = T]))
                for(em in eval.met)
                {
                  ind.em = which(tmp$model.name == mod & tmp$eval.metric == em)
                  out <- rbind(out, data.frame( Model.name = mod,
                                                Eval.metric = em,
                                                Testing.data = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Testing.data")), "value", drop = T]), 
                                                Evaluating.data = ifelse("Evaluating.data" %in% tmp$test
                                                                         , as.numeric(tmp[intersect(ind.em, which(tmp$test == "Evaluating.data")), "value", drop =  T]), NA ), 
                                                Cutoff = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Cutoff")), "value", drop = T]), 
                                                Sensitivity = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Sensitivity")), "value", drop = T]), 
                                                Specificity = as.numeric(tmp[intersect(ind.em, which(tmp$test == "Specificity")), "value", drop = T]))
                  )
                } # end loop on eval metric
              } # end loop on models names
            } # end as.data.frame == TRUE
            return(out)
          }
)

setMethod("get_variables_importance", "BIOMOD.EnsembleModeling.out",
          function(obj, ...) {
            vi <- NULL
            for (mod in get_built_models(obj)) {
              vi_tmp <- obj@em.models[[mod]]@model_variables_importance
              vi <- abind::abind(vi, vi_tmp, along = 3)
            }
            dimnames(vi)[[3]] <- get_built_models(obj)
            return(vi)
          })


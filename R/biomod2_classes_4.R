# getters.bm doc ----------------------------------------------------------
##' @name getters.bm
##' @aliases get_formal_model
##' @aliases get_scaling_model
##' @author Damien Georges
##' 
##' @title Functions to extract informations from \code{\link{biomod2_model}} objects
##' 
##' @description These functions allow the user to easily retrieve single models (formal or scaled) 
##' from \code{\link{biomod2_model}} objects from the modeling step.
##' 
##' 
##' @param object a \code{\link{biomod2_model}} object
##' 
##' 
##' @return 
##' 
##' \describe{
##'   \item{\code{get_formal_model}}{an object from the \code{model} slot of a 
##'   \code{\link{biomod2_model}} object}
##'   \item{\code{get_scaling_model}}{an object from the \code{scaling_model} slot of a 
##'   \code{\link{biomod2_model}} object}
##' }
##' 
##' 
##' @seealso \code{\link{biomod2_model}}
##' @family Toolbox functions
##' 
NULL

setGeneric("get_formal_model", def = function(object) { standardGeneric("get_formal_model") })

setGeneric("get_scaling_model", def = function(object) { standardGeneric("get_scaling_model") })


##' @name predict.bm
## @aliases predict
##' @aliases .predict.ANN_biomod2_model.RasterStack
##' @aliases .predict.ANN_biomod2_model.data.frame
##' @aliases .predict.CTA_biomod2_model.RasterStack
##' @aliases .predict.CTA_biomod2_model.data.frame
##' @aliases .predict.FDA_biomod2_model.RasterStack
##' @aliases .predict.FDA_biomod2_model.data.frame
##' @aliases .predict.GAM_biomod2_model.RasterStack
##' @aliases .predict.GAM_biomod2_model.data.frame
##' @aliases .predict.GBM_biomod2_model.RasterStack
##' @aliases .predict.GBM_biomod2_model.data.frame
##' @aliases .predict.GLM_biomod2_model.RasterStack
##' @aliases .predict.GLM_biomod2_model.data.frame
##' @aliases .predict.MARS_biomod2_model.RasterStack
##' @aliases .predict.MARS_biomod2_model.data.frame
##' @aliases .predict.MAXENT.Phillips_biomod2_model.data.frame
##' @aliases .predict.MAXENT.Phillips.2_biomod2_model.RasterStack
##' @aliases .predict.MAXENT.Phillips.2_biomod2_model.data.frame
##' @aliases .predict.RF_biomod2_model.RasterStack
##' @aliases .predict.RF_biomod2_model.data.frame
##' @aliases .predict.SRE_biomod2_model.RasterStack
##' @aliases .predict.SRE_biomod2_model.data.frame
##' @author Damien Georges
##' 
##' @title Functions to get predictions from \code{\link{biomod2_model}} objects
##' 
##' @description This function allows the user to predict single models from 
##' \code{\link{biomod2_model}} on (new) explanatory variables.
##' 
##' 
##' @param object a \code{\link{biomod2_model}} object
##' @param newdata a \code{data.frame} or \code{\link[raster:stack]{RasterStack}} object 
##' containing data for new predictions
##' @param \ldots (\emph{optional)}) 
##' 
##' 
##' @seealso \code{\link{biomod2_model}}
##' @family Toolbox functions
##' 
##' 
##' @importFrom raster raster as.matrix is.factor subset calc writeRaster readAll
##' predict reclassify rasterToPoints cellFromXY inMemory
##' @importFrom sp read.asciigrid
##' @importFrom gbm predict.gbm
##' 
NULL

#setGeneric("predict", def = function(object, ...) { standardGeneric("predict") })


#----------------------------------------------------------------------------- #
# 7. biomod2_model -----------------------------------------------------------
#----------------------------------------------------------------------------- #

##' @name biomod2_model
##' @aliases biomod2_model-class
##' @aliases ANN_biomod2_model
##' @aliases CTA_biomod2_model
##' @aliases FDA_biomod2_model
##' @aliases GAM_biomod2_model
##' @aliases GBM_biomod2_model
##' @aliases GLM_biomod2_model
##' @aliases MARS_biomod2_model
##' @aliases MAXENT.Phillips_biomod2_model
##' @aliases MAXENT.Phillips.2_biomod2_model
##' @aliases RF_biomod2_model
##' @aliases SRE_biomod2_model
##' @author Damien Georges
##' 
##' @title Single model output object class (when running \code{BIOMOD_Modeling()})
##' 
##' @description Class created by \code{\link{BIOMOD_Modeling}} and \code{\link{bm_RunModel}}
##' 
##' 
##' @slot model_name a \code{character} corresponding to the model name
##' @slot model_class a \code{character} corresponding to the model class
##' @slot model_options a \code{list} containing the model options
##' @slot model the corresponding model object
##' @slot scaling_model the corresponding scaled model object
##' @slot dir_name a \code{character} corresponding to the modeling folder
##' @slot resp_name a \code{character} corresponding to the species name
##' @slot expl_var_names a \code{vector} containing names of explanatory variables
##' @slot expl_var_type a \code{vector} containing classes of explanatory variables
##' @slot expl_var_range a \code{list} containing ranges of explanatory variables
##' @slot model_evaluation a \code{matrix} containing the model evaluations
##' @slot model_variables_importance a \code{matrix} containing the model variables importance
##' 
##' @param object a \code{\link{biomod2_model}} object
##' 
##' @details 
##' 
##' \code{biomod2_model} is the basic object for \pkg{biomod2} single species distribution models. 
##' \cr All listed classes below are derived from \code{biomod2_model}, and have a 
##' \code{model_class} slot specific value :
##' 
##' \itemize{
##'   \item{\code{ANN_biomod2_model} : }{\code{model_class} is \code{ANN}}
##'   \item{\code{CTA_biomod2_model} : }{\code{model_class} is \code{CTA}}
##'   \item{\code{FDA_biomod2_model} : }{\code{model_class} is \code{FDA}}
##'   \item{\code{GBM_biomod2_model} : }{\code{model_class} is \code{GBM}}
##'   \item{\code{GLM_biomod2_model} : }{\code{model_class} is \code{GLM}}
##'   \item{\code{MARS_biomod2_model} : }{\code{model_class} is \code{MARS}}
##'   \item{\code{MAXENT.Phillips_biomod2_model} : }{\code{model_class} is \code{MAXENT.Phillips}}
##'   \item{\code{MAXENT.Phillips.2_biomod2_model} : }{\code{model_class} is 
##'   \code{MAXENT.Phillips.2}}
##'   \item{\code{RF_biomod2_model} : }{\code{model_class} is \code{RF}}
##'   \item{\code{SRE_biomod2_model} : }{\code{model_class} is \code{SRE}}
##' }
##' 
##' 
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{bm_RunModel}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("biomod2_model")
##' showClass("ANN_biomod2_model")
##' showClass("CTA_biomod2_model")
##' showClass("FDA_biomod2_model")
##' showClass("GAM_biomod2_model")
##' showClass("GBM_biomod2_model")
##' showClass("GLM_biomod2_model")
##' showClass("MARS_biomod2_model")
##' showClass("MAXENT.Phillips_biomod2_model")
##' showClass("MAXENT.Phillips.2_biomod2_model")
##' showClass("RF_biomod2_model")
##' showClass("SRE_biomod2_model")
##' 
NULL

##' @name biomod2_model-class
##' @rdname biomod2_model
##' @export
##' 

# 7.1 Class Definition ----------------------------------------------------------------------------
setClass('biomod2_model',
         representation(model_name = 'character',
                        model_class = 'character',
                        model_options = 'list',
                        model = 'ANY',
                        scaling_model = 'ANY',
                        dir_name = 'character',
                        resp_name = 'character',
                        expl_var_names = 'character',
                        expl_var_type = 'character',
                        expl_var_range = 'list',
                        model_evaluation = 'matrix',
                        model_variables_importance = 'matrix'),
         prototype = list(model_name = 'mySpecies_DataSet_RunName_myModelClass',
                          model_class = 'myModelClass',
                          model_options = list(),
                          model = list(),
                          scaling_model = list(),
                          dir_name = '.',
                          resp_name = 'mySpecies',
                          expl_var_names = 'myRespVar',
                          expl_var_type = 'unknown',
                          expl_var_range = list(),
                          model_evaluation = matrix(),
                          model_variables_importance = matrix()),
         validity = function(object) { # check that scaler is a glm if it is defined
           if (length(object@scaling_model) > 0 && 
               !inherits(object@scaling_model, c("glm", "lm"))) {
             return(FALSE)
           } else { 
             return(TRUE) 
           }
         })

# 7.2 Getters -------------------------------------------------------------------------------------

##' 
##' @rdname getters.bm
##' @export
##' 

setMethod('get_formal_model', signature('biomod2_model'), 
          function(object) { return(object@model) })

##' 
##' @rdname getters.bm
##' @export
##' 

setMethod('get_scaling_model', signature('biomod2_model'), 
          function(object) { return(object@scaling_model) })


# 7.3 Other Functions -----------------------------------------------------------------------------
##' 
##' @rdname biomod2_model
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('biomod2_model'),
          function(object) {
            .bm_cat("biomod2_model")
            cat("\n\t model name :", object@model_name, fill = .Options$width)
            cat("\n\t model class :", object@model_class, fill = .Options$width)
            cat("\n\t This model", ifelse(length(object@scaling_model), "has", "doesn't have"), "its own scale", fill = .Options$width)
            
            cat("\n")
            cat("\n\t modeling folder:", object@dir_name, fill = .Options$width)
            cat("\n\t response modelled:", object@resp_name, fill = .Options$width)
            cat("\n\n\t explanatory variables used:", fill = .Options$width)
            cat("\n\t", "name", "\t", "type", "\t", "range", fill = .Options$width)
            for (i in 1:length(object@expl_var_names)) {
              cat("\n\t", object@expl_var_names[i],"\t", object@expl_var_type[i], "\t", object@expl_var_range[[i]], fill = .Options$width)
            }
            
            cat("\n")
            cat("\n\t NOTE : ")
            cat("\n\t\t You can access 'formal' model with get_formal_model function")
            cat(ifelse(length(object@scaling_model), "\n\t\t You can access scaling model with get_scaling_model function\n", "\n"))
            
            .bm_cat()
          })


## 7.4 biomod2_model predict method ----------------------------------------

# this method just dispatch onto method predict2 so that it may handle 
# a second argument dispatch (newdata = data.frame or RasterStack)

setMethod('predict', signature(object = 'biomod2_model'),
          function(object, newdata, ...) {
            predict2(object, newdata, ...)
          }
)

## 7.5 biomod2_model predict2 method ----------------------------------------

### generic method definition -----------------------------------------------

setGeneric("predict2", function(object, newdata, ...) { standardGeneric("predict2") }) 

### biomod2_model + Raster  -------------------------------------------------
setMethod('predict2', signature(object = 'biomod2_model', newdata = "RasterStack"),
          function(object, newdata, predfun, seedval = NULL, ...) {
            args <- list(...)
            namefile <- args$namefile
            overwrite <- args$overwrite
            on_0_1000 <- args$on_0_1000
            
            if (is.null(overwrite)) { overwrite <- TRUE }
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            
            set.seed(seedval)
            # eval(parse(text = paste0("proj <- ", predcommand)))
            proj <- predfun(object, newdata)
            
            if (length(get_scaling_model(object))) {
              names(proj) <- "pred"
              proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
            }
            if (on_0_1000) { proj <- round(proj * 1000) }
            
            # save raster on hard drive ?
            if (!is.null(namefile)) {
              cat("\n\t\tWriting projection on hard drive...")
              if (on_0_1000) { ## projections are stored as positive integer
                writeRaster(proj, filename = namefile, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
              } else { ## keep default data format for saved raster
                writeRaster(proj, filename = namefile, overwrite = overwrite)
              }
              proj <- raster(namefile, RAT = FALSE)
            }
            return(proj)
          }
)

### biomod2_model + data.frame ---------------------------------------------

setMethod('predict2', signature(object = 'biomod2_model', newdata = "data.frame"),
          function(object, newdata, predfun, seedval = NULL, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            omit.na <- args$omit.na
            
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            if (is.null(omit.na)) { omit.na <- FALSE }
            
            ## check if na occurs in newdata cause they are not well supported
            if (omit.na) {
              not_na_rows <- apply(newdata, 1, function(x) { sum(is.na(x)) == 0 })
            } else {
              not_na_rows <- rep(T, nrow(newdata))
            }
            
            set.seed(seedval)
            # eval(parse(text = paste0("proj <- ", predcommand)))
            proj <- predfun(object, newdata, not_na_rows)
            
            ## add original NA from formal dataset
            if (sum(!not_na_rows) > 0) {
              tmp <- rep(NA, length(not_na_rows))
              tmp[not_na_rows] <- proj
              proj <- tmp
              rm('tmp')
            }
            
            if (length(get_scaling_model(object))) {
              proj <- data.frame(pred = proj)
              proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5, dat = proj)
            }
            if (on_0_1000) { proj <- round(proj * 1000) }
            
            return(proj)
          }
)

#----------------------------------------------------------------------------- #
# 8 biomod2_model subclass ---------------------------------------------------
#----------------------------------------------------------------------------- #

#----------------------------------------------------------------------------- #
## 8.1 ANN_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('ANN_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'ANN'),
         validity = function(object) { 
           if (!inherits(object@model, "nnet")) {
             return(FALSE) } 
           else {
             return(TRUE)
           }})

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict2', signature(object = 'ANN_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata){
              predict(newdata, get_formal_model(object), type = 'raw')
            }
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

setMethod('predict2', signature(object = 'ANN_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), newdata[not_na_rows, , drop = FALSE], type = 'raw'))            
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

#----------------------------------------------------------------------------- #
## 8.2 CTA_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('CTA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'CTA'),
         validity = function(object) { # check model class
           if (!inherits(object@model, "rpart")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict2', signature(object = 'CTA_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            if (any(object@expl_var_type == "factor")) {
              # CTA with factors cannot yet be predicted on raster 
              stop("\n\t! CTA raster prediction not possible with categorical variables !")
            } 
            predfun <- function(object, newdata){
              predict(newdata, model = get_formal_model(object), type = 'prob', index = 2) 
            }
            
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

setMethod('predict2', signature(object = 'CTA_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'prob')[, 2])
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)



#----------------------------------------------------------------------------- #
## 8.3 FDA_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('FDA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'FDA'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "fda")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict2', signature(object = 'FDA_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata){
              predict(newdata, model = get_formal_model(object), type = 'posterior', index = 2)            
            }
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

setMethod('predict2', signature(object = 'FDA_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'posterior')[, 2])
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, omit.na = TRUE, ...)
          }
)

#----------------------------------------------------------------------------- #
## 8.4 GAM_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('GAM_biomod2_model',
         representation(model_subclass = 'character'), 
         contains = 'biomod2_model',
         prototype = list(model_class = 'GAM', model_subclass = 'GAM_mgcv'),
         validity = function(object) { ## check model class
           if ((!(object@model_subclass %in% c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv'))) ||
               (object@model_subclass %in% c('GAM_mgcv', 'GAM_gam') && !inherits(object@model, c("gam", "Gam"))) ||
               (object@model_subclass == 'BAM_mgcv' && !inherits(object@model, c("bam")))) {
             return(FALSE)
           } else {
             return(TRUE)
           }
         })

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict2', signature(object = 'GAM_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            .load_gam_namespace(object@model_subclass)
            
            predfun <- function(object, newdata){
              .run_pred(object = get_formal_model(object), Prev = 0.5 , dat = newdata)
            }
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

setMethod('predict2', signature(object = 'GAM_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            .load_gam_namespace(object@model_subclass)
            
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(.run_pred(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows, , drop = FALSE])))
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


#----------------------------------------------------------------------------- #
## 8.5 GBM_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('GBM_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GBM', n.trees_optim = 1000),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "gbm")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict2', signature(object = 'GBM_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata){
              predict(newdata, model = get_formal_model(object), fun = predict.gbm, n.trees = object@n.trees_optim, type = 'response')
            }
            
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

setMethod('predict2', signature(object = 'GBM_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), n.trees = object@n.trees_optim, type = 'response'))
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


#----------------------------------------------------------------------------- #
## 8.6 GLM_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('GLM_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GLM'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "glm")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict.bm
##' @export
##' 


setMethod('predict2', signature(object = 'GLM_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata){
              .run_pred(object = get_formal_model(object), Prev = 0.5 , dat = newdata)   
            }
            
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

setMethod('predict2', signature(object = 'GLM_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(.run_pred(object = get_formal_model(object), Prev = 0.5 , dat = as.data.frame(newdata[not_na_rows, , drop = FALSE])))
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


#----------------------------------------------------------------------------- #
## 8.7 MARS_biomod2_model ----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('MARS_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MARS'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, c('earth', 'MARS', 'mars'))) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict.bm
##' @export
##' 

# predcommand <- ".run_pred(object = get_formal_model(object), Prev = 0.5 , dat = newdata)"
# # redirect to predict2.biomod2_model.RasterStack
# callNextMethod(object, newdata, predfun = predfun, ...)

setMethod('predict2', signature(object = 'MARS_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            
            args <- list(...)
            namefile <- args$namefile
            overwrite <- args$overwrite
            on_0_1000 <- args$on_0_1000
            seedval <- args$seedval
            
            if (is.null(overwrite)) { overwrite <- TRUE }
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            
            ## handle separately rasterstack depending on the presence or not of factorial variable
            fact.var <- is.factor(newdata)
            set.seed(seedval)
            if (any(fact.var)) {
              ## get factor levels
              fact.var.levels <- subset(levels(newdata), fact.var)
              proj <- calc(newdata, function(x)
              {
                xx <- data.frame(x)
                ## ensure that the data.frame has the right set of levels
                for (i in which(fact.var)) {
                  xx[[i]] <- factor(xx[[i]], levels = unlist(fact.var.levels[[i]]))
                }
                ## do the projection
                proj.out <- as.numeric(predict(get_formal_model(object), xx, type = 'response'))
                return(proj.out)
              })
            } else {
              proj <- predict(newdata, model = get_formal_model(object), type = 'response')
            }
            
            if (length(get_scaling_model(object))) {
              names(proj) <- "pred"
              proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
            }
            
            if (on_0_1000) { proj <- round(proj * 1000)}
            
            # save raster on hard drive ?
            if (!is.null(namefile)) {
              cat("\n\t\tWriting projection on hard drive...")
              if (on_0_1000) { ## projections are stored as positive integer
                writeRaster(proj, filename = namefile, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
              } else { ## keep default data format for saved raster
                writeRaster(proj, filename = namefile, overwrite = overwrite)
              }
              proj <- raster(namefile, RAT = FALSE)
            }
            return(proj)
          }
)


setMethod('predict2', signature(object = 'MARS_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'response'))
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

#----------------------------------------------------------------------------- #
## 8.8 MAXENT.Phillips_biomod2_model -----------------------------------------
#----------------------------------------------------------------------------- #

setClass('MAXENT.Phillips_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXENT.Phillips'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict', signature(object = 'MAXENT.Phillips_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "MAXENT.Phillips", object, newdata, ...))
          })

.predict.MAXENT.Phillips_biomod2_model.data.frame <- function(object, newdata, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  temp_workdir <- args$temp_workdir
  do_raster <- args$do_raster
  newraster <- args$newraster
  
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(do_raster)) { do_raster <- FALSE }
  
  ## check if na occurs in newdata cause they are not well supported
  not_na_rows <- apply(newdata, 1, function(x){ sum(is.na(x)) == 0 })
  newdata = as.data.frame(newdata[not_na_rows, , drop = FALSE])
  
  ## Prediction data
  Pred_swd <- read.csv(file.path(temp_workdir, "Predictions/Pred_swd.csv"))
  if (nrow(Pred_swd) != nrow(newdata)) {
    tmp = newdata[, 1:3]
    colnames(tmp) = c("predict", "x", "y")
    Pred_swd = cbind(tmp, newdata)
  } else {
    Pred_swd <- cbind(Pred_swd[, 1:3], newdata)
  }
  # m_predictFile <- file.path(temp_workdir, "Predictions/Pred_swdBis.csv")
  m_predictFile <- file.path(temp_workdir, "Predictions", paste0("Pred_swdBis_", sample(1:100000, 1), ".csv"))
  while (file.exists(m_predictFile)) {
    m_predictFile <- file.path(temp_workdir, "Predictions", paste0("Pred_swdBis_", sample(1:100000, 1), ".csv"))
  }
  write.table(Pred_swd, file = m_predictFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  
  # checking maxent.jar is present
  path_to_maxent.jar <- file.path(object@model_options$path_to_maxent.jar, "maxent.jar")
  if (!file.exists(path_to_maxent.jar)) {
    path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
  }
  
  # cat("\n\t\tRunning Maxent...")
  maxent.command <- paste0("java ", ifelse(is.null(object@model_options$memory_allocated), "", paste0("-mx", object@model_options$memory_allocated, "m")),
                           " -cp ", "\"", path_to_maxent.jar, "\"",
                           " density.Project ",
                           "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                           "\"", m_predictFile, "\" ",
                           "\"", file.path(temp_workdir, "projMaxent.asc") , "\" ",
                           "doclamp=false visible=false autorun nowarnings notooltips")
  system(command = maxent.command, wait = TRUE, intern = TRUE)
  
  # cat("\n\t\tReading Maxent outputs...")
  proj <- as.numeric(read.asciigrid(file.path(temp_workdir, "projMaxent.asc"))@data[, 1])
  
  if (do_raster) {
    newraster[which(newraster[] == 1)] = proj
    proj <- newraster
    
    if (!inMemory(proj)) {
      proj <- readAll(proj) # to prevent from tmp files removing
    }
  } else {
    ## add original NA from formal dataset
    if (sum(!not_na_rows) > 0) {
      tmp <- rep(NA, length(not_na_rows))
      tmp[not_na_rows] <- proj
      proj <- tmp
      rm('tmp')
    }
  }
  
  if (on_0_1000) { proj <- round(proj * 1000) }
  return(proj)
}


#----------------------------------------------------------------------------- #
## 8.9 MAXENT.Phillips.2_biomod2_model ---------------------------------------
#----------------------------------------------------------------------------- #

setClass('MAXENT.Phillips.2_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXENT.Phillips.2'),
         validity = function(object) { 
           if (!inherits(object@model, "maxnet")) {
             return(FALSE)
           } else { 
             return(TRUE) 
           }
         })

##' 
##' @rdname predict.bm
##' @export
##' 

setMethod('predict', signature(object = 'MAXENT.Phillips.2_biomod2_model'),
          function(object, newdata, ...)
          {
            return(.template_predict(mod = "MAXENT.Phillips.2", object, newdata, ...))
          })

.predict.MAXENT.Phillips.2_biomod2_model.RasterStack <- function(object, newdata, ...)
{
  newdata.df <- newdata %>% as.matrix()
  
  args <- list(...)
  namefile <- args$namefile
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  seedval <- args$seedval
  
  if (is.null(overwrite)) { overwrite <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  set.seed(seedval)
  proj <- predict(object = get_formal_model(object), newdata = newdata.df, clamp = FALSE, type = 'logistic')[, 1]
  
  if (length(get_scaling_model(object))) {
    proj.to.scale <- data.frame(pred = proj)
    proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5 , dat = proj.to.scale)
  }
  
  if (on_0_1000) { proj <- round(proj * 1000) }
  
  ## convert back to raster file
  proj.ras <- raster(newdata)
  proj.ras[apply(newdata.df, 1, function(.x) { all(!is.na(.x)) })] <- proj
  proj <- proj.ras
  
  # save raster on hard drive ?
  if (!is.null(namefile)) {
    cat("\n\t\tWriting projection on hard drive...")
    if (on_0_1000) { ## projections are stored as positive integer
      writeRaster(proj, filename = namefile, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename = namefile, overwrite = overwrite)
    }
    proj <- raster(namefile, RAT = FALSE)
  }
  
  return(proj)
}

.predict.MAXENT.Phillips.2_biomod2_model.data.frame <- function(object, newdata, ...)
{
  return(.template_predict.data.frame(predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'logistic')[, 1])", object, newdata, ...))
}


#----------------------------------------------------------------------------- #
## 8.9 MAXENT.Tsuruoka_biomod2_model -----------------------------------------
#----------------------------------------------------------------------------- #

# setClass('MAXENT.Tsuruoka_biomod2_model',
#          representation(),
#          contains = 'biomod2_model',
#          prototype = list(model_class = 'MAXENT.Tsuruoka'),
#          validity = function(object) { ## check model class
#            if( sum(!(c("maxent") %in% class(object@model))) > 0) { return(FALSE) } else { return(TRUE) }
#          })
# 
# setMethod('predict', signature(object = 'MAXENT.Tsuruoka_biomod2_model'),
#           function(object, newdata, ...)
#           {
#             return(.template_predict(mod = "MAXENT.Tsuruoka", object, newdata, ...))
#           })
# 
# .predict.MAXENT.Tsuruoka_biomod2_model.RasterStack <- function(object, newdata, ...)*
# {
#   args <- list(...)
#   namefile <- args$namefile
#   overwrite <- args$overwrite
#   on_0_1000 <- args$on_0_1000
# 
#   if (is.null(overwrite)) { overwrite <- TRUE }
#   if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
# 
#   proj <- calc(newdata, function(x) {
#     proj.out <- rep(NA, nrow(x))
#     x.no.na <- na.omit(x)
#     if(nrow(x.no.na)){
#       proj.not.na <- as.numeric(predict.maxent(get_formal_model(object), x.no.na)[, '1'])
#       proj.out[-attr(x.no.na, "na.action")] <- proj.not.na
#     }
#     return(proj.out)
#     })
# 
#   if (length(get_scaling_model(object))) {
#     names(proj) <- "pred"
#     proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
#   }
# 
#   if (on_0_1000) { proj <- round(proj * 1000) }
# 
#   # save raster on hard drive ?
#   if (!is.null(namefile)) {
#     cat("\n\t\tWriting projection on hard drive...")
#     if (on_0_1000) { ## projections are stored as positive integer
#       writeRaster(proj, filename = namefile, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
#     } else { ## keep default data format for saved raster
#       writeRaster(proj, filename = namefile, overwrite = overwrite)
#     }
#     proj <- raster(namefile, RAT = FALSE)
#   }
#   return(proj)
# }
# 
# .predict.MAXENT.Tsuruoka_biomod2_model.data.frame <- function(object, newdata, ...)
# {
#   return(.template_predict.data.frame(seedval = NULL, predcommand = "as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]))[,'1'])", object, newdata, ...))
# }


#----------------------------------------------------------------------------- #
## 8.10 RF_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('RF_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'RF'),
         validity = function(object) { # check model class
           if (!inherits(object@model, "randomForest")) { return(FALSE)} else { return(TRUE) }
         })

##' 
##' @rdname predict.bm
##' @export
##' 


setMethod('predict2', signature(object = 'RF_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata){
              predict(newdata, model = get_formal_model(object), type = 'prob', index = 2)            
              }
            
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

setMethod('predict2', signature(object = 'RF_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'prob')[, '1'])        
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


#----------------------------------------------------------------------------- #
## 8.11 SRE_biomod2_model ----------------------------------------------------
#----------------------------------------------------------------------------- #

setClass('SRE_biomod2_model',
         representation(extremal_conditions = 'data.frame'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'SRE'),
         validity = function(object){ return(TRUE) })

##' 
##' @rdname predict.bm
##' @export
##' 


setMethod('predict2', signature(object = 'SRE_biomod2_model', newdata = "RasterStack"),
          function(object, newdata, ...) {

            predfun <- function(object, newdata){
              .sre_projection(new.env = newdata, extrem.cond = object@extremal_conditions)
              }
            
            # redirect to predict2.biomod2_model.RasterStack
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

setMethod('predict2', signature(object = 'SRE_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            seedval <- args$seedval
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            
            set.seed(seedval)
            proj <- .sre_projection(new.env = newdata, extrem.cond = object@extremal_conditions)
            if (on_0_1000) { proj <- round(proj * 1000) }
            return(proj)
          }
)



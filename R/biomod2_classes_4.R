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


# predict.bm doc ----------------------------------------------------------
##' @name predict.bm
## @aliases predict
##' @aliases predict.biomod2_model
##' @author Damien Georges
##' 
##' @title Functions to get predictions from \code{\link{biomod2_model}} objects
##' 
##' @description This function allows the user to predict single models from 
##' \code{\link{biomod2_model}} on (new) explanatory variables.
##' 
##' 
##' @param object a \code{\link{biomod2_model}} object
##' @param newdata a \code{data.frame} or
##' \code{\link[terra:SpatRaster]{SpatRaster}} object 
##' containing data for new predictions
##' @param \ldots (\emph{optional}) 
##' 
##' 
##' @seealso \code{\link{biomod2_model}}
##' @family Toolbox functions
##' 
NULL

#setGeneric("predict", def = function(object, ...) { standardGeneric("predict") })


#----------------------------------------------------------------------------- #
# 7. biomod2_model class doc -------------------------------------------------
#----------------------------------------------------------------------------- #

##' @name biomod2_model
##' @aliases biomod2_model-class
##' @aliases ANN_biomod2_model-class
##' @aliases CTA_biomod2_model-class
##' @aliases FDA_biomod2_model-class
##' @aliases GAM_biomod2_model-class
##' @aliases GBM_biomod2_model-class
##' @aliases GLM_biomod2_model-class
##' @aliases MARS_biomod2_model-class
##' @aliases MAXENT_biomod2_model-class
##' @aliases MAXNET_biomod2_model-class
##' @aliases RF_biomod2_model-class
##' @aliases SRE_biomod2_model-class
##' @aliases XGBOOST_biomod2_model-class
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
##' @slot model_evaluation a \code{data.frame} containing the model evaluations
##' @slot model_variables_importance a \code{data.frame} containing the model variables importance
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
##'   \item \code{ANN_biomod2_model} : \code{model_class} is \code{ANN}
##'   \item \code{CTA_biomod2_model} : \code{model_class} is \code{CTA}
##'   \item \code{FDA_biomod2_model} : \code{model_class} is \code{FDA}
##'   \item \code{GBM_biomod2_model} : \code{model_class} is \code{GBM}
##'   \item \code{GLM_biomod2_model} : \code{model_class} is \code{GLM}
##'   \item \code{MARS_biomod2_model} : \code{model_class} is \code{MARS}
##'   \item \code{MAXENT_biomod2_model} : \code{model_class} is \code{MAXENT}
##'   \item \code{MAXNET_biomod2_model} : \code{model_class} is 
##'   \code{MAXNET}
##'   \item \code{RF_biomod2_model} : \code{model_class} is \code{RF}
##'   \item \code{SRE_biomod2_model} : \code{model_class} is \code{SRE}
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
##' showClass("MAXENT_biomod2_model")
##' showClass("MAXNET_biomod2_model")
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
                        model_options = 'BIOMOD.options.dataset',
                        model = 'ANY',
                        scaling_model = 'ANY',
                        dir_name = 'character',
                        resp_name = 'character',
                        expl_var_names = 'character',
                        expl_var_type = 'character',
                        expl_var_range = 'list',
                        model_evaluation = 'data.frame',
                        model_variables_importance = 'data.frame'),
         prototype = list(model_name = 'mySpecies_DataSet_RunName_myModelClass', ## REMOVE prototype ??
                          model_class = 'myModelClass',
                          model_options = new('BIOMOD.options.dataset'),
                          model = list(),
                          scaling_model = list(),
                          dir_name = '.',
                          resp_name = 'mySpecies',
                          expl_var_names = 'myRespVar',
                          expl_var_type = 'unknown',
                          expl_var_range = list(),
                          model_evaluation = data.frame(),
                          model_variables_importance = data.frame()),
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


# 7.3 Show method  -------------------------------------------------------------
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
# a second argument dispatch (newdata = data.frame or SpatRaster)
# the methods is used both for simple model (class biomod2_model) as well as
# ensemble models (class biomod2_ensemble_model)
##' @rdname predict.bm
##' @export

setMethod('predict', signature(object = 'biomod2_model'),
          function(object, newdata, ...) {
            predict2(object, newdata, ...)
          }
)

## 7.5 biomod2_model predict2 method ----------------------------------------

### generic method definition -----------------------------------------------
### method used both for biomod2_model and biomod2_ensemble_model
### predict2.bm doc ----------------------------------------------------------
##' @name predict2.bm
##' @aliases predict2
##' @aliases predict2.biomod2_model.SpatRaster
##' @aliases predict2.biomod2_model.data.frame
##' @aliases predict2.ANN_biomod2_model.SpatRaster
##' @aliases predict2.ANN_biomod2_model.data.frame
##' @aliases predict2.CTA_biomod2_model.SpatRaster
##' @aliases predict2.CTA_biomod2_model.data.frame
##' @aliases predict2.FDA_biomod2_model.SpatRaster
##' @aliases predict2.FDA_biomod2_model.data.frame
##' @aliases predict2.GAM_biomod2_model.SpatRaster
##' @aliases predict2.GAM_biomod2_model.data.frame
##' @aliases predict2.GBM_biomod2_model.SpatRaster
##' @aliases predict2.GBM_biomod2_model.data.frame
##' @aliases predict2.GLM_biomod2_model.SpatRaster
##' @aliases predict2.GLM_biomod2_model.data.frame
##' @aliases predict2.MARS_biomod2_model.SpatRaster
##' @aliases predict2.MARS_biomod2_model.data.frame
##' @aliases predict2.MAXENT_biomod2_model.data.frame
##' @aliases predict2.MAXNET_biomod2_model.SpatRaster
##' @aliases predict2.MAXNET_biomod2_model.data.frame
##' @aliases predict2.RF_biomod2_model.SpatRaster
##' @aliases predict2.RF_biomod2_model.data.frame
##' @aliases predict2.SRE_biomod2_model.SpatRaster
##' @aliases predict2.SRE_biomod2_model.data.frame
##' @author Remi Patin
##' 
##' @title Functions to get predictions from \code{\link{biomod2_model}} objects
##' 
##' @description Internal S4 method used to predict single models from 
##' \code{\link{biomod2_model}} on (new) explanatory variables. \code{predict2} 
##' was introduced to allow a signature with two arguments : \code{object}, 
##' a type of \code{\link{biomod2_model}} and \code{newdata}, either a
##' \code{\link[terra:SpatRaster]{SpatRaster}} or a \code{data.frame}.
##' 
##' 
##' @param object a \code{\link{biomod2_model}} object
##' @param newdata a \code{data.frame} or
##'   \code{\link[terra:SpatRaster]{SpatRaster}} object containing data for new
##'   predictions
##' @param predfun a \code{function}, generated by the \code{predict2} method
##'   specific to each \code{\link{biomod2_model}} subclass and used within the
##'   generic \code{predict2.biomod2_model.SpatRaster} or
##' \code{predict2.biomod2_model.data.frame} to do the prediction
##' @param seedval (\emph{optional, default} \code{NULL}) \cr An \code{integer}
##'   value corresponding to the new seed value to be set
##' @param \ldots (\emph{optional)}) 
##' 
##' 
##' @seealso \code{\link{biomod2_model}}
##' @family Toolbox functions
##' 
##' 
##' @importFrom terra rast as.matrix is.factor subset writeRaster predict cellFromXY inMemory classify
##' @importFrom gbm predict.gbm
##' @importFrom methods callNextMethod
##' @keywords internal

setGeneric("predict2", function(object, newdata, ...) {
  standardGeneric("predict2") 
}) 

### biomod2_model + SpatRaster  -------------------------------------------------
##' @rdname predict2.bm

setMethod('predict2', signature(object = 'biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, predfun, seedval = NULL, ...) {
            args <- list(...)
            filename <- args$filename
            overwrite <- args$overwrite
            on_0_1000 <- args$on_0_1000
            mod.name <- args$mod.name
            
            if (is.null(overwrite)) { overwrite <- TRUE }
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            set.seed(seedval)
            proj <- predfun(object, newdata, mod.name)
            
            if (length(get_scaling_model(object)) > 0) {
              names(proj) <- "pred"
              proj <- .run_pred(object = get_scaling_model(object), 
                                Prev = 0.5 , 
                                dat = proj, 
                                mod.name = mod.name)
            }
            if (on_0_1000) { 
              proj <- round(proj * 1000) 
            }
            
            # save raster on hard drive ?
            if (!is.null(filename)) {
              cat("\n\t\tWriting projection on hard drive...")
              if (on_0_1000) { ## projections are stored as positive integer
                writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
              } else { ## keep default data format for saved raster
                writeRaster(proj, filename = filename, overwrite = overwrite)
              }
              proj <- rast(filename)
            }
            return(proj)
          }
)

### biomod2_model + data.frame ---------------------------------------------
##' @rdname predict2.bm

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
              not_na_rows <- rep(TRUE, nrow(newdata))
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
            
            if (length(get_scaling_model(object)) > 0) {
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
##' @name ANN_biomod2_model-class
##' @rdname biomod2_model
##' @export

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
##' @rdname predict2.bm
##' 

setMethod('predict2', signature(object = 'ANN_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              predict(newdata, get_formal_model(object), type = 'raw', wopt = list(names = mod.name))
            }
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.bm
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
##' @name CTA_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('CTA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'CTA'),
         validity = function(object) { # check model class
           if (!inherits(object@model, "rpart")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict2.bm
##' 

setMethod('predict2', signature(object = 'CTA_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              proj <- 
                subset(predict(newdata,
                               model = get_formal_model(object), 
                               type = 'prob', na.rm = TRUE,
                               wopt = list(names = rep(mod.name,2))), 
                       2)    
              proj
            }
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.bm
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
##' @name FDA_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('FDA_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'FDA'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "fda")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict2.bm
##' 

setMethod('predict2', signature(object = 'FDA_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              # new predict command used with terra
              proj <- 
                subset(
                  predict(newdata, 
                          model = get_formal_model(object), 
                          type = 'posterior',
                          na.rm = TRUE,
                          wopt = list(names = mod.name)), 
                  2)   
              # datamask <- classify(any(is.na(newdata)), 
              # matrix(c(0,0,1,NA),ncol = 2, byrow = TRUE))
              # mask(proj, datamask)
              
              # old predict function used with raster
              # predict(newdata, model = get_formal_model(object), type = 'posterior', index = 2) 
              proj
            }
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.bm
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
##' @name GAM_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('GAM_biomod2_model',
         representation(model_subclass = 'character'), 
         contains = 'biomod2_model',
         prototype = list(model_class = 'GAM', model_subclass = 'GAM_mgcv_gam'),
         validity = function(object) { ## check model class
           if ((!(object@model_subclass %in% c('GAM_mgcv_gam', 'GAM_gam_gam', 'GAM_mgcv_bam'))) ||
               (object@model_subclass %in% c('GAM_mgcv_gam', 'GAM_gam_gam') && !inherits(object@model, c("gam", "Gam"))) ||
               (object@model_subclass == 'GAM_mgcv_bam' && !inherits(object@model, c("bam")))) {
             return(FALSE)
           } else {
             return(TRUE)
           }
         })

##' 
##' @rdname predict2.bm
##' 

setMethod('predict2', signature(object = 'GAM_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            .load_gam_namespace(object@model_subclass)
            
            predfun <- function(object, newdata, mod.name){
              .run_pred(object = get_formal_model(object), Prev = 0.5 , dat = newdata, mod.name = mod.name)
            }
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

##' @rdname predict2.bm
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
##' @name GBM_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('GBM_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GBM', n.trees_optim = 1000),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "gbm")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict2.bm
##' 

setMethod('predict2', signature(object = 'GBM_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              proj <- predict(newdata,
                              model = get_formal_model(object),
                              fun = predict.gbm,
                              n.trees = object@n.trees_optim, 
                              type = 'response',
                              na.rm = TRUE,
                              wopt = list(names = mod.name))
              # datamask <- classify(any(is.na(newdata)),
              #                      matrix(c(0,0,1,NA),ncol = 2, byrow = TRUE))
              # mask(proj, datamask)
              proj
            }
            
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

##' @rdname predict2.bm
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
##' @name GLM_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('GLM_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'GLM'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "glm")) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict2.bm
##' 


setMethod('predict2', signature(object = 'GLM_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              .run_pred(object = get_formal_model(object), Prev = 0.5 , dat = newdata, mod.name = mod.name)   
            }
            
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.bm
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
##' @name MARS_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('MARS_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MARS'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, c('earth', 'MARS', 'mars'))) { return(FALSE) } else { return(TRUE) }
         })

##' 
##' @rdname predict2.bm
##' 

# predcommand <- ".run_pred(object = get_formal_model(object), Prev = 0.5 , dat = newdata)"
# # redirect to predict2.biomod2_model.SpatRaster
# callNextMethod(object, newdata, predfun = predfun, ...)

setMethod('predict2', signature(object = 'MARS_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              predict(newdata, model = get_formal_model(object),
                      type = 'response', wopt = list(names = mod.name))
            }
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)


##' @rdname predict2.bm
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
## 8.8 MAXENT_biomod2_model -----------------------------------------
#----------------------------------------------------------------------------- #
##' @name MAXENT_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('MAXENT_biomod2_model',
         representation(model_output_dir = 'character'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXENT'),
         validity = function(object) { return(TRUE) })

##' 
##' @rdname predict2.bm
##' @importFrom terra rast as.points crds values
##' 

setMethod('predict2', signature(object = 'MAXENT_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            temp_workdir <- args$temp_workdir
            filename <- args$filename
            overwrite <- args$overwrite
            mod.name <- args$mod.name
            
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            
            pa <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "PA")
            run <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "run")
            dataset <- paste0("_", pa, "_", run)
            
            # Proj Data
            vec_data_filename <- foreach(thislayername = names(newdata), .combine = 'c') %do% {
              current_data_filename <-
                file.path(temp_workdir, paste0(thislayername,'.asc'))
              writeRaster(subset(newdata,thislayername), 
                          filename = current_data_filename,
                          overwrite = overwrite,
                          NAflag = -9999)
              return(current_data_filename)
            }
            
            # checking maxent.jar is present
            path_to_maxent.jar <- file.path(object@model_options@args.values[[dataset]]$path_to_maxent.jar, "maxent.jar")
            if (!file.exists(path_to_maxent.jar)) {
              path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
            }
            
            maxent.command <- 
              paste0("java ",
                     ifelse(is.null(object@model_options@args.values[[dataset]]$memory_allocated), "",
                            paste0("-mx", object@model_options@args.values[[dataset]]$memory_allocated, "m")),
                     ifelse(is.null(object@model_options@args.values[[dataset]]$initial_heap_size), "",
                            paste0(" -Xms", object@model_options@args.values[[dataset]]$initial_heap_size)),
                     ifelse(is.null(object@model_options@args.values[[dataset]]$max_heap_size), "",
                            paste0(" -Xmx", object@model_options@args.values[[dataset]]$max_heap_size)),
                     " -cp ", "\"", path_to_maxent.jar, "\"",
                     " density.Project ",
                     "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                     "\"", temp_workdir, "\" ",
                     "\"", file.path(temp_workdir, "projMaxent.asc"), "\" ",
                     " doclamp=false visible=false autorun nowarnings notooltips")
            system(command = maxent.command, wait = TRUE, intern = TRUE)
            
            file.remove(vec_data_filename)
            proj.spdf <- read.asciigrid(file.path(temp_workdir, "projMaxent.asc"))
            names(proj.spdf) <- mod.name
            proj <- rast(proj.spdf)
            
            # Remi 11/2022 Not sure the following lines are necessary
            if (!inMemory(proj)) {
              if (!isNamespaceLoaded("raster")) {
                if(!requireNamespace('raster', quietly = TRUE)) stop("Package 'raster' not found")
              }
              
              proj <- raster::readAll(proj@ptr) # to prevent from tmp files removing
              x <- message(proj, "readAll") # to have message if need be ?
            }
            
            if (on_0_1000) { 
              proj <- round(proj * 1000) 
            }
            
            # save raster on hard drive ?
            if (!is.null(filename)) {
              cat("\n\t\tWriting projection on hard drive...")
              if (on_0_1000) { ## projections are stored as positive integer
                writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
              } else { ## keep default data format for saved raster
                writeRaster(proj, filename = filename, overwrite = overwrite)
              }
              proj <- rast(filename)
            }
            proj
          }
)


##' 
##' @rdname predict2.bm
##' @importFrom sp read.asciigrid
##' 

setMethod('predict2', signature(object = 'MAXENT_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            
            args <- list(...)
            on_0_1000 <- args$on_0_1000
            temp_workdir <- args$temp_workdir
            
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            
            pa <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "PA")
            run <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "run")
            dataset <- paste0("_", pa, "_", run)
            
            ## check if na occurs in newdata cause they are not well supported
            not_na_rows <- apply(newdata, 1, function(x){ all(!is.na(x)) })
            if(!all(not_na_rows)){
              newdata  <-  newdata[not_na_rows, , drop = FALSE]
            }
            
            
            # get categorical variables and transform them into numeric
            categorical_var <- .get_categorical_names(newdata)
            newdata <- .categorical2numeric(newdata, categorical_var)
            ## Prediction data
            ## when newdata have as many rows as calibration prediction
            ## the function re-use the first three columns of Pred_swd (predict, 
            ## x and y). Otherwise it uses column from newdata to generate the
            ## given columns. Those columns are likely ignored by maxent so 
            ## there is no need to provide meaningful values.
            Pred_swd <- read.csv(file.path(temp_workdir, "Predictions/Pred_swd.csv"))
            if (nrow(Pred_swd) != nrow(newdata)) {
              tmp = cbind("predict", newdata[, 1], newdata[, 1])
              colnames(tmp) = c("predict", "x", "y")
              Pred_swd = cbind(tmp, newdata)
            } else {
              Pred_swd <- cbind(Pred_swd[, 1:3], newdata)
            }
            m_predictFile <- file.path(temp_workdir, "Predictions", paste0("Pred_swdBis_", sample(1:100000, 1), ".csv"))
            while (file.exists(m_predictFile)) {
              m_predictFile <- file.path(temp_workdir, "Predictions", paste0("Pred_swdBis_", sample(1:100000, 1), ".csv"))
            }
            write.table(Pred_swd, file = m_predictFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
            
            # checking maxent.jar is present
            path_to_maxent.jar <- file.path(object@model_options@args.values[[dataset]]$path_to_maxent.jar, "maxent.jar")
            if (!file.exists(path_to_maxent.jar)) {
              path_to_maxent.jar <-  file.path(getwd(), "maxent.jar")
            }
            
            maxent.command <- 
              paste0("java ",
                     ifelse(is.null(object@model_options@args.values[[dataset]]$memory_allocated), "",
                            paste0("-mx", object@model_options@args.values[[dataset]]$memory_allocated, "m")),
                     ifelse(is.null(object@model_options@args.values[[dataset]]$initial_heap_size), "",
                            paste0(" -Xms", object@model_options@args.values[[dataset]]$initial_heap_size)),
                     ifelse(is.null(object@model_options@args.values[[dataset]]$max_heap_size), "",
                            paste0(" -Xmx", object@model_options@args.values[[dataset]]$max_heap_size)),
                     " -cp ", "\"", path_to_maxent.jar, "\"",
                     " density.Project ",
                     "\"", list.files(path = object@model_output_dir, pattern = ".lambdas$", full.names = TRUE), "\" ",
                     "\"", m_predictFile, "\" ",
                     "\"", file.path(temp_workdir, "projMaxent.asc") , "\" ",
                     "doclamp=false visible=false autorun nowarnings notooltips")
            system(command = maxent.command, wait = TRUE, intern = TRUE)
            
            
            # remove maxent data files
            file.remove(m_predictFile)
            
            
            # cat("\n\t\tReading Maxent outputs...")
            # As of 23/11/2022 rast does not seem to work for large asciigrid. 
            # So dependence to sp was added again.
            # proj <- as.numeric(values(rast(file.path(temp_workdir, "projMaxent.asc")), mat = FALSE))
            proj <- as.numeric(read.asciigrid(file.path(temp_workdir, "projMaxent.asc"))@data[, 1])
            ## add original NA from formal dataset
            if (sum(!not_na_rows) > 0) {
              tmp <- rep(NA, length(not_na_rows))
              tmp[not_na_rows] <- proj
              proj <- tmp
              rm('tmp')
            }
            if (on_0_1000) { proj <- round(proj * 1000) }
            return(proj)
          }
)

#----------------------------------------------------------------------------- #
## 8.9 MAXNET_biomod2_model ---------------------------------------
#----------------------------------------------------------------------------- #
##' @name MAXNET_biomod2_model-class
##' @rdname biomod2_model
##' @export
setClass('MAXNET_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'MAXNET'),
         validity = function(object) { 
           if (!inherits(object@model, "maxnet")) {
             return(FALSE)
           } else { 
             return(TRUE) 
           }
         })

##' 
##' @rdname predict2.bm
##' @importFrom terra as.points rasterize crds cats subset
##' 


setMethod('predict2', signature(object = 'MAXNET_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            args <- list(...)
            filename <- args$filename
            overwrite <- args$overwrite
            on_0_1000 <- args$on_0_1000
            seedval <- args$seedval
            mod.name <- args$mod.name 
            
            if (is.null(overwrite)) { overwrite <- TRUE }
            if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
            
            
            
            newdata.points <- as.points(newdata)
            newdata.df <- as.data.frame(newdata.points)
            categorical_var <- names(newdata)[is.factor(newdata)]
            
            if(length(categorical_var) > 0){
              for(this_var in categorical_var){
                this_levels <- cats(newdata[[categorical_var]])[[1]][,2]
                newdata.df[,this_var] <- factor(newdata.df[,this_var], levels = this_levels)
              }
            }
            set.seed(seedval)
            proj <- predict(object = get_formal_model(object), newdata = newdata.df, clamp = FALSE, type = 'logistic')[, 1]
            
            if (length(get_scaling_model(object)) > 0) {
              proj.to.scale <- data.frame(pred = proj)
              proj <- .run_pred(object = get_scaling_model(object), Prev = 0.5 , dat = proj.to.scale)
            }
            
            if (on_0_1000) { 
              proj <- round(proj * 1000) 
            }
            
            ## convert back to raster file
            proj <- rasterize(crds(newdata.points), 
                              y = newdata, 
                              values = proj,
                              wopt = list(names = mod.name))
            
            # save raster on hard drive ?
            if (!is.null(filename)) {
              cat("\n\t\tWriting projection on hard drive...")
              if (on_0_1000) { ## projections are stored as positive integer
                writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
              } else { ## keep default data format for saved raster
                writeRaster(proj, filename = filename, overwrite = overwrite)
              }
              proj <- rast(filename)
            }
            
            return(proj)
            
          }
)


##' @rdname predict2.bm
setMethod('predict2', signature(object = 'MAXNET_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'logistic')[, 1])
            }
            
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

#----------------------------------------------------------------------------- #
## 8.10 RF_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #
##' @name RF_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('RF_biomod2_model',
         representation(),
         contains = 'biomod2_model',
         prototype = list(model_class = 'RF'),
         validity = function(object) { # check model class
           if (!inherits(object@model, "randomForest")) { return(FALSE)} else { return(TRUE) }
         })

##' 
##' @rdname predict2.bm
##' 


setMethod('predict2', signature(object = 'RF_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            pa <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "PA")
            run <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "run")
            dataset <- paste0("_", pa, "_", run)
            
            if (!is.null(object@model_options@args.values[[dataset]]$type) && object@model_options@args.values[[dataset]]$type == "classification") {
              predfun <- function(object, newdata, mod.name){
                # new predict command used with terra
                subset(predict(newdata, model = get_formal_model(object),
                               type = 'prob',
                               wopt = list(names = rep(mod.name,2))), 
                       2)   
                # old predict function used with raster
                # predict(newdata, model = get_formal_model(object), type = 'prob', index = 2)
              }
            } else { #regression case
              predfun <- function(object, newdata, mod.name){
                predict(newdata, model = get_formal_model(object),
                        type = 'response',
                        wopt = list(names = rep(mod.name,2))) 
              }
            }
            
            # old predict function used with raster
            # predict(newdata, model = get_formal_model(object), type = 'prob', index = 2)            
            # redirect to predict2.biomod2_model.SpatRaster
            
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

##' @rdname predict2.bm
setMethod('predict2', signature(object = 'RF_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            pa <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "PA")
            run <- .extract_modelNamesInfo(object@model_name, obj.type = "mod", info = "run")
            dataset <- paste0("_", pa, "_", run)
            
            if (!is.null(object@model_options@args.values[[dataset]]$type) && object@model_options@args.values[[dataset]]$type == "classification") {
              predfun <- function(object, newdata, not_na_rows) {
                as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'prob')[, '1'])        
              }
            } else { # regression case
              predfun <- function(object, newdata, not_na_rows) {
                as.numeric(predict(get_formal_model(object), as.data.frame(newdata[not_na_rows, , drop = FALSE]), type = 'response'))        
              }
            }
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)


#----------------------------------------------------------------------------- #
## 8.11 SRE_biomod2_model ----------------------------------------------------
#----------------------------------------------------------------------------- #
##' @name SRE_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('SRE_biomod2_model',
         representation(extremal_conditions = 'data.frame'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'SRE'),
         validity = function(object){ return(TRUE) })

##' 
##' @rdname predict2.bm
##' 


setMethod('predict2', signature(object = 'SRE_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            
            predfun <- function(object, newdata, mod.name){
              .sre_projection(new.env = newdata, extrem.cond = object@extremal_conditions, mod.name = mod.name)
            }
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

##' @rdname predict2.bm
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



#----------------------------------------------------------------------------- #
## 8.12 XGBOOST_biomod2_model -----------------------------------------------------
#----------------------------------------------------------------------------- #
##' @name XGBOOST_biomod2_model-class
##' @rdname biomod2_model
##' @export

setClass('XGBOOST_biomod2_model',
         representation(n.trees_optim = 'numeric'),
         contains = 'biomod2_model',
         prototype = list(model_class = 'XGBOOST'),
         validity = function(object) { ## check model class
           if (!inherits(object@model, "xgb.Booster")) { 
             return(FALSE) 
           } else { 
             return(TRUE) 
           }
         })

##' 
##' @rdname predict2.bm

setMethod('predict2', signature(object = 'XGBOOST_biomod2_model', newdata = "SpatRaster"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, mod.name){
              proj <- predict(newdata,
                              model = get_formal_model(object),
                              fun = xgbpred,
                              na.rm = TRUE,
                              wopt = list(names = mod.name))
              proj
            }
            
            # redirect to predict2.biomod2_model.SpatRaster
            callNextMethod(object, newdata, predfun = predfun, ...)
            
          }
)

##' @rdname predict2.bm
setMethod('predict2', signature(object = 'XGBOOST_biomod2_model', newdata = "data.frame"),
          function(object, newdata, ...) {
            predfun <- function(object, newdata, not_na_rows){
              as.numeric(predict(get_formal_model(object), as.matrix(newdata[not_na_rows, , drop = FALSE])))
            }
            # redirect to predict2.biomod2_model.data.frame
            callNextMethod(object, newdata, predfun = predfun, ...)
          }
)

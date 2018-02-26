# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Compilation of deprecated function that will be removed from 
# the package one day or another
# Damien G. - 26/06/2013
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

setGeneric("getModelsPrediction",
           function(obj,...){
             standardGeneric("getModelsPrediction")
           })

setMethod("getModelsPrediction", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_predictions(obj, eval_data=FALSE)")
            if(!as.data.frame){
              if(obj@models.prediction@inMemory ){
                return(obj@models.prediction@val)
              } else{
                if(obj@models.prediction@link != ''){
                  #                   load(obj@models.prediction@link)
                  #                   return(models.prediction)
                  
                  return(get(load(obj@models.prediction@link)))
                } else{ return(NULL) }
              }              
            } else {
              if(obj@models.prediction@inMemory ){
                mod.pred <- as.data.frame(obj@models.prediction@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                 }))
                return(mod.pred)
              } else{
                if(obj@models.prediction@link != ''){
                  #                   load(obj@models.prediction@link)
                  #                   mod.pred <- as.data.frame(models.prediction)
                  mod.pred <- as.data.frame(get(load(obj@models.prediction@link)))                  
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                   function(x){
                                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                  return(mod.pred)
                } else{ return(NULL) }
              }
              
            }
          }
)


setGeneric("getModelsPredictionEval",
           function(obj,...){
             standardGeneric("getModelsPredictionEval")
           })

setMethod("getModelsPredictionEval", "BIOMOD.models.out",
          function(obj, as.data.frame = FALSE){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_predictions(obj, eval_data=TRUE)")
            if(!as.data.frame){
              if(obj@models.prediction.eval@inMemory ){
                return(obj@models.prediction.eval@val)
              } else{
                if(obj@models.prediction.eval@link != ''){
                  models.prediction.eval <- get(load(obj@models.prediction.eval@link))
                  return(models.prediction.eval)
                } else{ return(NULL) }
              }              
            } else {
              if(obj@models.prediction.eval@inMemory ){
                mod.pred <- as.data.frame(obj@models.prediction.eval@val)
                names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                 function(x){
                                                   return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                 }))
                return(mod.pred)
              } else{
                if(obj@models.prediction.eval@link != ''){
                  load(obj@models.prediction.eval@link)
                  mod.pred <- as.data.frame(models.prediction.eval)
                  names(mod.pred) <- unlist(lapply(strsplit(names(mod.pred),".", fixed=TRUE), 
                                                   function(x){
                                                     return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                   }))
                  return(mod.pred)
                } else{ return(NULL) }
              }
              
            }
          }
)


setGeneric("getModelsEvaluations",
           function(obj,...){
             standardGeneric("getModelsEvaluations")
           })

setMethod("getModelsEvaluations", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_evaluations(obj)")
            return(get_evaluations(obj=obj))
          }
)


setGeneric("getModelsVarImport",
           function(obj,...){
             standardGeneric("getModelsVarImport")
           })

setMethod("getModelsVarImport", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_variables_importance(obj)")
            if(obj@variables.importances@inMemory ){
              return(obj@variables.importances@val)
            } else{
              if(obj@variables.importances@link != ''){
                #                 load(obj@variables.importances@link)
                #                 return(variables.importances)
                return(get(load(obj@variables.importances@link)))
              } else{ return(NA) }
            }
          }
)

setGeneric("getModelsOptions",
           function(obj,...){
             standardGeneric("getModelsOptions")
           })

setMethod("getModelsOptions", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_options(obj)")
            if(obj@models.options@inMemory ){
              return(obj@models.options@val)
            } else{
              if(obj@models.options@link != ''){
                #                 load(obj@models.options@link)
                #                 return(models.options)
                return(get(load(obj@models.options@link)))                
              } else{ return(NA) }
            }
          }
)

setGeneric("getModelsInputData",
           function(obj, ...){
             standardGeneric("getModelsInputData")
           })

setMethod("getModelsInputData", "BIOMOD.models.out",
          function(obj, subinfo = NULL){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_predictions(obj)")
            if(is.null(subinfo)){
              if(obj@formated.input.data@inMemory ){
                return(obj@formated.input.data@val)
              } else{
                if(obj@formated.input.data@link != ''){
                  data <- get(load(obj@formated.input.data@link))
                  return(data)
                } else{ cat("\n***"); return(NA) }
              }              
            } else if(subinfo == 'MinMax'){
              return(apply(getModelsInputData(obj)@data.env.var,2, function(x){
                if(is.numeric(x)){
                  return( list(min = min(x,na.rm=T), max = max(x, na.rm=T) ) )
                } else if(is.factor(x)){
                  return(list(levels = levels(x)))
                }
              }) )
            } else if(subinfo == 'expl.var'){
              return(as.data.frame(getModelsInputData(obj)@data.env.var))
            } else if(subinfo == 'expl.var.names'){
              return(obj@expl.var.names)
            } else if(subinfo == 'resp.var'){
              return(as.numeric(getModelsInputData(obj)@data.species))
            } else if(subinfo == 'eval.resp.var'){
              return(as.numeric(getModelsInputData(obj)@eval.data.species))
            } else if(subinfo == 'eval.expl.var'){
              return(as.data.frame(getModelsInputData(obj)@eval.data.env.var))
            } else{
              stop("Unknow subinfo tag")
            }
            
          }
)

setGeneric("getModelsBuiltModels",
           function(obj,...){
             standardGeneric("getModelsBuiltModels")
           })

setMethod("getModelsBuiltModels", "BIOMOD.models.out",
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_built_models(obj)")
            return(obj@models.computed)
          }
)

setGeneric("getProjection",
           function(obj, ...){
             standardGeneric("getProjection")
           })

setMethod("getProjection", "BIOMOD.projection.out",
          function(obj, model = NULL, as.data.frame = FALSE){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_predictions(obj)")
            if(!as.data.frame & is.null(model)){
              if(obj@proj@inMemory ){
                return(obj@proj@val)
              } else {
                if( grepl(".RData", obj@proj@link) ){
                  return(get(load(obj@proj@link)))
                } else if(grepl(".grd", obj@proj@link) | grepl(".img", obj@proj@link)){
                  return(raster::stack(obj@proj@link, RAT=FALSE))
                } else {
                  filesToLoad <- list.files(path=sub("/individual_projections","", obj@proj@link), full.names=T)
                  toMatch <- c('.grd$','.img$')
                  filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)  
                  if(length(filesToLoad)){
                    return(raster::stack(filesToLoad[1], RAT=FALSE))
                  } else {
                    filesToLoad <- list.files(path=obj@proj@link, full.names=T)
                    toMatch <- c('.grd$','.img$')
                    filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
                    toMatch <- obj@models.projected
                    filesToLoad <- grep(pattern=paste(toMatch,collapse="|"), filesToLoad, value=T)
                    proj <- raster::stack(filesToLoad, RAT=FALSE)
                    toMatch <- c(obj@proj@link,".img$",'.grd$', .Platform$file.sep)
                    names(proj) <- gsub(pattern=paste(toMatch,collapse="|"), "", filesToLoad)
                    return(proj)
                  }   
                }
              } 
            } else if(as.data.frame){
              if(obj@proj@inMemory ){
                proj <- as.data.frame(obj@proj@val)
                names(proj) <- unlist(lapply(strsplit(names(proj),".", fixed=TRUE), 
                                             function(x){
                                               return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                             }))
                return(proj)
              } else{
                if(obj@proj@link != ''){
                  load(obj@proj@link)
                  project <- as.data.frame(proj)
                  names(project) <- unlist(lapply(strsplit(names(project),".", fixed=TRUE), 
                                                  function(x){
                                                    return(paste(obj@sp.name, x[3], x[2], x[1],sep="_"))
                                                  }))
                  return(project)
                } else{ return(NA) }
              }
            }
            
          }
)


setGeneric("getEMalgos",
           function(obj,...){
             standardGeneric("getEMalgos")
           })

setMethod("getEMalgos", "BIOMOD.EnsembleModeling.out",
          function(obj, model){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_built_models(obj)")
            return(get_built_models(obj))
            
          }
)

setGeneric("getEM_needed_models",
           function(obj, ...){
             standardGeneric("getEM_needed_models")
           })

setMethod("getEM_needed_models", "BIOMOD.EnsembleModeling.out",
          function(obj, subset='all', ...){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_needed_models(obj)")
            add.args <- list(...)
            needed_models <- lapply(obj@em.models, function(x){
              return(x@model)
            })
            needed_models <- unique(unlist(needed_models))
            return(needed_models)
          }
)


setGeneric("getEMkeptModels",
           function(obj,...){
             standardGeneric("getEMkeptModels")
           })

setMethod("getEMkeptModels", "BIOMOD.EnsembleModeling.out",
          function(obj, model){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_kept_models(obj)")
            if(is.character(model) | is.numeric(model)){
              return(obj@em.res[[model]]$em.models.kept)
            } else{
              return(NULL)
            }
            
          }
)

setGeneric("getEMeval",
           function(obj, ...){
             standardGeneric("getEMeval")
           })

setMethod("getEMeval", "BIOMOD.EnsembleModeling.out",
          function(obj, model=NULL, met=NULL){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_evaluations(obj)")
            return(get_evaluations(obj=obj, model=model, met=met))

            
          }
)

setGeneric("getEMbuiltModels",
           function(obj,...){
             standardGeneric("getEMbuiltModels")
           })

setMethod("getEMbuiltModels", "BIOMOD.EnsembleModeling.out",
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_built_models(obj)")
            return(obj@em.computed)
          })

setGeneric( "getFormalModel", 
            def = function(obj,...){
              standardGeneric( "getFormalModel" )
            } )

setGeneric( "getScalingModel", 
            def = function(obj,...){
              standardGeneric( "getScalingModel" )
            } )

setMethod('getFormalModel', signature('biomod2_model'),
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_formal_model(obj)")
            return(obj@model)
          })

setMethod('getScalingModel', signature('biomod2_model'),
          function(obj){
            cat("\n!! deprecated function that will be remove in next package update")
            cat("\n please prefere to use get_scaling_model(obj)")
            return(obj@scaling_model)
          })



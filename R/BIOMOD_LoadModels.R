BIOMOD_LoadModels <- function(bm.out, ... ){
  
  ####
  # ... can be models, run.eval, data.set to make a models subselection
  add.args <- list(...)
  args <- .BIOMOD_LoadModels.check.args(bm.out, add.args)
  add.args <- args$add.args
  rm(args)
  
  models.to.load <- get_built_models(bm.out)
  
  envir <- parent.frame()
  
  
  if(!is.null(add.args$full.name)){
    models.to.load <- add.args$full.name
  } else{ ## make a subselection
    if(!is.null(add.args$models)){
      model.to.load.tmp <- c()
      for(mod in add.args$models){
        if(sum(grepl(mod,  models.to.load)) > 0){
          model.to.load.tmp <- c(model.to.load.tmp, grep(mod, models.to.load, value=TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }
    
    if(!is.null(add.args$run.eval)){
      model.to.load.tmp <- c()
      for(re in add.args$run.eval){
        if(sum(grepl(re,  models.to.load)) > 0){
          model.to.load.tmp <- c(model.to.load.tmp, grep(re, models.to.load, value=TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }
    
    if(!is.null(add.args$data.set)){
      model.to.load.tmp <- c()
      for(ds in add.args$data.set){
        if(sum(grepl(ds,  models.to.load)) > 0){
          model.to.load.tmp <- c(model.to.load.tmp, grep(ds, models.to.load, value=TRUE))
        }
      }
      models.to.load <-  model.to.load.tmp
    }    
    
  }

  if(length(models.to.load) == 0){
    cat("\n   ! No models computed matched, No models loaded !")
    return(NULL)
  }
  
  if(!is.null(add.args$as) & length(models.to.load)==1){
    assign(x=add.args$as, 
           value=get(load(file=file.path(bm.out@sp.name, "models", bm.out@modeling.id, models.to.load))),
           envir=envir)
    invisible(TRUE)
  } else {
    for(mtl in models.to.load){
      load(file=file.path(bm.out@sp.name, "models", bm.out@modeling.id, mtl),envir=envir)
    }
    
    return(models.to.load)
  }

  
}

.BIOMOD_LoadModels.check.args <- function(bm.out, add.args){
  if(! (inherits(bm.out, 'BIOMOD.models.out') | inherits(bm.out, 'BIOMOD.EnsembleModeling.out'))){
    stop("bm.out arg must be a BIOMOD.models.out or a BIOMOD.EnsembleModeling.out object")
  }
  
  available.args <- c("models", "run.eval", "data.set", "path", "as", "full.name")
  given.args <- names(add.args)
  
  ## is all additional args are known ?
  if(sum(given.args %in% available.args) != length(given.args)){
    cat("\n   !", toString( given.args[which(!(given.args %in% available.args))] ), "arguments are unknown. Please reffer to function help file to see available ones.", fill=.Options$width)
    ## remove unknown args
    for(ga in given.args[which(!(given.args %in% available.args))]){
      add.args[[ga]] <- NULL
    }
  }
  
  ## get all available model names
  avail_models <- get_built_models(bm.out)
  
  ## check additional args values
  ### models names
  if(!is.null(add.args$models)){
    if(sum(add.args$models %in% .extractModelNamesInfo(model.names=avail_models, info='models') ) != length(add.args$models) ){
      stop(paste("models argument must be one of ", toString(.extractModelNamesInfo(model.names=avail_models, info='models')), sep="") )
    } else{
      add.args$models = paste("_", add.args$models, sep="")
    }
  }
  
  ### run.eval names
  if(!is.null(add.args$run.eval)){
    if(sum(add.args$run.eval %in% .extractModelNamesInfo(model.names=avail_models, info='run.eval') != length(add.args$run.eval)) ){
      stop(paste("run.eval argument must be one of ", toString(.extractModelNamesInfo(model.names=avail_models), info='run.eval'), sep="") )
    } else{
      add.args$run.eval = paste("_", add.args$run.eval, "_", sep="")
    }
  }
  
  ### data.set names
  if(!is.null(add.args$data.set)){
    if(sum(add.args$data.set %in% .extractModelNamesInfo(model.names=avail_models, info='data.set') != length(add.args$data.set)) ){
      stop(paste("data.set argument must be one of ", toString(.extractModelNamesInfo(model.names=avail_models), info='data.set'), sep="") )
    } else{
      add.args$data.set = paste("_", add.args$data.set, "_", sep="")
    }
  }
  
  ### path to sim
  if(!is.null(add.args$path)){
    if(!(bm.out@sp.name %in% list.dirs(path = add.args$path))){
      stop("invalid path given")
    }
  } else{
    add.args$path = "."
  }
  
  ### full.name
  if(!is.null(add.args$full.name)){
    if(sum(!(add.args$full.name %in% avail_models)) > 0){
      stop("full.name arg must be one of : ", toString(avail_models))
    }
  }
  
  
  return(list(add.args = add.args))
  
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
                
'.extractModelNamesInfo' <- function(model.names, info = 'species'){
  if(!is.character(model.names)){
    stop("model.names must be a character vector")
  }
  if(!is.character(info) | length(info) != 1 | !(info %in% c('species', 'data.set', 'models', 'run.eval')) ){
    stop("info must be 'species', 'data.set', 'models' or 'run.eval'")
  }
                
  info.tmp <- as.data.frame(strsplit(model.names, "_"))
  
  return( switch(info,
                 species = paste(unique(unlist(info.tmp[-c(nrow(info.tmp), nrow(info.tmp)-1, nrow(info.tmp)-2),])), collapse="_"),
                 data.set = paste(unique(unlist(info.tmp[(nrow(info.tmp)-2),]))),
                 run.eval = paste(unique(unlist(info.tmp[(nrow(info.tmp)-1),]))),
                 models = paste(unique(unlist(info.tmp[(nrow(info.tmp)),]))) ) )
              
}
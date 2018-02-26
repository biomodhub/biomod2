'BIOMOD_EnsembleForecasting' <- function( EM.output,
                                          projection.output = NULL,
                                          new.env = NULL,
                                          xy.new.env = NULL,
                                          selected.models = 'all',
                                          proj.name = NULL,
                                          binary.meth = NULL,
                                          filtered.meth = NULL,
                                          compress = TRUE,
                                          ...){
  .bmCat("Do Ensemble Models Projections")
  
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  args <- list(...)
  
  args_checked <- .BIOMOD_EnsembleForecasting.check.args( EM.output,
                                                          projection.output,
                                                          new.env,
                                                          selected.models,
                                                          proj.name,
                                                          total.consensus=FALSE,
                                                          binary.meth,
                                                          filtered.meth )
  
#   total.consensus <- args_checked$total.consensus
  proj.name <- args_checked$proj.name
  selected.models <- args_checked$selected.models
    
  output.format <- args$output.format # raster output format
  compress <- args$compress # compress or not output
  do.stack <- args$do.stack # save raster as stack or layers
  keep.in.memory <- args$keep.in.memory # store results on memory or only on hard drive
  on_0_1000 <- args$on_0_1000 # convert 0-1 predictions on 0-1000 scale to limit memory consuming


  
  
  if(is.null(output.format)){
    if(length(projection.output)){
      if(projection.output@type != 'RasterStack')
        output.format <- ".RData"
      else
        output.format <- ".grd"
    } else{
      if(! inherits(new.env,'Raster'))
        output.format <- ".RData"
      else
        output.format <- ".grd"
    }

  }
  
  if(is.null(compress)) compress <- FALSE
  
  if(is.null(do.stack)) {
    do.stack <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than projection.output ones
    if(!is.null(projection.output)){
      if(all(grepl("individual_projections", projection.output@proj@link))){
        do.stack <- FALSE
      }
    }
  }
  
  if(is.null(keep.in.memory)){
    keep.in.memory <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than projection.output ones
    if(!is.null(projection.output)){
      keep.in.memory <- projection.output@proj@inMemory
    }
  } 
  
  if(is.null(xy.new.env)) { 
    if(!is.null(projection.output)){
      xy.new.env <- projection.output@xy.coord
    } else{
      xy.new.env <- matrix()
    }
  } 

  if(is.null(on_0_1000)) on_0_1000 = TRUE # by default outputs are return on a 0 - 1000 scale 
  
  rm(list=c('args_checked','args'))
  
  # get argument from projection
  proj_out <- new('BIOMOD.projection.out',
                  proj.names = proj.name,
                  sp.name =  EM.output@sp.name,
                  expl.var.names = EM.output@expl.var.names,
                  models.projected = selected.models,
#                   scaled.models = modeling.output@rescal.all.models,
                  xy.coord = xy.new.env,
                  modeling.object.id = EM.output@modeling.id)
  
  proj_out@modeling.object@link = EM.output@link
  
  proj_is_raster <- FALSE
  if(inherits(new.env, 'Raster')){
    proj_is_raster <- TRUE
  } else if( length(projection.output) ){
    if(inherits(projection.output@proj, 'BIOMOD.stored.raster.stack')){
      proj_is_raster <- TRUE
    }
  }
  
  
  if(proj_is_raster){
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj <- new('BIOMOD.stored.array')
    do.stack = TRUE
  }
  
  # 1.c creating output directory
  dir.create(file.path(EM.output@sp.name,paste("proj_", proj.name, sep="")), 
             showWarnings = FALSE, recursive = TRUE, mode = "777")
  
  indiv_proj_dir <- file.path(EM.output@sp.name,paste("proj_", proj.name, sep=""), "individual_projections")
#   if(!do.stack){
#     rasterOptions(todisk=T)
    dir.create(indiv_proj_dir, 
               showWarnings = FALSE, recursive = TRUE, mode = "777")
#   }
  
  saved.files <- c()
                                                  
  
  ### get needed models prediction ###
  needed_predictions <- get_needed_models(EM.output, selected.models=selected.models)

  if (length(projection.output)){
    formal_pred <- get_predictions(projection.output, full.name=needed_predictions, as.data.frame=ifelse(projection.output@type=='array',T,F) )
  } else{
    # make prediction according to given environment
    tmp_dir <- paste('Tmp', format(Sys.time(), "%s"), sep="")
    formal_pred <- BIOMOD_Projection( modeling.output = load_stored_object(EM.output@models.out.obj),
                                      new.env = new.env,
                                      proj.name = tmp_dir,
                                      xy.new.env = NULL,
                                      selected.models = needed_predictions,
                                      compress = TRUE,
                                      build.clamping.mask = F,
                                      do.stack=T, silent = T, on_0_1000 = on_0_1000 )
    # getting the results
    formal_pred <- get_predictions(formal_pred, full.name=needed_predictions, as.data.frame=ifelse(inherits(new.env,'Raster'),F,T))
    
    # remove tmp directory
    unlink(file.path(EM.output@sp.name,paste("proj_",tmp_dir,sep="")),recursive = TRUE, force = TRUE)
  }
  
  ef.out <- NULL
  # 2. Do the ensemble modeling
  for( em.comp in EM.output@em.computed[which(EM.output@em.computed %in% selected.models)]){
    cat("\n\n\t> Projecting", em.comp, "...")
    model.tmp <- NULL
    file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp,output.format))
    BIOMOD_LoadModels(EM.output, full.name=em.comp, as='model.tmp')
    if(inherits(formal_pred,'Raster')){
      ef.tmp <- predict(model.tmp, 
                        formal_predictions = raster::subset(formal_pred, subset=model.tmp@model, drop=FALSE),
                        on_0_1000 = on_0_1000, 
                        filename = ifelse(output.format == '.RData', '', file_name_tmp))
    } else {
      ef.tmp <- predict(model.tmp, formal_predictions = formal_pred[,model.tmp@model, drop=FALSE], on_0_1000 = on_0_1000)
    }
        
    if(inherits(ef.tmp,'Raster')){
      if(do.stack){
        if(length(ef.out)) ef.out <- stack(ef.out,ef.tmp) else ef.out <- raster::stack(ef.tmp)
      } else {
        file_name_tmp <- file.path(indiv_proj_dir,paste(em.comp,output.format,sep=""))
        if(output.format== '.RData'){
          save(ef.tmp, file=file_name_tmp, compress=compress)
        } 
        saved.files <- c(saved.files, file_name_tmp)
      }
    } else{
      ef.out <- cbind(ef.out,ef.tmp)
    } 
  }
  
  proj_out@models.projected <- EM.output@em.computed[which(EM.output@em.computed %in% selected.models)]
  
  if(do.stack){
    if( inherits(ef.out, "Raster") ) {
      names(ef.out) <- proj_out@models.projected
    } else {
      colnames(ef.out) <- proj_out@models.projected
    }
    # save object
    file_name_tmp <- file.path(EM.output@sp.name,paste("proj_", proj.name, sep=""),paste("proj_", proj.name,"_",EM.output@sp.name,"_ensemble",output.format,sep=""))
    if(output.format== '.RData'){
      save(ef.out, file=file_name_tmp, compress=compress)
    } else if( inherits(ef.out, "Raster") ){
      ## TODO : define the raster dataformat (depends if em.cv has been computed)
      writeRaster(ef.out,filename=file_name_tmp, overwrite=TRUE)
    }
    saved.files <- c(saved.files, file_name_tmp)
    proj_out@proj@link <- file_name_tmp
  } else {
    proj_out@proj@link <- saved.files #EM.output@em.computed
  }
  
  if(!is.null(ef.out)){
    proj_out@proj@val <- ef.out
    proj_out@proj@inMemory <- TRUE
  }
  
  
  
  ## binary && filtering transformations #############################################################
  if(length(binary.meth) | length(filtered.meth)){
    cat("\n")
    eval.meth <- unique(c(binary.meth,filtered.meth))
    
    ## get all treshold
    thresholds <- sapply(selected.models, function(x){
      get_evaluations(EM.output)[[x]][eval.meth, "Cutoff"]  
    })
    
    ## convert thresholds deopending on the chosen scale
    if(! on_0_1000 ) { thresholds <- thresholds / 1000 }
    
    ## do binary transformation
    for(eval.meth in binary.meth){
      cat("\n\t> Building", eval.meth,"binaries")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          proj_bin <- BinaryTransformation(raster(file.tmp, RAT=FALSE),thres.tmp)
          writeRaster(x = proj_bin,
                      filename = sub(output.format, paste("_",eval.meth,"bin", output.format, sep=""), file.tmp), 
                      overwrite=TRUE ) ## , datatype = "INT2S",NAflag=-9999)
        }
      } else {
        assign(x = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep=""),
               value = BinaryTransformation(ef.out,thresholds))
        
        if(output.format == '.RData'){
          save(list = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep=""), 
               file = file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", output.format ,sep="")), compress=compress)   
        } else {
          writeRaster(x=get(paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep="")),
                      filename=file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", output.format ,sep="")), 
                      overwrite=TRUE) ## , datatype = "INT2S", NAflag=-9999)
                      
        }
        
        rm(list=paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"bin", sep=""))
      }
    }
    
    ## do filtered transformation
    for(eval.meth in filtered.meth){
      cat("\n\t> Building", eval.meth,"filtered")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- thresholds[i]
          ## TODO : define the raster dataformat (depends if em.cv has been computed)
          filt_proj <- FilteringTransformation(raster(file.tmp, RAT=FALSE),thres.tmp)
          writeRaster(x = filt_proj,
                      filename = sub(output.format, paste("_",eval.meth,"filt", output.format, sep=""), file.tmp), 
                      overwrite=TRUE)
        }
      } else {
        assign(x = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep=""),
               value = FilteringTransformation(ef.out,thresholds))
        
        if(output.format == '.RData'){
          save(list = paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep=""), 
               file = file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", output.format ,sep="")), compress=compress)   
        } else {
          ## TODO : define the raster dataformat (depends if em.cv has been computed)
          writeRaster(x=get(paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep="")),
                      filename=file.path(EM.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", output.format ,sep="")), overwrite=TRUE)
        }
        
        rm(list=paste("proj_",proj.name, "_", EM.output@sp.name,"_ensemble_",eval.meth,"filt", sep=""))
      }
    }
  }
  
  
  # save object
  if(!keep.in.memory){
    proj_out <- free(proj_out)
  }
  
  assign(paste(EM.output@sp.name,".", proj.name, ".ensemble.projection.out", sep=""), proj_out)
  save(list = paste(EM.output@sp.name,".", proj.name, ".ensemble.projection.out", sep=""),
       file = file.path(EM.output@sp.name, paste("proj_", proj.name, sep=""), paste(EM.output@sp.name,".", proj.name, ".ensemble.projection.out", sep="")))
  
  cat("\n")
  .bmCat("Done")
  return(proj_out)
}







# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
.BIOMOD_EnsembleForecasting.check.args <- function( EM.output,
                                                    projection.output,
                                                    new.env,
                                                    selected.models,
                                                    proj.name,
                                                    total.consensus,
                                                    binary.meth,
                                                    filtered.meth ){

  # check needed data are referred
  if( (is.null(projection.output) & is.null(new.env)) | (!is.null(projection.output) & !is.null(new.env))){
    stop("You have to refer at one of 'projection.output' or 'new.env' argument")
  }
  
  if(!inherits(EM.output,"BIOMOD.EnsembleModeling.out")){
    stop("EM.output must be a 'BIOMOD.EnsembleModeling.out' object")
  }
  
  # check all needed predictions are available
  needed_pred <- get_needed_models(EM.output, selected.models=selected.models)  
  
  if(!is.null(projection.output)){
    if(!inherits(projection.output, "BIOMOD.projection.out")){
      stop("projection.output must be a 'BIOMOD.projection.out' object")
    }
    missing_pred <- needed_pred[! (needed_pred %in% projection.output@models.projected)]
    if( length(missing_pred) ){
      stop("Some models prediction missing :", toString(missing_pred))
    }
  }
  
  ## selected.models
  if(selected.models[1] == 'all'){
    selected.models <- get_built_models(EM.output)
  } else{
    selected.models <- intersect(selected.models, get_built_models(EM.output))
  }
  if(length(selected.models) < 1){
    stop('No models selected')
  }
  
  ## projection name
  if(!length(proj.name) & !length(projection.output)){
    stop("You have to give a valid 'proj.name' if you don't work with projection.output")
  } else if(!length(proj.name)){
    proj.name <- projection.output@proj.names
  }
  
  
                              
  if(total.consensus){
    if(length(EM.output@em.computed) < 2){
      cat("\n      ! Total consensus projection was switched off because only one Ensemble modeling was done")
      total.consensus <- FALSE
    }
  }
    
  if(!is.null(binary.meth)){
#     if(sum(!(binary.meth %in% EM.output@eval.metric))){
#       stop(paste("binary methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
#                  toString(EM.output@eval.metric)," )", sep=""))
#     }
  }
  
  if(!is.null(filtered.meth)){
#     if(sum(!(filtered.meth %in% EM.output@eval.metric))){
#       stop(paste("filtering methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
#                  toString(EM.output@eval.metric)," )", sep=""))
#     }
  }
                                                    
  return(list(projection.output = projection.output,
              EM.output = EM.output,
              selected.models = selected.models,
              proj.name = proj.name,
              total.consensus = total.consensus,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth))
  
}
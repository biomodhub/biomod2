
BIOMOD_EnsembleForecasting <- function( EM.output,
                                        projection.output = NULL,
                                        new.env = NULL,
                                        xy.new.env = NULL,
                                        selected.models = 'all',
                                        proj.name = NULL,
                                        binary.meth = NULL,
                                        filtered.meth = NULL,
                                        compress = TRUE,
                                        ...)
{
  .bmCat("Do Ensemble Models Projections")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- list(...)
  output.format <- args$output.format # raster output format
  compress <- args$compress # compress or not output
  do.stack <- args$do.stack # save raster as stack or layers
  keep.in.memory <- args$keep.in.memory # store results on memory or only on hard drive
  on_0_1000 <- args$on_0_1000 # convert 0-1 predictions on 0-1000 scale to limit memory consuming
  
  if (is.null(compress)) { compress <- FALSE }
  if (is.null(on_0_1000)) { on_0_1000 <- TRUE } # by default outputs are return on a 0 - 1000 scale 
  
  args <- .BIOMOD_EnsembleForecasting.check.args(EM.output, projection.output, new.env,
                                                 selected.models, proj.name, binary.meth, filtered.meth)
  projection.output <- args$projection.output
  EM.output <- args$EM.output
  selected.models <- args$selected.models
  proj.name <- args$proj.name
  # total.consensus <- args$total.consensus
  binary.meth <- args$binary.meth
  filtered.meth <- args$filtered.meth
  
  if (is.null(output.format)) {
    if (length(projection.output)) {
      output.format = ifelse(projection.output@type != 'RasterStack', ".RData", ".grd")
    } else {
      output.format = ifelse(!inherits(new.env, 'Raster'), ".RData", ".grd")
    }
  }
  
  if (is.null(do.stack)) {
    do.stack <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than projection.output ones
    if (!is.null(projection.output) &&
        all(grepl("individual_projections", projection.output@proj@link))) {
      do.stack <- FALSE
    }
  }
  
  if (is.null(keep.in.memory)) {
    keep.in.memory <- TRUE # if no info at all set it TRUE
    # if not explicitly defined apply same rules than projection.output ones
    if (!is.null(projection.output)) {
      keep.in.memory <- projection.output@proj@inMemory
    }
  } 
  
  if (is.null(xy.new.env)) {
    if (!is.null(projection.output)) {
      xy.new.env <- projection.output@xy.coord
    } else {
      xy.new.env <- matrix()
    }
  } 
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  proj_out <- new('BIOMOD.projection.out',
                  proj.names = proj.name,
                  sp.name =  EM.output@sp.name,
                  expl.var.names = EM.output@expl.var.names,
                  models.projected = selected.models,
                  xy.coord = xy.new.env,
                  modeling.object.id = EM.output@modeling.id)
  proj_out@modeling.object@link = EM.output@link
  
  proj_is_raster <- FALSE
  if (inherits(new.env, 'Raster') || (length(projection.output) && inherits(projection.output@proj, 'BIOMOD.stored.raster.stack'))) {
    proj_is_raster <- TRUE
  }
  if (proj_is_raster) {
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj <- new('BIOMOD.stored.array')
    do.stack = TRUE
  }
  
  ## 2. Create simulation directories -------------------------------------------------------------
  nameProj <- paste0("proj_", proj.name)
  nameProjSp <- paste0(nameProj, "_", EM.output@sp.name, "_ensemble")
  nameEval <- paste0(nameProjSp, "_", eval.meth)
  namePath <- file.path(EM.output@sp.name, nameProj)
  .BIOMOD_EnsembleForecasting.prepare.workdir(sp.name = EM.output@sp.name, proj.folder = nameProj)
  
  
  ## 3. Get needed projections --------------------------------------------------------------------
  needed_predictions <- get_needed_models(EM.output, selected.models=selected.models)
  if (length(projection.output)) {
    formal_pred <- get_predictions(projection.output,
                                   full.name = needed_predictions,
                                   as.data.frame = ifelse(projection.output@type == 'array', TRUE, FALSE))
  } else {
    # make prediction according to given environment
    tmp_dir <- paste0('Tmp', format(Sys.time(), "%s"))
    formal_pred <- BIOMOD_Projection(modeling.output = load_stored_object(EM.output@models.out.obj),
                                     new.env = new.env,
                                     proj.name = tmp_dir,
                                     xy.new.env = NULL,
                                     selected.models = needed_predictions,
                                     compress = TRUE,
                                     build.clamping.mask = FALSE,
                                     do.stack = TRUE,
                                     silent = TRUE,
                                     on_0_1000 = on_0_1000)
    # getting the results
    formal_pred <- get_predictions(formal_pred,
                                   full.name = needed_predictions,
                                   as.data.frame = ifelse(inherits(new.env, 'Raster'), FALSE, TRUE))
    # remove tmp directory
    unlink(file.path(EM.output@sp.name, paste0("proj_", tmp_dir)), recursive = TRUE, force = TRUE)
  }
  
  ## 4. MAKING PROJECTIONS ------------------------------------------------------------------------
  ef.out <- NULL
  saved.files <- proj_names <- vector()
  for (em.comp in EM.output@em.computed[which(EM.output@em.computed %in% selected.models)]) {
    cat("\n\n\t> Projecting", em.comp, "...")
    file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp,output.format))
    
    model.tmp <- NULL
    BIOMOD_LoadModels(EM.output, full.name = em.comp, as = 'model.tmp')
    if (inherits(formal_pred, 'Raster')) {
      ef.tmp <- predict(model.tmp,
                        formal_predictions = raster::subset(formal_pred, subset = model.tmp@model, drop = FALSE),
                        on_0_1000 = on_0_1000,
                        filename = ifelse(output.format == '.RData', '', file_name_tmp))
    } else {
      ef.tmp <- predict(model.tmp,
                        formal_predictions = formal_pred[, model.tmp@model, drop = FALSE],
                        on_0_1000 = on_0_1000)
    }
    
    if (!is.null(ef.tmp)) {
      proj_names <- c(proj_names, em.comp)
      if (inherits(ef.tmp, 'Raster')) {
        if (do.stack) {
          if (length(ef.out)) {
            ef.out <- stack(ef.out, ef.tmp)
          } else {
            ef.out <- raster::stack(ef.tmp)
          }
        } else {
          file_name_tmp <- file.path(indiv_proj_dir, paste0(em.comp, output.format))
          if (output.format == '.RData') {
            save(ef.tmp, file = file_name_tmp, compress = compress)
          }
          saved.files <- c(saved.files, file_name_tmp)
        }
      } else {
        ef.out <- cbind(ef.out, ef.tmp)
      }
    }
  }
  
  # proj_out@models.projected <- EM.output@em.computed[which(EM.output@em.computed %in% selected.models)]
  proj_out@models.projected <- proj_names
  
  if(do.stack) {
    if (inherits(ef.out, "Raster")) {
      names(ef.out) <- proj_out@models.projected
    } else {
      colnames(ef.out) <- proj_out@models.projected
    }
    # save object
    file_name_tmp <- file.path(namePath, paste0(nameProjSp, output.format))
    if (output.format == '.RData') {
      save(ef.out, file = file_name_tmp, compress = compress)
    } else if (inherits(ef.out, "Raster")) {
      ## TODO : define the raster dataformat (depends if em.cv has been computed)
      writeRaster(ef.out, filename = file_name_tmp, overwrite = TRUE)
    }
    saved.files <- c(saved.files, file_name_tmp)
  } 
  proj_out@proj@link <- saved.files #EM.output@em.computed
  
  if (!is.null(ef.out)) {
    proj_out@proj@val <- ef.out
    proj_out@proj@inMemory <- TRUE
  }
  
  
  ## 5. Compute binary and/or filtered transformation ---------------------------------------------
  if (length(binary.meth) | length(filtered.meth))
  {
    cat("\n")
    eval.meth <- unique(c(binary.meth, filtered.meth))
    
    ## Get all evaluation thresholds
    thresholds <- sapply(selected.models, function(x) {
      get_evaluations(EM.output)[[x]][eval.meth, "Cutoff"]
    })
    if (!on_0_1000) { thresholds <- thresholds / 1000 }
    
    ## TODO : define the raster dataformat (depends if em.cv has been computed)
        ## Do binary transformation
    for (eval.meth in binary.meth) {
      cat("\n\t> Building", eval.meth, "binaries")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          writeRaster(x = bm_BinaryTransformation(raster(file.tmp, RAT = FALSE), thresholds[i]),
                      filename = sub(output.format,
                                     paste0("_", eval.meth, "bin", output.format),
                                     file.tmp),
                      overwrite = TRUE)
        }
      } else {
        nameBin = paste0(nameEval, "bin")
        assign(x = nameBin, value = bm_BinaryTransformation(ef.out, thresholds))
        
        if (output.format == '.RData') {
          save(list = nameBin,
               file = file.path(namePath, paste0(nameBin, output.format)),
               compress = compress)
        } else {
          writeRaster(x = get(nameBin),
                      filename = file.path(namePath, paste0(nameBin, output.format)),
                      overwrite = TRUE,)
        }
      }
    }
    
    ## Do filtered transformation
    for (eval.meth in filtered.meth) {
      cat("\n\t> Building", eval.meth, "filtered")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          writeRaster(x = bm_BinaryTransformation(raster(file.tmp, RAT = FALSE), thresholds[i], doFiltering = TRUE),
                      filename = sub(output.format,
                                     paste0("_", eval.meth, "filt", output.format),
                                     file.tmp),
                      overwrite = TRUE)
        }
      } else {
        nameFilt = paste0(nameEval, "filt")
        assign(x = nameFilt, value = bm_BinaryTransformation(ef.out, thresholds, doFiltering = TRUE))
        if (output.format == '.RData') {
          save(list = nameFilt,
               file = file.path(namePath, paste0(nameFilt, output.format)),
               compress = compress)
        } else {
          writeRaster(x = get(nameFilt),
                      filename = file.path(namePath, paste0(nameFilt, output.format)),
                      overwrite = TRUE)
        }
      }
    }
  }
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  ## save a copy of output object without value to be lighter
  nameOut <- paste0(EM.output@sp.name, ".", proj.name, ".ensemble.projection.out")
  if (!keep.in.memory) { proj_out <- free(proj_out) }
  assign(nameOut, proj_out)
  save(list = nameOut, file = file.path(namePath, nameOut))
  
  .bmCat("Done")
  return(proj_out)
}



###################################################################################################

.BIOMOD_EnsembleForecasting.prepare.workdir <- function(sp.name, proj.folder)
{
  cat("\nCreating suitable Workdir...\n")
  dir.create(file.path(sp.name, proj.folder), showWarnings = FALSE, recursive = TRUE, mode = "777")
  indiv_proj_dir <- file.path(sp.name, proj.folder, "individual_projections")
  dir.create( indiv_proj_dir, showWarnings = FALSE, recursive = TRUE, mode = "777")
}

###################################################################################################

.BIOMOD_EnsembleForecasting.check.args <- function(EM.output,
                                                   projection.output,
                                                   new.env,
                                                   selected.models,
                                                   proj.name,
                                                   # total.consensus,
                                                   binary.meth,
                                                   filtered.meth)
{
  ## 1. Check EM.output -------------------------------------------------------
  .fun_testIfInherits(TRUE, "EM.output", EM.output, "BIOMOD.EnsembleModeling.out")
  
  ## 2. Check needed data and predictions -------------------------------------
  if ((is.null(projection.output) && is.null(new.env)) ||
      (!is.null(projection.output) & !is.null(new.env))) {
    stop("You have to refer at one of 'projection.output' or 'new.env' argument")
  }
  
  if (!is.null(projection.output)) {
    .fun_testIfInherits(TRUE, "projection.output", projection.output, "BIOMOD.projection.out")
    
    ## check all needed predictions are available
    needed_pred <- get_needed_models(EM.output, selected.models = selected.models)  
    missing_pred <- needed_pred[!(needed_pred %in% projection.output@models.projected)]
    if (length(missing_pred)) {
      stop("Some models predictions missing :", toString(missing_pred))
    }
  }
  
  ## 3. Check selected.models -------------------------------------------------
  if (selected.models[1] == 'all') {
    selected.models <- get_built_models(EM.output)
  } else {
    selected.models <- intersect(selected.models, get_built_models(EM.output))
  }
  if (length(selected.models) < 1) {
    stop('No models selected')
  }
  
  ## 4. Check proj.name -------------------------------------------------------
  if (!length(proj.name) && !length(projection.output)) {
    stop("You have to give a valid 'proj.name' if you don't work with projection.output")
  } else if (!length(proj.name)) {
    proj.name <- projection.output@proj.names
  }
  
  ## 5. Check total.consensus -------------------------------------------------
  # if (total.consensus && length(EM.output@em.computed) < 2) {
  #   cat("\n      ! Total consensus projection was switched off because only one Ensemble modeling was done")
  #   total.consensus <- FALSE
  # }

  ## 6. Check binary.meth and filtered.meth -----------------------------------
  if (!is.null(binary.meth)) {
    #     if(sum(!(binary.meth %in% EM.output@eval.metric))){
    #       stop(paste0("binary methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
    #                  toString(EM.output@eval.metric)," )"))
    #     }
  }
  
  if (!is.null(filtered.meth)) {
    #     if(sum(!(filtered.meth %in% EM.output@eval.metric))){
    #       stop(paste0("filtering methods must be compatible with Ensemble Modeling evaluation metrics (e.g. ",
    #                  toString(EM.output@eval.metric)," )"))
    #     }
  }
  
  
  return(list(projection.output = projection.output,
              EM.output = EM.output,
              selected.models = selected.models,
              proj.name = proj.name,
              # total.consensus = total.consensus,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth))
}

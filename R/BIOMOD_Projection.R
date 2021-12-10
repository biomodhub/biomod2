####################################################################################################
# BIOMOD_Projection
# Damien.G
# feb 2012
####################################################################################################

# AIM :
#   Project models from BIOMOD_Modeling with different explanatory variables

# INPUT :


# OUTPUT :



# NOTE :
#   It would be nice to add done projections to input Biomod.models.object
#   .BIOMOD_Projection.check.args <- may be reorder variables if necessary

####################################################################################################

BIOMOD_Projection <- function(modeling.output,
                              new.env,
                              proj.name,
                              xy.new.env = NULL,
                              selected.models = 'all',
                              binary.meth = NULL,
                              filtered.meth = NULL,
                              compress = TRUE,
                              build.clamping.mask = TRUE,
                              ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- list(...)
  omit.na <- args$omit.na # omit all non full filled environmental cells (always TRUE if env is a raster object)
  silent <- args$silent # echo advancement or not
  do.stack <- args$do.stack # store output in a lone stack
  clamping.level <- args$clamping.levels # remove all cells where at least clamping.level variables are out of their calibrating range
  output.format <- args$output.format # raster output format
  keep.in.memory <- args$keep.in.memory # store results on memory or only on hard drive
  on_0_1000 <- args$on_0_1000 # transform projections on a 0 - 1000 scale to limit memory consumption
  split.proj <- args$split.proj
  
  if (is.null(omit.na)) { omit.na <- TRUE }
  if (is.null(silent)) { silent <- FALSE }
  if (is.null(do.stack)) { do.stack <- TRUE }
  if (is.null(keep.in.memory)) { keep.in.memory <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- TRUE } # by default we return projections on a 0 -  1000 scale.
  if (is.null(split.proj)) { split.proj <- 1 } # by default we return projections on a 0 -  1000 scale.
  if (!silent) { .bmCat("Do Models Projections") }
  
  args <- .BIOMOD_Projection.check.args(modeling.output, new.env, proj.name, xy.new.env, selected.models,
                                        binary.meth, filtered.meth, compress, do.stack, output.format)
  
  proj.name <- args$proj.name
  selected.models <- args$selected.models
  binary.meth <- args$binary.meth
  filtered.meth <- args$filtered.meth
  compress <- args$compress
  do.stack <- args$do.stack
  xy.new.env <- args$xy.new.env
  output.format <- args$output.format
  rm(args)
  
  ## 1. Create output object ----------------------------------------------------------------------
  proj_out <- new('BIOMOD.projection.out',
                  proj.names = proj.name,
                  sp.name =  modeling.output@sp.name,
                  expl.var.names = modeling.output@expl.var.names,
                  models.projected = selected.models,
                  scaled.models = modeling.output@rescal.all.models,
                  xy.coord = xy.new.env,
                  modeling.object.id = modeling.output@modeling.id)
  proj_out@modeling.object@link = modeling.output@link
  if (inherits(new.env, 'Raster')) {
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else {
    proj_out@proj <- new('BIOMOD.stored.array')
  }
  
  ## 2. Create simulation directories -------------------------------------------------------------
  nameProj <- paste0("proj_", proj.name)
  nameProjSp <- paste0(nameProj, "_", modeling.output@sp.name)
  namePath <- file.path(modeling.output@sp.name, nameProj)
  dir.create(namePath, showWarnings = FALSE, recursive = TRUE, mode = "777")
  if (!do.stack) {
    dir.create(file.path(namePath, "individual_projections"),
               showWarnings = FALSE, recursive = TRUE, mode = "777")
  }
  
  ## 3. Define the clamping mask ------------------------------------------------------------------
  if (build.clamping.mask) {
    if (!silent) { cat("\n\t> Building clamping mask\n") }
    nameMask <- paste0(nameProjSp, "_ClampingMask")
    MinMax <- get_formal_data(modeling.output, 'MinMax')
    assign(x = nameMask, value = .build.clamping.mask(new.env, MinMax))
    
    if (output.format == '.RData') {
      save(list = nameMask,
           file = file.path(namePath, paste0(nameProj, "_ClampingMask", output.format)),
           compress = compress)
    } else {
      writeRaster(x = get(nameMask),
                  filename = file.path(namePath, paste0(nameProj, "_ClampingMask", output.format)),
                  datatype = "INT2S",
                  NAflag = -9999,
                  overwrite = TRUE)
    }
  }
  
  ## 4. MAKING PROJECTIONS ------------------------------------------------------------------------
  if (!do.stack) {
    proj <- sapply(selected.models, function(mod.name)
    {
      cat("\n\t> Projecting", mod.name, "...")
      filename <- file.path(namePath, "individual_projections"
                            , paste0(nameProj, "_", mod.name, ifelse(output.format == ".RData", ".grd", output.format)))
      return(filename)
    })
    
  } else {
    proj <- lapply(selected.models, function(mod.name)
    {
      cat("\n\t> Projecting", mod.name, "...")
      BIOMOD_LoadModels(modeling.output, full.name = mod.name, as = "mod")
      pred.tmp <- predict(mod, new.env, on_0_1000 = on_0_1000, filename = filename, omit.na = omit.na, split.proj = split.proj)
      return(pred.tmp)
    })
    ## Putting predictions into the right format
    if (inherits(new.env, "Raster")) {
      proj <- stack(proj)
      names(proj) <- selected.models
    } else {
      proj <- as.data.frame(proj)
      names(proj) <- selected.models
      proj <- DF_to_ARRAY(proj)
    }
    if (keep.in.memory) {
      proj_out@proj@val <- proj
      proj_out@proj@inMemory <- TRUE
    }
  }

  ## save projections
  assign(x = nameProjSp, value = proj)
  
  saved.files <- file.path(namePath, paste0(nameProjSp, output.format))
  if (output.format == '.RData') {
    save(list = nameProjSp, file = saved.files, compress = compress)
  } else if(do.stack) {
    writeRaster(x = get(nameProjSp), filename = saved.files,
                overwrite = TRUE, datatype = "INT2S", NAflag = -9999)
  } else {
    saved.files = unlist(proj)
  }
  proj_out@type <- class(proj_out@proj@val)
  proj_out@proj@link <- saved.files
  
  
  ## 5. Compute binary and/or filtered transformation ---------------------------------------------
  if (!is.null(binary.meth) | !is.null(filtered.meth))
  {
    cat("\n")
    eval.meth <- unique(c(binary.meth, filtered.meth))
    
    ## Get all evaluation thresholds
    if (inherits(new.env, "Raster")) {
      thresholds <- matrix(0, nrow = length(eval.meth), ncol = length(selected.models), dimnames = list(eval.meth, selected.models))
      for (mod in selected.models) {
        PA.run <- .extractModelNamesInfo(model.names = mod, info = 'data.set')
        eval.run <- .extractModelNamesInfo(model.names = mod, info = 'run.eval')
        algo.run <- .extractModelNamesInfo(model.names = mod, info = 'models')
        thresholds[eval.meth, mod] <- get_evaluations(modeling.output)[eval.meth, "Cutoff", algo.run, eval.run, PA.run]
        if (!on_0_1000) { thresholds[eval.meth, mod]  <- thresholds[eval.meth, mod] / 1000 }
      }
    } else {
      thresholds <- array(0, dim = c(length(eval.meth), dim(proj)[-1]), dimnames = c(list(eval.meth), dimnames(proj)[-1]))
      for (mod in selected.models) {
        PA.run <- .extractModelNamesInfo(model.names = mod, info = 'data.set')
        eval.run <- .extractModelNamesInfo(model.names = mod, info = 'run.eval')
        algo.run <- .extractModelNamesInfo(model.names = mod, info = 'models')
        thresholds[eval.meth, algo.run, eval.run, PA.run] <- get_evaluations(modeling.output)[eval.meth, "Cutoff", algo.run, eval.run, PA.run]
        if (!on_0_1000) {
          thresholds[eval.meth, algo.run, eval.run, PA.run]  <- thresholds[eval.meth, algo.run, eval.run, PA.run] / 1000
        }
      }
    }
    
    ## Do binary transformation
    for (eval.meth in binary.meth) {
      cat("\n\t> Building", eval.meth, "binaries")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- asub(thresholds, eval.meth[drop = FALSE], 1, drop = FALSE)[, i]
          writeRaster(x = bm_BinaryTransformation(raster(file.tmp, RAT = FALSE), thres.tmp),
                      filename = sub(output.format,
                                     paste0("_", eval.meth, "bin", output.format),
                                     file.tmp),
                      overwrite = TRUE,
                      datatype = "INT2S",
                      NAflag = -9999)
        }
      } else {
        nameBin <- paste0(nameProjSp, "_", eval.meth, "bin")
        assign(x = nameBin, value = bm_BinaryTransformation(proj, asub(thresholds, eval.meth[drop = FALSE], 1, drop = FALSE)))
        
        if (output.format == '.RData') {
          save(list = nameBin,
               file = file.path(namePath, paste0(nameBin, output.format)),
               compress = compress)
        } else {
          writeRaster(x = get(nameBin),
                      filename = file.path(namePath, paste0(nameBin, output.format)),
                      overwrite = TRUE,
                      datatype = "INT2S",
                      NAflag = -9999)
        }
      }
    }
    
    ## Do filtered transformation
    for (eval.meth in filtered.meth) {
      cat("\n\t> Building", eval.meth, "filtered")
      if (!do.stack) {
        for (i in 1:length(proj_out@proj@link)) {
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- asub(thresholds, eval.meth[drop = FALSE], 1, drop = FALSE)[, i]
          writeRaster(x = bm_BinaryTransformation(raster(file.tmp, RAT = FALSE), thres.tmp, doFiltering = TRUE),
                      filename = sub(output.format,
                                     paste0("_", eval.meth, "filt", output.format),
                                     file.tmp),
                      overwrite = TRUE,
                      datatype = ifelse(on_0_1000, "INT2S", "FLT4S"),
                      NAflag = -9999)
        }
      } else {
        nameFilt <- paste0(nameProjSp, "_", eval.meth, "filt")
        assign(x = nameFilt, value = bm_BinaryTransformation(proj, asub(thresholds, eval.meth[drop = FALSE], 1, drop = FALSE), doFiltering = TRUE))
        
        if (output.format == '.RData') {
          save(list = nameFilt,
               file = file.path(namePath, paste0(nameFilt, output.format)),
               compress = compress)
        } else {
          writeRaster(x = get(nameFilt),
                      filename = file.path(namePath, paste0(nameFilt, output.format)),
                      overwrite = TRUE ,
                      datatype = ifelse(on_0_1000, "INT2S", "FLT4S"),
                      NAflag = -9999)
        }
      }
    }
  }
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  ## save a copy of output object without value to be lighter
  nameOut <- paste0(modeling.output@sp.name, ".", proj.name, ".projection.out")
  assign(nameOut, free(proj_out))
  save(list = nameOut, file = file.path(namePath, nameOut))
  
  if (!silent) { .bmCat("Done") }
  return(proj_out)
}

###################################################################################################

.BIOMOD_Projection.check.args <- function(modeling.output, new.env, proj.name, xy.new.env,
                                          selected.models, binary.meth, filtered.meth,
                                          compress, do.stack, output.format)
{
  ## 1. Check modeling.output -------------------------------------------------
  .fun_testIfInherits(TRUE, "modeling.output", modeling.output, "BIOMOD.models.out")
  
  ## 2. Check proj.name -------------------------------------------------------
  if (is.null(proj.name)) {
    stop("\nYou must define a name for Projection Outputs")
  } else {
    dir.create(paste0(modeling.output@sp.name, '/proj_', proj.name, '/'), showWarnings = FALSE)
  }
  
  ## 3. Check new.env ---------------------------------------------------------
  .fun_testIfInherits(TRUE, "new.env", new.env, c('matrix', 'data.frame', 'RasterStack'))
  if (inherits(new.env, 'RasterStack')) {
    .fun_testIfIn(TRUE, "names(new.env)", names(new.env), modeling.output@expl.var.names))
  } else {
    .fun_testIfIn(TRUE, "colnames(new.env)", colnames(new.env), modeling.output@expl.var.names))
  }
  
  ## 4. Check xy.new.env ------------------------------------------------------
  if (!is.null(xy.new.env)  & !inherits(new.env, 'Raster')) {
    xy.new.env = data.matrix(xy.new.env)
    if (ncol(xy.new.env) != 2 || nrow(xy.new.env) != nrow(new.env)) {
      stop("invalid xy coordinates argument given -- dimensions mismatch !")
    }
  } else {
    xy.new.env = matrix()
  }
  
  ## 5. Check selected.models -------------------------------------------------
  if (selected.models[1] == 'all') {
    selected.models <- modeling.output@models.computed
  } else {
    selected.models <- intersect(selected.models, modeling.output@models.computed)
  }
  if (length(selected.models) < 1) {
    stop('No models selected')
  }
  
  ## check that given models exist
  files.check <- paste0(modeling.output@sp.name, '/models/', modeling.output@modeling.id, "/", selected.models)
  not.checked.files <- grep('MAXENT.Phillips|SRE', files.check)
  if (length(not.checked.files) > 0) {
    files.check <- files.check[-not.checked.files]
  }
  missing.files <- files.check[!file.exists(files.check)]
  if (length(missing.files) > 0) {
    stop(paste0("Projection files missing : ", toString(missing.files)))
    if (length(missing.files) == length(files.check)) {
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }
  
  ## 6. Check binary.meth & filtered.meth -------------------------------------
  if (!is.null(binary.meth) | !is.null(filtered.meth)) {
    models.evaluation <- get_evaluations(modeling.output)
    if (is.null(models.evaluation)) {
      warning("Binary and/or Filtered transformations of projection not ran because of models evaluation information missing")
    } else {
      available.evaluation <- unique(unlist(dimnames(models.evaluation)[1]))
      if (!is.null(binary.meth) && sum(!(binary.meth %in% available.evaluation)) > 0) {
        warning(paste0(toString(binary.meth[!(binary.meth %in% available.evaluation)]),
                       " Binary Transformation were switched off because no corresponding evaluation method found"))
        binary.meth <- binary.meth[binary.meth %in% available.evaluation]
      }
      
      if (!is.null(filtered.meth) && sum(!(filtered.meth %in% available.evaluation)) > 0) {
        warning(paste0(toString(filtered.meth[!(filtered.meth %in% available.evaluation)]),
                       " Filtered Transformation were switched off because no corresponding evaluation method found"))
        filtered.meth <- filtered.meth[filtered.meth %in% available.evaluation]
      }
    }
  }
  
  ## 7. Check compress --------------------------------------------------------
  if (compress == 'xz') {
    compress <- ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }
  
  ## 8. Check do.stack --------------------------------------------------------
  if (!inherits(new.env, 'RasterStack')) {
    if (!do.stack) { cat("\n\t\t! 'do.stack' arg is always set as TRUE for data.frame/matrix dataset") }
    do.stack <- TRUE
  } else if (do.stack) { # test if there is enough memory to work with RasterStack
    test = canProcessInMemory(raster::subset(new.env, 1), 2 * length(selected.models) + nlayers(new.env))
    if (!test) { rasterOptions(todisk = TRUE) }
  }
  
  ## 9. Check output.format ---------------------------------------------------
  if (!is.null(output.format)) {
    if (!output.format %in% c(".img", ".grd", ".RData")) {
      stop(paste0("output.format argument should be one of '.img','.grd' or '.RData'\n"
                  , "Note : '.img','.grd' are only available if you give environmental condition as a rasterStack object"))
    }
    if (output.format %in% c(".img", ".grd") && !inherits(new.env, "Raster")) {
      warning("output.format was automatically set to '.RData' because environmental conditions are not given as a raster object")
    }
  }
  
  ## set default values
  if (is.null(output.format)) {
    output.format <- ifelse(!inherits(new.env, "Raster"), ".RData", ".grd")
  }
  
  
  return(list(proj.name = proj.name,
              xy.new.env = xy.new.env,
              selected.models = selected.models,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth,
              compress = compress,
              do.stack = do.stack,
              output.format = output.format))
}


###################################################################################################

.build.clamping.mask <- function(env, MinMax)
{
  if (inherits(env, 'Raster'))
  { ## raster case ------------------------------------------------------------
    env <- raster::stack(env)
    env <- raster::subset(env, names(MinMax))
    
    ## create an empty mask
    clamp.mask <- ref.mask <- raster::subset(env, 1, drop = TRUE)
    clamp.mask[!is.na(clamp.mask[])] <- 0
    ref.mask[!is.na(clamp.mask[])] <- 1
    
    for (e.v in names(MinMax)) {
      if (!is.null(MinMax[[e.v]]$min)) { # numeric variable
        clamp.mask <- clamp.mask + 
          bm_BinaryTransformation(raster::subset(env, e.v, drop = TRUE), MinMax[[e.v]]$max) + ## pix with values outside of calib range
          (ref.mask - bm_BinaryTransformation(raster::subset(env, e.v, drop = TRUE), MinMax[[e.v]]$min)) ## pix with no values (within env[[1]] area)
      } else if (!is.null(MinMax[[e.v]]$levels)) { # factorial variable
        clamp.mask <- clamp.mask + 
          (raster::subset(env, e.v, drop = TRUE) %in% MinMax[[e.v]]$levels) ## pix with values outside of calib range
      }
    }
    
  } else if (is.data.frame(env) | is.matrix(env) | is.numeric(env))
  { ## matrix and data.frame case ---------------------------------------------
    env <- as.data.frame(env)
    
    # create an empty mask
    clamp.mask <- rep(0, nrow(env))
    
    for (e.v in names(MinMax)) {
      if (!is.null(MinMax[[e.v]]$min)) { # numeric variable
        clamp.mask <- clamp.mask + 
          bm_BinaryTransformation(env[, e.v], MinMax[[e.v]]$max) + ## pix with values outside of calib range
          (1 - bm_BinaryTransformation(env[, e.v], MinMax[[e.v]]$min)) ## pix with no values (within env[[1]] area)
      } else if (!is.null(MinMax[[e.v]]$levels)) { # factorial variable
        clamp.mask <- clamp.mask + (env[, e.v] %in% MinMax[[e.v]]$levels)
      }
    }
  } else { stop("Unsupported env arg") }
  return(clamp.mask)
}

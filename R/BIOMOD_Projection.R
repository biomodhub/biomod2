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
'BIOMOD_Projection' <- function(modeling.output,
                                new.env,
                                proj.name,
                                xy.new.env = NULL,
                                selected.models = 'all',
                                binary.meth = NULL,
                                filtered.meth = NULL,
                                compress = TRUE,
                                build.clamping.mask = TRUE,
                                ...){

  # 0. get additional args =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- list(...)
  omit.na <- args$omit.na # omit all non full filled environemental cells (always TRUE if env is a raster object)
  silent <- args$silent # echo advancement or not
  do.stack <- args$do.stack # store output in a lone stack
  clamping.level <- args$clamping.levels # remove all cells where at least clamping.level variables are out of their calibrating range
#   clamping.value <- args$clamping.value # reference value for clamped cells
  output.format <- args$output.format # raster output format
  keep.in.memory <- args$keep.in.memory # store results on memory or only on hard drive
  on_0_1000 <- args$on_0_1000 # transform projections on a 0 - 1000 scale to limit memory consumption
  split.proj <- args$split.proj

  if(is.null(omit.na)) omit.na <- TRUE
  if(is.null(silent)) silent <- FALSE
  if(is.null(do.stack)) do.stack <- TRUE
#   if(is.null(clamping.value)) clamping.value <- -1
  if(is.null(keep.in.memory)) keep.in.memory <- TRUE
  if(is.null(on_0_1000)) on_0_1000 <- TRUE # by default we return projections on a 0 -  1000 scale.
  if(is.null(split.proj)) split.proj <- 1 # by default we return projections on a 0 -  1000 scale.

#   if(!do.stack | !keep.in.memory) rasterOptions(todisk=TRUE)

#   if(is.null(output.format)){
#     if(!inherits(new.env,"Raster"))
#       output.format <- ".RData"
#     else
#       output.format <- ".grd"
#   }

  if(!silent){
    .bmCat("Do Models Projections")
  }

  # 1. Some Inputs args checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  args <- .BIOMOD_Projection.check.args(modeling.output,
                                        new.env,
                                        proj.name,
                                        xy.new.env,
                                        selected.models,
                                        binary.meth,
                                        filtered.meth,
                                        compress,
                                        do.stack,
                                        output.format)#, clamping.level)

  proj.name <- args$proj.name
  selected.models <- args$selected.models
  binary.meth <- args$binary.meth
  filtered.meth <- args$filtered.meth
  compress <- args$compress
  do.stack <- args$do.stack
  xy.new.env <- args$xy.new.env
  output.format <- args$output.format
#   clamping.level <- args$clamping.level

  rm(args)

  # 1b. Creating the outpput object
  proj_out <- new('BIOMOD.projection.out',
              proj.names = proj.name,
              sp.name =  modeling.output@sp.name,
              expl.var.names = modeling.output@expl.var.names,
              models.projected = selected.models,
              scaled.models = modeling.output@rescal.all.models,
              xy.coord = xy.new.env,
              modeling.object.id = modeling.output@modeling.id)

  # add link to biomod2 modeling object
  proj_out@modeling.object@link = modeling.output@link


  # adapting the proj slot to projection data type (e.g. rasterStack, array)
#   if(!do.stack){
#     proj_out@proj <- new('BIOMOD.stored.files')
#   } else if(inherits(new.env, 'Raster')){
#     proj_out@proj <- new('BIOMOD.stored.raster.stack')
#   } else{
#     proj_out@proj <- new('BIOMOD.stored.array')
#   }
  if(inherits(new.env, 'Raster')){
    proj_out@proj <- new('BIOMOD.stored.raster.stack')
  } else{
    proj_out@proj <- new('BIOMOD.stored.array')
  }

  # 1.c creating output directory
  dir.create(file.path(modeling.output@sp.name,paste("proj_", proj.name, sep="")),
             showWarnings = FALSE, recursive = TRUE, mode = "777")

  # 1.c Define the clamping mask
  if(build.clamping.mask){
    if(!silent) cat("\n\t> Building clamping mask\n")
    MinMax <- get_formal_data(modeling.output,'MinMax')

    assign(x = paste("proj_",proj.name,"_",modeling.output@sp.name,"_ClampingMask",sep=""),
           value = .build.clamping.mask(new.env, MinMax) )

    if(output.format == '.RData'){
      save(list=paste("proj_",proj.name,"_",modeling.output@sp.name,"_ClampingMask",sep=""),
           file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_ClampingMask", output.format ,sep="") ), compress=compress )
    } else {
      writeRaster(x=get(paste("proj_",proj.name,"_",modeling.output@sp.name,"_ClampingMask",sep="")),
                  filename=file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_ClampingMask", output.format ,sep="")),
                  datatype = "INT2S", NAflag=-9999,
                  overwrite=TRUE)
    }


  }

  # 2. Making projections -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!do.stack){
    dir.create(file.path(modeling.output@sp.name,paste("proj_", proj.name, sep=""), "individual_projections"),
               showWarnings = FALSE, recursive = TRUE, mode = "777")
  }

  saved.files <- c()

  proj <- lapply(selected.models, function(mod.name){
    cat("\n\t> Projecting",mod.name,"...")
    ## load biomod model
    BIOMOD_LoadModels(modeling.output, full.name=mod.name, as="mod")

    filename <- NULL
    if(!do.stack){
      filename <- file.path(modeling.output@sp.name, paste("proj_", proj.name, sep=""), "individual_projections", paste("proj_", proj.name, "_", mod.name,ifelse(output.format==".RData",".grd",output.format), sep="") )
    }

    pred.tmp <- predict(mod, new.env, on_0_1000=on_0_1000, filename=filename, omit.na=omit.na, split.proj = split.proj)

    if(do.stack){ # return the prediction only if stack has to be build
      return(pred.tmp)
    } else {
      return(filename)
    }

  })

  # 2b. Puting outputs in the right format =-=-=-=-=-=-=-=-=-=-=-= #
  if(!do.stack){
    saved.files <- unlist(proj)
  } else {
    if(inherits(new.env, "Raster")){
      proj <- stack(proj)
      names(proj) <- selected.models
    } else {
      proj <- as.data.frame(proj)
      names(proj) <- selected.models
      proj <- DF_to_ARRAY(proj)
    }
    if( keep.in.memory ){
      proj_out@proj@val <- proj
      proj_out@proj@inMemory <- TRUE
    }
  }

  ## save projections
  assign(x = paste("proj_",proj.name, "_", modeling.output@sp.name, sep=""),
         value = proj)

  if(output.format == '.RData'){
    saved.files <- file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name, output.format ,sep=""))
    save(list = paste("proj_",proj.name, "_", modeling.output@sp.name, sep=""),
         file = saved.files, compress=compress)
  } else {
    if(do.stack){
      saved.files <- file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name, output.format ,sep=""))
      writeRaster(x=get(paste("proj_",proj.name, "_", modeling.output@sp.name, sep="")),
                  filename=saved.files, overwrite=TRUE,
                  datatype = "INT2S", NAflag=-9999
                  # datatype = "INT2U"
                  )

    }
  }

  proj_out@type <- class(proj_out@proj@val)
  proj_out@proj@link <- saved.files


  # 3. Compute binary and filtered transformation =-=-=-=-=-=-=-=- #
  if(!is.null(binary.meth) | !is.null(filtered.meth)){
    cat("\n")
    eval.meth <- unique(c(binary.meth,filtered.meth))

    ## get all treshold
    if(inherits(new.env, "Raster")){
      thresholds <- matrix(0,nrow=length(eval.meth), ncol=length(selected.models),
                           dimnames=list(eval.meth,selected.models))
      for(mod in selected.models){
        PA.run   <- .extractModelNamesInfo(model.names=mod, info='data.set')
        eval.run <- .extractModelNamesInfo(model.names=mod, info='run.eval')
        algo.run <- .extractModelNamesInfo(model.names=mod, info='models')
        thresholds[eval.meth,mod] <- get_evaluations(modeling.output)[eval.meth,"Cutoff",algo.run,eval.run,PA.run]
        if(! on_0_1000) thresholds[eval.meth,mod]  <- thresholds[eval.meth,mod] / 1000
      }
    } else{
      thresholds <- array(0,dim=c(length(eval.meth),dim(proj)[-1]) ,
                         dimnames=c(list(eval.meth), dimnames(proj)[-1]) )
      for(mod in selected.models){
        PA.run   <- .extractModelNamesInfo(model.names=mod, info='data.set')
        eval.run <- .extractModelNamesInfo(model.names=mod, info='run.eval')
        algo.run <- .extractModelNamesInfo(model.names=mod, info='models')
        thresholds[eval.meth,algo.run,eval.run,PA.run] <- get_evaluations(modeling.output)[eval.meth,"Cutoff",algo.run,eval.run,PA.run]
        if(! on_0_1000) thresholds[eval.meth,algo.run,eval.run,PA.run]  <- thresholds[eval.meth,algo.run,eval.run,PA.run] / 1000
      }
    }

    ## do binary transformation
    for(eval.meth in binary.meth){
      cat("\n\t> Building", eval.meth,"binaries")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)[,i]
          writeRaster(x = BinaryTransformation(raster(file.tmp, RAT=FALSE),thres.tmp),
                      filename = sub(output.format, paste("_",eval.meth,"bin", output.format, sep=""), file.tmp),
                      overwrite=TRUE,
                      datatype = "INT2S",NAflag=-9999)
        }
      } else {
      assign(x = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep=""),
             value = BinaryTransformation(proj,asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)))

      if(output.format == '.RData'){
        save(list = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep=""),
             file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"bin", output.format ,sep="")), compress=compress)
      } else {
        writeRaster(x=get(paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep="")),
                    filename=file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"bin", output.format ,sep="")),
                    overwrite=TRUE,
                    datatype = "INT2S", NAflag=-9999)
      }

      rm(list=paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"bin", sep=""))
      }
    }

    ## do filtered transformation
    for(eval.meth in filtered.meth){
      cat("\n\t> Building", eval.meth,"filtered")
      if(!do.stack){
        for(i in 1:length(proj_out@proj@link)){
          file.tmp <- proj_out@proj@link[i]
          thres.tmp <- asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)[,i]
          writeRaster(x = FilteringTransformation(raster(file.tmp, RAT=FALSE),thres.tmp),
                      filename = sub(output.format, paste("_",eval.meth,"filt", output.format, sep=""), file.tmp),
                      overwrite=TRUE, datatype=ifelse(on_0_1000,"INT2S","FLT4S"), NAflag=-9999
                      )
        }
      } else {
        assign(x = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep=""),
               value = FilteringTransformation(proj,asub(thresholds, eval.meth[drop=FALSE], 1, drop=FALSE)))

        if(output.format == '.RData'){
          save(list = paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep=""),
               file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"filt", output.format ,sep="")), compress=compress)
        } else {
          writeRaster(x=get(paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep="")),
                      filename=file.path(modeling.output@sp.name, paste("proj_", proj.name, sep= ""), paste("proj_",proj.name,"_", modeling.output@sp.name,"_",eval.meth,"filt", output.format ,sep="")),
                      overwrite=TRUE , datatype=ifelse(on_0_1000,"INT2S","FLT4S"), NAflag=-9999
                      )
        }

        rm(list=paste("proj_",proj.name, "_", modeling.output@sp.name,"_",eval.meth,"filt", sep=""))
      }
    }

  }
  # End compute binary and filtered transformation -=-=-=-=-=-=-=- #

#   # 2.b a posteriori clamping...
#   #### TO DO : make an a priori claming
#   if(!is.null(clamping.level)){
#     cat("\n   > clamping projections...")
#     if(inherits(proj_out@proj@val, 'Raster')){
#       proj_out@proj@val <- raster::stack(reclassify(proj_out@proj@val * reclassify( (- 1 * (clampMask)), c(-0.5,0.5,1)), c(-Inf,0,clamped.value) ))
#     } else{
#       proj_out@proj@val[which(clampMask > 0),,,] <- clamped.value
#     }
#
#   }

  # save a copy of output object without value to be lighter
  assign(paste(modeling.output@sp.name,".", proj.name, ".projection.out", sep=""), free(proj_out))
  save(list = paste(modeling.output@sp.name,".", proj.name, ".projection.out", sep=""),
       file = file.path(modeling.output@sp.name, paste("proj_", proj.name, sep=""), paste(modeling.output@sp.name,".", proj.name, ".projection.out", sep="")))

  if(!silent) .bmCat("Done")

  # 4. Returning output
  return(proj_out)
}

####################################################################################################
### Utilities Fuctions #############################################################################
####################################################################################################

.BIOMOD_Projection.check.args <- function(modeling.output, new.env, proj.name, xy.new.env,
                                          selected.models, binary.meth, filtered.meth,
                                          compress, do.stack, output.format){#, clamping.level){
  ## modeling.output
  if(!inherits(modeling.output, 'BIOMOD.models.out')){
    stop("'modeling.output' must be the result of BIOMOD_Modeling() computation")
  }

  ## new.env
  # NOTE : may be reorder variables if necessary
  if(!inherits(new.env, c('matrix', 'data.frame', 'RasterStack'))){
    stop("'new.env' must be a matrix, a data.frame or a RasterStack")
  }
  if(inherits(new.env, 'RasterStack')){
    if(sum(!(names(new.env) %in% modeling.output@expl.var.names)) > 0 ){
      stop("'new.env' layer names don't match with explanatory variables used for buiding models")
    }
  } else{
    if(sum(!(colnames(new.env) %in% modeling.output@expl.var.names)) > 0 ){
      stop("'new.env' colnames don't match with explanatory variables used for buiding models")
    }
  }

  ## proj.name
  # The projection Name
  if(is.null(proj.name)){
    stop("\nYou must define a name for Projection Outputs")
  } else{
    dir.create(paste(modeling.output@sp.name,'/proj_',proj.name,'/',sep=''),
               showWarnings=FALSE)
  }

  ## xy.new.env
  if(!is.null(xy.new.env)  & !inherits(new.env,'Raster')){
    xy.new.env = data.matrix(xy.new.env)
    if(ncol(xy.new.env) != 2 | nrow(xy.new.env) != nrow(new.env)) stop("invalid xy coordinates argument given -- dimensions mismatch !")
  } else {
    xy.new.env = matrix()
  }

  ## selected.models
  if(selected.models[1] == 'all'){
    selected.models <- modeling.output@models.computed
  } else{
    selected.models <- intersect(selected.models, modeling.output@models.computed)
  }
  if(length(selected.models) < 1){
    stop('No models selected')
  }

  # check that given models exits
  files.check <- paste(modeling.output@sp.name,'/models/',modeling.output@modeling.id,"/",selected.models,sep='')
  not.checked.files <- c(grep('MAXENT.Phillips', files.check), grep('SRE', files.check))
  if(length(not.checked.files) > 0){files.check <- files.check[-not.checked.files]}
  missing.files <- files.check[!file.exists(files.check)]
  if( length(missing.files) > 0 ){
    stop(paste("Projection files missing : ", toString(missing.files), sep=''))
    if(length(missing.files) == length(files.check)){
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }

  # The binaries  and filtering transformations
  if(!is.null(binary.meth) | !is.null(filtered.meth)){
    models.evaluation <- get_evaluations(modeling.output)
    if(is.null(models.evaluation)){
      warning("Binary and/or Filtered transformations of projection not ran because of models
              evaluation information missing")
    } else{
      available.evaluation <- unique(unlist(dimnames(models.evaluation)[1]))
      if(!is.null(binary.meth)){
        if(sum(!(binary.meth %in% available.evaluation)) > 0){
          warning(paste(toString(binary.meth[!(binary.meth %in% available.evaluation)]),
                        " Binary Transformation were switched off because no corresponding",
                        " evaluation method found ", sep=""))
          binary.meth <- binary.meth[binary.meth %in% available.evaluation]
        }
      }

      if(!is.null(filtered.meth)){
        if(sum(!(filtered.meth %in% available.evaluation)) > 0){
          warning(paste(toString(filtered.meth[!(filtered.meth %in% available.evaluation)]),
                        " Filtered Transformation were switched off because no corresponding",
                        " evaluation method found ", sep=""))
          filtered.meth <- filtered.meth[filtered.meth %in% available.evaluation]
        }
      }
    }
  }

  ## compress
  if(compress == 'xz'){
    compress <- ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }

  ## do.stack
  if(!inherits(new.env, 'RasterStack')){
    if(!do.stack) cat("\n\t\t! 'do.stack' arg is always set as TRUE for data.frame/matrix dataset")
    do.stack <- TRUE
  } else{
    if(do.stack){
      # test if there is memory enough to work with RasterStack
      test = canProcessInMemory( raster::subset(new.env,1), 2*length(selected.models) + nlayers(new.env) )
      if (!test) rasterOptions(todisk=T)
#       if (!do.stack){
#         cat("\n   ! Results will be saved as individual RasterLayers because of a lack of memory !")
#       }
    }
  }

  ## output.format ##
  if(!is.null(output.format)){
    if(! output.format %in% c(".img",".grd",".RData")){
      stop("output.format argument should be one of '.img','.grd' or '.RData'\n Note : '.img','.grd' are only available if you give environmental condition as a rasterStack object")
    }
    if( output.format %in% c(".img",".grd") & !inherits(new.env,"Raster") ){
      warning("output.format was automatically set to '.RData' because environmental conditions are not given as a raster object")
    }
  }
  ## set default values
  if(is.null(output.format)){
    if(!inherits(new.env,"Raster"))
      output.format <- ".RData"
    else
      output.format <- ".grd"
  }


#   ## clamping checking
#   if(!is.null(clamping.level)){
#     if(!is.numeric(clamping.level)){
#       stop("clamping.level must be NULL or numeric")
#     }
#
#     # limit clamping level
#     if( clamping.level > length(modeling.output@expl.var.names)){
#       cat("\n   ! clamping.level was down to", length(modeling.output@expl.var.names))
#       clamping.level <- length(modeling.output@expl.var.names)
#     }
#
#     if( clamping.level < 1){
#       cat("\n   ! clamping was swich off")
#       clamping.level <- NULL
#     }
#   }


  return(list(#modeling.output = modeling.output,
              #new.env = new.env,
              proj.name = proj.name,
              xy.new.env = xy.new.env,
              selected.models = selected.models,
              binary.meth = binary.meth,
              filtered.meth = filtered.meth,
              compress = compress,
              do.stack = do.stack,
              output.format = output.format))#, clamping.level = clamping.level))

}

####################################################################################################
####################################################################################################

.build.clamping.mask <- function(env, MinMax){
  if(inherits(env,'Raster')){ # raster case
    env <- raster::stack(env)
    env <- raster::subset(env,names(MinMax))

    # create an empty mask
#     clamp.mask <- reclassify( raster::subset(env,1, drop=TRUE), c(-Inf,Inf,0) )
    clamp.mask <- ref.mask <- raster::subset(env,1, drop=TRUE)
    clamp.mask[!is.na(clamp.mask[])] <- 0
    ref.mask[!is.na(clamp.mask[])] <- 1

    for(e.v in names(MinMax)){

      if(!is.null(MinMax[[e.v]]$min)){ # numeric variable
        clamp.mask <- clamp.mask + ( BinaryTransformation(raster::subset(env, e.v, drop=TRUE), MinMax[[e.v]]$max ) +  (ref.mask - BinaryTransformation(raster::subset(env, e.v, drop=TRUE), MinMax[[e.v]]$min )) )

      } else if(!is.null(MinMax[[e.v]]$levels)){ # factorial variable
        clamp.mask <- clamp.mask + (raster::subset(env, e.v, drop=TRUE) %in% MinMax[[e.v]]$levels)
      }
    }

    ## fix projection system

  } else if(is.data.frame(env) | is.matrix(env) | is.numeric(env)){ # matrix and data.frame case
    env <- as.data.frame(env)

    # create an empty mask
    clamp.mask <- rep(0,nrow(env))

    for(e.v in names(MinMax)){
      if(!is.null(MinMax[[e.v]]$min)){ # numeric variable
        clamp.mask <- clamp.mask + ( BinaryTransformation(env[,e.v], MinMax[[e.v]]$max ) +
          (1 - BinaryTransformation(env[,e.v], MinMax[[e.v]]$min )) )
      } else if(!is.null(MinMax[[e.v]]$levels)){ # factorial variable
        clamp.mask <- clamp.mask + (env[,e.v] %in% MinMax[[e.v]]$levels)
      }
    }

  } else{
    stop("Unsupported env arg")
  }

  return(clamp.mask)

}

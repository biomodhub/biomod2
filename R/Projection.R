setGeneric( "Projection", 
            def = function(models.name,
                           modeling.work.dir = getwd(),
                           new.env.data ,
                           ...){
                            standardGeneric( "Projection" )
                            } )

setMethod( 'Projection', signature(new.env.data = 'data.frame'),
  function(models.name,
           modeling.work.dir = getwd(),
           new.env.data ,
           xy = NULL,
           proj.name = NULL,
           binary.proj = NULL,
           filtred.proj = NULL,
           models.evaluation = NULL,
           models.options = NULL,
           compress = TRUE,
           scaled.models=TRUE,
           do.stack = FALSE){
        
#     # 1. loading resuired libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#     .Models.dependencies(silent=TRUE, models.options=models.options )
  
    # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    sp.name  <- .extractModelNamesInfo(model.names=models.name, info='species')
    PA.run   <- .extractModelNamesInfo(model.names=models.name, info='data.set')
    eval.run <- .extractModelNamesInfo(model.names=models.name, info='run.eval')
    algo.run <- .extractModelNamesInfo(model.names=models.name, info='models')
    
    # 3. Printing Projection Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
#     cat('\nDoing Models Projections...')
#     if(length(grep('EF.',models.name)) > 0 ){
#       kept.models.name <- models.name[-grep('EF.',models.name)] 
#       kept.algo.run <- algo.run[-grep('EF.',algo.run)]
#     } else {
#       kept.models.name <- models.name
#       kept.algo.run <- algo.run
#     }
    
    proj.array <- lapply(models.name, .Projection.do.proj, env=new.env.data, xy=xy, scaled.models=scaled.models, proj.name=paste("proj_",proj.name, sep=""), models.options=models.options)
    proj.array <- as.data.frame(proj.array)
    names(proj.array) <- models.name
    
    # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(binary.proj)>0){
      cat("\nBinary transformations...")
      lapply(binary.proj, function(bin.proj){
        
        cuts <- unlist(lapply(names(proj.array), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
          }))

        proj.bin.array <- BinaryTransformation(proj.array, cuts)
        proj.bin.array <- DF_to_ARRAY(proj.bin.array)

        eval(parse(text = paste(proj.name,"_",sp.name,"_bin_",bin.proj, "_array <- proj.bin.array", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_bin_",bin.proj,
                                "_array, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_bin_",bin.proj,"_array' )",sep="")))
        
        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_bin_",bin.proj,"_array , proj.bin.array, cuts)", sep="" )))
      })
    }

    
    # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    if(length(filtred.proj)>0){
      cat("\nFiltred transformations...")
      lapply(filtred.proj, function(filt.proj){
        
        cuts <- unlist(lapply(names(proj.array), function(x){
          mod <- tail(unlist(strsplit(x,"_")), 3)[3]
          run <- tail(unlist(strsplit(x,"_")), 3)[2]
          dat <- tail(unlist(strsplit(x,"_")), 3)[1]
          return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
        }))
        
        proj.filt.array <- FilteringTransformation(proj.array, cuts)
        proj.filt.array <- DF_to_ARRAY(proj.filt.array)
        
        eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_array <- proj.filt.array", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                "_array, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                proj.name,"_",sp.name,"_filt_",filt.proj,"_array' )",sep="")))
        
        eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_array , proj.filt.array, cuts)", sep="" )))
      })
    }
    
    proj.array <- DF_to_ARRAY(proj.array)
    
    # 7. Saving projection on hard disk
    eval(parse(text = paste(proj.name,"_",sp.name, " <- proj.array", sep="")))
    eval(parse(text = paste("save(",proj.name,"_",sp.name, ", file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                            proj.name,"_",sp.name,"' )",sep="")))
    eval(parse(text = paste("rm(",proj.name,"_",sp.name,")", sep="" )))
    gc(reset=TRUE)
    
    return(invisible(proj.array))   
  })








setMethod( 'Projection', signature(new.env.data = 'RasterStack'),
  function(models.name,
           modeling.work.dir = getwd(),
           new.env.data ,
           xy = NULL,           
           proj.name = NULL,
           binary.proj = NULL,
           filtred.proj = NULL,
           models.evaluation = NULL,
           models.options = NULL,
           stack = TRUE,
           compress = TRUE,
           scaled.models=TRUE,
           do.stack = FALSE){
        
    # 1. loading resuired libraries =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
    .Models.dependencies(silent=TRUE, models.options=models.options)
    
    # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    sp.name  <- .extractModelNamesInfo(model.names=models.name, info='species')
    PA.run   <- .extractModelNamesInfo(model.names=models.name, info='data.set')
    eval.run <- .extractModelNamesInfo(model.names=models.name, info='run.eval')
    algo.run <- .extractModelNamesInfo(model.names=models.name, info='models')
       
    
    # 3. Printing Projection Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    
    # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
    cat('\nDoing Models Projections...')
#     if(length(grep('EF.',models.name)) > 0 ){
#       kept.models.name <- models.name[-grep('EF.',models.name)] 
#       kept.algo.run <- algo.run[-grep('EF.',algo.run)]
#     } else {
#       kept.models.name <- models.name
#       kept.algo.run <- algo.run
#     }
    
    if(do.stack){
      proj.ras <- lapply(models.name, .Projection.do.proj, env=new.env.data, scaled.models=scaled.models, proj.name=paste("proj_",proj.name, sep=""), models.options=models.options)
  
      # transform list of rasterLayers into a rasterStack
      proj.stack <- stack(x = proj.ras)
  
      names(proj.stack) <- models.name #names(proj.ras.mod)
      rm(proj.ras)

      # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(length(binary.proj)>0){
        cat("\nBinary transformations...")
        lapply(binary.proj, function(bin.proj){
          
          cuts <- unlist(lapply(names(proj.stack), function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
            }))
  
          proj.bin.stack <- BinaryTransformation(proj.stack, cuts)
          names(proj.bin.stack) <- paste(names(proj.stack), ".bin", sep="")
  
          eval(parse(text = paste(proj.name,"_",sp.name,"_bin_",bin.proj, "_RasterStack <- proj.bin.stack", sep="")))
          eval(parse(text = paste("save(",proj.name,"_",sp.name,"_bin_",bin.proj,
                                  "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                  proj.name,"_",sp.name,"_bin_",bin.proj,"_RasterStack' )",sep="")))
          
          eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_bin_",bin.proj,"_RasterStack , proj.bin.stack, cuts)", sep="" )))
        })
      }
      
      # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
      if(length(filtred.proj)>0){
        cat("\nFiltered transformations...")
        lapply(filtred.proj, function(filt.proj){
          
          cuts <- unlist(lapply(names(proj.stack), function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
          }))
          
          proj.filt.stack <- FilteringTransformation(proj.stack, cuts)
          names(proj.filt.stack) <- paste(names(proj.stack), ".filt", sep="")
          
          eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_RasterStack <- proj.filt.stack", sep="")))
          eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                  "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                  proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterStack' )",sep="")))
          
          eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterStack , proj.filt.stack, cuts)", sep="" )))
        })
      }
      
      # 7. Saving projection on hard disk
       eval(parse(text = paste(proj.name,"_",sp.name, "_RasterStack <- proj.stack", sep="")))
       eval(parse(text = paste("save(",proj.name,"_",sp.name, "_RasterStack, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                               proj.name,"_",sp.name,"_RasterStack' )",sep="")))
       eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_RasterStack)", sep="" )))
       gc(reset=TRUE)
    } else{
      
      # all models will be saved separatly
      proj.stack <- c() # list of saved files
      for(m.n in models.name){
        
        proj.ras <- .Projection.do.proj(m.n, env=new.env.data, scaled.models=scaled.models, proj.name=paste("proj_",proj.name, sep=""), models.options=models.options)
        names(proj.ras) <- m.n #names(proj.ras.mod)

        # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(length(binary.proj)>0){
          cat("\nBinary transformations...")
          lapply(binary.proj, function(bin.proj){
            
            cuts <- unlist(lapply(names(proj.ras), function(x){
              mod <- tail(unlist(strsplit(x,"_")), 3)[3]
              run <- tail(unlist(strsplit(x,"_")), 3)[2]
              dat <- tail(unlist(strsplit(x,"_")), 3)[1]
              return(models.evaluation[bin.proj,"Cutoff", mod, run, dat])
              }))
    
            proj.bin.ras <- BinaryTransformation(proj.ras, cuts)
            names(proj.bin.ras) <- paste(names(proj.ras), ".bin", sep="")
    
            eval(parse(text = paste(proj.name,"_",m.n,"_bin_",bin.proj, "_RasterLayer <- proj.bin.ras", sep="")))
            eval(parse(text = paste("save(",proj.name,"_",m.n,"_bin_",bin.proj,
                                    "_RasterLayer, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                    proj.name,"_",m.n,"_bin_",bin.proj,"_RasterLayer' )",sep="")))
            
            eval(parse(text = paste("rm(",proj.name,"_",m.n,"_bin_",bin.proj,"_RasterLayer , proj.bin.ras, cuts)", sep="" )))
          })
        }
        
        # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(length(filtred.proj)>0){
          cat("\nFiltered transformations...")
          lapply(filtred.proj, function(filt.proj){
            
            cuts <- unlist(lapply(names(proj.ras), function(x){
              mod <- tail(unlist(strsplit(x,"_")), 3)[3]
              run <- tail(unlist(strsplit(x,"_")), 3)[2]
              dat <- tail(unlist(strsplit(x,"_")), 3)[1]
              return(models.evaluation[filt.proj,"Cutoff", mod, run, dat])
            }))
            
            proj.filt.ras <- FilteringTransformation(proj.ras, cuts)
            names(proj.filt.ras) <- paste(names(proj.ras), ".filt", sep="")
            
            eval(parse(text = paste(proj.name,"_",sp.name,"_filt_",filt.proj, "_RasterLayer <- proj.filt.ras", sep="")))
            eval(parse(text = paste("save(",proj.name,"_",sp.name,"_filt_",filt.proj,
                                    "_RasterLayer, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                    proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterLayer' )",sep="")))
            
            eval(parse(text = paste("rm(",proj.name,"_",sp.name,"_filt_",filt.proj,"_RasterLayer , proj.filt.ras, cuts)", sep="" )))
          })
        }
        
        # 7. Saving projection on hard disk
        eval(parse(text = paste(proj.name,"_",m.n, "_RasterLayer <- proj.ras", sep="")))
        eval(parse(text = paste("save(",proj.name,"_",m.n, "_RasterLayer, file = '",modeling.work.dir,"/",sp.name,"/proj_",proj.name,"/",
                                 proj.name,"_",m.n,"_RasterLayer' )",sep="")))
        proj.stack <- c( proj.stack, paste(proj.name,"_",m.n,"_RasterLayer",sep="") )
        eval(parse(text = paste("rm(",proj.name,"_",m.n,"_RasterLayer)", sep="" )))
        gc(reset=TRUE)
      }
      
    }
    
    ## remove MAXENT.Phillips tmp dir if exists
    if(file.exists(file.path(sp.name, proj.name, 'MaxentTmpData'))){
      .Delete.Maxent.WorkDir( file.path(sp.name, proj.name) )
    }
    
    return(invisible(proj.stack))   
  })


setGeneric( ".Projection.do.proj", 
            def = function(model.name, env, model.dir = NULL,...){
                    standardGeneric( ".Projection.do.proj" )
                    } )

setMethod('.Projection.do.proj', signature(env='data.frame'),
  function(model.name, env, xy = NULL, model.dir = NULL, scaled.models=TRUE, proj.name=NULL, models.options=NULL){
    cat('\n\t>', model.name)
    # automaticly fill model.dir if not given
    if(is.null(model.dir)){
      model.dir <- paste(getwd(),'/',.extractModelNamesInfo(model.name, info='species'),'/models', sep="")
    }
    
    # loading model
    if(length(c(grep('SRE',model.name) )) == 0){
      model.sp = eval(parse(text = load(paste(model.dir,'/',model.name, sep=""))) )
      eval(parse(text=paste("rm(",model.name,")",sep="")))
    }
  
    # check model.type
    model.type <- tail(unlist(strsplit(model.name, split="_")),1)
    if(!( model.type %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF', 'MAXENT.Phillips') )){
      if(!grep('EF.',model.type))
        stop('Unknow model type')
    }
    
    if(model.type == 'ANN'){
      set.seed(555) # to be able to refind our trees MAY BE BAD
      # proj automaticly scaled
      return(data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env, type = "raw")),
                                              name = model.name ) * 1000)))
    }
    
    if(model.type == 'CTA'){
      proj <- as.integer(as.numeric(predict(model.sp, env,type="prob")[,2]) * 1000)
    
      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
      
      return( data.frame( proj = proj) )
    }
    
    if(model.type == 'FDA'){
      # proj automaticly scaled
      return( data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env, type = "posterior")[, 2]),
                                              name = model.name) * 1000)))
    }
    
    if(model.type == 'GBM'){
      best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
      proj <- as.integer(predict.gbm(model.sp, env, best.iter, type = "response") * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }

      return(data.frame( proj = proj))
    }
    
    if(model.type == 'GLM'){
      proj <- as.integer(.testnull(model.sp, Prev=0.5, env) * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
      
      return( data.frame(proj = proj) )
    }
    
    if(model.type == 'GAM'){
      proj <- as.integer(.testnull(model.sp, Prev=0.5, env) * 1000)

      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
 
      return( data.frame( proj = proj ) )
    }
    
#     if(model.type == 'MARS'){
#       # proj automaticly scaled
#       return(data.frame( proj = as.integer(.Rescaler5(as.numeric(predict(model.sp, env)), 
#                                               name = model.name) * 1000)))
#     }
    
    if(model.type == 'RF'){
      proj <- as.integer(as.numeric(predict(model.sp,env, type='prob')[,'1']) *1000)
      
      if(scaled.models){
        proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
      }
      
      return( data.frame( proj = proj ))
    }
    
    if(model.type == 'SRE'){
      # loading data of the correspunding run
      load(paste(model.dir,'/',model.name,'/Data_',model.name, sep=""))                
      return(eval(parse(text=paste("sre(Data_",model.name,"$Response, Data_",
                                   model.name,"$Explanatory, env, Data_",model.name,
                                   "$Quant)*1000", sep=""))))
    }
    
    if(model.type == 'MAXENT.Phillips'){
      if(!is.null(xy)){
        proj <- as.integer(predict( object=model.sp, newdata=env, proj_name=proj.name, xy=xy) * 1000)
        
#         if(scaled.models){
          proj <- as.integer(.Rescaler5(proj/1000, name = model.name ) * 1000)
#         }        
        
        return( data.frame( proj = proj ))
#         
#         .Prepare.Maxent.Proj.WorkDir(env, xy, proj.name=file.path(.extractModelNamesInfo(model.name, info='species'), proj.name ))
#         
#         cat("\t Running Maxent...")
#         
#         system(command=paste("java -cp ", file.path(models.options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"), " density.Project \"", model.dir,.Platform$file.sep,
#                              model.name, .Platform$file.sep ,sub("_MAXENT.Phillips","",model.name),
#                              ".lambdas\" ", .extractModelNamesInfo(model.name, info='species'), "/", proj.name, "/MaxentTmpData/Proj_swd.csv ", .extractModelNamesInfo(model.name, info='species'), "/", proj.name, "/MaxentTmpData/projMaxent", sep=""), wait = TRUE)
#         
#         maxent.proj <- read.asciigrid(paste(species.name=.extractModelNamesInfo(model.name, info='species'), "/", proj.name , "/MaxentTmpData/projMaxent.asc", sep=""))@data
#         .Delete.Maxent.WorkDir(species.name=paste(species.name=.extractModelNamesInfo(model.name, info='species'), "/", proj.name,sep=""))
#         return(proj = as.integer(.Rescaler5(as.numeric(maxent.proj[,1]), 
#                                               name = model.name) * 1000))
      } else {
        cat('\n MAXENT.Phillips need coordinates to run! NA returned ')
        return(data.frame(rep(NA,nrow(env))))
      }
    }
          
  })








setMethod('.Projection.do.proj', signature(env='RasterStack'),
  function(model.name, env, model.dir = NULL, scaled.models=TRUE, proj.name=NULL, models.options=NULL){
    cat('\n\t>', model.name)
    
    # automaticly fill model.dir if not given
    if(is.null(model.dir)){
      model.dir <- paste(getwd(),'/',.extractModelNamesInfo(model.name, info='species'),'/models', sep="")
    }
    
    # loading model
    if(length(c(grep('SRE',model.name))) == 0){
      model.sp = eval(parse(text = load(paste(model.dir,'/',model.name, sep=""))) )
      eval(parse(text=paste("rm(",model.name,")",sep="")))
    }
  
    # check model.type
    model.type <- tail(unlist(strsplit(model.name, split="_")),1)
    if(!( model.type %in% c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF', 'MAXENT.Phillips') )){
      if(!grep('EF.',model.type))
        stop('Unknow model type')
    }
    
    if(model.type == 'ANN'){
      set.seed(555) # to be able to refind our trees MAY BE BAD
      proj.ras <- predict(env, model.sp, type="raw")
      proj.ras[!is.na(proj.ras[])] <- .Rescaler5(proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(proj.ras * 1000))
    }
    
    if(model.type == 'CTA'){
      set.seed(123) # to be able to refind our trees MAY BE BAD
      proj.ras <- predict(env, model=model.sp, type='prob', index=2)
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'FDA'){
      pred.ras <- predict(env, model.sp, type="post", index=2)
      pred.ras[!is.na(pred.ras[])] <- .Rescaler5(pred.ras[!is.na(pred.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(pred.ras * 1000))
    }
    
    if(model.type == 'GBM'){
      if(file.exists(paste(model.dir,'/',model.name,'_best.iter'))){
        load(paste(model.dir,'/',model.name,'_best.iter'))
      } else{
        best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
      }
        
      proj.ras <- predict(env, model.sp, n.trees=best.iter, type='response')
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'GLM'){
      proj.ras <- predict(env, model=model.sp, type='response')
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'GAM'){
      proj.ras <- predict(env, model=model.sp, type='response')
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )
    }
    
    if(model.type == 'MARS'){
      pred.ras <- predict(env, model.sp)
      pred.ras[!is.na(pred.ras[])] <- .Rescaler5(pred.ras[!is.na(pred.ras[])], ref=NULL,
                                                 name=model.name, original=FALSE)
      return( round(pred.ras * 1000) )
    }
    
    if(model.type == 'RF'){
      proj.ras <- predict(env, model=model.sp, type='prob', index=2)
      if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5( proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                    name=model.name, original=FALSE)
      }
      return( round(proj.ras*1000) )      
    }
    
    if(model.type == 'SRE'){
#       cat('\n SRE prediction not supported yet ! ')
      load(paste(model.dir,'/',model.name,'/Data_',model.name, sep=""))
      data.sre <- get(paste('Data_',model.name, sep=""))
      rm(list=paste('Data_',model.name, sep=""))
#       sre.out <- eval(parse(text=paste("sre(Data_",model.name,"$Response, Data_",
#                                    model.name,"$Explanatory, env, Data_",model.name,
#                                    "$Quant)*1000", sep="")))
      sre.out <- raster::subset(sre(data.sre$Response, data.sre$Explanatory, env, data.sre$Quant), 1, drop=TRUE) * 1000
      
      return(sre.out)
    }
    
    if(model.type == 'MAXENT.Phillips'){
      proj.ras <- predict( object=model.sp, newdata=env, proj_name=proj.name)
      
#       if(scaled.models){
        proj.ras[!is.na(proj.ras[])] <- .Rescaler5(proj.ras[!is.na(proj.ras[])], ref=NULL,
                                                   name=model.name, original=FALSE)
#       }        
      
      return(round(proj.ras*1000))
      
      
         
    }   
  })


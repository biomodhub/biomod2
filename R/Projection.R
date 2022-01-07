
.getBinary = function(i.proj, names.proj, obj.eval, obj.data, doFiltering = FALSE, doStack = FALSE)
{
  cuts <- unlist(lapply(names.proj, function(x) {
    mod <- tail(unlist(strsplit(x, "_")), 3)[3]
    run <- tail(unlist(strsplit(x, "_")), 3)[2]
    dat <- tail(unlist(strsplit(x, "_")), 3)[1]
    return(obj.eval[i.proj,"Cutoff", mod, run, dat])
  }))
  
  res <- bm_BinaryTransformation(obj.data, cuts, doFiltering)
  if (doStack) {
    if (doFiltering) {
      names(res) = paste0(names(obj.data), ".filt")
    } else {
      names(res) = paste0(names(obj.data), ".bin")
    }
  } else {
    res = DF_to_ARRAY(res)
  }
  return(res)
}

.saveBinary = function(i.proj, names.proj, obj.eval, obj.data, doFiltering = FALSE, doStack = FALSE
                       , name.obj, name.file)
{
  val.obj = .getBinary(i.proj, names.proj, obj.eval, obj.data, doFiltering, doStack)
  assign(x = name.obj, value = val.obj)
  eval(parse(text = paste0("save(", name.obj, ", file = name.file)")))
}



setGeneric("Projection",
           def = function(models.name, modeling.work.dir = getwd(), new.env.data, ...)
           {
             standardGeneric( "Projection" )
           })

setMethod('Projection', signature(new.env.data = 'data.frame'),
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
                   do.stack = FALSE)
          {
            # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
            sp.name  <- .extractModelNamesInfo(model.names = models.name, info = 'species')
            name.folder = paste0("proj_", proj.name)
            name.sp = paste0(proj.name, "_", sp.name)
            
            # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
            proj.array <- lapply(models.name,
                                 .Projection.do.proj,
                                 env = new.env.data,
                                 xy = xy,
                                 scaled.models = scaled.models,
                                 proj.name = name.folder,
                                 models.options = models.options)
            proj.array <- as.data.frame(proj.array)
            names(proj.array) <- models.name
            
            # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
            if (length(binary.proj) > 0) {
              cat("\nBinary transformations...")
              lapply(binary.proj, function(bin.proj)
              {
                name.obj = paste0(name.sp, "_bin_", bin.proj, "_array")
                name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.obj)
                .saveBinary(i.proj = bin.proj,
                            names.proj = names(proj.array),
                            obj.eval = models.evaluation,
                            obj.data = proj.array,
                            doFiltering = FALSE,
                            doStack = FALSE,
                            name.obj = name.obj,
                            name.file = name.file)
              })
            }
            
            if (length(filtred.proj) > 0) {
              cat("\nFiltered transformations...")
              lapply(filtred.proj, function(filt.proj)
              {
                name.obj = paste0(name.sp, "_filt_", bin.proj, "_array")
                name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.obj)
                .saveBinary(i.proj = filt.proj,
                            names.proj = names(proj.array),
                            obj.eval = models.evaluation,
                            obj.data = proj.array,
                            doFiltering = TRUE,
                            doStack = FALSE,
                            name.obj = name.obj,
                            name.file = name.file)
              })
            }
            
            # 7. Saving projection on hard disk
            proj.array <- DF_to_ARRAY(proj.array)
            assign(x = name.sp, value = proj.array)
            
            name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.sp)
            eval(parse(text = paste0("save(", name.sp, ", file = name.file)")))
            gc(reset = TRUE)
            
            return(invisible(proj.array))
          })



setMethod('Projection', signature(new.env.data = 'RasterStack'),
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
                   do.stack = FALSE)
          {
            # 2. extract model info  =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
            sp.name  <- .extractModelNamesInfo(model.names = models.name, info = 'species')
            name.folder = paste0("proj_", proj.name)
            name.sp = paste0(proj.name, "_", sp.name)
            
            # 4. Computing Projections =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
            cat('\nDoing Models Projections...')
            
            if (do.stack) {
              proj.stack <- lapply(models.name,
                                   .Projection.do.proj,
                                   env = new.env.data,
                                   scaled.models = scaled.models,
                                   proj.name = name.folder,
                                   models.options = models.options)
              proj.stack <- stack(proj.stack)
              names(proj.stack) <- models.name
              
              # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
              if (length(binary.proj) > 0) {
                cat("\nBinary transformations...")
                lapply(binary.proj, function(bin.proj)
                {
                  name.obj = paste0(name.sp, "_bin_", bin.proj, "_RasterStack")
                  name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.obj)
                  .saveBinary(i.proj = filt.proj,
                              names.proj = names(proj.stack),
                              obj.eval = models.evaluation,
                              obj.data = proj.stack,
                              doFiltering = FALSE,
                              doStack = TRUE,
                              name.obj = name.obj,
                              name.file = name.file)
                })
              }
              
              # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
              if (length(filtred.proj) > 0) {
                cat("\nFiltered transformations...")
                lapply(filtred.proj, function(filt.proj)
                {
                  name.obj = paste0(name.sp, "_filt_", bin.proj, "_RasterStack")
                  name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.obj)
                  .saveBinary(i.proj = filt.proj,
                              names.proj = names(proj.stack),
                              obj.eval = models.evaluation,
                              obj.data = proj.stack,
                              doFiltering = TRUE,
                              doStack = TRUE,
                              name.obj = name.obj,
                              name.file = name.file)
                })
              }
              
              # 7. Saving projection on hard disk
              assign(x = paste0(name.sp, "_RasterStack"), value = proj.stack)
              
              name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.sp, "_RasterStack")
              eval(parse(text = paste0("save(", name.sp, "_RasterStack", ", file = name.file)")))
              gc(reset = TRUE)
            } else {
              
              # all models will be saved separatly
              proj.stack <- c() # list of saved files
              for (m.n in models.name)
              {
                proj.ras <- .Projection.do.proj(m.n,
                                                env = new.env.data,
                                                scaled.models = scaled.models,
                                                proj.name = name.folder,
                                                models.options = models.options)
                names(proj.ras) <- m.n
                
                # 5. Computing Binary transformation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
                if (length(binary.proj) > 0) {
                  cat("\nBinary transformations...")
                  lapply(binary.proj, function(bin.proj)
                  {
                    name.obj = paste0(name.sp, "_bin_", bin.proj, "_RasterLayer")
                    name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.obj)
                    .saveBinary(i.proj = filt.proj,
                                names.proj = names(proj.ras),
                                obj.eval = models.evaluation,
                                obj.data = proj.ras,
                                doFiltering = FALSE,
                                doStack = TRUE,
                                name.obj = name.obj,
                                name.file = name.file)
                  })
                }
                
                # 6. Computing Filtering transformation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
                if (length(filtred.proj) > 0) {
                  cat("\nFiltered transformations...")
                  lapply(filtred.proj, function(filt.proj)
                  {
                    name.obj = paste0(name.sp, "_filt_", bin.proj, "_RasterLayer")
                    name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.obj)
                    .saveBinary(i.proj = filt.proj,
                                names.proj = names(proj.ras),
                                obj.eval = models.evaluation,
                                obj.data = proj.ras,
                                doFiltering = TRUE,
                                doStack = TRUE,
                                name.obj = name.obj,
                                name.file = name.file)
                  })
                }
                
                # 7. Saving projection on hard disk
                assign(x = paste0(name.sp, "_RasterLayer"), value = proj.ras)
                
                name.file = paste0(modeling.work.dir, "/", sp.name, "/", name.folder, "/", name.sp, "_RasterLayer")
                eval(parse(text = paste0("save(", name.sp, "_RasterLayer", ", file = name.file)")))
                gc(reset = TRUE)
                
                proj.stack <- c(proj.stack, paste0(proj.name, "_", m.n, "_RasterLayer"))
              }
            }
            
            ## remove MAXENT.Phillips tmp dir if exists
            if (file.exists(file.path(sp.name, proj.name, 'MaxentTmpData'))) {
              bm_MAXENTdeleteWorkdir(file.path(sp.name, proj.name))
            }
            
            return(invisible(proj.stack))
          })


setGeneric(".Projection.do.proj",
           def = function(model.name, env, model.dir = NULL, ...) {
             standardGeneric(".Projection.do.proj")
           })

setMethod('.Projection.do.proj', signature(env = 'data.frame'),
          function(model.name, env, xy = NULL, model.dir = NULL
                   , scaled.models = TRUE, proj.name = NULL, models.options = NULL)
          {
            cat('\n\t>', model.name)
            
            if (is.null(model.dir)) {
              model.dir <- paste0(getwd(), '/',
                                  .extractModelNamesInfo(model.name, info = 'species'),
                                  '/models')
            }
            
            # loading model
            if(length(c(grep('SRE', model.name))) == 0) {
              model.sp = get(load(paste0(model.dir, '/', model.name)))
              rm(list = model.name)
            }
            
            # check model.type
            model.type <- tail(unlist(strsplit(model.name, split = "_")), 1)
            avail.models.list <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE'
                                   , 'FDA', 'MARS', 'RF', 'MAXENT.Phillips')
            if (!(model.type %in% avail.models.list) && !grep('EF.',model.type)) {
              stop('Unknown model type')
            }
            
            proj = switch(model.type
                          , 'ANN' = {
                            set.seed(555) # to be able to refind our trees MAY BE BAD
                            return(as.numeric(predict(model.sp, env, type = "raw")))
                          }
                          , 'CTA' = {
                            return(as.integer(as.numeric(predict(model.sp, env, type = "prob")[, 2]) * 1000))
                          }
                          , 'FDA' = {
                            return(as.numeric(predict(model.sp, env, type = "posterior")[, 2]))
                          }
                          , 'GBM' = {
                            best.iter = gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
                            return(as.integer(predict.gbm(model.sp, env, best.iter, type = "response") * 1000))
                          }
                          , 'GLM' = {
                            return(as.integer(.run_pred(model.sp, Prev = 0.5, env) * 1000))
                          }
                          , 'GAM' = {
                            return(as.integer(.run_pred(model.sp, Prev = 0.5, env) * 1000))
                          }
                          , 'RF' = {
                            return(as.integer(as.numeric(predict(model.sp, env, type = 'prob')[, '1']) * 1000))
                          }
                          , 'SRE' = {
                            # loading data of the corresponding run
                            load(paste0(model.dir, '/', model.name, '/Data_', model.name))
                            return(eval(parse(text = paste0("bm_SRE(Data_", model.name, "$Response, Data_",
                                                            model.name, "$Explanatory, env, Data_",
                                                            model.name, "$Quant) * 1000"))))
                          }
                          , 'MAXENT.Phillips' = {
                            if (!is.null(xy)) {
                              return(as.integer(predict(object = model.sp, newdata = env, proj_name = proj.name, xy = xy) * 1000))
                            } else {
                              cat('\n MAXENT.Phillips need coordinates to run! NA returned ')
                              return(data.frame(rep(NA, nrow(env))))
                            }
                          }
            )
            
            if (model.name != "SRE" ||
                model.name != "MAXENT.Phillips" ||
                (model.name == "MAXENT.Phillips" && !is.null(xy))) {
              if (model.name %in% c('ANN', 'FDA')) {
                proj = as.integer(bm_Rescaler(proj, name = model.name) * 1000)
              }
              if (scaled.models && model.name %in% c('CTA', 'GBM', 'GLM', 'GAM', 'RF', 'MAXENT.Phillips')) {
                proj = as.integer(bm_Rescaler(proj / 1000, name = model.name) * 1000)
              }
              return(data.frame(proj = proj))
            } else {
              return(proj)
            }
          })


setMethod('.Projection.do.proj', signature(env = 'RasterStack'),
          function(model.name, env, model.dir = NULL
                   , scaled.models = TRUE, proj.name = NULL, models.options = NULL)
          {
            cat('\n\t>', model.name)
            
            if (is.null(model.dir)) {
              model.dir <- paste0(getwd(), '/',
                                  .extractModelNamesInfo(model.name, info = 'species'),
                                  '/models')
            }
            
            # loading model
            if(length(c(grep('SRE', model.name))) == 0) {
              model.sp = get(load(paste0(model.dir, '/', model.name)))
              rm(list = model.name)
            }
            
            # check model.type
            model.type <- tail(unlist(strsplit(model.name, split = "_")), 1)
            avail.models.list <- c('GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE'
                                   , 'FDA', 'MARS', 'RF', 'MAXENT.Phillips')
            if (!(model.type %in% avail.models.list) && !grep('EF.',model.type)) {
              stop('Unknown model type')
            }
            
            proj.ras = switch(model.type
                              , 'ANN' = {
                                set.seed(555) # to be able to refind our trees MAY BE BAD
                                return(predict(env, model.sp, type = "raw"))
                              }
                              , 'CTA' = {
                                return(predict(env, model.sp, type = 'prob', index = 2))
                              }
                              , 'FDA' = {
                                return(predict(env, model.sp, type = 'post', index = 2))
                              }
                              , 'GBM' = {
                                if (file.exists(paste(model.dir, '/', model.name, '_best.iter'))) {
                                  load(paste(model.dir, '/', model.name, '_best.iter'))
                                } else {
                                  best.iter <- gbm.perf(model.sp, method = "cv", plot.it = FALSE) # may be better to load it
                                }
                                return(predict(env, model.sp, n.trees = best.iter, type = 'response'))
                              }
                              , 'GLM' = {
                                return(predict(env, model.sp, type = 'response'))
                              }
                              , 'GAM' = {
                                return(predict(env, model.sp, type = 'response'))
                              }
                              , 'MARS' = {
                                return(predict(env, model.sp))
                              }
                              , 'RF' = {
                                return(predict(env, model.sp, type = 'prob', index = 2))
                              }
                              , 'SRE' = {
                                # loading data of the corresponding run
                                data.sre = get(load(paste0(model.dir, '/', model.name, '/Data_', model.name)))
                                return(raster::subset(bm_SRE(data.sre$Response, data.sre$Explanatory, env, data.sre$Quant), 1, drop = TRUE) * 1000)
                              }
                              , 'MAXENT.Phillips' = {
                                return(predict(object = model.sp, newdata = env, proj_name = proj.name))
                              }
            )
            
            if (model.name != "SRE") {
              if (model.name %in% c('ANN', 'FDA', 'MARS') ||
                  (model.name %in% c('CTA', 'GBM', 'GLM', 'GAM', 'RF') && scaled.models)) {
                proj.ras[!is.na(proj.ras[])] <- bm_Rescaler(proj.ras[!is.na(proj.ras[])]
                                                            , ref = NULL, name = model.name, original = FALSE)
              }
              return(round(proj.ras * 1000))
            } else {
              return(proj.ras)
            }
          })

.Biomod.Models.loop <- function(X,
                                modeling.id,
                                Model,
                                Options,
                                VarImport,
                                mod.eval.method,
                                SavePred,
                                xy=NULL,
                                scal.models = TRUE){
  cat("\n\n-=-=-=- Run : ",X$name, '\n')
  res.sp.run <- list()

  for(i in 1:ncol(X$calibLines)){ # loop on RunEval
    cat('\n\n-=-=-=--=-=-=-',paste(X$name,dimnames(X$calibLines)[[2]][i],sep=""),'\n')

    res.sp.run[[dimnames(X$calibLines)[[2]][i]]] <- lapply(Model, .Biomod.Models,
                                                      Data = X$dataBM,
                                                      Options = Options,
                                                      calibLines = na.omit(X$calibLines[,i,]), ## transform 3D calibLines obj into a 1D vector
                                                      Yweights = na.omit(X$Yweights),
                                                      nam = paste(X$name,dimnames(X$calibLines)[[2]][i], sep=""),
                                                      VarImport = VarImport,
                                                      mod.eval.method = mod.eval.method,
                                                      evalData = X$evalDataBM,
                                                      SavePred = T,#SavePred,
                                                      xy = X$xy,
                                                      eval.xy = X$eval.xy,
                                                      scal.models = scal.models,
                                                      modeling.id = modeling.id)

    names(res.sp.run[[dimnames(X$calibLines)[[2]][i]]]) <- Model

  }

  return(res.sp.run)
}


.Biomod.Models <- function (Model, Data, Options, calibLines, Yweights, nam, VarImport = 0,
                            mod.eval.method = c('ROC','TSS','KAPPA'), evalData = NULL,
                            SavePred = FALSE,
                            xy = NULL, eval.xy = NULL, scal.models = TRUE, modeling.id = ''){

  ################################################################################################
  # 1. Print model running and getting model options =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  # check and get modified args if nececary
  args <- .Biomod.Models.check(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData, scal.models)

  if(is.null(args)){ # trouble in input data -> Not Run
    return(0)
  } else {
    Data <- args$Data
    Yweights <- args$Yweights
    evalLines <- args$evalLines
    Type <- args$Type
    criteria <- args$criteria
    Prev <- args$Prev
    mod.eval.method <- args$mod.eval.method
    evalData <- args$evalData
    scal.models <- args$scal.models
    resp_name <- args$resp_name
    expl_var_names <- args$expl_var_names
    compress.arg <- TRUE # ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  }

  categorial_var <- unlist(sapply(expl_var_names, function(x){if(is.factor(Data[,x])) return(x) else return(NULL)} ))

  model_name <- paste(nam,'_',Model,sep="")


  # defining the function outputs
  ListOut <- list(evaluation = NULL,
                  var.import = NULL,
                  pred = NULL,
                  pred.eval = NULL,
                  calib.failure = NULL)


  # CTA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "CTA") {

    # converting cost argument
    if(is.null(Options@CTA$cost)){
      cost.tmp <- rep(1,(ncol(Data)-2))
    } else{
      cost.tmp <- Options@CTA$cost
    }
    if(Options@CTA$parms == 'default'){
      model.sp <- try( rpart(makeFormula(colnames(Data)[1],
                                         head(Data[,-c(1,ncol(Data)), drop=FALSE]),
                                         'simple', 0),
                             data = Data[calibLines,],
                             weights = Yweights,
                             method = Options@CTA$method,
                             cost = cost.tmp,
                             control = eval(Options@CTA$control)) )
    } else{
      model.sp <- try( rpart(makeFormula(colnames(Data)[1],
                                         head(Data[,-c(1,ncol(Data)), drop=FALSE]),
                                         'simple', 0),
                             data = Data[calibLines,],
                             weights = Yweights,
                             method = Options@CTA$method,
                             parms = Options@CTA$parms,
                             cost = cost.tmp,
                             control = eval(Options@CTA$control)) )
    }



    if( !inherits(model.sp,"try-error") ){
      # select best trees --------------- May be done otherway
      tr <- as.data.frame(model.sp$cptable)
      tr$xsum <- tr$xerror + tr$xstd
      tr <- tr[tr$nsplit > 0, ]
      Cp <- tr[tr$xsum == min(tr$xsum), "CP"]

      model.sp <- prune(model.sp, cp = Cp[length(Cp)])

      # creation of biomod2 model object
      model.bm <- new("CTA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'CTA',
                      model_options = Options@CTA,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end CTA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #




  # GAM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "GAM"){

    # NOTE : To be able to take into account GAM options and weights we have to do a eval(parse(...))
    # it's due to GAM implementation ( using of match.call() troubles)

#     cat("\n\tUser defined control args building..")
#     user.control.list <- Options@GAM$control
#
#     if(Options@GAM$algo == 'GAM_gam'){
#       default.control.list <- gam::gam.control()
#     } else{
#       default.control.list <- mgcv::gam.control()
#     }
#
#     control.list <- lapply(names(default.control.list), function(x){
#       if(x %in% names(user.control.list)){
#         return(user.control.list[[x]])
#       } else {
#         return(default.control.list[[x]])
#       }
#     })
#     names(control.list) <- names(default.control.list)

    ### Old version
    if(Options@GAM$algo == 'GAM_gam'){ ## gam package
      # package loading
      if(isNamespaceLoaded("mgcv")){
        if(isNamespaceLoaded("caret")){unloadNamespace("caret")} ## need to unload caret before car
        if(isNamespaceLoaded("car")){unloadNamespace("car")} ## need to unload car before mgcv
        unloadNamespace("mgcv")
      }
      # if(!isNamespaceLoaded("gam")){requireNamespace("gam", quietly = TRUE)}
      requireNamespace("gam", quietly = TRUE)

      cat('\n\t> GAM (gam) modelling...')

      gamStart <- eval(parse(text=paste("gam::gam(",colnames(Data)[1] ,"~1 ," ,
                                        " data = Data[calibLines,,drop=FALSE], family = ", Options@GAM$family$family,"(link = '",Options@GAM$family$link,"')",#eval(Options@GAM$family),
                                        ", weights = Yweights[calibLines])" ,sep="")))
      model.sp <- try( gam::step.Gam(gamStart, .scope(Data[1:3,-c(1,ncol(Data))], "gam::s", Options@GAM$k),
                                     data = Data[calibLines,,drop=FALSE],
      #                                keep = .functionkeep,
                                     direction = "both",
                                     trace= Options@GAM$control$trace,
                                     control = Options@GAM$control))#eval(control.list)) )
    } else { ## mgcv package
      # package loading
#       if( ("package:gam" %in% search()) ){ detach("package:gam", unload=TRUE)}
#       if( ! ("package:mgcv" %in% search()) ){ require("mgcv",quietly=TRUE) }
      if(isNamespaceLoaded("gam")){unloadNamespace("gam")}
      if(!isNamespaceLoaded("mgcv")){requireNamespace("mgcv", quietly = TRUE)}


      if(is.null(Options@GAM$myFormula)){
        cat("\n\tAutomatic formula generation...")
        gam.formula <- makeFormula(resp_name,head(Data[,expl_var_names,drop=FALSE]),Options@GAM$type, Options@GAM$interaction.level, k=Options@GAM$k)
      } else{
        gam.formula <- Options@GAM$myFormula
      }

      if (Options@GAM$algo == 'GAM_mgcv'){
        cat('\n\t> GAM (mgcv) modelling...')
        model.sp <- try(mgcv::gam(gam.formula,
                                   data= Data[calibLines,,drop=FALSE],
                                   family= Options@GAM$family,
                                   weights = Yweights,
                                   control = Options@GAM$control))

      } else if (Options@GAM$algo == 'BAM_mgcv'){ ## big data.frame gam version
        cat('\n\t> BAM (mgcv) modelling...')
        model.sp <- try(mgcv::bam(gam.formula,
                                   data=Data[calibLines,,drop=FALSE],
                                   family=Options@GAM$family,
                                   weights = Yweights,
                                   control = Options@GAM$control))
      }
    }


    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("GAM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GAM',
                      model_subclass = Options@GAM$algo,
                      model_options = Options@GAM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end GAM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # GBM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "GBM") {

    model.sp <- try(gbm(formula = makeFormula(colnames(Data)[1],head(Data)[,expl_var_names,drop=FALSE], 'simple',0),
                        data = Data[calibLines,,drop=FALSE],
                        distribution = Options@GBM$distribution,
                        var.monotone = rep(0, length = ncol(Data)-2), # -2 because of removing of sp and weights
                        weights = Yweights,
                        interaction.depth = Options@GBM$interaction.depth,
                        n.minobsinnode = Options@GBM$n.minobsinnode,
                        shrinkage = Options@GBM$shrinkage,
                        bag.fraction = Options@GBM$bag.fraction,
                        train.fraction = Options@GBM$train.fraction,
                        n.trees = Options@GBM$n.trees,
                        verbose = Options@GBM$verbose,
                        #class.stratify.cv = Options@GBM$class.stratify.cv,
                        cv.folds = Options@GBM$cv.folds, 
                        n.cores = Options@GBM$n.cores)) ## to prevent from parallel issues

    if( !inherits(model.sp,"try-error") ){
      best.iter <- try(gbm.perf(model.sp, method = Options@GBM$perf.method , plot.it = FALSE))

      model.bm <- new("GBM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GBM',
                      n.trees_optim = best.iter,
                      model_options = Options@GBM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))

    }
  }
  # end GBM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # GLM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "GLM"){

    ## build the most complete model formula
    if(is.null(Options@GLM$myFormula)){
      glm.formula <- makeFormula(colnames(Data)[1],head(Data),Options@GLM$type, Options@GLM$interaction.level)
    } else{
      glm.formula <- Options@GLM$myFormula
    }

    if(Options@GLM$test != 'none'){
      ## make the model selection
      glmStart <- glm(eval(parse(text=paste(colnames(Data)[1],"~1",sep=""))),
                      data = Data[calibLines,,drop=FALSE],
                      family = Options@GLM$family,
                      control = eval(Options@GLM$control),
                      weights = Yweights[calibLines],
                      mustart = rep(Options@GLM$mustart, sum(calibLines)),
                      model = TRUE)

      ## remove warnings
      warn <- options('warn')
      options(warn=-1)
      model.sp <- try( stepAIC(glmStart,
                               glm.formula,
                               data = Data[calibLines,,drop=FALSE],
                               direction = "both", trace = FALSE,
                               k = criteria,
                               weights = Yweights[calibLines],
                               steps = 10000,
                               mustart = rep(Options@GLM$mustart, sum(calibLines))) )

      ## reexec warnings
      options(warn)

    } else {
      ## keep the total model
      model.sp <- try( glm(glm.formula,
                           data = cbind(Data[calibLines,,drop=FALSE],matrix(Yweights[calibLines], ncol=1, dimnames=list(NULL, "Yweights"))),
                           family = Options@GLM$family,
                           control = eval(Options@GLM$control),
                           weights = Yweights,
#                            mustart = rep(Options@GLM$mustart, sum(calibLines)),
                           model = TRUE) )
    }

    if( !inherits(model.sp,"try-error") ){
      # print the selected formula
      cat("\n\tselected formula : ")
      print(model.sp$formula, useSource=FALSE)



      model.bm <- new("GLM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GLM',
                      model_options = Options@GLM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end GLM models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



#   # MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#   if (Model == "MARS"){
#     ## deal with nk argument
#     ## if not defined, it will be setted up to default mars value i.e max(21, 2 * ncol(x) + 1)
#     nk <- Options@MARS$nk
#     if(is.null(nk)){
#       nk <- max(21, 2 * length(expl_var_names) + 1)
#     }
#
#     model.sp <- try( mars(x = Data[calibLines,expl_var_names,drop=FALSE],
#                           y = Data[calibLines,1],
#                           degree = Options@MARS$degree,
#                           nk = nk,
#                           penalty = Options@MARS$penalty,
#                           thresh = Options@MARS$thresh,
#                           prune = Options@MARS$prune,
#                           w = Yweights[calibLines]) )
#
#     if( !inherits(model.sp,"try-error") ){
#
#       model.bm <- new("MARS_biomod2_model",
#                       model = model.sp,
#                       model_name = model_name,
#                       model_class = 'MARS',
#                       model_options = Options@MARS,
#                       resp_name = resp_name,
#                       expl_var_names = expl_var_names,
#                       expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
#                       expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
#     }
#   }
#   # end MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  # MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "MARS"){

    ## build the most complete model formula
    if(is.null(Options@MARS$myFormula)){
      mars.formula <- makeFormula(colnames(Data)[1],head(Data)[, -ncol(Data), drop = FALSE],Options@MARS$type, Options@MARS$interaction.level)
    } else{
      mars.formula <- Options@MARS$myFormula
    }

    ## deal with nk argument
    ## if not defined, it will be setted up to default mars value i.e max(21, 2 * ncol(x) + 1)
    nk <- Options@MARS$nk
    if(is.null(nk)){
      # nk <- max(21, 2 * length(expl_var_names) + 1)
      nk <- min(200, max(20, 2 * length(expl_var_names))) + 1
    }

    model.sp <- try(earth(formula = mars.formula,
                          data = Data[calibLines, , drop=FALSE],
                          weights = Yweights,
                          glm = list(family = binomial),
                          ncross = 0,
                          keepxy = FALSE,
                          # degree = Options@MARS$degree,
                          pmethod = Options@MARS$pmethod,
                          nprune = Options@MARS$nprune,
                          nk = nk,
                          penalty = Options@MARS$penalty,
                          thresh = Options@MARS$thresh))

    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("MARS_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MARS',
                      model_options = Options@MARS,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end MARS models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # FDA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "FDA") {
    model.sp <- try( do.call(fda,
                             c( list( formula = makeFormula(colnames(Data)[1],head(Data)[,expl_var_names,drop=FALSE], 'simple',0),
                                      data = Data[calibLines,,drop=FALSE],
                                      method = eval(parse(text=call(Options@FDA$method))),
                                      weights = Yweights[calibLines] ),
                                Options@FDA$add_args) ) )

    if( !inherits(model.sp,"try-error") ){

      model.bm <- new("FDA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'FDA',
                      model_options = Options@FDA,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end FDA models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # ANN models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "ANN") {
    size = Options@ANN$size
    decay = Options@ANN$decay

    if(is.null(size) | is.null(decay) | length(size)>1 | length(decay)>1 ){
      ## define the size and decay to test
      if(is.null(size)) size <- c(2, 4, 6, 8)
      if(is.null(decay)) decay <- c(0.001, 0.01, 0.05, 0.1)

      ## do cross validation test to find the optimal values of size and decay parameters (prevent from overfitting)
      CV_nnet <- 
        .CV.nnet(
          Input = Data[,expl_var_names,drop=FALSE],
          Target = Data[calibLines,1],
          size = size,
          decay = decay,
          maxit = Options@ANN$maxit,
          nbCV = Options@ANN$NbCV,
          W = Yweights[calibLines]
        )

      ## get the optimised parameters values
      decay <- CV_nnet[1, 2]
      size <- CV_nnet[1,1]
    }

    model.sp <- try(nnet(formula = makeFormula(resp_name,head(Data[,expl_var_names,drop=FALSE]), 'simple',0),
                         data = Data[calibLines,,drop=FALSE],
                         size = size,
                         rang = Options@ANN$rang,
                         decay = decay,
                         weights=Yweights,
                         maxit = Options@ANN$maxit,
                         trace = FALSE))

    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("ANN_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'ANN',
                      model_options = Options@ANN,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # ANN models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # RF models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "RF") {
    if(Options@RF$do.classif){
      # defining occurences as factor for doing classification and not regression in RF
      Data <- Data %>% mutate_at(resp_name, factor)
    }

    if(Options@RF$mtry == 'default'){
      model.sp <- try(randomForest(formula = makeFormula(resp_name,head(Data), 'simple',0),
                                   data = Data[calibLines,],
                                   ntree = Options@RF$ntree,
                                   #mtry = ifelse(Options@RF$ntree == 'default', round((ncol(Data)-1)/2), Options@RF$ntree ),
                                   importance = FALSE,
                                   norm.votes = TRUE,
                                   strata = factor(c(0,1)),
                                   nodesize = Options@RF$nodesize,
                                   maxnodes = Options@RF$maxnodes) )
    } else {
      model.sp <- try(randomForest(formula = makeFormula(resp_name,head(Data), 'simple',0),
                                   data = Data[calibLines,],
                                   ntree = Options@RF$ntree,
                                   mtry = Options@RF$mtry,
                                   importance = FALSE,
                                   norm.votes = TRUE,
                                   strata = factor(c(0,1)),
                                   nodesize = Options@RF$nodesize,
                                   maxnodes = Options@RF$maxnodes) )
    }


    if(Options@RF$do.classif){
      # canceling occurences class modifications
      Data <- Data %>% mutate_at(resp_name, function(.x) .x %>% as.character() %>% as.numeric())
    }

    if( !inherits(model.sp,"try-error") ){

      model.bm <- new("RF_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'RF',
                      model_options = Options@RF,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end RF models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #





  # SRE models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "SRE"){
    model.sp <- try(sre(Response = Data[calibLines,1],
                        Explanatory = Data[calibLines,expl_var_names,drop=FALSE],
                        NewData = NULL,
                        Quant = Options@SRE$quant,
                        return_extremcond=TRUE))

    if( !inherits(model.sp,"try-error") ){
      model.bm <- new("SRE_biomod2_model",
                      extremal_conditions = model.sp,
                      model_name = model_name,
                      model_class = 'SRE',
                      model_options = Options@SRE,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                      expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
    }
  }
  # end SRE models creation =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #





  # MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (Model == "MAXENT.Phillips"){
    MWD <- .Prepare.Maxent.WorkDir(Data, xy, calibLines, nam, VarImport = 0,
									evalData, eval.xy, species.name=resp_name,
									modeling.id = modeling.id,
                                   	background_data_dir = Options@MAXENT.Phillips$background_data_dir)

    # run MaxEnt:
    cat("\n Running Maxent...")
    maxent.cmd <- paste0("java ",
                        ifelse(is.null(Options@MAXENT.Phillips$memory_allocated),"",paste("-mx",Options@MAXENT.Phillips$memory_allocated,"m",sep="")),
                        " -jar ", file.path(Options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"),
                        " environmentallayers=\"", MWD$m_backgroundFile,
                        "\" samplesfile=\"", MWD$m_speciesFile,
                        "\" projectionlayers=\"", gsub(", ",",",toString(MWD$m_predictFile)),
                        "\" outputdirectory=\"", MWD$m_outdir, "\"",
                        " outputformat=logistic ",
                        #                            "jackknife maximumiterations=",Options@MAXENT.Phillips$maximumiterations,
                        ifelse(length(categorial_var),
                               paste(" togglelayertype=",categorial_var, collapse=" ",sep=""),
                               ""),
                        " redoifexists",
                        " visible=", Options@MAXENT.Phillips$visible,
                        " linear=", Options@MAXENT.Phillips$linear,
                        " quadratic=", Options@MAXENT.Phillips$quadratic,
                        " product=", Options@MAXENT.Phillips$product,
                        " threshold=", Options@MAXENT.Phillips$threshold,
                        " hinge=", Options@MAXENT.Phillips$hinge,
                        " lq2lqptthreshold=", Options@MAXENT.Phillips$lq2lqptthreshold,
                        " l2lqthreshold=", Options@MAXENT.Phillips$l2lqthreshold,
                        " hingethreshold=", Options@MAXENT.Phillips$hingethreshold,
                        " beta_threshold=", Options@MAXENT.Phillips$beta_threshold,
                        " beta_categorical=", Options@MAXENT.Phillips$beta_categorical,
                        " beta_lqp=", Options@MAXENT.Phillips$beta_lqp,
                        " beta_hinge=", Options@MAXENT.Phillips$beta_hinge,
                        " betamultiplier=", Options@MAXENT.Phillips$betamultiplier,
                        " defaultprevalence=", Options@MAXENT.Phillips$defaultprevalence,
                        " autorun nowarnings notooltips noaddsamplestobackground")

    system(command = maxent.cmd, wait = TRUE, intern = TRUE,
           ignore.stdout = FALSE, ignore.stderr = FALSE)

    model.bm <- new("MAXENT.Phillips_biomod2_model",
                    model_output_dir = MWD$m_outdir,
                    model_name = model_name,
                    model_class = 'MAXENT.Phillips',
                    model_options = Options@MAXENT.Phillips,
                    resp_name = resp_name,
                    expl_var_names = expl_var_names,
                    expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
                    expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))

    # for MAXENT.Phillips predicitons are calculated in the same time than models building to save time.
    cat("\n Getting predictions...")
    g.pred <- try(round(as.numeric(read.csv(MWD$m_outputFile)[,3]) * 1000))

    # remove tmp dir
    .Delete.Maxent.WorkDir(MWD)
  }
  # end MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-= #
  
  # MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(Model == "MAXENT.Phillips.2")
  {
    # browser()
    model.sp <- 
      try(
        maxnet::maxnet(
          p = Data %>% filter(calibLines) %>% pull(resp_name), 
          data = Data %>% filter(calibLines) %>% select_at(expl_var_names)
          # f = if(!is.null(Options@MAXENT.Phillips.2@))
        )
      )
    
    
    if( !inherits(model.sp,"try-error") )
    {
      model.bm <- 
        new(
          "MAXENT.Phillips.2_biomod2_model",
          model = model.sp,
          model_name = model_name,
          model_class = 'MAXENT.Phillips.2',
          model_options = Options@MAXENT.Phillips.2,
          resp_name = resp_name,
          expl_var_names = expl_var_names,
          expl_var_type = get_var_type(Data %>% filter(calibLines) %>% select_at(expl_var_names)),
          expl_var_range = get_var_range(Data %>% filter(calibLines) %>% select_at(expl_var_names))
        )
    }
  }
  # end MAXENT.Phillips models creation -=-=-=-=-=-=-=-=-=-=-=-=-= #

  # # MAXENT.Tsuruoka models creation -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  # if(Model == "MAXENT.Tsuruoka"){
  #   model.sp <- try(stop('MAXENT.Tsuruoka is depreacated(because maxent package is not maintained anymore)'))
  #   # model.sp <- try(maxent::maxent(feature_matrix = Data[calibLines, expl_var_names, drop = FALSE],
  #   #                                code_vector = as.factor(Data[calibLines, 1]),
  #   #                                l1_regularizer = Options@MAXENT.Tsuruoka$l1_regularizer,
  #   #                                l2_regularizer = Options@MAXENT.Tsuruoka$l2_regularizer,
  #   #                                use_sgd = Options@MAXENT.Tsuruoka$use_sgd,
  #   #                                set_heldout = Options@MAXENT.Tsuruoka$set_heldout,
  #   #                                verbose = Options@MAXENT.Tsuruoka$verbose))
  # 
  #   if( !inherits(model.sp,"try-error") ){
  #     model.bm <- new("MAXENT.Tsuruoka_biomod2_model",
  #                     model = model.sp,
  #                     model_name = model_name,
  #                     model_class = 'MAXENT.Tsuruoka',
  #                     model_options = Options@MAXENT.Tsuruoka,
  #                     resp_name = resp_name,
  #                     expl_var_names = expl_var_names,
  #                     expl_var_type = get_var_type(Data[calibLines,expl_var_names,drop=F]),
  #                     expl_var_range = get_var_range(Data[calibLines,expl_var_names,drop=F]))
  #   }
  # }
  # # end of MAXENT.Tsuruoka models creation -=-=-=-=-=-=-=-=-=-=-=- #

  # make prediction =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if((Model != "MAXENT.Phillips")){
    g.pred <- try(predict(model.bm, Data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE))
  }


  # scale or not predictions =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(scal.models & !inherits(g.pred,'try-error')){
    cat("\n\tModel scaling...")
    #     model.bm@scaling_model <- try(.scaling_model(g.pred/1000, Data[, 1])) ## without weigths
    model.bm@scaling_model <- try( .scaling_model(g.pred/1000, Data[, 1, drop = TRUE], weights= Yweights )) ## with weights
    g.pred <- try(predict(model.bm, Data[,expl_var_names,drop=FALSE], on_0_1000=TRUE))
  }

  # check predictions existance and stop execution if not ok -=-=- #
  test_pred_ok <- TRUE
  if (inherits(g.pred,"try-error")) { # model calibration or prdiction failed
    test_pred_ok <- FALSE
    cat("\n*** inherits(g.pred,'try-error')")
  } else if (sum(!is.na(g.pred))<=1){ # only NA predicted
    test_pred_ok <- FALSE
    cat("\n*** only NA predicted")
  } else if(length(unique(na.omit(g.pred))) <=1){ # single value predicted
    test_pred_ok <- FALSE
    cat("\n*** single value predicted")
  }

  if(test_pred_ok){
    # keep the model name
    ListOut$ModelName <- model_name
  } else{
    # keep the name of uncompleted modelisations
    cat("\n   ! Note : ", model_name, "failed!\n")
    ListOut$calib.failure = model_name
    return(ListOut) ## end of function.
  }

  # make prediction on evaluation data =-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(!is.null(evalData)){
    g.pred.eval <- try(predict(model.bm, evalData[,expl_var_names,drop=FALSE], on_0_1000=TRUE))
  }

  # save predictions -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if(SavePred){
    ListOut$pred <- g.pred
    if(exists("g.pred.eval"))
      ListOut$pred.eval <- g.pred.eval
  }


  # Model evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
  if(length(mod.eval.method) > 0){
    cat("\n\tEvaluating Model stuff...")

    ## Check no NA in g.pred to avoid evaluation failures
    na_cell_id <- which(is.na(g.pred))
    if(length(na_cell_id)){
#       g.pred.without.na <- g.pred[-na_cell_id]
      evalLines <- evalLines[!(evalLines %in% na_cell_id)]
      cat('\n\tNote : some NA occurs in predictions')
    } #else {
#       g.pred.without.na <- g.pred
#     }

    cross.validation <- 
      sapply(
        mod.eval.method,
        function(.x){
          Find.Optim.Stat(
            Stat = .x,
            Fit = g.pred[evalLines],
            Obs = Data %>% filter(evalLines) %>% pull(1)
          )
        }
      )

    rownames(cross.validation) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")

    if(exists('g.pred.eval')){

      ## Check no NA in g.pred to avoid evaluation failures
      na_cell_id <- which(is.na(g.pred.eval))
      if(length(na_cell_id)){
        g.pred.eval.without.na <- g.pred.eval[-na_cell_id]
        evalData <- evalData[-na_cell_id,]
        cat('\n\tNote : some NA occurs in evaluation predictions')
      } else {
        g.pred.eval.without.na <- g.pred.eval
      }

      true.evaluation <- sapply(mod.eval.method,
                                function(x){
                                  return( Find.Optim.Stat(Stat = x,
                                                          Fit = g.pred.eval.without.na,
                                                          Obs = evalData[,1],
                                                          Fixed.thresh = cross.validation["Cutoff",x]) )
                                })


      cross.validation <- rbind(cross.validation["Testing.data",], true.evaluation)

      rownames(cross.validation) <- c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity")
    }

    ListOut$evaluation <- t(round(cross.validation,digits=3))

    ## store results
    cross.validation <- t(round(cross.validation,digits=3))
    ListOut$evaluation <- cross.validation
    model.bm@model_evaluation <- cross.validation

    ## remove useless objects
    rm(list=c('cross.validation') )
  }
  # End evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #




  # Variables Importance -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  if (VarImport > 0){ # do Varimp stuff
    cat("\n\tEvaluating Predictor Contributions...", "\n")
    variables.importance <- variables_importance(model.bm, Data[, expl_var_names,drop=FALSE], nb_rand=VarImport)
    model.bm@model_variables_importance <- variables.importance$mat
    ## we stored only the mean of variables importance run
    ListOut$var.import <- round(rowMeans(variables.importance$mat, na.rm=T),digits=3)
    ## remove useless objects
    rm(list=c('variables.importance') )
  }
#   if (VarImport > 0){ # do Varimp stuff
#     # Create data.frame vhere corelation between predictions with and without permutation will be
#     # stored
#
#     cat("\n\tEvaluating Predictor Contributions...", "\n")
#     VarImpTable <- matrix(data = 0, nrow = VarImport, ncol = length(expl_var_names))
#     dimnames(VarImpTable) <- list(paste('rand', 1:VarImport, sep=""), expl_var_names)
#
#     for(vari in expl_var_names){
#       for (run in 1:VarImport) {
#         ## create a new dataset with interest variable suffled
#         TempDS <- Data[, expl_var_names,drop=FALSE]
#         TempDS[, vari] <- sample(TempDS[, vari])
#
#         if(Model != "MAXENT.Phillips"){
#           ## make projection on suffled dataset
#           shuffled.pred <- try(predict(model.bm, TempDS, on_0_1000=TRUE))
#         } else{
#           ## for MAXENT.Phillips, we have created all the permutation at model building step
#           shuffled.pred <- try(round(as.numeric(read.csv(file.path(model.bm@model_output_dir, paste(nam, vari, run, "swd.csv", sep="_")))[,3])*1000) )
#           ## scal suffled.pred if necessary
#           if(length(getScalingModel(model.bm))){
#             shuffled.pred <- try( round(.testnull(object = getScalingModel(model.bm), Prev = 0.5 , dat = data.frame(pred = shuffled.pred/1000) ) *1000) )
#             #               shuffled.pred <- round(as.numeric(predict(getScalingModel(model.bm), shuffled.pred/1000))*1000)
#           }
#           ## remove useless files on hard drive
#           file.remove(list.files(path=model.bm@model_output_dir,
#                                  pattern=paste(nam, vari, run, "swd", sep="_"),
#                                  full.names=TRUE))
#         }
#
#         ## test if differences exist between the 2 vectors
#         # check predictions existance and stop execution if not ok -=-=- #
#         test_shuffled.pred_ok <- TRUE
#         if (inherits(shuffled.pred,"try-error")) { # model calibration or prdiction failed
#           test_shuffled.pred_ok <- FALSE
#         } else if (sum(!is.na(shuffled.pred))<=1){ # only NA predicted
#           test_shuffled.pred_ok <- FALSE
#         } else if(length(unique(na.omit(shuffled.pred))) <=1){ # single value predicted
#           test_shuffled.pred_ok <- FALSE
#         } else if(length(shuffled.pred)!= length(g.pred)){
#           test_shuffled.pred_ok <- FALSE
#         }
#
#         if(!test_shuffled.pred_ok){
#           cat("\n   ! Note : ", model_name, "variable importance for",vari,run,"failed!\n")
#           VarImpTable[run,vari] <- 0
#         } else{
#           if(sum( g.pred != shuffled.pred, na.rm=T) == 0){
#             VarImpTable[run,vari] <- 0
#           } else {
#             ## calculate correlation between vectors as proxy for variables importance
#             VarImpTable[run,vari] <- 1 - max(round(cor(x=g.pred, y=shuffled.pred, use="pairwise.complete.obs", method="pearson"),digits=3),0,na.rm=T)
#           }
#         }
#
#       }
#     }
#
#     ## store results
#     model.bm@model_variables_importance <- VarImpTable
#     ## we stored only the mean of variables importance run
#     ListOut$var.import <- round(apply(VarImpTable, 2, mean, na.rm=T),digits=3)
#   }
  # End Variables Importance -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



  # Model saving step =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  assign(x=paste( nam, Model, sep = "_"),
         value= model.bm)
  save(list=paste( nam, Model, sep = "_"),
       file=file.path(resp_name, "models", modeling.id, paste( nam, Model, sep = "_")),
       compress=compress.arg)


  # End model saving step =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

  return(ListOut)
}


.Biomod.Models.check <- function(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData, scal.models, criteria=NULL, Prev=NULL){
  # get species and expanatory variables names
  resp_name <- colnames(Data)[1]
  expl_var_names <- colnames(Data)[-1]

  # replace Pseudo absences selected (NA) into true absences (0).. for model computing purpose

  if(sum(is.na(Data[,1])))
    Data[which(is.na(Data[,1])),1] <- 0

  # Calib Lines checking
  # & Test if there is absences AND presences in data given
  if (sum(!calibLines)>0){ # data are splited into 2 set one for evaluation and an other for evaluation stuff
    evalLines <- !calibLines
    if(sum(Data[calibLines,1] == 0 ) == 0 || sum(Data[calibLines,1] == 0 ) == sum(calibLines) ||
      sum(Data[evalLines,1] == 0) == 0 || sum(Data[evalLines,1] == 0) == sum(evalLines)){
      warning(paste(colnames(Data)[1], " ", Model," was switched off because of no both
                   presences and absences data given",sep=""), immediate.=T)
     return(NULL)
    }
  } else { # all values are taken for cali and valid stuff -----> Not so good but usefull for little data set
    evalLines <- calibLines
    if(sum(Data[,1] == 0 ) == 0 || sum(Data[,1] == 0 ) == nrow(Data)){
      warning(paste(colnames(Data)[1], " ", Model," was switched off because of no both
                   presences and absences data given (full model)",sep=""), immediate.=T)
      return(NULL)
    }
  }

  # weights checking
  if(is.null(Yweights)){
    Yweights <- rep(1,nrow(Data))
  }

  if(Model %in% c('GBM', 'CTA', 'ANN', 'FDA', 'GAM', 'MARS')){ # this models required data and weights to be in a same datdaset
    Data <- cbind(Data,Yweights)
  }

  # scaling parameter checking
  # never scal SRE
  if(Model == "SRE") scal.models <- FALSE
  # always scal ANN, FDA
  if(Model %in% c("ANN", "FDA") ) scal.models <- TRUE


  # models options checking and printing
  if (Model == "GLM"){
    cat('\nModel=GLM')
    if(!is.null(Options@GLM$myFormula)){
      cat('\n\tformula = ', paste(Options@GLM$myFormula[2],Options@GLM$myFormula[1],Options@GLM$myFormula[3]))
    } else{
      cat(' (',Options@GLM$type,'with', ifelse(Options@GLM$interaction.level == 0, 'no interaction )', paste('order',Options@GLM$interaction.level,'interaction level )')))
    }

    if(Options@GLM$test == "AIC"){
      criteria <- 2
      cat("\n\tStepwise procedure using AIC criteria")
    } else if(Options@GLM$test == "BIC"){
      criteria <- log(ncol(Data))
      cat("\n\tStepwise procedure using BIC criteria")
    } else if(Options@GLM$test == "none"){
      criteria <- 0
      cat("\n\tNo stepwise procedure")
      cat("\n\t! You might be confronted to models convergence issues !")
    }

  }

  if (Model == "GBM") {
    cat("\nModel=Generalised Boosting Regression \n")
    cat("\t", Options@GBM$n.trees, "maximum different trees and ", Options@GBM$cv.folds,
        " Fold Cross-Validation")
    set.seed(456) # to be able to refind our trees MAY BE BAD
  }

  if (Model == "GAM") {
    cat("\nModel=GAM")
    cat("\n\t",Options@GAM$algo,"algorithm chosen")
    #         cat("\t", Options@GAM$spline, " Degrees of smoothing")
    #         if(is.null(Yweights)) Yweights <- rep(1,nrow(Data))
    #         Data <- cbind(Data,Yweights)
  }

  if (Model == "CTA") {
    cat("\nModel=Classification tree \n")
    cat("\t", Options@CTA$control$xval, "Fold Cross-Validation")
    set.seed(123) # to be able to refind our trees MAY BE BAD
  }

  if (Model == "ANN") {
    cat("\nModel=Artificial Neural Network \n")
    cat("\t", Options@ANN$NbCV, "Fold Cross Validation + 3 Repetitions")
    #         cat("\tCalibration and evaluation phase: Nb of cross-validations: ",
    #             ncol(Ids), "\n")
    set.seed(555) # to be able to refind our trees MAY BE BAD
  }

  if (Model == "SRE")
    cat("\nModel=Surface Range Envelop")

  if (Model == "FDA"){
    cat("\nModel=Flexible Discriminant Analysis")
  }

  if (Model == "MARS"){
    cat("\nModel=Multiple Adaptive Regression Splines")
    if(!is.null(Options@MARS$myFormula)){
      cat('\n\tformula = ', paste(Options@MARS$myFormula[2],Options@MARS$myFormula[1],Options@MARS$myFormula[3]))
    } else{
      cat(' (',Options@MARS$type,'with', ifelse(Options@MARS$interaction.level == 0, 'no interaction )', paste('order',Options@MARS$interaction.level,'interaction level )')))
    }
    cat("\n")
  }

  if (Model == "RF"){
    cat("\nModel=Breiman and Cutler's random forests for classification and regression")
    set.seed(71)
  }

  if(Model == 'MAXENT.Phillips'){
    cat('\nModel=MAXENT.Phillips')
  }
  
  if(Model == 'MAXENT.Phillips.2'){
    cat('\nModel=MAXENT.Phillips (maxnet)')
  }

  # if(Model == 'MAXENT.Tsuruoka'){
  #   cat('\nModel=MAXENT.Tsuruoka')
  # }

  #     if (Model == "GLM" | Model == "GAM")
  #         Prev <- sum(DataBIOMOD[, i + Biomod.material$NbVar])/nrow(DataBIOMOD)
  ## not exactly same as before
  if (Model == "GLM" | Model == "GAM"){
    Prev <- sum(Data[,1], na.rm=T)/length(Data[,1])
  }

  # Evaluation Check
  available.eval.meth <- c('ROC','KAPPA','TSS','ACCURACY','BIAS','POD','FAR','POFD','SR','CSI',
                           'ETS','HK','HSS','OR','ORSS')

  #   if( Model %in% c('SRE') ) available.eval.meth <- available.eval.meth[which(available.eval.meth!='ROC')]
  if(sum(!(mod.eval.method %in% available.eval.meth)) > 0 ){
    warnings(paste(toString(mod.eval.method[!which(mod.eval.method %in% available.eval.meth)]),
                   ' were switched off !', sep='' ), imediate = TRUE)
  }
  mod.eval.method <- mod.eval.method[which(mod.eval.method %in% available.eval.meth)]

  return(list(Data=Data,
              Yweights=Yweights,
              evalLines=evalLines,
              criteria=criteria,
              Prev=Prev,
              mod.eval.method=mod.eval.method,
              evalData=evalData,
              scal.models=scal.models,
              resp_name=resp_name,
              expl_var_names=expl_var_names))

}

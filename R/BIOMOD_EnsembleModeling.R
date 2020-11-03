'BIOMOD_EnsembleModeling' <- function( modeling.output,
                                       chosen.models = 'all',
                                       em.by = 'PA_dataset+repet',
                                       eval.metric = 'all',
                                       eval.metric.quality.threshold = NULL,
                                       models.eval.meth = c('KAPPA','TSS','ROC'),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = FALSE,
                                       prob.mean.weight.decay = 'proportional',
                                       VarImport=0){
  .bmCat("Build Ensemble Models")
  # 1. args checking -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  args <- .BIOMOD_EnsembleModeling.check.args( modeling.output,
                                               chosen.models,
                                               eval.metric,
                                               eval.metric.quality.threshold,
                                               models.eval.meth,
                                               prob.mean,
                                               prob.cv,
                                               prob.ci,
                                               prob.ci.alpha,
                                               prob.median,
                                               committee.averaging,
                                               prob.mean.weight,
                                               prob.mean.weight.decay,
                                               em.by)

  modeling.output <- args$modeling.output
  chosen.models <- args$chosen.models
  eval.metric <- args$eval.metric
  eval.metric.quality.threshold <- args$eval.metric.quality.threshold
  models.eval.meth <- args$models.eval.meth
  prob.mean <- args$prob.mean
  prob.cv <- args$prob.cv
  prob.ci <- args$prob.ci
  prob.ci.alpha <- args$prob.ci.alpha
  prob.median <- args$prob.median
  committee.averaging <- args$committee.averaging
  prob.mean.weight <- args$prob.mean.weight
  prob.mean.weight.decay  <- args$prob.mean.weight.decay
  em.by <- args$em.by

  rm('args')

  em.avail <- c('prob.mean', 'prob.cv', 'prob.ci.inf', 'prob.ci.sup', 'prob.median', 'committee.averaging', 'prob.mean.weight')
  em.algo <- em.avail[c(prob.mean, prob.cv, prob.ci, prob.ci, prob.median, committee.averaging, prob.mean.weight)]

  # create a EM option list
  Options <- list(em.by=em.by)
  expl_var_type = get_var_type(get_formal_data(modeling.output,'expl.var'))
  expl_var_range = get_var_range(get_formal_data(modeling.output,'expl.var'))


  # 1b. creating output object and begin to fill it
#   EM <- list()
  EM <- new('BIOMOD.EnsembleModeling.out',
            sp.name = modeling.output@sp.name,
            expl.var.names = modeling.output@expl.var.names,
            em.by = em.by,
            modeling.id = modeling.output@modeling.id
#             models.out.obj = new('BIOMOD.stored.models.out',
#                                  inMemory = FALSE,
#                                  link = paste(modeling.output@sp.name,"/",modeling.output@sp.name,".models.out",sep="")),
#             eval.metric = eval.metric,
#             eval.metric.quality.threshold = eval.metric.quality.threshold#,
#             em.ci.alpha = prob.ci.alpha
            )

  EM@models.out.obj@link <- file.path(modeling.output@sp.name,paste(modeling.output@sp.name,".", modeling.output@modeling.id,".models.out",sep="") )

  # 2. doing Ensemble modeling

  ## 2.1 make a list of models names that will be combined together according to by argument.
  em.mod.assemb <- .em.models.assembling(chosen.models, em.by)

  for(assemb in names(em.mod.assemb) ){
    cat("\n\n  >", assemb, "ensemble modeling")
    models.kept <- em.mod.assemb[[assemb]]

    #### defined data that will be used for models performances calculation ####
    if(modeling.output@has.evaluation.data){
      eval.obs <- get_formal_data(modeling.output,'eval.resp.var')
      eval.expl <- get_formal_data(modeling.output,'eval.expl.var')
    }

    ##### !!!!!! TO DO -> select appropriate part of dataset according to em.by
    obs <-  get_formal_data(modeling.output,'resp.var')
    expl <- get_formal_data(modeling.output,'expl.var')

    ## subselection of observations according to dataset used to produce ensemble models
    if(em.by %in% c("PA_dataset",'PA_dataset+algo','PA_dataset+repet')){
      if(unlist(strsplit(assemb,"_"))[3] != 'AllData'){
        if(inherits(get_formal_data(modeling.output), "BIOMOD.formated.data.PA")){
          kept_cells <- get_formal_data(modeling.output)@PA[, unlist(strsplit(assemb,"_"))[3]]
        } else {
          kept_cells <- rep(T, length(obs))
        }

        obs <- obs[kept_cells]
        expl <- expl[kept_cells, ,drop=F]
      }
    }

    ## subselection of observations according to dataset used to produce ensemble models is done at evaluation step


#     if(em.by %in% c("algo","all") ){
#       ## we need to take all data even if it is much better to have
#       obs <-  get_formal_data(modeling.output,'resp.var')
#       expl <- get_formal_data(modeling.output,'expl.var')
#     }
    # remove na
    obs[is.na(obs)] <- 0




    #### get needed models predictions ####
    needed_predictions <- .get_needed_predictions(modeling.output, em.by, models.kept, eval.metric, eval.metric.quality.threshold)
    # if no prediction selected => swith to next model
    if(!length(needed_predictions)) next

    ## loop on evaluation metrics ##
    for(eval.m in eval.metric){

      # define model name
#       base_model_name <- paste(modeling.output@sp.name,"_",assemb,"_",eval.m,'_' ,sep="")
      models.kept <- needed_predictions$models.kept[[eval.m]]
      models.kept.scores <- needed_predictions$models.kept.scores[[eval.m]]

      ## Loop over em.algo ##

      for(algo in em.algo){
        #### Models building ####

        # 1. Mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.mean'){
          cat("\n   > Mean of probabilities...")

#           model_name <- paste(base_model_name,"EMmean",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMmeanBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMmean_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMmean',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id)
        }

        # 2. CV of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.cv'){
          cat("\n   > Coef of variation of probabilities...")
#           model_name <- paste(base_model_name,"EMcv",sep="")
          model_name <-paste(modeling.output@sp.name,"_","EMcvBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMcv_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMcv',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id)
        }

        # 3. Median of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.median'){
          cat("\n   > Median of probabilities...")
#           model_name <- paste(base_model_name,"EMmedian",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMmedianBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMmedian_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMmedian',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id)
        }

        # 4. CI of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
        if(algo == 'prob.ci.inf'){
          cat("\n   > Confidence Interval...")

          ## Quantile inferior
#           model_name <- paste(base_model_name,"EMciInf",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMciInfBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMci_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMci',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           alpha = prob.ci.alpha,
                           side = 'inferior')
        }

        if(algo == 'prob.ci.sup'){
          ## Quantile superior
#           model_name <- paste(base_model_name,"EMciSup",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMciSupBy",eval.m, "_", assemb ,sep="")

          model.bm <- new("EMci_biomod2_model",
                           model = models.kept,
                           model_name = model_name,
                           model_class = 'EMci',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           alpha = prob.ci.alpha,
                           side = 'superior')

        }

        # 5. Comitee averaging of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(algo == 'committee.averaging'){
          cat("\n   >  Committee averaging...")
#           model_name <- paste(base_model_name,"EMca",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMcaBy",eval.m, "_", assemb ,sep="")

          models.kept.tresh <- unlist(lapply(models.kept, function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(get_evaluations(modeling.output)[eval.m, "Cutoff", mod, run, dat])
          }))
          names(models.kept.tresh) <- models.kept

          ## remove models if some thresholds are undefined
          to_keep <- is.finite(models.kept.tresh)

          model.bm <- new("EMca_biomod2_model",
                           model = models.kept[to_keep],
                           model_name = model_name,
                           model_class = 'EMca',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           tresholds = models.kept.tresh[to_keep])

        }

        # 6. weighted mean of probabilities -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if(algo == 'prob.mean.weight'){
          cat("\n   > Probabilities weighting mean...")
#           model_name <- paste(base_model_name,"EMwmean",sep="")
          model_name <- paste(modeling.output@sp.name,"_","EMwmeanBy",eval.m, "_", assemb ,sep="")

          # remove SRE models if ROC
          models.kept.tmp <- models.kept
          models.kept.scores.tmp <- models.kept.scores
#           prediction.kept.tmp <- prediction.kept

          if(eval.m == 'ROC'){
            sre.id <- grep("_SRE", models.kept)
            if(length(sre.id)>0){
              cat("\n      ! SRE modeling were switched off")
              models.kept.tmp <- models.kept[-sre.id]
              models.kept.scores.tmp <- models.kept.scores[-sre.id]
#               prediction.kept.tmp <- prediction.kept[,models.kept]
            }
          }

          ## remove models if score is not defined
          models.kept.tmp <- models.kept.tmp[is.finite(models.kept.scores.tmp)]
          models.kept.scores.tmp <- models.kept.scores.tmp[is.finite(models.kept.scores.tmp)]

          # weights are "decay" times decreased for each subsequent model in model quality order.
          models.kept.scores.tmp <- round(models.kept.scores.tmp, 3) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.

          # dealing with numerical decay
          cat("\n\t\t", " original models scores = ", models.kept.scores.tmp)
          if(is.numeric(prob.mean.weight.decay)){
            DecayCount <- sum(models.kept.scores.tmp>0)
            WOrder <- order(models.kept.scores.tmp, decreasing=T)
            Dweights <- models.kept.scores.tmp
            ## old version
            # for(J in 1:DecayCount) Dweights[WOrder[J]] <- (DecayCount - J + 1) * prob.mean.weight.decay
            ## end old version
            for(J in 1:DecayCount) Dweights[WOrder[J]] <- I(prob.mean.weight.decay^(DecayCount - J + 1))
            #If 2 or more score are identical -> make a mean weight between the ones concerned
            for(J in 1:length(models.kept.scores.tmp)){
              if(sum(models.kept.scores.tmp[J]==models.kept.scores.tmp)>1) Dweights[which(models.kept.scores.tmp[J]==models.kept.scores.tmp)] <- mean(Dweights[which(models.kept.scores.tmp[J]==models.kept.scores.tmp)])
            }
            models.kept.scores.tmp <- round(Dweights, digits=3)
            rm(list=c('Dweights','DecayCount','WOrder'))
          } else if ( is.function(prob.mean.weight.decay) ){ # dealing with function decay
            models.kept.scores.tmp <- sapply(models.kept.scores.tmp, prob.mean.weight.decay)
          }

          ### Standardise model weights
          models.kept.scores.tmp <- round(models.kept.scores.tmp/sum(models.kept.scores.tmp, na.rm=T), digits=3)

          cat("\n\t\t", " final models weights = ", models.kept.scores.tmp)

          model.bm <- new("EMwmean_biomod2_model",
                           model = models.kept.tmp,
                           model_name = model_name,
                           model_class = 'EMwmean',
                           model_options = Options,
                           resp_name = modeling.output@sp.name,
                           expl_var_names = modeling.output@expl.var.names,
                           expl_var_type = expl_var_type,
                           expl_var_range = expl_var_range,
                           modeling.id = modeling.output@modeling.id,
                           penalization_scores = models.kept.scores.tmp)
        }

        #### Models Evaluation ####
        pred.bm <- predict(model.bm, expl, formal_predictions=needed_predictions$predictions[,model.bm@model, drop=F], on_0_1000 = T )

        ## store models prediction on the hard drive ---------------------------
        ## create the suitable directory architecture
        pred.bm.name <- paste0(model_name, ".predictions")
        pred.bm.outfile <- file.path(model.bm@resp_name, ".BIOMOD_DATA", model.bm@modeling.id,
                                     "ensemble.models", "ensemble.models.predictions",
                                     pred.bm.name)
        dir.create(dirname(pred.bm.outfile), showWarnings = FALSE, recursive = TRUE)
        ## save models predictions
        assign(pred.bm.name, pred.bm)
        save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
        rm(list = pred.bm.name)
        ## end strore models preciciton on the hard drive ----------------------

        if(exists('eval.obs') & exists('eval.expl')){
          eval_pred.bm <- predict(model.bm, eval.expl)
          ## store models prediction on the hard drive -------------------------
          pred.bm.name <- paste0(model_name, ".predictionsEval")
          ## save models predictions
          assign(pred.bm.name, eval_pred.bm)
          save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
          rm(list = pred.bm.name)
          ## end strore models preciciton on the hard drive --------------------
        }


        # Model evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
        if( length(models.eval.meth) ){
          cat("\n\t\t\tEvaluating Model stuff...")

          if(algo == 'prob.cv'){ ## switch of evalutaion process
            cross.validation <- matrix(NA,4,length(models.eval.meth),
                                       dimnames = list(c("Testing.data","Cutoff","Sensitivity", "Specificity"),
                                                       models.eval.meth))
          } else {
            if(em.by == "PA_dataset+repet"){ ## select the same evaluation data than formal models
              ## get formal models calib/eval lines
              calib_lines <- get_calib_lines(modeling.output)
              ## get info on wich dataset and which repet this ensemble model is based on
              pa_dataset_id <- paste("_", unlist(strsplit(assemb,"_"))[3], sep="")
              repet_id <- paste("_", unlist(strsplit(assemb,"_"))[2], sep="")
              ## define and extract the subset of points model will be evaluated on
              if (repet_id == "_Full"){
                eval_lines <- rep(T, length(pred.bm))
              } else {
                eval_lines <- ! na.omit(calib_lines[ , repet_id, pa_dataset_id])
                ## trick to detect when it is a full model but with a non common name
                if(all(!eval_lines)){ ## i.e. all lines used for calib => full model
                  eval_lines <- !eval_lines
                }
              }
            } else {
              eval_lines <- rep(T, length(pred.bm))
            }

            cross.validation <- sapply(models.eval.meth,
                                       Find.Optim.Stat,
                                       Fit = pred.bm[eval_lines],
                                       Obs = obs[eval_lines])
            rownames(cross.validation) <- c("Testing.data","Cutoff","Sensitivity", "Specificity")
          }



          if(exists('eval_pred.bm')){

            if(algo == 'prob.cv'){ ## switch of evalutaion process
              cross.validation <- matrix(NA,5,length(models.eval.meth),
                                         dimnames = list(c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity"),
                                                         models.eval.meth))
            } else {
              true.evaluation <- sapply(models.eval.meth, function(x){
                                        Find.Optim.Stat(
                                        Fit = eval_pred.bm * 1000,
                                        Obs = eval.obs,
                                        Fixed.thresh = cross.validation["Cutoff",x])})

              cross.validation <- rbind(cross.validation["Testing.data",], true.evaluation)
              rownames(cross.validation) <- c("Testing.data","Evaluating.data","Cutoff","Sensitivity", "Specificity")
            }
          }

          ## store results
          cross.validation <- t(round(cross.validation,digits=3))
          model.bm@model_evaluation <- cross.validation

          ## remove useless objects
          rm(list=c('cross.validation') )
        }
        # End evaluation stuff =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

        #### Var Importance calculation ####
        if (VarImport > 0){ # do Varimp stuff
          cat("\n\t\t\tEvaluating Predictor Contributions...", "\n")
          variables.importance <- variables_importance(model.bm, expl, nb_rand=VarImport)
          model.bm@model_variables_importance <- variables.importance$mat
          ## remove useless objects
          rm(list=c('variables.importance') )
        }

        #### Models saving #####
        assign(model_name,model.bm)
        save(list=model_name,file=file.path(modeling.output@sp.name,
                                            "models",
                                            modeling.output@modeling.id,
                                            model_name))

        #### Add to sumary objects ####
        EM@em.models <- c(EM@em.models, model.bm)
        EM@em.computed <- c(EM@em.computed, model_name)


      }
    }
  }

  ### fix models names ###
  names(EM@em.models) <- EM@em.computed

  model.name <- paste(EM@sp.name, '.', EM@modeling.id, 'ensemble.models.out', sep="")
  assign(x=model.name,
         value=EM)
  save(list=model.name,
       file=file.path(EM@sp.name,model.name))

  .bmCat("Done")
  return(EM)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

.BIOMOD_EnsembleModeling.check.args <- function(  modeling.output,
                                                   chosen.models,
                                                   eval.metric,
                                                   eval.metric.quality.threshold,
                                                   models.eval.meth,
                                                   prob.mean,
                                                   prob.cv,
                                                   prob.ci,
                                                   prob.ci.alpha,
                                                   prob.median,
                                                   committee.averaging,
                                                   prob.mean.weight,
                                                   prob.mean.weight.decay,
                                                   em.by ){
  # 1. modeling.output checking
  if(!(inherits(modeling.output, "BIOMOD.models.out"))){
    stop("Invalid modeling.output argument !\nIt must be a 'BIOMOD.models.out' object")
  }

  # 2. chosen.models checking
  if(!length(chosen.models) | (length(chosen.models)==1 & chosen.models[1] == 'all')){ # select all models
    cat("\n   ! all models available will be included in ensemble.modeling")
    chosen.models <- modeling.output@models.computed
  } else{
    chosen.models.check <- chosen.models %in% modeling.output@models.computed
    if(sum(!chosen.models.check) > 0){
      stop(paste("Some selected models do not exist: ", toString(chosen.models[!chosen.models.check]),
                 "\nPlease choose models in computed models ( ",
                 toString(modeling.output@models.computed), " )",sep=""))
    }
  }

  # 3. eval.metric checking
  if(!is.null(eval.metric)){
    if(!is.character(eval.metric)){
      stop("eval.metric must be a character vector or NULL")
    }
    if('all' %in% eval.metric){
      eval.metric <- dimnames(get_evaluations(modeling.output))[[1]]
    }
    eval.metric.check <- eval.metric %in% dimnames(get_evaluations(modeling.output))[[1]]
    if(sum(!eval.metric.check) > 0){
      stop(paste("Some selected evaluation metrics are not available: ", toString(eval.metric[!eval.metric.check]),
                 "\nPlease choose some in those computed yet ( ",
                 toString(dimnames(get_evaluations(modeling.output))[[1]]), " )",sep=""))
    }
  }

  # 4. eval.metric.quality.threshold
  if(!is.null(eval.metric)){
    if(!is.null(eval.metric.quality.threshold)){
      if(!is.numeric(eval.metric.quality.threshold)){
        stop("eval.metric.quality.threshold must be NULL or a numeric vector")
      }
      if(length(eval.metric) != length(eval.metric.quality.threshold)){
        stop("you must specify as many eval.metric.quality.threshold as eval.metric (if you specify some)")
      }
      cat("\n   > Evaluation & Weighting methods summary :\n")
      cat(paste(eval.metric, eval.metric.quality.threshold,  sep = " over ", collapse = "\n      "), fill=TRUE, labels = "     ")
    } else{
      cat("\n   ! No eval.metric.quality.threshold -> All models will be kept for Ensemble Modeling")
      eval.metric.quality.threshold <- rep(0, length(eval.metric))
    }
  }

  # 4b. model.eval.meth checking
  models.eval.meth <- unique(models.eval.meth)

  if(sum(models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS','KAPPA','ACCURACY','BIAS',
                                 'POD','PODFD','CSI','ETS','HK','ROC')) != length(models.eval.meth)){
    stop(paste(models.eval.meth[which( (models.eval.meth %in% c('FAR','SR','HSS','ORSS','TSS',
                                                                'KAPPA','ACCURACY','BIAS', 'POD',
                                                                'PODFD','CSI', 'ETS','HK','ROC'))
                                       == FALSE) ]," is not an availabe model evaluation metric !",sep=""))
  }

  # 5. check selected EM algo
  if( !is.logical(prob.mean) | !is.logical(prob.cv) | !is.logical(prob.ci) | !is.logical(prob.median) |
      !is.logical(committee.averaging) | !is.logical(prob.mean.weight) ){
    stop("prob.mean, prob.cv, prob.ci, prob.median, committee.averaging and prob.mean.weight arguments must be logical")
  }
  if(is.null(eval.metric)){
    if(committee.averaging | prob.mean.weight){
      stop("You must choose eval.metric if you want to compute Committee Averaging or Probability weighting mean algorithms")
    }
  }

  # 6. alpha for Confident interval
  if(prob.ci){
    if(!is.numeric(prob.ci.alpha)){
      stop("prob.ci.alpha must be numeric")
    }
    if(prob.ci.alpha <= 0 | prob.ci.alpha>= 0.5){
      stop("prob.ci.alpha must be a numeric between 0 and 0.5")
    }
  }

  # 7. decay checking
  if(prob.mean.weight){
    test.prob.mean.weight.decay <- TRUE
    ## check compatibility of prob.mean.weight.decay class
    if(!is.numeric(prob.mean.weight.decay) & !is.character(prob.mean.weight.decay) & !is.function(prob.mean.weight.decay)){
      test.prob.mean.weight.decay <- FALSE
    } else if(is.numeric(prob.mean.weight.decay)){ ## check numeric prob.mean.weight.decay
      if(prob.mean.weight.decay < 0){
        test.prob.mean.weight.decay <- FALSE
      }
    } else if(is.character(prob.mean.weight.decay)){ ## check character prob.mean.weight.decay
      if(prob.mean.weight.decay != 'proportional'){
        test.prob.mean.weight.decay <- FALSE
      }
    }

    if(!test.prob.mean.weight.decay){
      stop("'prob.mean.weight.decay' should be either 'proportional', a numeric value > 0 or a function")
    }
  }

  if(is.null(eval.metric)){
    eval.metric <- 'none'
  }

  # 8. by arg checking
  available.em.by <- c('PA_dataset', 'algo', 'all', 'PA_dataset+repet', 'PA_dataset+algo')
  if(!(em.by %in% available.em.by) ){
    stop("Invalid 'em.by' argument given. It must be one of : 'PA_dataset', 'algo', 'all', 'PA_dataset+repet' or 'PA_dataset+algo'")
  }


  return( list( modeling.output = modeling.output,
                chosen.models = chosen.models,
                eval.metric = eval.metric,
                eval.metric.quality.threshold = eval.metric.quality.threshold,
                models.eval.meth = models.eval.meth,
                prob.mean = prob.mean,
                prob.cv = prob.cv,
                prob.ci = prob.ci,
                prob.ci.alpha = prob.ci.alpha,
                prob.median = prob.median,
                committee.averaging = committee.averaging,
                prob.mean.weight = prob.mean.weight,
                prob.mean.weight.decay  = prob.mean.weight.decay,
                em.by = em.by))

}


# =-=-=-=-=-=-=-=- em.models.assembling function -=-=-=-=-=-=-=- #
.em.models.assembling <- function(chosen.models, em.by){
  assembl.list = list()

  if(em.by == 'PA_dataset'){
    for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
#       assembl.list[[paste(dat,"_AllRun", sep="")]] <- chosen.models[grep(paste("_",dat,"_",sep=""), chosen.models)]
      assembl.list[[paste("mergedAlgo_mergedRun_", dat, sep="")]] <- chosen.models[grep(paste("_",dat,"_",sep=""), chosen.models)]
    }
    return(assembl.list)
  }

  if(em.by == 'algo'){
    for(algo in .extractModelNamesInfo(chosen.models, info='models')){
#       assembl.list[[paste(algo,"_AllRun", sep="")]] <- chosen.models[grep(paste("_",algo,sep=""), chosen.models)]
      assembl.list[[paste(algo,"_mergedRun_mergedData", sep="")]] <- chosen.models[grep(paste("_",algo,sep=""), chosen.models)]
    }
    return(assembl.list)
  }

  if(em.by == 'all'){
#     assembl.list[["TotalConsensus"]] <- chosen.models
    assembl.list[[paste("mergedAlgo_mergedRun_mergedData", sep="")]] <- chosen.models
    return(assembl.list)
  }

  if(em.by == 'PA_dataset+repet'){
    for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
      for(repet in .extractModelNamesInfo(chosen.models, info='run.eval')){
        mod.tmp <- intersect(x=grep(paste("_",dat,"_",sep=""), chosen.models),
                             y=grep(paste("_",repet,"_",sep=""), chosen.models))
        if(length(mod.tmp)){
#           assembl.list[[paste(dat,"_",repet,'_AllAlgos', sep="")]] <- chosen.models[mod.tmp]
          assembl.list[[paste("mergedAlgo_",repet,"_",dat, sep="")]] <- chosen.models[mod.tmp]
        }
      }
    }
    return(assembl.list)
  }

  if(em.by == 'PA_dataset+algo'){
    for(dat in .extractModelNamesInfo(chosen.models, info='data.set')){
      for(algo in .extractModelNamesInfo(chosen.models, info='models')){
        mod.tmp <- intersect(x=grep(paste("_",dat,"_",sep=""), chosen.models),
                             y=grep(paste("_",algo,sep=""), chosen.models))
        if(length(mod.tmp)){
#           assembl.list[[paste(dat,"_AllRepet_",algo, sep="")]] <- chosen.models[mod.tmp]
          assembl.list[[paste(algo,"_mergedRun_", dat, sep="")]] <- chosen.models[mod.tmp]
        }
      }
    }
    return(assembl.list)
  }

}
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #



.get_needed_predictions <- function(modeling.output, em.by, models.kept, eval.metric, eval.metric.quality.threshold){
  out <- list(predictions = NULL,
              models.kept = NULL,
              models.kept.scores = NULL)
  for(eval.m in eval.metric){
    if( eval.m != 'none'){
      models.kept.scores <- unlist(lapply(models.kept, function(x){
        mod <- tail(unlist(strsplit(x,"_")), 3)[3]
        run <- tail(unlist(strsplit(x,"_")), 3)[2]
        dat <- tail(unlist(strsplit(x,"_")), 3)[1]
        # select evaluations scores obtained for Evaluation Data if exists or CV if not
        if(modeling.output@has.evaluation.data){
          return(get_evaluations(modeling.output)[eval.m, "Evaluating.data", mod, run, dat])
        } else{
          return(get_evaluations(modeling.output)[eval.m, "Testing.data", mod, run, dat])
        }

      }))
      ## set NA to -1
      if(!is.null(models.kept.scores)){
        models.kept.scores[is.na(models.kept.scores)] <- -1
      }
      out$models.kept[[eval.m]] <- models.kept[models.kept.scores > eval.metric.quality.threshold[which(eval.metric == eval.m)]]
      out$models.kept.scores[[eval.m]] <- models.kept.scores[models.kept.scores > eval.metric.quality.threshold[which(eval.metric == eval.m)]]
    } else{
      out$models.kept[[eval.m]] <- models.kept
    }
  }

  models.kept.union <- unique(unlist(out$models.kept))

  if(length(models.kept.union) ){
#     if(modeling.output@has.evaluation.data){
#       out$predictions <- as.data.frame(get_predictionsEval(modeling.output, as.data.frame = TRUE)[,models.kept.union, drop=F])
#     } else{
      ## load prediction on each PA dataset
      if(em.by %in% c("PA_dataset",'PA_dataset+algo','PA_dataset+repet')){
        out$predictions <- as.data.frame(get_predictions(modeling.output, as.data.frame = TRUE)[,models.kept.union, drop=F])
      } else{ ## redo prediction on full data.set
        cat("\n   ! Models projections for whole zonation required...")
        temp_name <- paste('tmp_',sub(".","",as.character(format(Sys.time(), "%OS6")), fixed=T),sep="")
        out$predictions <- BIOMOD_Projection(modeling.output = modeling.output,
                                             new.env = get_formal_data(modeling.output)@data.env.var,
                                             proj.name = temp_name,
                                             xy.new.env = get_formal_data(modeling.output)@coord,
                                             selected.models = models.kept.union,
                                             compress = TRUE,
                                             build.clamping.mask = F,
                                             do.stack=T, silent = T)@proj@val

        # transform array into data.frame
        out$predictions <- as.data.frame(out$predictions)
        names(out$predictions) <- unlist(lapply(strsplit(names(out$predictions),".", fixed=TRUE),
                                         function(x){
                                           x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
                                           data.set.id <- x.rev[1]
                                           cross.valid.id <- x.rev[2]
                                           algo.id <- paste(rev(x.rev[3:length(x.rev)]), collapse = ".", sep = "")
                                           model.id <- paste(modeling.output@sp.name,
                                                             data.set.id,
                                                             cross.valid.id,
                                                             algo.id, sep="_")
                                           return(model.id)
                                         }))
        # keep only wanted columns
        out$predictions <- out$predictions[,models.kept.union, drop=F]
        unlink(file.path(modeling.output@sp.name,paste("proj_", temp_name, sep="") ),recursive = TRUE, force = TRUE)
        cat("\n")
      }
#     }
    return(out)
  } else {
    cat("\n   ! No models kept due to threshold filtering... Ensemble Modeling was skipped!")
    return(NULL)
  }
}





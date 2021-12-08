'BIOMOD_EnsembleModeling' <- function(modeling.output,
                                      chosen.models = 'all',
                                      em.by = 'PA_dataset+repet',
                                      eval.metric = 'all',
                                      eval.metric.quality.threshold = NULL,
                                      eval.metric.user.data = NULL,
                                      models.eval.meth = c('KAPPA','TSS','ROC'),
                                      prob.mean = TRUE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      prob.median = FALSE,
                                      committee.averaging = FALSE,
                                      prob.mean.weight = FALSE,
                                      prob.mean.weight.decay = 'proportional',
                                      VarImport = 0)
{
  .bmCat("Build Ensemble Models")
  
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .BIOMOD_EnsembleModeling.check.args(modeling.output,
                                              chosen.models,
                                              eval.metric,
                                              eval.metric.quality.threshold,
                                              eval.metric.user.data,
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
  eval.metric.user <- args$eval.metric.user
  eval.metric.user.data <- args$eval.metric.user.data
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
  rm(args)
  
  ## Get selected options
  em.avail <- c('prob.mean', 'prob.cv', 'prob.ci.inf', 'prob.ci.sup',
                'prob.median', 'committee.averaging', 'prob.mean.weight')
  em.algo <- em.avail[c(prob.mean, prob.cv,  prob.ci,  prob.ci,
                        prob.median, committee.averaging, prob.mean.weight)]
  em.options <- list(em.by = em.by)
  expl_var_type = get_var_type(get_formal_data(modeling.output, 'expl.var'))
  expl_var_range = get_var_range(get_formal_data(modeling.output, 'expl.var'))
  
  
  ## 1. Create output object ----------------------------------------------------------------------
  EM <- new('BIOMOD.EnsembleModeling.out',
            sp.name = modeling.output@sp.name,
            expl.var.names = modeling.output@expl.var.names,
            em.by = em.by,
            modeling.id = modeling.output@modeling.id)
  EM@models.out.obj@link <- file.path(modeling.output@sp.name,
                                      paste0(modeling.output@sp.name, ".",
                                             modeling.output@modeling.id, ".models.out"))
  
  ## 2. Do Ensemble modeling ----------------------------------------------------------------------
  
  ## make a list of models names that will be combined together according to em.by argument
  em.mod.assemb <- .em.models.assembling(chosen.models, em.by)
  for (assemb in names(em.mod.assemb))
  {
    cat("\n\n  >", assemb, "ensemble modeling")
    models.kept <- em.mod.assemb[[assemb]]
    
    #### defined data that will be used for models performances calculation ####
    if (modeling.output@has.evaluation.data) {
      eval.obs <- get_formal_data(modeling.output, 'eval.resp.var')
      eval.expl <- get_formal_data(modeling.output, 'eval.expl.var')
    }
    
    ## subselection of observations according to dataset used to produce ensemble models
    obs <-  get_formal_data(modeling.output, 'resp.var')
    expl <- get_formal_data(modeling.output, 'expl.var')
    if (em.by %in% c("PA_dataset", 'PA_dataset+algo', 'PA_dataset+repet') &&
        unlist(strsplit(assemb, "_"))[3] != 'AllData') {
      if (inherits(get_formal_data(modeling.output), "BIOMOD.formated.data.PA")) {
        kept_cells <- get_formal_data(modeling.output)@PA[, unlist(strsplit(assemb, "_"))[3]]
      } else {
        kept_cells <- rep(TRUE, length(obs))
      }
      obs <- obs[kept_cells]
      expl <- expl[kept_cells, , drop = FALSE]
    }
    obs[is.na(obs)] <- 0
    
    
    #### get needed models predictions ############################################################
    ## if no prediction selected => swith to next model
    needed_predictions <- .get_needed_predictions(modeling.output, em.by, models.kept
                                                  , eval.metric, eval.metric.quality.threshold
                                                  , eval.metric.user, eval.metric.user.data)
    if (!length(needed_predictions)) next
    
    ## LOOP over evaluation metrics ##
    for (eval.m in eval.metric)
    {
      models.kept <- needed_predictions$models.kept[[eval.m]]
      models.kept.scores <- needed_predictions$models.kept.scores[[eval.m]]
      
      ## LOOP over em.algo ##
      for (algo in em.algo)
      {
        algo.1 = algo.2 = algo.3 = NULL
        models.kept.tmp = models.kept
        if (algo == 'prob.mean') {
          algo.1 = "Mean of probabilities"
          algo.2 = algo.3 = "EMmean"
        } else if (algo == 'prob.cv') {
          algo.1 = "Coef of variation of probabilities"
          algo.2 = algo.3 = "EMcv"
        } else if (algo == 'prob.median') {
          algo.1 = "Median of probabilities"
          algo.2 = algo.3 = "EMmedian"
        } else if (algo == 'prob.ci.inf') {
          algo.1 = "Confidence Interval"
          algo.2 = "EMciInf"
          algo.3 = "EMci"
        } else if (algo == 'prob.ci.sup') {
          algo.1 = "Confidence Interval"
          algo.2 = "EMciSup"
          algo.3 = "EMci"
        } else if (algo == 'committee.averaging') {
          algo.1 = "Committee averaging"
          algo.2 = algo.3 = "EMca"
          
          ## remove models if some thresholds are undefined
          models.kept.thresh <- unlist(lapply(models.kept.tmp, function(x){
            mod <- tail(unlist(strsplit(x,"_")), 3)[3]
            run <- tail(unlist(strsplit(x,"_")), 3)[2]
            dat <- tail(unlist(strsplit(x,"_")), 3)[1]
            return(get_evaluations(modeling.output)[eval.m, "Cutoff", mod, run, dat])
          }))
          names(models.kept.thresh) <- models.kept.tmp
          models.kept.tmp = models.kept.tmp[is.finite(models.kept.thresh)]
          models.kept.thresh.tmp = models.kept.thresh[is.finite(models.kept.thresh)]
        } else if (algo == 'prob.mean.weight') {
          algo.1 = "Probabilities weighting mean"
          algo.2 = algo.3 = "EMwmean"
          
          # remove SRE models if ROC
          models.kept.scores.tmp <- models.kept.scores
          if (eval.m == 'ROC') {
            sre.id <- grep("_SRE", models.kept.tmp)
            if (length(sre.id) > 0) {
              cat("\n      ! SRE modeling were switched off")
              models.kept.tmp <- models.kept.tmp[-sre.id]
              models.kept.scores.tmp <- models.kept.scores[-sre.id]
            }
          }
          
          ## remove models if score is not defined
          models.kept.tmp <- models.kept.tmp[is.finite(models.kept.scores.tmp)]
          models.kept.scores.tmp <- models.kept.scores.tmp[is.finite(models.kept.scores.tmp)]
          
          # weights are "decay" times decreased for each subsequent model in model quality order.
          models.kept.scores.tmp <- round(models.kept.scores.tmp, 3) # sometimes there can be a rounding issue in R, so here I make sure all values are rounded equally.
          
          # dealing with numerical decay
          cat("\n\t\t", " original models scores = ", models.kept.scores.tmp)
          if (is.numeric(prob.mean.weight.decay)) {
            DecayCount <- sum(models.kept.scores.tmp > 0)
            WOrder <- order(models.kept.scores.tmp, decreasing = TRUE)
            Dweights <- models.kept.scores.tmp
            for (J in 1:DecayCount) {
              Dweights[WOrder[J]] <- I(prob.mean.weight.decay ^ (DecayCount - J + 1))
            }
            # If 2 or more scores are identical -> make a mean weight between the ones concerned
            for (J in 1:length(models.kept.scores.tmp)) {
              comp = models.kept.scores.tmp[J] == models.kept.scores.tmp
              if (sum(comp) > 1) {
                Dweights[which(comp == TRUE)] <- mean(Dweights[which(comp == TRUE)])
              }
            }
            models.kept.scores.tmp <- round(Dweights, digits = 3)
            rm(list = c('Dweights', 'DecayCount', 'WOrder'))
          } else if (is.function(prob.mean.weight.decay)) { # dealing with function decay
            models.kept.scores.tmp <- sapply(models.kept.scores.tmp, prob.mean.weight.decay)
          }
          
          ### Standardise model weights
          models.kept.scores.tmp <- round(models.kept.scores.tmp / sum(models.kept.scores.tmp, na.rm = TRUE)
                                          , digits = 3)
          cat("\n\t\t", " final models weights = ", models.kept.scores.tmp)
        }
        
        
        #### Models building ##################################################
        cat("\n   >", algo.1, "...")
        model_name <- paste0(modeling.output@sp.name, "_", algo.2, "By", eval.m, "_", assemb)
        model.bm <- new(paste0(algo.3, "_biomod2_model"),
                        model = models.kept.tmp,
                        model_name = model_name,
                        model_class = algo.3,
                        model_options = em.options,
                        resp_name = modeling.output@sp.name,
                        expl_var_names = modeling.output@expl.var.names,
                        expl_var_type = expl_var_type,
                        expl_var_range = expl_var_range,
                        modeling.id = modeling.output@modeling.id)
        
        if (algo == 'prob.ci.inf') {
          model.bm@alpha <- prob.ci.alpha
          model.bm@side <- 'inferior'
        } else if (algo == 'prob.ci.sup') {
          model.bm@alpha <- prob.ci.alpha
          model.bm@side <- 'superior'
        } else if (algo == 'committee.averaging') {
          model.bm@thresholds <- models.kept.thresh.tmp
        } else if(algo == 'prob.mean.weight'){
          model.bm@penalization_scores <- models.kept.scores.tmp
        }
        
        #### Models Predictions ###############################################
        
        ## create the suitable directory architecture
        pred.bm.name <- paste0(model_name, ".predictions")
        pred.bm.outfile <- file.path(model.bm@resp_name, ".BIOMOD_DATA", model.bm@modeling.id,
                                     "ensemble.models", "ensemble.models.predictions",
                                     pred.bm.name)
        dir.create(dirname(pred.bm.outfile), showWarnings = FALSE, recursive = TRUE)
        
        ## store models prediction on the hard drive
        pred.bm <- predict(model.bm, expl
                           , formal_predictions = needed_predictions$predictions[, model.bm@model, drop = FALSE]
                           , on_0_1000 = TRUE)
        assign(pred.bm.name, pred.bm)
        save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
        rm(list = pred.bm.name)
        
        ## do the same for evaluation data
        if (exists('eval.obs') & exists('eval.expl')) {
          pred.bm.name <- paste0(model_name, ".predictionsEval")
          eval_pred.bm <- predict(model.bm, eval.expl)
          assign(pred.bm.name, eval_pred.bm)
          save(list = pred.bm.name, file = pred.bm.outfile, compress = TRUE)
          rm(list = pred.bm.name)
        }
        
        #### Models Evaluation ################################################
        
        if (length(models.eval.meth)) {
          cat("\n\t\t\tEvaluating Model stuff...")
          
          if (algo == 'prob.cv') { ## switch off evaluation process
            cross.validation <- matrix(NA, 4, length(models.eval.meth),
                                       dimnames = list(c("Testing.data", "Cutoff", "Sensitivity", "Specificity"),
                                                       models.eval.meth))
          } else {
            if (em.by == "PA_dataset+repet") {
              ## select the same evaluation data than formal models
              ## get formal models calib/eval lines
              calib_lines <- get_calib_lines(modeling.output)
              ## get info on wich dataset and which repet this ensemble model is based on
              pa_dataset_id <- paste0("_", unlist(strsplit(assemb, "_"))[3])
              repet_id <- paste0("_", unlist(strsplit(assemb, "_"))[2])
              ## define and extract the subset of points model will be evaluated on
              if (repet_id == "_Full") {
                eval_lines <- rep(TRUE, length(pred.bm))
              } else {
                ## trick to detect when it is a full model but with a non common name
                ## i.e. all lines used for calib => full model
                eval_lines <- !na.omit(calib_lines[, repet_id, pa_dataset_id])
                if (all(!eval_lines)) { eval_lines <- !eval_lines }
              }
            } else {
              eval_lines <- rep(TRUE, length(pred.bm))
            }
            
            cross.validation <- sapply(models.eval.meth,
                                       bm_FindOptimStat,
                                       Fit = pred.bm[eval_lines],
                                       Obs = obs[eval_lines])
            rownames(cross.validation) <- c("Testing.data", "Cutoff", "Sensitivity", "Specificity")
          }
          
          if (exists('eval_pred.bm')) {
            if (algo == 'prob.cv') { ## switch off evaluation process
              cross.validation <- matrix(NA, 5, length(models.eval.meth),
                                         dimnames = list(c("Testing.data", "Evaluating.data", "Cutoff"
                                                           , "Sensitivity", "Specificity"),
                                                         models.eval.meth))
            } else {
              true.evaluation <- sapply(models.eval.meth, function(x) {
                bm_FindOptimStat(Fit = eval_pred.bm * 1000,
                                 Obs = eval.obs,
                                 Fixed.thresh = cross.validation["Cutoff", x])
              })
              cross.validation <- rbind(cross.validation["Testing.data", ], true.evaluation)
              rownames(cross.validation) <- c("Testing.data", "Evaluating.data", "Cutoff", "Sensitivity", "Specificity")
            }
          }
          
          ## store results
          model.bm@model_evaluation <- t(round(cross.validation, digits = 3))
        }
        
        #### Models Variable Importance #######################################
        
        if (VarImport > 0) {
          cat("\n\t\t\tEvaluating Predictor Contributions...", "\n")
          model.bm@model_variables_importance <- bm_VariablesImportance(model.bm, expl, nb_rand = VarImport)
        }
        
        #### Models saving #####
        assign(model_name, model.bm)
        save(list = model_name, file = file.path(modeling.output@sp.name, "models",
                                                 modeling.output@modeling.id, model_name))
        
        #### Add to sumary objects ####
        EM@em.models <- c(EM@em.models, model.bm)
        EM@em.computed <- c(EM@em.computed, model_name)
      }
    }
  }
  
  ### fix models names ###
  names(EM@em.models) <- EM@em.computed
  model.name <- paste0(EM@sp.name, '.', EM@modeling.id, 'ensemble.models.out')
  assign(x = model.name, value = EM)
  save(list = model.name, file = file.path(EM@sp.name, model.name))
  
  .bmCat("Done")
  return(EM)
}

###################################################################################################

.BIOMOD_EnsembleModeling.check.args <- function(modeling.output,
                                                chosen.models,
                                                eval.metric,
                                                eval.metric.quality.threshold,
                                                eval.metric.user.data,
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
{
  ## 1. Check modeling.output -------------------------------------------------
  .fun_testIfInherits(TRUE, "modeling.output", modeling.output, "BIOMOD.models.out")
  
  ## 2. Check chosen.models ---------------------------------------------------
  if (!length(chosen.models) | (length(chosen.models) == 1 && chosen.models[1] == 'all')) {
    cat("\n   ! all models available will be included in ensemble.modeling")
    chosen.models <- modeling.output@models.computed
  } else {
    .fun_testIfIn(TRUE, "chosen.models", chosen.models, modeling.output@models.computed)
  }
  
  ## 3. Check eval.metric -----------------------------------------------------
  if (!is.null(eval.metric)) {
    if (!is.character(eval.metric)) {
      stop("eval.metric must be a character vector or NULL")
    }
    eval.metric.user = FALSE
    if ('user.defined' %in% eval.metric) {
      eval.metric.user = TRUE
      if (!is.null(eval.metric.user.data)) {
        .fun_testIfIn(TRUE, "chosen.models", chosen.models, colnames(eval.metric.user.data))
        # chosen.models.check <- chosen.models %in% colnames(eval.metric.user.data)
        # if (sum(!chosen.models.check) > 0) {
        #   stop(paste0("Some selected models do not have eval values in `eval.metric.user.data`: "
        #               , toString(chosen.models[!chosen.models.check])
        #               , "\nPlease provide eval values for all selected models ( "
        #               , toString(chosen.models), " )"))
        # }
        eval.metric.user.data <- eval.metric.user.data[, chosen.models, drop = FALSE]
        eval.metric <- rownames(eval.metric.user.data)
      } else {
        stop("eval.metric.user.data must be a data.frame or NULL")
      }
    } else {
      if ('all' %in% eval.metric) {
        eval.metric <- dimnames(get_evaluations(modeling.output))[[1]]
      }
      .fun_testIfIn(TRUE, "eval.metric", eval.metric, dimnames(get_evaluations(modeling.output))[[1]])
      # eval.metric.check <- eval.metric %in% dimnames(get_evaluations(modeling.output))[[1]]
      # if(sum(!eval.metric.check) > 0){
      #   stop(paste0("Some selected evaluation metrics are not available: "
      #               , toString(eval.metric[!eval.metric.check])
      #               , "\nPlease choose some in those computed yet ( "
      #               , toString(dimnames(get_evaluations(modeling.output))[[1]]), " )"))
      # }
    }
  }
  
  ## 4. Check eval.metric.quality.threshold -----------------------------------
  if (!is.null(eval.metric)) {
    if (!is.null(eval.metric.quality.threshold)) {
      if (!is.numeric(eval.metric.quality.threshold)) {
        stop("eval.metric.quality.threshold must be NULL or a numeric vector")
      }
      if (length(eval.metric) != length(eval.metric.quality.threshold)) {
        stop("you must specify as many eval.metric.quality.threshold as eval.metric (if you specify some)")
      }
      cat("\n   > Evaluation & Weighting methods summary :\n")
      cat(paste(eval.metric, eval.metric.quality.threshold, sep = " over ", collapse = "\n      ")
          , fill = TRUE, labels = "     ")
    } else {
      cat("\n   ! No eval.metric.quality.threshold -> All models will be kept for Ensemble Modeling")
      eval.metric.quality.threshold <- rep(0, length(eval.metric))
    }
  } else {
    eval.metric <- 'none'
  }
  
  ## 5. Check model.eval.meth -------------------------------------------------
  models.eval.meth <- unique(models.eval.meth)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC')
  .fun_testIfIn(TRUE, "models.eval.meth", models.eval.meth, avail.eval.meth.list)
  
  ## 6. Check selected EM algo ------------------------------------------------
  if (!is.logical(prob.mean) | !is.logical(prob.median) |
      !is.logical(prob.cv) | !is.logical(prob.ci) | 
      !is.logical(committee.averaging) | !is.logical(prob.mean.weight)) {
    stop("prob.mean, prob.cv, prob.ci, prob.median, committee.averaging and prob.mean.weight arguments must be logical")
  }
  if (is.null(eval.metric) && (committee.averaging | prob.mean.weight)) {
    stop("You must choose eval.metric if you want to compute Committee Averaging or Probability Weighted Mean algorithms")
  }
  
  ## 6.1 Check alpha for Confident interval
  if (prob.ci) {
    .fun_testIfPosNum(TRUE, "prob.ci.alpha", prob.ci.alpha)
    if (prob.ci.alpha <= 0 | prob.ci.alpha >= 0.5) {
      stop("prob.ci.alpha must be a numeric between 0 and 0.5")
    }
  }
  
  ## 6.2 Check decay for wmean
  if (prob.mean.weight) {
    if ((!is.numeric(prob.mean.weight.decay) &&
         !is.character(prob.mean.weight.decay) &&
         !is.function(prob.mean.weight.decay)) ||
        (is.numeric(prob.mean.weight.decay) && prob.mean.weight.decay < 0) ||
        (is.character(prob.mean.weight.decay) && prob.mean.weight.decay != 'proportional')) {
      stop("'prob.mean.weight.decay' should be either 'proportional', a numeric value > 0 or a function")
    }
  }
  
  ## 7. Check em.by -----------------------------------------------------------
  .fun_testIfIn(TRUE, "em.by", em.by, c('PA_dataset', 'algo', 'all', 'PA_dataset+repet', 'PA_dataset+algo'))
  
  
  return(list(modeling.output = modeling.output,
              chosen.models = chosen.models,
              eval.metric = eval.metric,
              eval.metric.quality.threshold = eval.metric.quality.threshold,
              eval.metric.user = eval.metric.user,
              eval.metric.user.data = eval.metric.user.data,
              models.eval.meth = models.eval.meth,
              prob.mean = prob.mean,
              prob.cv = prob.cv,
              prob.ci = prob.ci,
              prob.ci.alpha = prob.ci.alpha,
              prob.median = prob.median,
              committee.averaging = committee.averaging,
              prob.mean.weight = prob.mean.weight,
              prob.mean.weight.decay = prob.mean.weight.decay,
              em.by = em.by))
}

###################################################################################################

.em.models.assembling <- function(chosen.models, em.by)
{
  assembl.list = list()
  if (em.by == 'PA_dataset') {
    for (dat in .extractModelNamesInfo(chosen.models, info = 'data.set')) {
      assembl.list[[paste0("mergedAlgo_mergedRun_", dat)]] <- chosen.models[grep(paste0("_", dat, "_"), chosen.models)]
    }
  } else if (em.by == 'algo') {
    for (algo in .extractModelNamesInfo(chosen.models, info = 'models')) {
      assembl.list[[paste0(algo, "_mergedRun_mergedData")]] <- chosen.models[grep(paste0("_", algo), chosen.models)]
    }
  } else if (em.by == 'all') {
    assembl.list[["mergedAlgo_mergedRun_mergedData"]] <- chosen.models
  } else if (em.by == 'PA_dataset+repet') {
    for (dat in .extractModelNamesInfo(chosen.models, info = 'data.set')) {
      for (repet in .extractModelNamesInfo(chosen.models, info = 'run.eval')) {
        mod.tmp <- intersect(x = grep(paste0("_", dat, "_"), chosen.models)
                             , y = grep(paste0("_", repet, "_"), chosen.models))
        if (length(mod.tmp)) {
          assembl.list[[paste0("mergedAlgo_", repet, "_", dat)]] <- chosen.models[mod.tmp]
        }
      }
    }
  } else if (em.by == 'PA_dataset+algo') {
    for (dat in .extractModelNamesInfo(chosen.models, info = 'data.set')) {
      for (algo in .extractModelNamesInfo(chosen.models, info = 'models')) {
        mod.tmp <- intersect(x = grep(paste0("_", dat, "_"), chosen.models)
                             , y = grep(paste0("_", algo), chosen.models))
        if (length(mod.tmp)) {
          assembl.list[[paste0(algo, "_mergedRun_", dat)]] <- chosen.models[mod.tmp]
        }
      }
    }
  }
  return(assembl.list)
}


###################################################################################################

.get_needed_predictions <- function(modeling.output, em.by, models.kept, eval.metric
                                    , eval.metric.quality.threshold, eval.metric.user
                                    , eval.metric.user.data)
{
  out <- list(predictions = NULL, models.kept = NULL, models.kept.scores = NULL)
  for (eval.m in eval.metric) {
    if (eval.m != 'none') {
      if (eval.metric.user) {
        models.kept.scores <- eval.metric.user.data[eval.m, models.kept]
      } else {
        models.kept.scores <- unlist(lapply(models.kept, function(x) {
          mod <- tail(unlist(strsplit(x, "_")), 3)[3]
          run <- tail(unlist(strsplit(x, "_")), 3)[2]
          dat <- tail(unlist(strsplit(x, "_")), 3)[1]
          # select evaluations scores obtained for Evaluation Data if exists or CV if not
          if (modeling.output@has.evaluation.data) {
            return(get_evaluations(modeling.output)[eval.m, "Evaluating.data", mod, run, dat])
          } else {
            return(get_evaluations(modeling.output)[eval.m, "Testing.data", mod, run, dat])
          }
        }))
      }
      ## set NA to -1
      if (!is.null(models.kept.scores)) {
        models.kept.scores[is.na(models.kept.scores)] <- -1
      }
      thresh = eval.metric.quality.threshold[which(eval.metric == eval.m)]
      out$models.kept[[eval.m]] <- models.kept[models.kept.scores > thresh]
      out$models.kept.scores[[eval.m]] <- models.kept.scores[models.kept.scores > thresh]
    } else {
      out$models.kept[[eval.m]] <- models.kept
    }
  }
  
  models.kept.union <- unique(unlist(out$models.kept))
  
  if (length(models.kept.union)) {
    ## load prediction on each PA dataset
    if (em.by %in% c("PA_dataset", 'PA_dataset+algo', 'PA_dataset+repet')) {
      out$predictions <- as.data.frame(get_predictions(modeling.output, as.data.frame = TRUE)[, models.kept.union, drop = FALSE])
    } else{
      ## redo prediction on full data.set
      cat("\n   ! Models projections for whole zonation required...")
      temp_name <- paste0('tmp_', sub(".", "", as.character(format(Sys.time(), "%OS6")), fixed = TRUE))
      out$predictions <- BIOMOD_Projection(modeling.output = modeling.output,
                                           new.env = get_formal_data(modeling.output)@data.env.var,
                                           proj.name = temp_name,
                                           xy.new.env = get_formal_data(modeling.output)@coord,
                                           selected.models = models.kept.union,
                                           compress = TRUE,
                                           build.clamping.mask = FALSE,
                                           do.stack = TRUE,
                                           silent = TRUE
      )@proj@val
      
      # transform array into data.frame
      out$predictions <- as.data.frame(out$predictions)
      names(out$predictions) <- unlist(lapply(strsplit(names(out$predictions), ".", fixed = TRUE), function(x)
      {
        x.rev <- rev(x) ## we reverse the order of the splitted vector to have algo a t the end
        data.set.id <- x.rev[1]
        cross.valid.id <- x.rev[2]
        algo.id <- paste0(rev(x.rev[3:length(x.rev)]), collapse = ".")
        model.id <- paste(modeling.output@sp.name,
                          data.set.id,
                          cross.valid.id,
                          algo.id,
                          sep = "_")
        return(model.id)
      }))
      # keep only wanted columns
      out$predictions <- out$predictions[, models.kept.union, drop = FALSE]
      unlink(file.path(modeling.output@sp.name, paste0("proj_", temp_name))
             , recursive = TRUE, force = TRUE)
      cat("\n")
    }
    return(out)
  } else {
    cat("\n   ! No models kept due to threshold filtering... Ensemble Modeling was skipped!")
    return(NULL)
  }
}

###################################################################################################
##' @name bm_RunModelsLoop
##' @aliases bm_RunModelsLoop
##' @aliases bm_RunModel
##' @author Damien Georges
##' 
##' @title Loop to compute all single species distribution models
##' 
##' @description This internal \pkg{biomod2} function allows the user to compute all single 
##' species distribution models (asked by the \code{\link{BIOMOD_Modeling}} function).
##' 
##' 
##' @param bm.format a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param weights a \code{matrix} containing observation weights for each pseudo-absence (or 
##' \code{allData}) dataset
##' @param calib.lines a \code{matrix} containing calibration / validation lines for each 
##' pseudo-absence (or \code{allData}) x repetition (or \code{allRun}) combination that can be 
##' obtained with the \code{\link{bm_CrossValidation}} function
##' 
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param models a \code{vector} containing model names to be computed, must be among 
##' \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, 
##' \code{MARS}, \code{RF}, \code{MAXENT}, \code{MAXNET}
##' @param models.pa (\emph{optional, default} \code{NULL}) \cr 
##' A \code{list} containing for each model a \code{vector} defining which pseudo-absence datasets 
##' are to be used, must be among \code{colnames(bm.format@PA.table)}
##' @param bm.options a \code{\link{BIOMOD.models.options}} object returned by the  
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param scale.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions should be scaled with a 
##' binomial GLM or not
##' @param nb.cpu (\emph{optional, default} \code{1}) \cr 
##' An \code{integer} value corresponding to the number of computing resources to be used to 
##' parallelize the single models computation
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' @param do.progress (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the progress bar is to be rendered or not
##' 
##' 
##' @param model a \code{character} corresponding to the model name to be computed, must be either 
##' \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, 
##' \code{MARS}, \code{RF}, \code{MAXENT}, \code{MAXNET}
##' @param run.name a \code{character} corresponding to the model to be run (sp.name + pa.id + 
##' run.id)
##' @param dir.name (\emph{optional, default} \code{.}) \cr
##' A \code{character} corresponding to the modeling folder
##' @param Data a \code{data.frame} containing observations, coordinates and environmental 
##' variables that can be obtained with the \code{get_species_data} function
##' @param weights.vec a \code{vector} containing observation weights the concerned pseudo-absence 
##' (or \code{allData}) dataset
##' @param calib.lines.vec a \code{vector} containing calibration / validation lines for the 
##' concerned pseudo-absence (or \code{allData}) x repetition (or \code{allRun}) combination
##' @param eval.data (\emph{optional, default} \code{NULL}) \cr
##' A \code{data.frame} containing validation observations, coordinates and environmental 
##' variables that can be obtained with the \code{get_eval_data} function
##' 
##' 
##' 
##' @return  
##' 
##' A \code{list} containing for each model a \code{list} containing the following elements :
##' \itemize{
##'   \item{\code{model} : }{the name of correctly computed model}
##'   \item{\code{calib.failure} : }{the name of incorrectly computed model}
##'   \item{\code{pred} : }{the prediction outputs for calibration data}
##'   \item{\code{pred.eval} : }{the prediction outputs for evaluation data}
##'   \item{\code{evaluation} : }{the evaluation outputs returned by the 
##'   \code{\link{bm_FindOptimStat}} function}
##'   \item{\code{var.import} : }{the mean of variables importance returned by the 
##'   \code{\link{bm_VariablesImportance}} function}
##' }
##' 
##' 
##' @keywords models formula options CTA GLM GBM GAM RF ANN FDA SRE MARS MAXENT
##' 
##' 
##' @seealso \code{\link[rpart]{rpart}}, \code{\link[rpart]{prune}}, \code{\link[gbm]{gbm}}, 
##' \code{\link[MASS]{stepAIC}}, \code{\link[nnet]{nnet}}, \code{\link[earth]{earth}}, 
##' \code{\link[mda]{fda}}, \code{\link[mda]{mars}}, \code{\link[maxnet]{maxnet}}, 
##' \code{\link[randomForest]{randomForest}}, 
##' \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{bm_MakeFormula}}, \code{\link{bm_SampleFactorLevels}}, 
##' \code{\link{bm_FindOptimStat}}, \code{\link{bm_VariablesImportance}}
##' @family Secundary functions
##' 
##' 
##' 
##' @importFrom foreach foreach %dopar% 
## @importFrom doParallel registerDoParallel 
##' @importFrom rpart rpart prune
## @importFrom caret 
## @importFrom car 
## @importFrom gam gam step.Gam s
## @importFrom mgcv gam bam
##' @importFrom gbm gbm gbm.perf
##' @importFrom MASS stepAIC
##' @importFrom nnet nnet
##' @importFrom earth earth
##' @importFrom mda fda mars
##' @importFrom dplyr mutate_at select_at %>%
##' @importFrom maxnet maxnet
##' @importFrom randomForest randomForest
##' 
##' @export
##' 
##'
###################################################################################################

bm_RunModelsLoop <- function(bm.format,
                             weights,
                             calib.lines,
                             modeling.id,
                             models,
                             models.pa,
                             bm.options,
                             metric.eval,
                             var.import,
                             scale.models = TRUE,
                             nb.cpu = 1,
                             seed.val = NULL,
                             do.progress = TRUE)
{
  
  if (nb.cpu > 1) {
    if (.getOS() != "windows") {
      if (!isNamespaceLoaded("doParallel")) {
        if(!requireNamespace('doParallel', quietly = TRUE)) stop("Package 'doParallel' not found")
      }
      doParallel::registerDoParallel(cores = nb.cpu)
    } else {
      warning("Parallelisation with `foreach` is not available for Windows. Sorry.")
    }
  }
  
  ## PREPARE DATA ---------------------------------------------------------------------------------
  list.data <- list()
  pa.list = sapply(colnames(calib.lines), function(x) strsplit(x, "_")[[1]][2])
  for (pa.id in unique(pa.list)) { # loop on PA -------------------------------
    models.subset = models
    if (!is.null(models.pa)) {
      ## optional : subset of models associated to the concerned PA dataset
      models.subset = sapply(models.pa, function(x) pa.id %in% x)
      models.subset = names(models.pa)[which(models.subset == TRUE)]
    }
    
    data.all <- get_species_data(bm.format)
    if (pa.id %in% colnames(data.all)) {
      ## optional : subset of species data associated to the concerned PA dataset
      data.all <- data.all[which(data.all[, pa.id] == TRUE), ]
    }
    data.all <- data.all[, c(bm.format@sp.name, "x", "y", colnames(bm.format@data.env.var))]
    
    for (i in which(pa.list == pa.id)) { # loop on RUN ------------------------
      run.id = strsplit(colnames(calib.lines)[i], "_")[[1]][3]
      run.name = paste0(bm.format@sp.name, "_", pa.id, "_", run.id)
      
      for (modi in models.subset) { # loop on models --------------------------
        all.name = paste0(run.name, "_", modi)
        list.data[[all.name]] <- list(modi = modi,
                                      run.name = run.name,
                                      data.all = data.all,
                                      calib.lines.vec = na.omit(calib.lines[, i]), ## ATTENTION na.omit
                                      weights.vec = na.omit(weights[, pa.id])) ## ATTENTION na.omit
      }}}
  
  ## RUN models -----------------------------------------------------------------------------------
  out <- foreach(ii = 1:length(list.data)) %dopar%
    {
      cat('\n\n-=-=-=--=-=-=-', names(list.data)[ii], '\n')
      bm_RunModel(model = list.data[[ii]]$modi,
                  run.name = list.data[[ii]]$run.name,
                  dir.name = bm.format@dir.name,
                  modeling.id = modeling.id,
                  bm.options = bm.options,
                  Data = list.data[[ii]]$data.all,
                  weights.vec = list.data[[ii]]$weights.vec,
                  calib.lines.vec = list.data[[ii]]$calib.lines.vec,
                  eval.data = get_eval_data(bm.format),
                  metric.eval = metric.eval,
                  var.import = var.import,
                  scale.models = scale.models,
                  seed.val = seed.val,
                  do.progress = TRUE)
    }
  
  return(out)
}


# ---------------------------------------------------------------------------- #

##' 
##' @rdname bm_RunModelsLoop
##' @export
##' 

bm_RunModel <- function(model, run.name, dir.name = '.'
                        , modeling.id = '', bm.options
                        , Data, weights.vec, calib.lines.vec
                        , eval.data = NULL
                        , metric.eval = c('ROC','TSS','KAPPA'), var.import = 0
                        , scale.models = TRUE, nb.cpu = 1, seed.val = NULL, do.progress = TRUE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_RunModel.check.args(model, bm.options, Data, weights.vec, calib.lines.vec
                                  , eval.data, metric.eval, scale.models, seed.val, do.progress
                                  , criteria = NULL, Prev = NULL)
  if (is.null(args)) { return(NULL) }
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  ## get model name and names of categorical variables
  dir_name = dir.name
  model_name <- paste0(run.name, '_', model)
  
  ## 1. Create output object ----------------------------------------------------------------------
  ListOut <- list(model = NULL,
                  calib.failure = NULL,
                  pred = NULL,
                  pred.eval = NULL,
                  evaluation = NULL,
                  var.import = NULL)
  
  ## 2. CREATE MODELS -----------------------------------------------------------------------------
  set.seed(seed.val)
  
  if (model == "CTA") {
    ### 2.1 CTA model ----------------------------------------------------------
    cat('\n\t> CTA modeling...')
    
    # converting cost argument
    cost.tmp = bm.options@CTA$cost
    if (is.null(bm.options@CTA$cost)) { cost.tmp = rep(1, ncol(data_env)) }
    
    # defining rpart parameters for splitting function
    parms.tmp = bm.options@CTA$parms
    if (bm.options@CTA$parms == 'default') { parms.tmp = NULL }
    
    model.sp <- try(rpart(bm_MakeFormula(resp.name = resp_name
                                         , expl.var = head(data_env)
                                         , type = 'simple'
                                         , interaction.level = 0),
                          data = data_mod[calib.lines.vec, , drop = FALSE],
                          weights = weights,
                          method = bm.options@CTA$method,
                          parms = parms.tmp,
                          cost = cost.tmp,
                          control = eval(bm.options@CTA$control)
    ))
    
    if (!inherits(model.sp, "try-error")) {
      # select best trees --------------- May be done otherway
      tr <- as.data.frame(model.sp$cptable)
      tr$xsum <- tr$xerror + tr$xstd
      tr <- tr[tr$nsplit > 0,]
      Cp <- tr[tr$xsum == min(tr$xsum), "CP"]
      model.sp <- prune(model.sp, cp = Cp[length(Cp)])
      
      model.bm <- new("CTA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'CTA',
                      model_options = bm.options@CTA,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]),
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "GAM") {
    ### 2.2 GAM model ----------------------------------------------------------
    
    # package loading
    .load_gam_namespace(bm.options@GAM$algo)
    
    if (bm.options@GAM$algo == 'GAM_gam') { ## gam package
      
      
      # NOTE : To be able to take into account GAM options and weights we have to do a eval(parse(...))
      # it's due to GAM implementation (using of match.call() troubles)
      gamStart <- eval(parse(text = paste0("gam::gam(", resp_name, "~1 ,"
                                           , " data = data_mod[calib.lines.vec, , drop = FALSE], family = ", bm.options@GAM$family$family
                                           , "(link = '", bm.options@GAM$family$link, "')"
                                           , ", weights = weights.vec[calib.lines.vec])")))
      model.sp <- try(gam::step.Gam(gamStart,
                                    .scope(head(data_env), "gam::s", bm.options@GAM$k),
                                    data = data_mod[calib.lines.vec, , drop = FALSE],
                                    direction = "both",
                                    trace = bm.options@GAM$control$trace,
                                    control = bm.options@GAM$control))
    } else { ## mgcv package
      
      if (is.null(bm.options@GAM$myFormula)) {
        cat("\n\tAutomatic formula generation...")
        gam.formula <- bm_MakeFormula(resp.name = resp_name
                                      , expl.var = head(data_env)
                                      , type = bm.options@GAM$type
                                      , interaction.level = bm.options@GAM$interaction.level
                                      , k = bm.options@GAM$k)
        tmp = gsub("gam::", "", gam.formula)
        gam.formula = as.formula(paste0(tmp[c(2,1,3)], collapse = " "))
      } else {
        gam.formula <- bm.options@GAM$myFormula
      }
      
      if (bm.options@GAM$algo == 'GAM_mgcv') {
        cat('\n\t> GAM (mgcv) modeling...')
        
        model.sp <- try(mgcv::gam(gam.formula,
                                  data = data_mod[calib.lines.vec, , drop = FALSE],
                                  family = bm.options@GAM$family,
                                  weights = weights.vec[calib.lines.vec],
                                  control = bm.options@GAM$control))
        
      } else if (bm.options@GAM$algo == 'BAM_mgcv') { ## big data.frame gam version
        cat('\n\t> BAM (mgcv) modeling...')
        model.sp <- try(mgcv::bam(gam.formula,
                                  data = data_mod[calib.lines.vec, , drop = FALSE],
                                  family = bm.options@GAM$family,
                                  weights = weights.vec[calib.lines.vec],
                                  control = bm.options@GAM$control))
      }
    }
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("GAM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GAM',
                      model_subclass = bm.options@GAM$algo,
                      model_options = bm.options@GAM,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]),
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "GBM") {
    ### 2.3 GBM model ----------------------------------------------------------
    cat('\n\t> GBM modeling...')
    
    model.sp <- try(gbm(formula = bm_MakeFormula(resp.name = resp_name
                                                 , expl.var = head(data_env)
                                                 , type = 'simple'
                                                 , interaction.level = 0),
                        data = data_mod[calib.lines.vec, , drop = FALSE],
                        distribution = bm.options@GBM$distribution,
                        var.monotone = rep(0, length = ncol(data_env)),
                        weights = weights,
                        interaction.depth = bm.options@GBM$interaction.depth,
                        n.minobsinnode = bm.options@GBM$n.minobsinnode,
                        shrinkage = bm.options@GBM$shrinkage,
                        bag.fraction = bm.options@GBM$bag.fraction,
                        train.fraction = bm.options@GBM$train.fraction,
                        n.trees = bm.options@GBM$n.trees,
                        verbose = bm.options@GBM$verbose,
                        #class.stratify.cv = bm.options@GBM$class.stratify.cv,
                        cv.folds = bm.options@GBM$cv.folds,
                        n.cores = bm.options@GBM$n.cores
    ))
    
    if (!inherits(model.sp, "try-error")) {
      best.iter <- try(gbm.perf(model.sp, method = bm.options@GBM$perf.method , plot.it = FALSE))
      
      model.bm <- new("GBM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GBM',
                      n.trees_optim = best.iter,
                      model_options = bm.options@GBM,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "GLM"){
    ### 2.4 GLM model ----------------------------------------------------------
    cat('\n\t> GLM modeling...')
    if (is.null(bm.options@GLM$myFormula)) {
      cat("\n\tAutomatic formula generation...")
      glm.formula <- bm_MakeFormula(resp.name = resp_name
                                    , expl.var = head(data_env)
                                    , type = bm.options@GLM$type
                                    , interaction.level = bm.options@GLM$interaction.level)
    } else {
      glm.formula <- bm.options@GLM$myFormula
    }
    if (bm.options@GLM$test != 'none') {
      ## make the model selection
      glmStart <- glm(eval(parse(text = paste0(resp_name, "~1"))), 
                      data = data_mod[calib.lines.vec, , drop = FALSE], 
                      family = bm.options@GLM$family,
                      control = eval(bm.options@GLM$control),
                      weights = weights.vec[calib.lines.vec],
                      mustart = rep(bm.options@GLM$mustart, sum(calib.lines.vec)), 
                      model = TRUE)
      
      ## remove warnings
      warn <- options('warn')
      options(warn = -1)
      model.sp <- try(stepAIC(glmStart,
                              glm.formula,
                              data = data_mod[calib.lines.vec, , drop = FALSE],
                              direction = "both",
                              trace = FALSE,
                              k = criteria,
                              weights = weights.vec[calib.lines.vec], 
                              steps = 10000,
                              mustart = rep(bm.options@GLM$mustart, sum(calib.lines.vec))))
      ## reexec warnings
      options(warn)
      
    } else {
      ## keep the total model
      model.sp <- try(glm(glm.formula,
                          data = cbind(data_mod[calib.lines.vec, , drop = FALSE], 
                                       data.frame("weights" = weights.vec[calib.lines.vec])), 
                          family = bm.options@GLM$family,
                          control = eval(bm.options@GLM$control),
                          weights = weights,
                          model = TRUE))
    }
    if (!inherits(model.sp, "try-error")) {
      cat("\n\tselected formula : ")
      print(model.sp$formula, useSource = FALSE, showEnv = FALSE)
      model.bm <- new("GLM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GLM',
                      model_options = bm.options@GLM,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]),
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "MARS"){
    ### 2.5 MARS model ---------------------------------------------------------
    
    cat('\n\t> MARS modeling...')
    if (is.null(bm.options@MARS$myFormula)) {
      cat("\n\tAutomatic formula generation...")
      mars.formula <- bm_MakeFormula(resp.name = resp_name
                                     , expl.var = head(data_env)
                                     , type = bm.options@MARS$type
                                     , interaction.level = bm.options@MARS$interaction.level)
    } else {
      mars.formula <- bm.options@MARS$myFormula
    }
    
    ## deal with nk argument : if not defined, set up to default mars value i.e max(21, 2 * ncol(x) + 1)
    nk <- bm.options@MARS$nk
    if (is.null(nk)) {
      nk <- min(200, max(20, 2 * length(expl_var_names))) + 1
    }
    
    model.sp <- try(earth(formula = mars.formula,
                          data = data_mod[calib.lines.vec, , drop = FALSE], 
                          weights = weights,
                          glm = list(family = binomial),
                          ncross = 0,
                          keepxy = FALSE,
                          # degree = bm.options@MARS$degree,
                          pmethod = bm.options@MARS$pmethod,
                          nprune = bm.options@MARS$nprune,
                          nk = nk,
                          penalty = bm.options@MARS$penalty,
                          thresh = bm.options@MARS$thresh))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("MARS_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MARS',
                      model_options = bm.options@MARS,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "FDA") {
    ### 2.6 FDA model ----------------------------------------------------------
    
    cat('\n\t> FDA modeling...')
    model.sp <- try(do.call(fda, c(list(formula = bm_MakeFormula(resp.name = resp_name
                                                                 , expl.var = head(data_env)
                                                                 , type = 'simple'
                                                                 , interaction.level = 0),
                                        data = data_mod[calib.lines.vec, , drop = FALSE], 
                                        method = eval(parse(text = call(bm.options@FDA$method))),
                                        weights = weights.vec[calib.lines.vec]),
                                   bm.options@FDA$add_args)))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("FDA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'FDA',
                      model_options = bm.options@FDA,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]),
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "ANN") {
    ### 2.7 ANN model ----------------------------------------------------------
    
    cat('\n\t> ANN modeling...')
    size = bm.options@ANN$size
    decay = bm.options@ANN$decay
    if (is.null(size) | is.null(decay) | length(size) > 1 | length(decay) > 1) {
      
      ## define the size and decay to test
      if (is.null(size)) { size <- c(2, 4, 6, 8) }
      if (is.null(decay)) { decay <- c(0.001, 0.01, 0.05, 0.1) }
      
      ## do cross validation test to find the optimal values of size and decay parameters (prevent from overfitting)
      CV_nnet <- bm_CVnnet(Input = data_env[calib.lines.vec, , drop = FALSE],
                           Target = data_sp[calib.lines.vec], 
                           size = size,
                           decay = decay,
                           maxit = bm.options@ANN$maxit,
                           nbCV = bm.options@ANN$NbCV,
                           weights = weights.vec[calib.lines.vec],
                           seedval = seed.val)
      
      ## get the optimised parameters values
      decay <- CV_nnet[1, 2]
      size <- CV_nnet[1, 1]
    }
    
    model.sp <- try(nnet(formula = bm_MakeFormula(resp.name = resp_name
                                                  , expl.var = head(data_env)
                                                  , type = 'simple'
                                                  , interaction.level = 0),
                         data = data_mod[calib.lines.vec, , drop = FALSE], 
                         size = size,
                         rang = bm.options@ANN$rang,
                         decay = decay,
                         weights = weights,
                         maxit = bm.options@ANN$maxit,
                         trace = FALSE))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("ANN_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'ANN',
                      model_options = bm.options@ANN,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "RF") {
    ### 2.8 RF model -----------------------------------------------------------
    
    cat('\n\t> RF modeling...')
    if (bm.options@RF$do.classif) {
      # defining occurences as factor for doing classification and not regression in RF
      data_mod <- data_mod %>% mutate_at(resp_name, factor)
    }
    
    # mtry.tmp = bm.options@RF$mtry
    # if (bm.options@RF$mtry == 'default') { mtry.tmp = NULL }
    
    model.sp <- try(randomForest(formula = bm_MakeFormula(resp.name = resp_name
                                                          , expl.var = head(data_env)
                                                          , type = 'simple'
                                                          , interaction.level = 0),
                                 data = data_mod[calib.lines.vec, , drop = FALSE],
                                 ntree = bm.options@RF$ntree,
                                 # weights = weights.vec[calib.lines.vec],
                                 # mtry = mtry.tmp, 
                                 importance = FALSE,
                                 norm.votes = TRUE,
                                 strata = factor(c(0, 1)),
                                 sampsize = unlist(ifelse(!is.null(bm.options@RF$sampsize), list(bm.options@RF$sampsize), length(data_sp[calib.lines.vec]))),
                                 nodesize = bm.options@RF$nodesize,
                                 maxnodes = bm.options@RF$maxnodes))
    
    if (bm.options@RF$do.classif) {
      # canceling occurences class modifications
      data_mod <- data_mod %>% mutate_at(resp_name, function(.x) {
        .x %>% as.character() %>% as.numeric()
      })
    }
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("RF_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'RF',
                      model_options = bm.options@RF,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "SRE") {
    ### 2.9 SRE model ----------------------------------------------------------
    
    cat('\n\t> SRE modeling...')
    model.sp <- try(bm_SRE(resp.var = data_sp[calib.lines.vec],
                           expl.var = data_env[calib.lines.vec, , drop = FALSE],
                           new.env = NULL,
                           quant = bm.options@SRE$quant,
                           do.extrem = TRUE))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("SRE_biomod2_model",
                      extremal_conditions = model.sp,
                      model_name = model_name,
                      model_class = 'SRE',
                      model_options = bm.options@SRE,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  } else if (model == "MAXENT") {
    ### 2.10 MAXENT model ---------------------------------------------
    cat('\n\t> MAXENT modeling...')
    categorical_var <- .get_categorical_names(data_env)
    
    MWD <- .maxent.prepare.workdir(sp_name = resp_name
                                   , run_name = run.name
                                   , data_sp = data_sp
                                   , data_xy = data_xy
                                   , data_env = data_env
                                   , categorical_var = categorical_var
                                   , calib.lines.vec = calib.lines.vec
                                   , data_eval = eval.data
                                   , dir.name = dir_name
                                   , modeling.id = modeling.id
                                   , background_data_dir = bm.options@MAXENT$background_data_dir)
    
    # file to log potential errors
    maxent_stderr_file <- paste0(MWD$m_outdir, "/maxent.stderr")
    
    maxent.args <- 
      c(
        ifelse(is.null(bm.options@MAXENT$memory_allocated),"",
               paste0("-mx", bm.options@MAXENT$memory_allocated, "m")), 
        ifelse(is.null(bm.options@MAXENT$initial_heap_size), "",
               paste0(" -Xms", bm.options@MAXENT$initial_heap_size)),
        ifelse(is.null(bm.options@MAXENT$max_heap_size), "",
               paste0(" -Xmx", bm.options@MAXENT$max_heap_size)),
        paste0(" -jar ", 
               file.path(bm.options@MAXENT$path_to_maxent.jar, "maxent.jar")),
        paste0(" environmentallayers=\"", MWD$m_backgroundFile, "\""), 
        paste0(" samplesfile=\"", MWD$m_speciesFile, "\""),
        paste0(" projectionlayers=\"", gsub(", ", ",", toString(MWD$m_predictFile)), "\""),
        paste0(" outputdirectory=\"", MWD$m_outdir, "\""),
        paste0(" outputformat=logistic "), 
        ifelse(length(categorical_var),
               paste0(" togglelayertype=", categorical_var, collapse = " "),
               ""),
        " redoifexists",
        paste0(" visible=", bm.options@MAXENT$visible),
        paste0(" linear=", bm.options@MAXENT$linear),
        paste0(" quadratic=", bm.options@MAXENT$quadratic),
        paste0( " product=", bm.options@MAXENT$product),
        paste0(" threshold=", bm.options@MAXENT$threshold),
        paste0(" hinge=", bm.options@MAXENT$hinge),
        paste0(" lq2lqptthreshold=", bm.options@MAXENT$lq2lqptthreshold),
        paste0(" l2lqthreshold=", bm.options@MAXENT$l2lqthreshold),
        paste0(" hingethreshold=", bm.options@MAXENT$hingethreshold),
        paste0(" beta_threshold=", bm.options@MAXENT$beta_threshold),
        paste0(" beta_categorical=", bm.options@MAXENT$beta_categorical),
        paste0(" beta_lqp=", bm.options@MAXENT$beta_lqp),
        paste0(" beta_hinge=", bm.options@MAXENT$beta_hinge),
        paste0(" betamultiplier=", bm.options@MAXENT$betamultiplier),
        paste0(" defaultprevalence=", bm.options@MAXENT$defaultprevalence),
        " autorun ",
        " nowarnings ", 
        " notooltips ",
        " noaddsamplestobackground"
      )
    
    system2(command = "java", args = maxent.args,
            wait = TRUE,
            stdout = "", stderr = maxent_stderr_file)
    
    maxent_exec_output <- readLines(maxent_stderr_file)
    
    if(any(grepl(pattern = "Error", x = maxent_exec_output))) {
      g.pred <- NA
      class(g.pred) <- "try-error"
      cat( 
        paste0("\n*** Error in MAXENT, more info available in ",
               maxent_stderr_file)
      )
      
    } else {
      model.bm <- new("MAXENT_biomod2_model",
                      model_output_dir = MWD$m_outdir,
                      model_name = model_name,
                      model_class = 'MAXENT',
                      model_options = bm.options@MAXENT,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
      
      # for MAXENT predictions are calculated in the same time than models building to save time.
      cat("\n Getting predictions...")
      g.pred <- try(round(as.numeric(read.csv(MWD$m_outputFile)[, 3]) * 1000))
      
      if (var.import > 0) {
        cat("\n Getting predictor contributions...")
        variables.importance <- bm_VariablesImportance(bm.model = model.bm
                                                       , expl.var = data_env
                                                       , nb.rep = var.import
                                                       , temp_workdir = MWD$m_outdir
                                                       , seed.val = seed.val
                                                       , do.progress = do.progress)
      }
    }
  } else if (model == "MAXNET") {
    ### 2.11 MAXNET model -------------------------------------------
    
    cat('\n\t> MAXNET modeling...')
    model.sp <- try(maxnet(p = data_sp[calib.lines.vec], data = data_env[calib.lines.vec, , drop = FALSE]))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("MAXNET_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MAXNET',
                      model_options = bm.options@MAXNET,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(data_env[calib.lines.vec, , drop = FALSE]), 
                      expl_var_range = get_var_range(data_env[calib.lines.vec, , drop = FALSE]))
    }
  }
  
  ## 3. CREATE PREDICTIONS ------------------------------------------------------------------------
  temp_workdir = NULL
  
  if (model != "MAXENT") {
    g.pred <- try(predict(model.bm, data_env, on_0_1000 = TRUE, seedval = seed.val, temp_workdir = temp_workdir))
  }
  
  if (model == "MAXENT" & !inherits(g.pred, 'try-error')) {
    temp_workdir = model.bm@model_output_dir
  }
  
  ## scale or not predictions -------------------------------------------------
  if (scale.models & !inherits(g.pred, 'try-error')) {
    cat("\n\tModel scaling...")
    model.bm@scaling_model <- try(.scaling_model(g.pred / 1000, data_sp, weights = weights.vec))
    ## with weights
    g.pred <- try(predict(model.bm, data_env, on_0_1000 = TRUE, seedval = seed.val, temp_workdir = temp_workdir))
  }
  
  ## check predictions existence and stop execution if not ok -----------------
  test_pred_ok <- TRUE
  if (inherits(g.pred, "try-error")) { # model calibration or prediction failed
    test_pred_ok <- FALSE
    cat("\n*** inherits(g.pred,'try-error')")
  } else if (all(is.na(g.pred))) { # only NA predicted
    test_pred_ok <- FALSE
    cat("\n*** only NA predicted")
  } else if (length(unique(na.omit(g.pred))) <= 1) { # single value predicted
    test_pred_ok <- FALSE
    cat("\n*** single value predicted")
  }
  
  if (test_pred_ok) {
    # keep the model name
    ListOut$model <- model_name
  } else {
    # keep the name of uncompleted modelisations
    cat("\n   ! Note : ", model_name, "failed!\n")
    ListOut$calib.failure = model_name
    return(ListOut) ## end of function.
  }
  
  ## make prediction on evaluation data ---------------------------------------
  if (!is.null(eval.data)) {
    g.pred.eval <- try(
      predict(model.bm, 
              eval.data[, expl_var_names, drop = FALSE], 
              on_0_1000 = TRUE, 
              seedval = seed.val, 
              temp_workdir = temp_workdir)
    )
  }
  
  ## SAVE predictions ---------------------------------------------------------
  ListOut$pred <- g.pred
  if (exists("g.pred.eval")) { ListOut$pred.eval <- g.pred.eval }
  
  
  ## 4. EVALUATE MODEL ----------------------------------------------------------------------------
  if (length(metric.eval) > 0) {
    cat("\n\tEvaluating Model stuff...")
    
    ## Check no NA in g.pred to avoid evaluation failures
    na_cell_id <- which(is.na(g.pred))
    if (length(na_cell_id) > 0) {
      eval.lines.vec <- eval.lines.vec[-na_cell_id]
      cat('\n\tNote : some NA occurs in predictions')
    }
    
    if (length(which(eval.lines.vec == TRUE)) < length(g.pred)) {
      ## CALIBRATION & VALIDATION LINES -------------------------------------------------
      cross.validation <- foreach(xx = metric.eval, .combine = "rbind") %do% {
        bm_FindOptimStat(metric.eval = xx,
                         obs = data_sp[which(eval.lines.vec == FALSE)],
                         fit = g.pred[which(eval.lines.vec == FALSE)])
      }
      colnames(cross.validation)[which(colnames(cross.validation) == "best.stat")] <- "calibration"
      
      stat.validation <- foreach(xx = metric.eval, .combine = "rbind") %do% {
        bm_FindOptimStat(metric.eval = xx,
                         obs = data_sp[which(eval.lines.vec == TRUE)],
                         fit = g.pred[which(eval.lines.vec == TRUE)],
                         threshold = cross.validation$cutoff[
                           which(cross.validation$metric.eval == xx)
                         ])
      }
      cross.validation$validation <- stat.validation$best.stat
    } else {
      ## NO VALIDATION LINES -----------------------------------------------------
      cross.validation <- foreach(xx = metric.eval, .combine = "rbind") %do% {
        bm_FindOptimStat(metric.eval = xx,
                         obs = data_sp[which(eval.lines.vec == TRUE)],
                         fit = g.pred[which(eval.lines.vec == TRUE)])
      }
      colnames(cross.validation)[which(colnames(cross.validation) == "best.stat")] <- "calibration"
      cross.validation$validation <- NA
    }
    
    
    
    if (exists('g.pred.eval')) {
      
      ## Check no NA in g.pred.eval to avoid evaluation failures
      na_cell_id <- which(is.na(g.pred.eval))
      if (length(na_cell_id) > 0) {
        g.pred.eval.without.na <- g.pred.eval[-na_cell_id]
        eval.data <- eval.data[-na_cell_id, ]
        cat('\n\tNote : some NA occurs in evaluation predictions')
      } else {
        g.pred.eval.without.na <- g.pred.eval
      }
      
      stat.evaluation <- foreach(xx = metric.eval, .combine = "rbind") %do% {
        bm_FindOptimStat(metric.eval = xx,
                         obs = eval.data[, 1],
                         fit = g.pred.eval.without.na,
                         threshold = cross.validation["cutoff", xx])
      }
      cross.validation$evaluation <- stat.evaluation$best.stat
    } else {
      cross.validation$evaluation <- NA
    }
    
    ## store results
    for (col.i in 2:ncol(cross.validation)) {
      cross.validation[, col.i] <- round(cross.validation[, col.i], digits = 3)
    }
    ListOut$evaluation <- cross.validation
    model.bm@model_evaluation <- cross.validation
    rm(cross.validation)
  }
  
  
  ## 5. COMPUTE VARIABLES IMPORTANCE --------------------------------------------------------------
  if (var.import > 0) {
    cat("\n\tEvaluating Predictor Contributions...")
    if (model != "MAXENT") {
      variables.importance <- bm_VariablesImportance(bm.model = model.bm
                                                     , expl.var = data_env
                                                     , nb.rep = var.import
                                                     , seed.val = seed.val
                                                     , do.progress = do.progress)
    }
    ListOut$var.import <- variables.importance
    model.bm@model_variables_importance <- variables.importance
    rm(variables.importance)
    cat("\n")
  }
  
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  nameModel = paste(run.name, model, sep = "_") 
  assign(x = nameModel, value = model.bm)
  save(list = nameModel, file = file.path(dir_name, resp_name, "models", modeling.id, nameModel), compress = TRUE)
  
  
  return(ListOut)
}


###################################################################################################

.bm_RunModel.check.args <- function(model, bm.options, Data, weights.vec, calib.lines.vec
                                    , eval.data, metric.eval, scale.models, seed.val = NULL, do.progress = TRUE
                                    , criteria = NULL, Prev = NULL)
{
  ## 0. Do some cleaning over Data argument -----------------------------------
  data_sp <- Data[, 1]
  data_xy <- Data[, c("x", "y")]
  data_env <- Data[, -c(1, which(colnames(Data) %in% c("x", "y"))), drop = FALSE]
  resp_name <- colnames(Data)[1] ## species name
  expl_var_names <- colnames(data_env) ## explanatory variable names
  # replace Pseudo absences selected (NA) into true absences (0).. for model computing purpose
  if (sum(is.na(Data[, 1]))) { data_sp[which(is.na(data_sp))] <- 0 }
  
  ## 1. Check CalibLines argument ---------------------------------------------
  if (any(calib.lines.vec == FALSE)) ## if some lines for evaluation...
  {
    eval.lines.vec <- !calib.lines.vec
    # ...test if there is (pseudo)absences AND presences in evaluation and calibration datasets
    if (length(which(data_sp[calib.lines.vec] == 0)) == 0 ||
        length(which(data_sp[calib.lines.vec] == 0)) == length(calib.lines.vec) ||
        length(which(data_sp[eval.lines.vec] == 0)) == 0 ||
        length(which(data_sp[eval.lines.vec] == 0)) == length(eval.lines.vec)) {
      warning(paste0(resp_name, " ", model,
                     " was switched off because of no both presences and absences data given"),
              immediate. = TRUE)
      return(NULL)
    }
  } else { ## evaluation = calibration dataset
    eval.lines.vec <- calib.lines.vec
    # ...test if there is absences AND presences in whole dataset
    if (length(which(data_sp == 0)) == 0 ||
        length(which(data_sp == 0)) == length(data_sp)) {
      warning(paste0(resp_name, " ", model,
                     " was switched off because of no both presences and absences data given (full model)"),
              immediate. = TRUE)
      return(NULL)
    }
  }
  
  ## 2. Check weights argument ------------------------------------------------
  if (is.null(weights.vec)) { weights.vec <- rep(1, nrow(Data)) }
  ## These models require data and weights to be in the same dataset
  if (model %in% c('ANN', 'MARS', 'CTA', 'GBM')) {
    data_env_w <- cbind(data_env, weights.vec)
    colnames(data_env_w) <- c(colnames(data_env), "weights")
  } else {
    data_env_w <- data_env
  }
  
  ## 3. Check scale.models argument -------------------------------------------
  if (model == "SRE") { scale.models <- FALSE } else if (model %in% c("ANN", "FDA")) { scale.models <- TRUE }
  
  
  ## 4. Check bm.options argument ---------------------------------------------
  seedval = NULL
  if (model == "GLM") {
    cat('\nModel=GLM')
    if (!is.null(bm.options@GLM$myFormula)) {
      cat('\n\tformula = ', paste(bm.options@GLM$myFormula[2],
                                  bm.options@GLM$myFormula[1],
                                  bm.options@GLM$myFormula[3]))
    } else {
      cat(' (', bm.options@GLM$type, 'with',
          ifelse(bm.options@GLM$interaction.level == 0,
                 'no interaction )',
                 paste('order', bm.options@GLM$interaction.level, 'interaction level )')
          ))
    }
    if (bm.options@GLM$test == "AIC") {
      criteria <- 2
      cat("\n\tStepwise procedure using AIC criteria")
    } else if (bm.options@GLM$test == "BIC") {
      criteria <- log(ncol(data_env))
      cat("\n\tStepwise procedure using BIC criteria")
    } else if (bm.options@GLM$test == "none") {
      criteria <- 0
      cat("\n\tNo stepwise procedure")
      cat("\n\t! You might be confronted to model convergence issues !")
    }
  } else if (model == "GBM") {
    cat("\nModel=Generalised Boosting Regression \n")
    cat("\t", bm.options@GBM$n.trees, "maximum different trees and ", bm.options@GBM$cv.folds, " Fold Cross-Validation")
    seedval = 456 # to be able to refind our trees MAY BE BAD
  } else if (model == "GAM") {
    cat("\nModel=GAM")
    cat("\n\t", bm.options@GAM$algo, "algorithm chosen")
    seedval = 321 # to be able to refind our trees MAY BE BAD
  } else if (model == "CTA") {
    cat("\nModel=Classification tree \n")
    cat("\t", bm.options@CTA$control$xval, "Fold Cross-Validation")
    seedval = 123 # to be able to refind our trees MAY BE BAD
  } else if (model == "ANN") {
    cat("\nModel=Artificial Neural Network \n")
    cat("\t", bm.options@ANN$NbCV, "Fold Cross Validation + 3 Repetitions")
    seedval = 555 # to be able to refind our trees MAY BE BAD
  } else if (model == "SRE") {
    cat("\nModel=Surface Range Envelop")
  } else if (model == "FDA"){
    cat("\nModel=Flexible Discriminant Analysis")
  } else if (model == "MARS"){
    cat("\nModel=Multiple Adaptive Regression Splines")
    if (!is.null(bm.options@MARS$myFormula)) {
      cat('\n\tformula = ', paste(bm.options@MARS$myFormula[2],
                                  bm.options@MARS$myFormula[1],
                                  bm.options@MARS$myFormula[3]))
    } else {
      cat(' (', bm.options@MARS$type, 'with',
          ifelse(bm.options@MARS$interaction.level == 0,
                 'no interaction )',
                 paste('order', bm.options@MARS$interaction.level, 'interaction level )')
          ))
    }
    cat("\n")
  } else if (model == "RF") {
    cat("\nModel=Breiman and Cutler's random forests for classification and regression")
    seedval = 71
  } else if (model == 'MAXENT') {
    cat('\nModel=MAXENT')
  } else if (model == 'MAXNET') {
    cat('\nModel=MAXNET')
  }
  # else if (model == 'MAXENT.Tsuruoka') {
  #   cat('\nModel=MAXENT.Tsuruoka')
  # }
  if (!is.null(seed.val)) {
    seedval = seed.val
  }
  
  ## 5. Check Prev argument ---------------------------------------------------
  if (model == "GLM" | model == "GAM") {
    Prev <- sum(data_sp, na.rm = TRUE) / length(data_sp)
  }
  
  ## 6. Check models.eval.meth arguments --------------------------------------
  metric.eval <- unique(metric.eval)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC')
  # .fun_testIfIn(TRUE, "metric.eval", metric.eval, avail.eval.meth.list)
  if (sum(!(metric.eval %in% avail.eval.meth.list)) > 0) {
    tmp = which(metric.eval %in% avail.eval.meth.list)
    warnings(paste0(toString(metric.eval[!tmp]), ' were switched off !'), imediate = TRUE)
    metric.eval <- metric.eval[tmp]
  }  
  
  data_mod <- cbind(data_sp, data_env_w)
  colnames(data_mod) <- c(resp_name, colnames(data_env_w))
  
  return(list(data_sp = data_sp,
              data_xy = data_xy,
              data_env = data_env,
              data_mod = data_mod,
              weights.vec = weights.vec,
              eval.lines.vec = eval.lines.vec,
              criteria = criteria,
              Prev = Prev, 
              metric.eval = metric.eval,
              eval.data = eval.data,
              scale.models = scale.models,
              resp_name = resp_name,
              expl_var_names = expl_var_names,
              seed.val = seedval,
              do.progress = do.progress))
}


.maxent.prepare.workdir <- function(sp_name, run_name = NULL, data_sp, data_xy, data_env
                                    , categorical_var = NULL, calib.lines.vec = NULL, data_eval
                                    , dir.name = '.', modeling.id = '', background_data_dir = 'default')
{
  cat('\n\t\tCreating Maxent Temp Proj Data...')
  
  ## initialise output
  MWD <- list()
  class(MWD) <- "maxent_workdir_info"
  
  ## default parameters setting
  if (is.null(run_name)) { run_name <- sp_name }
  if (is.null(calib.lines.vec)) { calib.lines.vec <- rep(TRUE, nrow(data_env)) }
  
  ## define all paths to files needed by MAXENT
  nameFolder = file.path(dir.name, sp_name, 'models', modeling.id)
  m_outdir <- file.path(nameFolder, paste0(run_name, '_MAXENT_outputs'))
  m_predictDir <- file.path(m_outdir, "Predictions")
  MWD$m_outdir <- m_outdir
  MWD$m_outputFile <- file.path(m_outdir, paste0(run_name, '_Pred_swd.csv'))
  MWD$m_predictDir <- m_predictDir
  
  ## directories creation
  dir.create(m_outdir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  dir.create(m_predictDir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  
  ## transform categorical variables into numeric to avoid factors being saved 
  ## as characters, which are not readable by maxent
  data_env <- .categorical2numeric(data_env, categorical_var)
  
  ## Presence Data --------------------------------------------------------------------------------
  presLines <- which((data_sp == 1) & calib.lines.vec)
  absLines <- which((data_sp == 0) & calib.lines.vec)
  Sp_swd <- cbind(rep(run_name, length(presLines))
                  , data_xy[presLines, ]
                  , data_env[presLines, , drop = FALSE])
  colnames(Sp_swd) <- c('species', 'X', 'Y', colnames(data_env))
  
  m_speciesFile <- file.path(m_outdir, "Sp_swd.csv")
  write.table(Sp_swd, file = m_speciesFile, quote = FALSE, row.names = FALSE, sep = ",")
  MWD$m_speciesFile <- m_speciesFile
  
  ## Background Data (create background file only if needed) --------------------------------------
  if (background_data_dir == 'default')  {
    # keep only 0 of calib lines
    Back_swd <- cbind(rep("background", length(absLines))
                      , data_xy[absLines, ]
                      , data_env[absLines, , drop = FALSE])
    colnames(Back_swd) <- c("background", colnames(Back_swd)[-1])
    
    m_backgroundFile <- file.path(m_outdir, "Back_swd.csv")
    write.table(Back_swd, file = m_backgroundFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_backgroundFile <- m_backgroundFile
  } else { ## use background directory given as an option
    MWD$m_backgroundFile <- background_data_dir
  }
  
  ## Prediction Data ------------------------------------------------------------------------------
  Pred_swd <- cbind(rep("predict", nrow(data_xy)), data_xy, data_env)
  colnames(Pred_swd)  <- c("predict", "x", "y", colnames(data_env))
  
  m_predictFile <- file.path(m_predictDir, "Pred_swd.csv")
  write.table(Pred_swd, file = m_predictFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  MWD$m_predictFile <- m_predictFile
  
  ## dealing with independent evaluation data -----------------------------------------------------
  if (!is.null(data_eval)) {
    Pred_eval_swd <- cbind(rep("predictEval", nrow(data_eval))
                           , data_eval[, c("x", "y")]
                           , data_eval[, colnames(data_env), drop = FALSE])
    colnames(Pred_eval_swd) <- c("predict", colnames(Back_swd)[-1])
    
    m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
    write.table(Pred_eval_swd, file = m_predictEvalFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_predictEvalFile <- m_predictEvalFile
  }
  
  return(MWD)
}

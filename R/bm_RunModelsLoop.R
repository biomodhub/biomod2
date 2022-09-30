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
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param model a \code{character} corresponding to the model name to be computed, must be either 
##' \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, 
##' \code{MARS}, \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}
##' @param bm.options a \code{\link{BIOMOD.models.options}} object returned by the  
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param var.import (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param save.output (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether all outputs should be saved on hard drive or not 
##' (\emph{! strongly recommended !})
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
##' @param dir.name (\emph{optional, default} \code{.}) \cr
##' A \code{character} corresponding to the modeling folder
##' @param weights a \code{vector} of \code{numeric} values corresponding to observation weights 
##' (one per observation)
##' @param nam a \code{character} corresponding to the model to be run (name + run.id)
##' @param Data a \code{data.frame} containing \code{data.species} and \code{data.env.var} slots 
##' of \code{bm.format} parameter
##' @param calib.lines a \code{data.frame} containing \code{data.split.table} slot of 
##' \code{bm.format} parameter, or an extraction of \code{data.species} slot (for a specific PA 
##' dataset extracted from \code{PA.table} slot)
##' @param xy a \code{data.frame} containing \code{coord} slot of \code{bm.format} 
##' parameter (for a specific PA dataset extracted from \code{PA.table} slot of \code{bm.format} 
##' parameter)
##' @param eval.data a \code{data.frame} containing \code{eval.data.species} and 
##' \code{eval.data.env.var} slots of \code{bm.format} parameter
##' @param eval.xy a \code{data.frame} containing \code{eval.coord} slot of \code{bm.format} 
##' parameter
##' 
##' 
##' @return  
##' 
##' A \code{list} containing for each model a \code{list} containing the following elements :
##' \itemize{
##'   \item{\code{model} : }{the name of correctly computed model}
##'   \item{\code{calib.failure} : }{the name of incorrectly computed model}
##'   \item{\code{pred} : }{the prediction outputs for calibration data}
##'   \item{\code{pred.eval} : }{the prediction outputs for validation data}
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
##' @importFrom dplyr mutate_at filter select_at %>% pull
##' @importFrom maxnet maxnet
##' @importFrom randomForest randomForest
##' 
##' @export
##' 
##'
###################################################################################################


bm_RunModelsLoop <- function(bm.format,
                             modeling.id,
                             model,
                             bm.options,
                             metric.eval,
                             var.import,
                             save.output = TRUE,
                             scale.models = TRUE,
                             nb.cpu = 1,
                             seed.val = NULL,
                             do.progress = TRUE)
{
  cat("\n\n-=-=-=- Run : ", bm.format$name, '\n')
  res.sp.run <- list()
  
  for (i in 1:ncol(bm.format$calib.lines)) { # loop on RunEval
    run.id = dimnames(bm.format$calib.lines)[[2]][i]
    run.name = paste0(bm.format$name, run.id)
    cat('\n\n-=-=-=--=-=-=-', run.name, '\n')
    
    if (nb.cpu > 1) {
      if (.getOS() != "windows") {
        if (!isNamespaceLoaded("doParallel")) { requireNamespace("doParallel") }
        doParallel::registerDoParallel(cores = nb.cpu)
      } else {
        warning("Parallelisation with `foreach` is not available for Windows. Sorry.")
      }
    }
    res.sp.run[[run.id]] = foreach(modi = model) %dopar%
      {
        bm_RunModel(model = modi,
                    Data = bm.format$dataBM,
                    modeling.id = modeling.id,
                    bm.options = bm.options,
                    calib.lines = na.omit(bm.format$calib.lines[, i, ]), ## transform 3D calib.lines obj into a 1D vector
                    weights = na.omit(bm.format$weights),
                    nam = run.name,
                    dir.name = bm.format$dir.name,
                    xy = bm.format$xy,
                    eval.data = bm.format$eval.data,
                    eval.xy = bm.format$eval.xy,
                    metric.eval = metric.eval,
                    var.import = var.import,
                    save.output = TRUE, ## save.output
                    scale.models = scale.models,
                    seed.val = seed.val,
                    do.progress = TRUE)
      }
    names(res.sp.run[[run.id]]) <- model
  }
  
  return(res.sp.run)
}


###################################################################################################

##' 
##' @rdname bm_RunModelsLoop
##' @export
##' 

bm_RunModel <- function(model, Data, modeling.id = '', bm.options, calib.lines, weights, nam,
                        dir.name = '.', xy = NULL, eval.data = NULL, eval.xy = NULL, 
                        metric.eval = c('ROC','TSS','KAPPA'), var.import = 0,
                        save.output = FALSE, scale.models = TRUE, nb.cpu = 1, seed.val = NULL,
                        do.progress = TRUE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_RunModel.check.args(model, Data, bm.options, calib.lines, weights, eval.data
                                  , metric.eval, scale.models, seed.val, do.progress)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## get model name and names of categorical variables
  dir_name = dir.name
  model_name <- paste0(nam, '_', model)
  categorical_var <- unlist(sapply(expl_var_names, function(x) {
    if (is.factor(Data[, x])) { return(x) } else { return(NULL) }
  }))
  
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
    ## 2.1 CTA model ----------------------------------------------------------
    cat('\n\t> CTA modeling...')
    
    # converting cost argument
    cost.tmp = bm.options@CTA$cost
    if (is.null(bm.options@CTA$cost)) { cost.tmp = rep(1, (ncol(Data) - 2)) }
    
    # defining rpart parameters for splitting function
    parms.tmp = bm.options@CTA$parms
    if (bm.options@CTA$parms == 'default') { parms.tmp = NULL }
    
    model.sp <- try(rpart(
      bm_MakeFormula(resp.name = colnames(Data)[1]
                     , expl.var = head(Data[, -c(1, ncol(Data)), drop = FALSE])
                     , type = 'simple'
                     , interaction.level = 0),
      data = Data[calib.lines, ],
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "GAM") {
    ## 2.2 GAM model ----------------------------------------------------------
   
    # package loading
    .load_gam_namespace(bm.options@GAM$algo)
    
    if (bm.options@GAM$algo == 'GAM_gam') { ## gam package

      
      # NOTE : To be able to take into account GAM options and weights we have to do a eval(parse(...))
      # it's due to GAM implementation (using of match.call() troubles)
      gamStart <- eval(parse(text = paste0("gam::gam(", colnames(Data)[1], "~1 ,"
                                           , " data = Data[calib.lines,,drop=FALSE], family = ", bm.options@GAM$family$family
                                           , "(link = '", bm.options@GAM$family$link, "')"
                                           , ", weights = weights[calib.lines])")))
      model.sp <- try(gam::step.Gam(gamStart,
                                    .scope(Data[1:3, -c(1, ncol(Data))], "gam::s", bm.options@GAM$k),
                                    data = Data[calib.lines, , drop = FALSE],
                                    direction = "both",
                                    trace = bm.options@GAM$control$trace,
                                    control = bm.options@GAM$control))
    } else { ## mgcv package

      if (is.null(bm.options@GAM$myFormula)) {
        cat("\n\tAutomatic formula generation...")
        gam.formula <- bm_MakeFormula(resp.name = resp_name
                                      , expl.var = head(Data[, expl_var_names, drop = FALSE])
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
                                  data = Data[calib.lines, , drop = FALSE],
                                  family = bm.options@GAM$family,
                                  weights = weights[calib.lines],
                                  control = bm.options@GAM$control))
        
      } else if (bm.options@GAM$algo == 'BAM_mgcv') { ## big data.frame gam version
        cat('\n\t> BAM (mgcv) modeling...')
        model.sp <- try(mgcv::bam(gam.formula,
                                  data = Data[calib.lines, , drop = FALSE],
                                  family = bm.options@GAM$family,
                                  weights = weights[calib.lines],
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "GBM") {
    ## 2.3 GBM model ----------------------------------------------------------
    
    cat('\n\t> GBM modeling...')
    model.sp <- try(gbm(formula = bm_MakeFormula(resp.name = colnames(Data)[1]
                                                 , expl.var = head(Data)[, expl_var_names, drop = FALSE]
                                                 , type = 'simple'
                                                 , interaction.level = 0),
                        data = Data[calib.lines, , drop = FALSE],
                        distribution = bm.options@GBM$distribution,
                        var.monotone = rep(0, length = ncol(Data) - 2), # -2 because of removing of sp and weights
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "GLM"){
    ## 2.4 GLM model ----------------------------------------------------------
    
    cat('\n\t> GLM modeling...')
    if (is.null(bm.options@GLM$myFormula)) {
      cat("\n\tAutomatic formula generation...")
      glm.formula <- bm_MakeFormula(resp.name = colnames(Data)[1]
                                    , expl.var = head(Data)
                                    , type = bm.options@GLM$type
                                    , interaction.level = bm.options@GLM$interaction.level)
    } else {
      glm.formula <- bm.options@GLM$myFormula
    }
    
    if (bm.options@GLM$test != 'none') {
      ## make the model selection
      glmStart <- glm(eval(parse(text = paste0(colnames(Data)[1], "~1"))), 
                      data = Data[calib.lines, , drop = FALSE], 
                      family = bm.options@GLM$family,
                      control = eval(bm.options@GLM$control),
                      weights = weights[calib.lines],
                      mustart = rep(bm.options@GLM$mustart, sum(calib.lines)), 
                      model = TRUE)
      
      ## remove warnings
      warn <- options('warn')
      options(warn = -1)
      model.sp <- try(stepAIC(glmStart,
                              glm.formula,
                              data = Data[calib.lines, , drop = FALSE],
                              direction = "both",
                              trace = FALSE,
                              k = criteria,
                              weights = weights[calib.lines], 
                              steps = 10000,
                              mustart = rep(bm.options@GLM$mustart, sum(calib.lines))))
      ## reexec warnings
      options(warn)
      
    } else {
      ## keep the total model
      model.sp <- try(glm(glm.formula,
                          data = cbind(Data[calib.lines, , drop = FALSE], 
                                       matrix(weights[calib.lines], ncol = 1, dimnames = list(NULL, "weights"))), 
                          family = bm.options@GLM$family,
                          control = eval(bm.options@GLM$control),
                          weights = weights,
                          model = TRUE))
    }
    
    if (!inherits(model.sp, "try-error")) {
      cat("\n\tselected formula : ")
      print(model.sp$formula, useSource = FALSE)
      
      model.bm <- new("GLM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GLM',
                      model_options = bm.options@GLM,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "MARS"){
    ## 2.5 MARS model ---------------------------------------------------------
    
    cat('\n\t> MARS modeling...')
    if (is.null(bm.options@MARS$myFormula)) {
      cat("\n\tAutomatic formula generation...")
      mars.formula <- bm_MakeFormula(resp.name = colnames(Data)[1]
                                     , expl.var = head(Data)[, -ncol(Data), drop = FALSE]
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
                          data = Data[calib.lines, , drop = FALSE], 
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "FDA") {
    ## 2.6 FDA model ----------------------------------------------------------
    
    cat('\n\t> FDA modeling...')
    model.sp <- try(do.call(fda, c(list(formula = bm_MakeFormula(resp.name = colnames(Data)[1]
                                                                 , expl.var = head(Data)[, expl_var_names, drop = FALSE]
                                                                 , type = 'simple'
                                                                 , interaction.level = 0),
                                        data = Data[calib.lines, , drop = FALSE], 
                                        method = eval(parse(text = call(bm.options@FDA$method))),
                                        weights = weights[calib.lines]),
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "ANN") {
    ## 2.7 ANN model ----------------------------------------------------------
    
    cat('\n\t> ANN modeling...')
    size = bm.options@ANN$size
    decay = bm.options@ANN$decay
    if (is.null(size) | is.null(decay) | length(size) > 1 | length(decay) > 1) {
      
      ## define the size and decay to test
      if (is.null(size)) { size <- c(2, 4, 6, 8) }
      if (is.null(decay)) { decay <- c(0.001, 0.01, 0.05, 0.1) }
      
      ## do cross validation test to find the optimal values of size and decay parameters (prevent from overfitting)
      CV_nnet <- bm_CVnnet(Input = Data[, expl_var_names, drop = FALSE],
                           Target = Data[calib.lines, 1], 
                           size = size,
                           decay = decay,
                           maxit = bm.options@ANN$maxit,
                           nbCV = bm.options@ANN$NbCV,
                           weights = weights[calib.lines],
                           seedval = seed.val)
      
      ## get the optimised parameters values
      decay <- CV_nnet[1, 2]
      size <- CV_nnet[1, 1]
    }
    
    model.sp <- try(nnet(formula = bm_MakeFormula(resp.name = resp_name
                                                  , expl.var = head(Data[, expl_var_names, drop = FALSE])
                                                  , type = 'simple'
                                                  , interaction.level = 0),
                         data = Data[calib.lines, , drop = FALSE], 
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "RF") {
    ## 2.8 RF model -----------------------------------------------------------
    
    cat('\n\t> RF modeling...')
    if (bm.options@RF$do.classif) {
      # defining occurences as factor for doing classification and not regression in RF
      Data <- Data %>% mutate_at(resp_name, factor)
    }
    
    # mtry.tmp = bm.options@RF$mtry
    # if (bm.options@RF$mtry == 'default') { mtry.tmp = NULL }
    
    model.sp <- try(randomForest(formula = bm_MakeFormula(resp.name = resp_name
                                                          , expl.var = head(Data)
                                                          , type = 'simple'
                                                          , interaction.level = 0),
                                 data = Data[calib.lines, ],
                                 ntree = bm.options@RF$ntree,
                                 # weights = weights,
                                 # mtry = mtry.tmp, 
                                 importance = FALSE,
                                 norm.votes = TRUE,
                                 strata = factor(c(0, 1)),
                                 sampsize = ifelse(!is.null(bm.options@RF$sampsize), bm.options@RF$sampsize, nrow(Data[calib.lines, ])),
                                 nodesize = bm.options@RF$nodesize,
                                 maxnodes = bm.options@RF$maxnodes))
    
    if (bm.options@RF$do.classif) {
      # canceling occurences class modifications
      Data <- Data %>% mutate_at(resp_name, function(.x) {
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "SRE") {
    ## 2.9 SRE model ----------------------------------------------------------
    
    cat('\n\t> SRE modeling...')
    model.sp <- try(bm_SRE(resp.var = Data[calib.lines, 1],
                           expl.var = Data[calib.lines, expl_var_names, drop = FALSE],
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
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  } else if (model == "MAXENT.Phillips") {
    ## 2.10 MAXENT.Phillips model ---------------------------------------------
    
    cat('\n\t> MAXENT.Phillips modeling...')
    MWD <- .maxent.prepare.workdir(Data, xy, calib.lines, RunName = nam,
                                   eval.data, eval.xy, dir.name = dir_name, species.name = resp_name,
                                   modeling.id = modeling.id,
                                   background_data_dir = bm.options@MAXENT.Phillips$background_data_dir)
    
    maxent.cmd <- paste0("java ",
                         ifelse(is.null(bm.options@MAXENT.Phillips$memory_allocated),
                                "",
                                paste0("-mx", bm.options@MAXENT.Phillips$memory_allocated, "m")), 
                         " -jar ", file.path(bm.options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"),
                         " environmentallayers=\"", MWD$m_backgroundFile, 
                         "\" samplesfile=\"", MWD$m_speciesFile,
                         "\" projectionlayers=\"", gsub(", ", ",", toString(MWD$m_predictFile)),
                         "\" outputdirectory=\"", MWD$m_outdir,
                         "\" outputformat=logistic ", 
                         ifelse(length(categorical_var),
                                paste0(" togglelayertype=", categorical_var, collapse = " "),
                                ""),
                         " redoifexists",
                         " visible=", bm.options@MAXENT.Phillips$visible,
                         " linear=", bm.options@MAXENT.Phillips$linear,
                         " quadratic=", bm.options@MAXENT.Phillips$quadratic,
                         " product=", bm.options@MAXENT.Phillips$product,
                         " threshold=", bm.options@MAXENT.Phillips$threshold,
                         " hinge=", bm.options@MAXENT.Phillips$hinge,
                         " lq2lqptthreshold=", bm.options@MAXENT.Phillips$lq2lqptthreshold,
                         " l2lqthreshold=", bm.options@MAXENT.Phillips$l2lqthreshold,
                         " hingethreshold=", bm.options@MAXENT.Phillips$hingethreshold,
                         " beta_threshold=", bm.options@MAXENT.Phillips$beta_threshold,
                         " beta_categorical=", bm.options@MAXENT.Phillips$beta_categorical,
                         " beta_lqp=", bm.options@MAXENT.Phillips$beta_lqp,
                         " beta_hinge=", bm.options@MAXENT.Phillips$beta_hinge,
                         " betamultiplier=", bm.options@MAXENT.Phillips$betamultiplier,
                         " defaultprevalence=", bm.options@MAXENT.Phillips$defaultprevalence,
                         " autorun nowarnings notooltips noaddsamplestobackground")
    
    system(command = maxent.cmd, wait = TRUE, intern = TRUE,
           ignore.stdout = FALSE, ignore.stderr = FALSE)
    
    model.bm <- new("MAXENT.Phillips_biomod2_model",
                    model_output_dir = MWD$m_outdir,
                    model_name = model_name,
                    model_class = 'MAXENT.Phillips',
                    model_options = bm.options@MAXENT.Phillips,
                    dir_name = dir_name,
                    resp_name = resp_name,
                    expl_var_names = expl_var_names,
                    expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                    expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    
    # for MAXENT.Phillips predicitons are calculated in the same time than models building to save time.
    cat("\n Getting predictions...")
    g.pred <- try(round(as.numeric(read.csv(MWD$m_outputFile)[, 3]) * 1000))
    
    if (var.import > 0) {
      cat("\n Getting predictor contributions...")
      variables.importance <- bm_VariablesImportance(bm.model = model.bm
                                                     , expl.var = Data[, expl_var_names, drop = FALSE]
                                                     , nb.rep = var.import
                                                     , temp_workdir = MWD$m_outdir
                                                     , seed.val = seed.val
                                                     , do.progress = do.progress)
    }
  } else if(model == "MAXENT.Phillips.2")
  {
    ## 2.11 MAXENT.Phillips.2 model -------------------------------------------
    
    cat('\n\t> MAXENT.Phillips modeling...')
    model.sp <- try(maxnet(p = Data %>% filter(calib.lines) %>% pull(resp_name), 
                           data = Data %>% filter(calib.lines) %>% select_at(expl_var_names)
                           # f = if(!is.null(bm.options@MAXENT.Phillips.2@))
    ))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("MAXENT.Phillips.2_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MAXENT.Phillips.2',
                      model_options = bm.options@MAXENT.Phillips.2,
                      dir_name = dir_name,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calib.lines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calib.lines, expl_var_names, drop = FALSE]))
    }
  }
  
  ## 2.12 MAXENT.Tsuruoka model -----------------------------------------------
  # if(model == "MAXENT.Tsuruoka"){
  #   model.sp <- try(stop('MAXENT.Tsuruoka is depreacated(because maxent package is not maintained anymore)'))
  #   # model.sp <- try(maxent::maxent(feature_matrix = Data[calib.lines, expl_var_names, drop = FALSE],
  #   #                                code_vector = as.factor(Data[calib.lines, 1]),
  #   #                                l1_regularizer = bm.options@MAXENT.Tsuruoka$l1_regularizer,
  #   #                                l2_regularizer = bm.options@MAXENT.Tsuruoka$l2_regularizer,
  #   #                                use_sgd = bm.options@MAXENT.Tsuruoka$use_sgd,
  #   #                                set_heldout = bm.options@MAXENT.Tsuruoka$set_heldout,
  #   #                                verbose = bm.options@MAXENT.Tsuruoka$verbose))
  # 
  #   if( !inherits(model.sp,"try-error") ){
  #     model.bm <- new("MAXENT.Tsuruoka_biomod2_model",
  #                     model = model.sp,
  #                     model_name = model_name,
  #                     model_class = 'MAXENT.Tsuruoka',
  #                     model_options = bm.options@MAXENT.Tsuruoka,
  #                     resp_name = resp_name,
  #                     expl_var_names = expl_var_names,
  #                     expl_var_type = get_var_type(Data[calib.lines,expl_var_names,drop = FALSE]),
  #                     expl_var_range = get_var_range(Data[calib.lines,expl_var_names,drop = FALSE]))
  #   }
  # }
  
  
  ## 3. CREATE PREDICTIONS ------------------------------------------------------------------------
  temp_workdir = NULL
  if (model == "MAXENT.Phillips") {
    temp_workdir = model.bm@model_output_dir
  }
  
  if (model != "MAXENT.Phillips") {
    g.pred <- try(predict(model.bm, Data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE
                          , seedval = seed.val, temp_workdir = temp_workdir))
  }
  
  ## scale or not predictions -------------------------------------------------
  if (scale.models & !inherits(g.pred, 'try-error')) {
    cat("\n\tModel scaling...")
    model.bm@scaling_model <- try(.scaling_model(g.pred / 1000, Data[, 1, drop = TRUE], weights = weights))
    ## with weights
    g.pred <- try(predict(model.bm, Data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE
                          , seedval = seed.val, temp_workdir = temp_workdir))
  }
  
  ## check predictions existence and stop execution if not ok -----------------
  test_pred_ok <- TRUE
  if (inherits(g.pred, "try-error")) { # model calibration or prediction failed
    test_pred_ok <- FALSE
    cat("\n*** inherits(g.pred,'try-error')")
  } else if (sum(!is.na(g.pred)) <= 1) { # only NA predicted
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
    g.pred.eval <- try(predict(model.bm, eval.data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE
                               , seedval = seed.val, temp_workdir = temp_workdir))
  }
  
  ## SAVE predictions ---------------------------------------------------------
  if (save.output) {
    ListOut$pred <- g.pred
    if (exists("g.pred.eval")) { ListOut$pred.eval <- g.pred.eval }
  }
  
  
  ## 4. EVALUATE MODEL ----------------------------------------------------------------------------
  if (length(metric.eval) > 0) {
    cat("\n\tEvaluating Model stuff...")
    
    ## Check no NA in g.pred to avoid evaluation failures
    na_cell_id <- which(is.na(g.pred))
    if (length(na_cell_id)) {
      evalLines <- evalLines[!(evalLines %in% na_cell_id)]
      cat('\n\tNote : some NA occurs in predictions')
    }
    
    cross.validation <- sapply(metric.eval, function(.x) {
      bm_FindOptimStat(metric.eval = .x,
                       obs = Data %>% filter(evalLines) %>% pull(1),
                       fit = g.pred[evalLines])
    })
    rownames(cross.validation) <- c("Testing.data", "Cutoff", "Sensitivity", "Specificity")
    
    if (exists('g.pred.eval')) {
      
      ## Check no NA in g.pred.eval to avoid evaluation failures
      na_cell_id <- which(is.na(g.pred.eval))
      if (length(na_cell_id)) {
        g.pred.eval.without.na <- g.pred.eval[-na_cell_id]
        eval.data <- eval.data[-na_cell_id, ]
        cat('\n\tNote : some NA occurs in evaluation predictions')
      } else {
        g.pred.eval.without.na <- g.pred.eval
      }
      
      true.evaluation <- sapply(metric.eval, function(x) {
        bm_FindOptimStat(metric.eval = x,
                         obs = eval.data[, 1],
                         fit = g.pred.eval.without.na,
                         threshold = cross.validation["Cutoff", x])
      })
      
      cross.validation <- rbind(cross.validation["Testing.data", ], true.evaluation)
      rownames(cross.validation) <- c("Testing.data", "Evaluating.data", "Cutoff", "Sensitivity", "Specificity")
    }
    
    ## store results
    cross.validation <- t(round(cross.validation, digits = 3))
    ListOut$evaluation <- cross.validation
    model.bm@model_evaluation <- cross.validation
    rm(cross.validation)
  }
  
  
  ## 5. COMPUTE VARIABLES IMPORTANCE --------------------------------------------------------------
  if (var.import > 0) {
    cat("\n\tEvaluating Predictor Contributions...")
    if (model != "MAXENT.Phillips") {
      variables.importance <- bm_VariablesImportance(bm.model = model.bm
                                                     , expl.var = Data[, expl_var_names, drop = FALSE]
                                                     , nb.rep = var.import
                                                     , seed.val = seed.val
                                                     , do.progress = do.progress)
    }
    model.bm@model_variables_importance <- variables.importance
    
    ## only the mean of variables importance run is returned
    ListOut$var.import <- round(rowMeans(variables.importance, na.rm = TRUE), digits = 3)
    rm(variables.importance)
    cat("\n")
  }
  
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  nameModel = paste(nam, model, sep = "_") 
  assign(x = nameModel, value = model.bm)
  save(list = nameModel, file = file.path(dir_name, resp_name, "models", modeling.id, nameModel), compress = TRUE)
  
  
  return(ListOut)
}


###################################################################################################

.bm_RunModel.check.args <- function(model, Data, bm.options, calib.lines, weights, eval.data
                                    , metric.eval, scale.models, criteria = NULL, Prev = NULL
                                    , seed.val = NULL, do.progress = TRUE)
{
  ## 0. Do some cleaning over Data argument -----------------------------------
  resp_name <- colnames(Data)[1] ## species name
  expl_var_names <- colnames(Data)[-1] ## explanatory variable names
  # replace Pseudo absences selected (NA) into true absences (0).. for model computing purpose
  if (sum(is.na(Data[, 1]))) { Data[which(is.na(Data[, 1])), 1] < - 0 }
  
  ## 1. Check CalibLines argument ---------------------------------------------
  if (sum(!calib.lines) > 0) ## if some lines for evaluation...
  {
    evalLines <- !calib.lines
    # ...test if there is absences AND presences in evaluation and calibration datasets
    if (sum(Data[calib.lines, 1] == 0) == 0 ||
        sum(Data[calib.lines, 1] == 0) == sum(calib.lines) ||
        sum(Data[evalLines, 1] == 0) == 0 ||
        sum(Data[evalLines, 1] == 0) == sum(evalLines)) {
      warning(paste0(colnames(Data)[1], " ", model,
                     " was switched off because of no both presences and absences data given"),
              immediate. = TRUE)
      return(NULL)
    }
  } else { ## evaluation = calibration dataset
    evalLines <- calib.lines
    # ...test if there is absences AND presences in whole dataset
    if (sum(Data[, 1] == 0) == 0 ||
        sum(Data[, 1] == 0) == nrow(Data)) {
      warning(paste0(colnames(Data)[1], " ", model,
                     " was switched off because of no both presences and absences data given (full model)"),
              immediate. = TRUE)
      return(NULL)
    }
  }
  
  ## 2. Check weights argument ------------------------------------------------
  if (is.null(weights)) { weights <- rep(1, nrow(Data)) }
  ## These models require data and weights to be in the same dataset
  if (model %in% c('GBM', 'CTA', 'ANN', 'FDA', 'GAM', 'MARS')) {
    Data <- cbind(Data, weights)
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
      criteria <- log(ncol(Data))
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
  } else if (model == 'MAXENT.Phillips') {
    cat('\nModel=MAXENT.Phillips')
  } else if (model == 'MAXENT.Phillips.2') {
    cat('\nModel=MAXENT.Phillips (maxnet)')
  }
  # else if (model == 'MAXENT.Tsuruoka') {
  #   cat('\nModel=MAXENT.Tsuruoka')
  # }
  if (!is.null(seed.val)) {
    seedval = seed.val
  }
  
  ## 5. Check Prev argument ---------------------------------------------------
  if (model == "GLM" | model == "GAM") {
    Prev <- sum(Data[, 1], na.rm = TRUE) / length(Data[, 1])
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
  
  
  return(list(Data = Data,
              weights = weights,
              evalLines = evalLines,
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


.maxent.prepare.workdir <- function(Data, xy, calib.lines = NULL, RunName = NULL,
                                    eval.data = NULL, evalxy =  NULL,
                                    dir.name = '.', species.name = NULL, modeling.id = '',
                                    background_data_dir = 'default')
{
  cat('\n\t\tCreating Maxent Temp Proj Data...')
  
  ## initialise output
  MWD <- list()
  class(MWD) <- "maxent_workdir_info"
  
  ## default parameters setting
  if (is.null(RunName)) { RunName <- colnames(Data)[1] }
  if (is.null(species.name)) { species.name <- colnames(Data)[1] }
  if (is.null(calib.lines)) { calib.lines <- rep(TRUE, nrow(Data)) }
  
  ## define all paths to files needed by MAXENT.Phillips
  nameFolder = file.path(dir.name, species.name, 'models', modeling.id)
  m_outdir <- file.path(nameFolder, paste0(RunName, '_MAXENT.Phillips_outputs'))
  m_predictDir <- file.path(m_outdir, "Predictions")
  MWD$m_outdir <- m_outdir
  MWD$m_outputFile <- file.path(m_outdir, paste0(RunName, '_Pred_swd.csv'))
  MWD$m_predictDir <- m_predictDir
  
  ## directories creation
  dir.create(m_outdir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  dir.create(m_predictDir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  
  ## Presence Data --------------------------------------------------------------------------------
  presLines <- which((Data[, 1] == 1) & calib.lines)
  absLines <- which((Data[, 1] == 0) & calib.lines)
  Sp_swd <- cbind(rep(RunName, length(presLines))
                  , xy[presLines, ]
                  , Data[presLines, 2:ncol(Data), drop = FALSE])
  colnames(Sp_swd) <- c('species', 'X', 'Y', colnames(Data)[2:ncol(Data)])
  
  m_speciesFile <- file.path(m_outdir, "Sp_swd.csv")
  write.table(Sp_swd, file = m_speciesFile, quote = FALSE, row.names = FALSE, sep = ",")
  MWD$m_speciesFile <- m_speciesFile
  
  ## Background Data (create background file only if needed) --------------------------------------
  if (background_data_dir == 'default')
  {
    # keep only 0 of calib lines
    Back_swd <- cbind(rep("background", length(absLines))
                      , xy[absLines, ]
                      , Data[absLines, 2:ncol(Data), drop = FALSE])
    colnames(Back_swd) <- c("background", colnames(Back_swd)[-1])
    
    m_backgroundFile <- file.path(m_outdir, "Back_swd.csv")
    write.table(Back_swd, file = m_backgroundFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_backgroundFile <- m_backgroundFile
  } else { ## use background directory given as an option
    MWD$m_backgroundFile <- background_data_dir
  }
  
  ## Prediction Data ------------------------------------------------------------------------------
  Pred_swd <- cbind(rep("predict", nrow(xy))
                    , xy
                    , Data[, 2:ncol(Data), drop = FALSE])
  colnames(Pred_swd)  <- c("predict", colnames(xy), colnames(Data)[-1])
  
  m_predictFile <- file.path(m_predictDir, "Pred_swd.csv")
  write.table(Pred_swd, file = m_predictFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
  MWD$m_predictFile <- m_predictFile
  
  ## dealing with independent evaluation data -----------------------------------------------------
  if (!is.null(eval.data)) {
    Pred_eval_swd <- cbind(rep("predictEval", nrow(evalxy))
                           , evalxy
                           , eval.data[, 2:ncol(eval.data), drop = FALSE])
    colnames(Pred_eval_swd) <- c("predict", colnames(Back_swd)[-1])
    
    m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
    write.table(Pred_eval_swd, file = m_predictEvalFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_predictEvalFile <- m_predictEvalFile
  }
  
  return(MWD)
}

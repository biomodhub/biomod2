###################################################################################################
##' @name bm_RunModelsLoop
##' @aliases bm_RunModelsLoop
##' @aliases bm_RunModel
##' @author Damien Georges
##' 
##' @title Loop to compute all single species distribution models
##' 
##' @description
##' 
##' This internal \pkg{biomod2} function allows the user to compute all single species 
##' distribution models (asked by the \code{\link{BIOMOD_Modeling}} function).
##' 
##' @param X a \code{BIOMOD.formated.data} or \code{BIOMOD.formated.data.PA} object returned by the 
##' \code{\link{BIOMOD_FormatingData}} function
##' @param modeling.id a \code{character} corresponding to the name (ID) of the simulation set 
##' (\emph{a random number by default})
##' @param Model a \code{character} corresponding to the model name to be computed, must be either 
##' \code{GLM}, \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, 
##' \code{MARS}, \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}
##' @param Options a \code{\link{BIOMOD.models.options}} object returned by the 
##' \code{\link{BIOMOD_ModelingOptions}} function
##' @param VarImport (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} corresponding to the number of permutations to be done for each variable to 
##' estimate variable importance
##' @param mod.eval.method a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param SavePred (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether all results and outputs must be saved on hard drive 
##' or not (\emph{! strongly recommended !})
##' @param scal.models (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether all models predictions must be scaled with a binomial 
##' GLM or not
##' 
##' 
##' @return  
##' 
##' A \code{list} containing for each model a \code{list} containing the following elements :
##' \itemize{
##'   \item{\code{ModelName} : }{the name of correctly computed model}
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
##' 
##' @seealso \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{bm_MakeFormula}}, \code{\link{bm_Rescaler}}, \code{\link{bm_FindOptimStat}}, 
##' \code{\link{bm_VariablesImportance}}
##' 
##' 
##' @keywords models, formula, options, CTA, GLM, GBM, GAM, RF, ANN, FDA, SRE, MARS, MAXENT
##' 
##' 
##' @importFrom rpart rpart prune
## @importFrom caret 
## @importFrom car 
## @importFrom gam gam step.Gam s
## @importFrom mgcv gam bam
##' @importFrom gbm gbm
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


bm_RunModelsLoop <- function(X,
                             modeling.id,
                             Model,
                             Options,
                             VarImport,
                             mod.eval.method,
                             SavePred = TRUE,
                             scal.models = TRUE)
{
  cat("\n\n-=-=-=- Run : ", X$name, '\n')
  res.sp.run <- list()
  
  for (i in 1:ncol(X$calibLines)) { # loop on RunEval
    run.id = dimnames(X$calibLines)[[2]][i]
    run.name = paste0(X$name, run.id)
    cat('\n\n-=-=-=--=-=-=-', run.name, '\n')
    
    res.sp.run[[run.id]] <- lapply(Model,
                                   bm_RunModel,
                                   Data = X$dataBM,
                                   Options = Options,
                                   calibLines = na.omit(X$calibLines[, i, ]), ## transform 3D calibLines obj into a 1D vector
                                   Yweights = na.omit(X$Yweights), 
                                   nam = run.name,
                                   VarImport = VarImport,
                                   mod.eval.method = mod.eval.method,
                                   evalData = X$evalDataBM,
                                   SavePred = TRUE, ## SavePred
                                   xy = X$xy,
                                   eval.xy = X$eval.xy,
                                   scal.models = scal.models,
                                   modeling.id = modeling.id)
    names(res.sp.run[[run.id]]) <- Model
  }
  
  return(res.sp.run)
}


###################################################################################################

##' 
##' @rdname bm_RunModelsLoop
##' @export
##' 

bm_RunModel <- function(Model, Data, Options, calibLines, Yweights, nam, VarImport = 0,
                        mod.eval.method = c('ROC','TSS','KAPPA'), evalData = NULL,
                        SavePred = FALSE,
                        xy = NULL, eval.xy = NULL, scal.models = TRUE, modeling.id = '')
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_RunModel.check.args(Model, Data, Options, calibLines, Yweights, mod.eval.method, evalData, scal.models)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## get model name and names of categorical variables
  model_name <- paste0(nam, '_', Model)
  categorial_var <- unlist(sapply(expl_var_names, function(x) {
    if (is.factor(Data[, x])) { return(x) } else { return(NULL) }
  }))
  
  ## 1. Create output object ----------------------------------------------------------------------
  ListOut <- list(ModelName = NULL,
                  calib.failure = NULL,
                  pred = NULL,
                  pred.eval = NULL,
                  evaluation = NULL,
                  var.import = NULL)
  
  ## 2. CREATE MODELS -----------------------------------------------------------------------------
  
  if (Model == "CTA") {
    ## 2.1 CTA model ------------------------------------------------------------
    cat('\n\t> CTA modeling...')
    
    # converting cost argument
    cost.tmp = Options@CTA$cost
    if (is.null(Options@CTA$cost)) { cost.tmp = rep(1, (ncol(Data) - 2)) }

    # defining rpart parameters for splitting function
    parms.tmp = Options@CTA$parms
    if (Options@CTA$parms == 'default') { parms.tmp = NULL }
    
    model.sp <- try(rpart(
      bm_MakeFormula(respName = colnames(Data)[1]
                     , explVar = head(Data[, -c(1, ncol(Data)), drop = FALSE])
                     , type = 'simple'
                     , interaction.level = 0),
      data = Data[calibLines, ],
      weights = Yweights,
      method = Options@CTA$method,
      parms = parms.tmp,
      cost = cost.tmp,
      control = eval(Options@CTA$control)
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
                      model_options = Options@CTA,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "GAM"){
    ## 2.2 GAM model ------------------------------------------------------------
    
    if (Options@GAM$algo == 'GAM_gam') { ## gam package
      # package loading
      if (isNamespaceLoaded("mgcv")) {
        if (isNamespaceLoaded("caret")) { unloadNamespace("caret") } ## need to unload caret before car
        if (isNamespaceLoaded("car")) { unloadNamespace("car") } ## need to unload car before mgcv
        unloadNamespace("mgcv")
      }
      requireNamespace("gam", quietly = TRUE)
      
      # NOTE : To be able to take into account GAM options and weights we have to do a eval(parse(...))
      # it's due to GAM implementation (using of match.call() troubles)
      gamStart <- eval(parse(text = paste0("gam::gam(", colnames(Data)[1], "~1 ,"
                                           , " data = Data[calibLines,,drop=FALSE], family = ", Options@GAM$family$family
                                           , "(link = '", Options@GAM$family$link, "')"
                                           , ", weights = Yweights[calibLines])")))
      model.sp <- try(gam::step.Gam(gamStart,
                                    .scope(Data[1:3, -c(1, ncol(Data))], "gam::s", Options@GAM$k),
                                    data = Data[calibLines, , drop = FALSE],
                                    direction = "both",
                                    trace = Options@GAM$control$trace,
                                    control = Options@GAM$control))
    } else { ## mgcv package
      # package loading
      if (isNamespaceLoaded("gam")) { unloadNamespace("gam") }
      if (!isNamespaceLoaded("mgcv")) { requireNamespace("mgcv", quietly = TRUE) }
      
      if (is.null(Options@GAM$myFormula)) {
        cat("\n\tAutomatic formula generation...")
        gam.formula <- bm_MakeFormula(respName = resp_name
                                      , explVar = head(Data[, expl_var_names, drop = FALSE])
                                      , type = Options@GAM$type
                                      , interaction.level = Options@GAM$interaction.level
                                      , k = Options@GAM$k)
      } else {
        gam.formula <- Options@GAM$myFormula
      }
      
      if (Options@GAM$algo == 'GAM_mgcv') {
        cat('\n\t> GAM (mgcv) modeling...')
        model.sp <- try(mgcv::gam(gam.formula,
                                  data = Data[calibLines, , drop = FALSE],
                                  family = Options@GAM$family,
                                  weights = Yweights,
                                  control = Options@GAM$control))
        
      } else if (Options@GAM$algo == 'BAM_mgcv') { ## big data.frame gam version
        cat('\n\t> BAM (mgcv) modeling...')
        model.sp <- try(mgcv::bam(gam.formula,
                                  data = Data[calibLines, , drop = FALSE],
                                  family = Options@GAM$family,
                                  weights = Yweights,
                                  control = Options@GAM$control))
      }
    }
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("GAM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GAM',
                      model_subclass = Options@GAM$algo,
                      model_options = Options@GAM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "GBM") {
    ## 2.3 GBM model ------------------------------------------------------------
    
    cat('\n\t> GBM modeling...')
    model.sp <- try(gbm(formula = bm_MakeFormula(respName = colnames(Data)[1]
                                                 , explVar = head(Data)[, expl_var_names, drop = FALSE]
                                                 , type = 'simple'
                                                 , interaction.level = 0),
                        data = Data[calibLines, , drop = FALSE],
                        distribution = Options@GBM$distribution,
                        var.monotone = rep(0, length = ncol(Data) - 2), # -2 because of removing of sp and weights
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
                        n.cores = Options@GBM$n.cores
    ))
    
    if (!inherits(model.sp, "try-error")) {
      best.iter <- try(gbm.perf(model.sp, method = Options@GBM$perf.method , plot.it = FALSE))
      
      model.bm <- new("GBM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GBM',
                      n.trees_optim = best.iter,
                      model_options = Options@GBM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "GLM"){
    ## 2.4 GLM model ------------------------------------------------------------
    
    cat('\n\t> GLM modeling...')
    if (is.null(Options@GLM$myFormula)) {
      cat("\n\tAutomatic formula generation...")
      glm.formula <- bm_MakeFormula(respName = colnames(Data)[1]
                                    , explVar = head(Data)
                                    , type = Options@GLM$type
                                    , interaction.level = Options@GLM$interaction.level)
    } else {
      glm.formula <- Options@GLM$myFormula
    }
    
    if (Options@GLM$test != 'none') {
      ## make the model selection
      glmStart <- glm(eval(parse(text = paste0(colnames(Data)[1], "~1"))), 
                      data = Data[calibLines, , drop = FALSE], 
                      family = Options@GLM$family,
                      control = eval(Options@GLM$control),
                      weights = Yweights[calibLines],
                      mustart = rep(Options@GLM$mustart, sum(calibLines)), 
                      model = TRUE)
      
      ## remove warnings
      warn <- options('warn')
      options(warn = -1)
      model.sp <- try(stepAIC(glmStart,
                              glm.formula,
                              data = Data[calibLines, , drop = FALSE],
                              direction = "both",
                              trace = FALSE,
                              k = criteria,
                              weights = Yweights[calibLines], 
                              steps = 10000,
                              mustart = rep(Options@GLM$mustart, sum(calibLines))))
      ## reexec warnings
      options(warn)
      
    } else {
      ## keep the total model
      model.sp <- try(glm(glm.formula,
                          data = cbind(Data[calibLines, , drop = FALSE], 
                                       matrix(Yweights[calibLines], ncol = 1, dimnames = list(NULL, "Yweights"))), 
                          family = Options@GLM$family,
                          control = eval(Options@GLM$control),
                          weights = Yweights,
                          model = TRUE))
    }
    
    if (!inherits(model.sp, "try-error")) {
      cat("\n\tselected formula : ")
      print(model.sp$formula, useSource = FALSE)
      
      model.bm <- new("GLM_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'GLM',
                      model_options = Options@GLM,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "MARS"){
    ## 2.5 MARS model -----------------------------------------------------------
    
    cat('\n\t> MARS modeling...')
    if (is.null(Options@MARS$myFormula)) {
      cat("\n\tAutomatic formula generation...")
      mars.formula <- bm_MakeFormula(respName = colnames(Data)[1]
                                     , explVar = head(Data)[, -ncol(Data), drop = FALSE]
                                     , type = Options@MARS$type
                                     , interaction.level = Options@MARS$interaction.level)
    } else {
      mars.formula <- Options@MARS$myFormula
    }
    
    ## deal with nk argument : if not defined, set up to default mars value i.e max(21, 2 * ncol(x) + 1)
    nk <- Options@MARS$nk
    if (is.null(nk)) {
      nk <- min(200, max(20, 2 * length(expl_var_names))) + 1
    }
    
    model.sp <- try(earth(formula = mars.formula,
                          data = Data[calibLines, , drop = FALSE], 
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
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("MARS_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MARS',
                      model_options = Options@MARS,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "FDA") {
    ## 2.6 FDA model ------------------------------------------------------------
    
    cat('\n\t> FDA modeling...')
    model.sp <- try(do.call(fda, c(list(formula = bm_MakeFormula(respName = colnames(Data)[1]
                                                                 , explVar = head(Data)[, expl_var_names, drop = FALSE]
                                                                 , type = 'simple'
                                                                 , interaction.level = 0),
                                        data = Data[calibLines, , drop = FALSE], 
                                        method = eval(parse(text = call(Options@FDA$method))),
                                        weights = Yweights[calibLines]),
                                   Options@FDA$add_args)))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("FDA_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'FDA',
                      model_options = Options@FDA,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]),
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "ANN") {
    ## 2.7 ANN model ------------------------------------------------------------
    
    cat('\n\t> ANN modeling...')
    size = Options@ANN$size
    decay = Options@ANN$decay
    if (is.null(size) | is.null(decay) | length(size) > 1 | length(decay) > 1) {
      
      ## define the size and decay to test
      if (is.null(size)) { size <- c(2, 4, 6, 8) }
      if (is.null(decay)) { decay <- c(0.001, 0.01, 0.05, 0.1) }
      
      ## do cross validation test to find the optimal values of size and decay parameters (prevent from overfitting)
      CV_nnet <- bm_CVnnet(Input = Data[, expl_var_names, drop = FALSE],
                           Target = Data[calibLines, 1], 
                           size = size,
                           decay = decay,
                           maxit = Options@ANN$maxit,
                           nbCV = Options@ANN$NbCV,
                           W = Yweights[calibLines])
      
      ## get the optimised parameters values
      decay <- CV_nnet[1, 2]
      size <- CV_nnet[1, 1]
    }
    
    model.sp <- try(nnet(formula = bm_MakeFormula(respName = resp_name
                                                  , explVar = head(Data[, expl_var_names, drop = FALSE])
                                                  , type = 'simple'
                                                  , interaction.level = 0),
                         data = Data[calibLines, , drop = FALSE], 
                         size = size,
                         rang = Options@ANN$rang,
                         decay = decay,
                         weights = Yweights,
                         maxit = Options@ANN$maxit,
                         trace = FALSE))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("ANN_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'ANN',
                      model_options = Options@ANN,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "RF") {
    ## 2.8 RF model -------------------------------------------------------------
    
    cat('\n\t> RF modeling...')
    if (Options@RF$do.classif) {
      # defining occurences as factor for doing classification and not regression in RF
      Data <- Data %>% mutate_at(resp_name, factor)
    }
    
    # mtry.tmp = Options@RF$mtry
    # if (Options@RF$mtry == 'default') { mtry.tmp = NULL }

    model.sp <- try(randomForest(formula = bm_MakeFormula(respName = resp_name
                                                          , explVar = head(Data)
                                                          , type = 'simple'
                                                          , interaction.level = 0),
                                 data = Data[calibLines, ],
                                 ntree = Options@RF$ntree,
                                 # mtry = mtry.tmp, 
                                 importance = FALSE,
                                 norm.votes = TRUE,
                                 strata = factor(c(0, 1)),
                                 nodesize = Options@RF$nodesize,
                                 maxnodes = Options@RF$maxnodes))

    if (Options@RF$do.classif) {
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
                      model_options = Options@RF,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "SRE") {
    ## 2.9 SRE model ------------------------------------------------------------
    
    cat('\n\t> SRE modeling...')
    model.sp <- try(bm_SRE(Response = Data[calibLines, 1],
                           Explanatory = Data[calibLines, expl_var_names, drop = FALSE],
                           NewData = NULL,
                           Quant = Options@SRE$quant,
                           return_extremcond = TRUE))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("SRE_biomod2_model",
                      extremal_conditions = model.sp,
                      model_name = model_name,
                      model_class = 'SRE',
                      model_options = Options@SRE,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  } else if (Model == "MAXENT.Phillips") {
    ## 2.10 MAXENT.Phillips model -----------------------------------------------
    
    cat('\n\t> MAXENT.Phillips modeling...')
    MWD <- .maxent.prepare.workdir(Data, xy, calibLines, RunName = nam,
                                   evalData, eval.xy, species.name = resp_name,
                                   modeling.id = modeling.id,
                                   background_data_dir = Options@MAXENT.Phillips$background_data_dir)
    
    maxent.cmd <- paste0("java ",
                         ifelse(is.null(Options@MAXENT.Phillips$memory_allocated),
                                "",
                                paste0("-mx", Options@MAXENT.Phillips$memory_allocated, "m")), 
                         " -jar ", file.path(Options@MAXENT.Phillips$path_to_maxent.jar, "maxent.jar"),
                         " environmentallayers=\"", MWD$m_backgroundFile, 
                         "\" samplesfile=\"", MWD$m_speciesFile,
                         "\" projectionlayers=\"", gsub(", ", ",", toString(MWD$m_predictFile)),
                         "\" outputdirectory=\"", MWD$m_outdir, "\"",
                         " outputformat=logistic ", 
                         ifelse(length(categorial_var),
                                paste0(" togglelayertype=", categorial_var, collapse = " "),
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
                    expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                    expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    
    # for MAXENT.Phillips predicitons are calculated in the same time than models building to save time.
    cat("\n Getting predictions...")
    g.pred <- try(round(as.numeric(read.csv(MWD$m_outputFile)[, 3]) * 1000))
    
    cat("\n Getting predictor contributions...")
    variables.importance <- bm_VariablesImportance(model = model.bm
                                                   , data = Data[, expl_var_names, drop = FALSE]
                                                   , nb_rand = VarImport
                                                   , temp_workdir = MWD$m_outdir)
  } else if(Model == "MAXENT.Phillips.2")
  {
    ## 2.11 MAXENT.Phillips.2 model ---------------------------------------------
    
    cat('\n\t> MAXENT.Phillips modeling...')
    model.sp <- try(maxnet(p = Data %>% filter(calibLines) %>% pull(resp_name), 
                           data = Data %>% filter(calibLines) %>% select_at(expl_var_names)
                           # f = if(!is.null(Options@MAXENT.Phillips.2@))
    ))
    
    if (!inherits(model.sp, "try-error")) {
      model.bm <- new("MAXENT.Phillips.2_biomod2_model",
                      model = model.sp,
                      model_name = model_name,
                      model_class = 'MAXENT.Phillips.2',
                      model_options = Options@MAXENT.Phillips.2,
                      resp_name = resp_name,
                      expl_var_names = expl_var_names,
                      expl_var_type = get_var_type(Data[calibLines, expl_var_names, drop = FALSE]), 
                      expl_var_range = get_var_range(Data[calibLines, expl_var_names, drop = FALSE]))
    }
  }
  
  ## 2.12 MAXENT.Tsuruoka model -----------------------------------------------
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
  
  
  ## 3. CREATE PREDICTIONS ------------------------------------------------------------------------
  if (Model != "MAXENT.Phillips") {
    g.pred <- try(predict(model.bm, Data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE))
  }
  
  ## scale or not predictions -------------------------------------------------
  if (scal.models & !inherits(g.pred, 'try-error')) {
    cat("\n\tModel scaling...")
    model.bm@scaling_model <- try(.scaling_model(g.pred / 1000, Data[, 1, drop = TRUE], weights = Yweights))
    ## with weights
    g.pred <- try(predict(model.bm, Data[, expl_var_names, drop = FALSE], on_0_1000 = TRUE))
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
    ListOut$ModelName <- model_name
  } else {
    # keep the name of uncompleted modelisations
    cat("\n   ! Note : ", model_name, "failed!\n")
    ListOut$calib.failure = model_name
    return(ListOut) ## end of function.
  }
  
  ## make prediction on evaluation data ---------------------------------------
  if (!is.null(evalData)) {
    g.pred.eval <- try(predict(model.bm, evalData[, expl_var_names, drop = FALSE], on_0_1000 = TRUE))
  }
  
  ## SAVE predictions ---------------------------------------------------------
  if (SavePred) {
    ListOut$pred <- g.pred
    if (exists("g.pred.eval")) { ListOut$pred.eval <- g.pred.eval }
  }
  
  
  ## 4. EVALUATE MODEL ----------------------------------------------------------------------------
  if (length(mod.eval.method) > 0) {
    cat("\n\tEvaluating Model stuff...")
    
    ## Check no NA in g.pred to avoid evaluation failures
    na_cell_id <- which(is.na(g.pred))
    if (length(na_cell_id)) {
      evalLines <- evalLines[!(evalLines %in% na_cell_id)]
      cat('\n\tNote : some NA occurs in predictions')
    }
    
    cross.validation <- sapply(mod.eval.method, function(.x) {
      bm_FindOptimStat(Stat = .x,
                       Fit = g.pred[evalLines],
                       Obs = Data %>% filter(evalLines) %>% pull(1))
    })
    rownames(cross.validation) <- c("Testing.data", "Cutoff", "Sensitivity", "Specificity")
    
    if (exists('g.pred.eval')) {
      
      ## Check no NA in g.pred.eval to avoid evaluation failures
      na_cell_id <- which(is.na(g.pred.eval))
      if (length(na_cell_id)) {
        g.pred.eval.without.na <- g.pred.eval[-na_cell_id]
        evalData <- evalData[-na_cell_id, ]
        cat('\n\tNote : some NA occurs in evaluation predictions')
      } else {
        g.pred.eval.without.na <- g.pred.eval
      }
      
      true.evaluation <- sapply(mod.eval.method, function(x) {
        bm_FindOptimStat(Stat = x,
                         Fit = g.pred.eval.without.na,
                         Obs = evalData[, 1],
                         Fixed.thresh = cross.validation["Cutoff", x])
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
  if (VarImport > 0) {
    cat("\n\tEvaluating Predictor Contributions...")
    if (Model != "MAXENT.Phillips") {
      variables.importance <- bm_VariablesImportance(model = model.bm
                                                     , data = Data[, expl_var_names, drop = FALSE]
                                                     , nb_rand = VarImport)
    }
    model.bm@model_variables_importance <- variables.importance
    
    ## only the mean of variables importance run is returned
    ListOut$var.import <- round(rowMeans(variables.importance, na.rm = TRUE), digits = 3)
    rm(variables.importance)
    cat("\n")
  }
  
  
  ## 6. SAVE MODEL OBJECT ON HARD DRIVE -----------------------------------------------------------
  nameModel = paste(nam, Model, sep = "_") 
  assign(x = nameModel, value = model.bm)
  save(list = nameModel, file = file.path(resp_name, "models", modeling.id, nameModel), compress = TRUE)
  
  
  return(ListOut)
}


###################################################################################################

.bm_RunModel.check.args <- function(Model, Data, Options, calibLines, Yweights, mod.eval.method
                                    , evalData, scal.models, criteria = NULL, Prev = NULL)
{
  ## 0. Do some cleaning over Data argument -----------------------------------
  resp_name <- colnames(Data)[1] ## species name
  expl_var_names <- colnames(Data)[-1] ## explanatory variable names
  # replace Pseudo absences selected (NA) into true absences (0).. for model computing purpose
  if (sum(is.na(Data[, 1]))) { Data[which(is.na(Data[, 1])), 1] < - 0 }
  
  ## 1. Check CalibLines argument ---------------------------------------------
  if (sum(!calibLines) > 0) ## if some lines for evaluation...
  {
    evalLines <- !calibLines
    # ...test if there is absences AND presences in evaluation and calibration datasets
    if (sum(Data[calibLines, 1] == 0) == 0 ||
        sum(Data[calibLines, 1] == 0) == sum(calibLines) ||
        sum(Data[evalLines, 1] == 0) == 0 ||
        sum(Data[evalLines, 1] == 0) == sum(evalLines)) {
      warning(paste0(colnames(Data)[1], " ", Model,
                     " was switched off because of no both presences and absences data given"),
              immediate. = TRUE)
      return(NULL)
    }
  } else { ## evaluation = calibration dataset
    evalLines <- calibLines
    # ...test if there is absences AND presences in whole dataset
    if (sum(Data[, 1] == 0) == 0 ||
        sum(Data[, 1] == 0) == nrow(Data)) {
      warning(paste0(colnames(Data)[1], " ", Model,
                     " was switched off because of no both presences and absences data given (full model)"),
              immediate. = TRUE)
      return(NULL)
    }
  }
  
  ## 2. Check Yweights argument -----------------------------------------------
  if (is.null(Yweights)) { Yweights <- rep(1, nrow(Data)) }
  ## These models require data and weights to be in the same dataset
  if (Model %in% c('GBM', 'CTA', 'ANN', 'FDA', 'GAM', 'MARS')) {
    Data <- cbind(Data, Yweights)
  }
  
  ## 3. Check scal.models argument --------------------------------------------
  if (Model == "SRE") { scal.models <- FALSE } else if (Model %in% c("ANN", "FDA")) { scal.models <- TRUE }
  
  
  ## 4. Check Options argument ------------------------------------------------
  if (Model == "GLM") {
    cat('\nModel=GLM')
    if (!is.null(Options@GLM$myFormula)) {
      cat('\n\tformula = ', paste(Options@GLM$myFormula[2],
                                  Options@GLM$myFormula[1],
                                  Options@GLM$myFormula[3]))
    } else {
      cat(' (', Options@GLM$type, 'with',
          ifelse(Options@GLM$interaction.level == 0,
                 'no interaction )',
                 paste('order', Options@GLM$interaction.level, 'interaction level )')
          ))
    }
    if (Options@GLM$test == "AIC") {
      criteria <- 2
      cat("\n\tStepwise procedure using AIC criteria")
    } else if (Options@GLM$test == "BIC") {
      criteria <- log(ncol(Data))
      cat("\n\tStepwise procedure using BIC criteria")
    } else if (Options@GLM$test == "none") {
      criteria <- 0
      cat("\n\tNo stepwise procedure")
      cat("\n\t! You might be confronted to model convergence issues !")
    }
  } else if (Model == "GBM") {
    cat("\nModel=Generalised Boosting Regression \n")
    cat("\t", Options@GBM$n.trees, "maximum different trees and ", Options@GBM$cv.folds, " Fold Cross-Validation")
    set.seed(456) # to be able to refind our trees MAY BE BAD
  } else if (Model == "GAM") {
    cat("\nModel=GAM")
    cat("\n\t", Options@GAM$algo, "algorithm chosen")
  } else if (Model == "CTA") {
    cat("\nModel=Classification tree \n")
    cat("\t", Options@CTA$control$xval, "Fold Cross-Validation")
    set.seed(123) # to be able to refind our trees MAY BE BAD
  } else if (Model == "ANN") {
    cat("\nModel=Artificial Neural Network \n")
    cat("\t", Options@ANN$NbCV, "Fold Cross Validation + 3 Repetitions")
    set.seed(555) # to be able to refind our trees MAY BE BAD
  } else if (Model == "SRE") {
    cat("\nModel=Surface Range Envelop")
  } else if (Model == "FDA"){
    cat("\nModel=Flexible Discriminant Analysis")
  } else if (Model == "MARS"){
    cat("\nModel=Multiple Adaptive Regression Splines")
    if (!is.null(Options@MARS$myFormula)) {
      cat('\n\tformula = ', paste(Options@MARS$myFormula[2],
                                  Options@MARS$myFormula[1],
                                  Options@MARS$myFormula[3]))
    } else {
      cat(' (', Options@MARS$type, 'with',
          ifelse(Options@MARS$interaction.level == 0,
                 'no interaction )',
                 paste('order', Options@MARS$interaction.level, 'interaction level )')
          ))
    }
    cat("\n")
  } else if (Model == "RF") {
    cat("\nModel=Breiman and Cutler's random forests for classification and regression")
    set.seed(71)
  } else if (Model == 'MAXENT.Phillips') {
    cat('\nModel=MAXENT.Phillips')
  } else if (Model == 'MAXENT.Phillips.2') {
    cat('\nModel=MAXENT.Phillips (maxnet)')
  }
  # else if (Model == 'MAXENT.Tsuruoka') {
  #   cat('\nModel=MAXENT.Tsuruoka')
  # }
  
  ## 5. Check Prev argument ---------------------------------------------------
  if (Model == "GLM" | Model == "GAM") {
    Prev <- sum(Data[, 1], na.rm = T) / length(Data[, 1])
  }
  
  ## 6. Check models.eval.meth arguments --------------------------------------
  mod.eval.method <- unique(mod.eval.method)
  avail.eval.meth.list <- c('TSS', 'KAPPA', 'ACCURACY', 'BIAS', 'POD', 'FAR', 'POFD'
                            , 'SR', 'CSI', 'ETS', 'HK', 'HSS', 'OR', 'ORSS', 'ROC')
  # .fun_testIfIn(TRUE, "mod.eval.method", mod.eval.method, avail.eval.meth.list)
  if (sum(!(mod.eval.method %in% avail.eval.meth.list)) > 0) {
    tmp = which(mod.eval.method %in% avail.eval.meth.list)
    warnings(paste0(toString(mod.eval.method[!tmp]), ' were switched off !'), imediate = TRUE)
    mod.eval.method <- mod.eval.method[tmp]
  }
  
  
  return(list(Data = Data,
              Yweights = Yweights,
              evalLines = evalLines,
              criteria = criteria,
              Prev = Prev, 
              mod.eval.method = mod.eval.method,
              evalData = evalData,
              scal.models = scal.models,
              resp_name = resp_name,
              expl_var_names = expl_var_names))
}


.maxent.prepare.workdir <- function(Data, xy, calibLines = NULL, RunName = NULL,
                                    evalData = NULL, evalxy =  NULL,
                                    species.name = NULL, modeling.id = '',
                                    background_data_dir = 'default')
{
  cat('\n\t\tCreating Maxent Temp Proj Data...')
  
  ## initialise output
  MWD <- list()
  class(MWD) <- "maxent_workdir_info"
  
  ## default parameters setting
  if (is.null(RunName)) { RunName <- colnames(Data)[1] }
  if (is.null(species.name)) { species.name <- colnames(Data)[1] }
  if (is.null(calibLines)) { calibLines <- rep(TRUE, nrow(Data)) }
  
  ## define all paths to files needed by MAXENT.Phillips
  nameFolder = file.path(species.name, 'models', modeling.id)
  m_outdir <- file.path(nameFolder, paste0(RunName, '_MAXENT.Phillips_outputs'))
  m_predictDir <- file.path(m_outdir, "Predictions")
  MWD$m_outdir <- m_outdir
  MWD$m_outputFile <- file.path(m_outdir, paste0(RunName, '_Pred_swd.csv'))
  MWD$m_predictDir <- m_predictDir
  
  ## directories creation
  dir.create(m_outdir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  dir.create(m_predictDir, showWarnings = FALSE, recursive = TRUE, mode = '777')
  
  ## Presence Data --------------------------------------------------------------------------------
  presLines <- which((Data[, 1] == 1) & calibLines)
  absLines <- which((Data[, 1] == 0) & calibLines)
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
  if (!is.null(evalData)) {
    Pred_eval_swd <- cbind(rep("predictEval", nrow(evalxy))
                           , evalxy
                           , evalData[, 2:ncol(evalData), drop = FALSE])
    colnames(Pred_eval_swd) <- c("predict", colnames(Back_swd)[-1])
    
    m_predictEvalFile <- file.path(m_predictDir, "PredEval_swd.csv")
    write.table(Pred_eval_swd, file = m_predictEvalFile, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
    MWD$m_predictEvalFile <- m_predictEvalFile
  }
  
  return(MWD)
}

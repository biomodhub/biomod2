##' @name BIOMOD_tuning
##' @aliases BIOMOD_tuning
##' 
##' @title Tune models parameters
##' @description Function to tune biomod single models parameters
##'
##' @param data            BIOMOD.formated.data object returned by BIOMOD_FormatingData
##' @param models          vector of models names choosen among 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips' 
##' @param models.options  BIOMOD.models.options object returned by BIOMOD_ModelingOptions. Default: BIOMOD_ModelingOptions()
##' @param trControl       global control parameters for runing (default trainControl(method="cv",summaryFunction=twoClassSummary,classProbs=T),returnData = FALSE). for details see trainControl
##' @param ctrl.GAM        specify control parameters only for GAM (default trControl)
##' @param ctrl.GBM        specify control parameters only for GBM (default trControl)
##' @param ctrl.GLM        specify control parameters only for GLM (default trControl)
##' @param ctrl.CTA        specify control parameters only for CTA (default trControl)
##' @param ctrl.RF         specify control parameters only for RF (default trControl)
##' @param ctrl.ANN        specify control parameters only for ANN (default trControl)
##' @param ctrl.MARS       specify control parameters only for MARS (default trControl)
##' @param ctrl.FDA        specify control parameters only for FDA (default trControl)
##' @param metric          metric to select the optimal model (Default ROC). TSS (maximizing Sensitivity and Specificity) is also possible. see ?train
##' @param metric.ME       metric to select the optimal model for MAXENT.Phillips (Default: ROC). One out of Mean.AUC (or ROC), Mean.AUC.DIFF, Mean.ORmin, Mean.OR10 and AICc. see ?ENMevaluate and Muscarella et al. 2014
##' @param tuneLength      see ?train (default 30)
##' @param method.RF       which classification or regression model to use for randomForest (default: "rf"). see http://topepo.github.io/caret/Random_Forest.html
##' @param method.ANN      which classification or regression model to use for artificial neural networks (default: "avNNet"). see http://topepo.github.io/caret/Neural_Network.html
##' @param method.MARS     which classification or regression model to use for mars (default: "earth"). see http://topepo.github.io/caret/Multivariate_Adaptive_Regression_Splines.html
##' @param method.GAM      which classification or regression model to use for GAM (default: "gam"). see http://topepo.github.io/caret/Generalized_Additive_Model.html
##' @param method.GLM      which classification or regression model to use for GLM: (default: 'glmStepAIC'). see http://topepo.github.io/caret/Generalized_Linear_Model.html
##' @param type.GLM        vector of modeling types choosen among 'simple', 'quadratic', 'polynomial' or 's_smoother' (default c('simple','quadratic','polynomial','s_smoother'))
##' @param interaction.GLM vector of interaction type choosen among 0, 1. Default c(0,1)
##' @param cvmethod.ME     method used for data partitioning for MAXENT.Phillips (default: 'randomkfold')
##' @param kfolds.ME       number of bins to use for k-fold cross-validation used for MAXENT.Phillips (Default: 10).
##' @param overlap.ME      logical; Calculates pairwise metric of niche overlap if TRUE (Default: FALSE). (see ?calc.niche.overlap)
##' @param clamp.ME        logical; If TRUE (Default) "Features are constrained to remain within the range of values in the training data" (Elith et al. 2011)
##' @param n.bg.ME         Number of Background points used to run MAXENT.Phillips (Default: 10000)
##' @param env.ME          RasterStack of model predictor variables
##' @param size.tune.ANN   size parameters (number of units in the hidden layer) for ANN used for tuning (default: c(2,4,6,8)).  Will be optimised using the method specified in ctrl.ANN (if not available trControl).
##' @param decay.tune.ANN  weight decay parameters used for tuning for ANN (default: c(0.001, 0.01, 0.05, 0.1))  Will be optimised by method specified in ctrl.ANN (if not available trControl).
##' @param maxit.ANN       maximum number of iterations for ANN (default 500) 
##' @param MaxNWts.ANN     The maximum allowable number of weights for ANN (default 10 * (ncol(myBiomodData'at'data.env.var) + 1) + 10 + 1). 
##' @param parallel.ME     logical. If TRUE, the parallel computing is enabled for MAXENT.Phillips 
##' @param numCores.ME     number of cores used to train MAXENT.Phillips 
##' @param Yweights        response points weights. This argument will only affect models that allow case weights. 
##' 
##' @return
##' BIOMOD.models.options object with optimized parameters
##' 
##' @author Frank Breiner \email{frank.breiner@wsl.ch}
##' 
##' @references 
##' Kuhn, Max. 2008. Building predictive models in R using the caret package. \emph{Journal of Statistical Software} \bold{28}, 1-26.
##' Kuhn, Max, and Kjell Johnson. 2013. Applied predictive modeling. New York: Springer.
##' Muscarella, Robert, et al. 2014. ENMeval: An R package for conducting spatially independent evaluations and estimating optimal model complexity for Maxent ecological niche models. \emph{Methods in Ecology and Evolution}, \bold{5}, 1198-1205.
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_ModelingOptions}}, \code{\link[caret]{train}}, \code{\link[ENMeval]{ENMevaluate}}, 
##' 
##' @examples
##' \dontrun{
##' # species occurrences
##' DataSpecies <- read.csv(system.file("external/species/mammals_table.csv",
##'                                     package="biomod2"))
##' head(DataSpecies)
##' 
##' # the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' # the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' # the XY coordinates of species data
##' myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]
##' 
##' # Environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
##'                              package="biomod2"),
##'                 system.file( "external/bioclim/current/bio4.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio7.grd", 
##'                              package="biomod2"),  
##'                 system.file( "external/bioclim/current/bio11.grd", 
##'                              package="biomod2"), 
##'                 system.file( "external/bioclim/current/bio12.grd", 
##'                              package="biomod2"))
##' # 1. Formatting Data
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # 2. Defining Models Options using default options.
##' ### Duration for turing all models sequential with default settings 
##' ### on 3.4 GHz processor: approx. 45 min tuning all models in parallel
##' ### (on 8 cores) using foreach loops runs much faster: approx. 14 min
##' 
##' #library(doParallel);cl<-makeCluster(8);doParallel::registerDoParallel(cl) 
##' 
##' 
##' time.seq<-system.time(Biomod.tuning <- BIOMOD_tuning(myBiomodData,
##'                                                              env.ME = myExpl,
##'                                                              n.bg.ME = ncell(myExpl)))
##' #stopCluster(cl)
##' 
##' myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
##'                                      models = c('RF','CTA'), 
##'                                      models.options = Biomod.tuning$models.options, 
##'                                      NbRunEval=1, 
##'                                      DataSplit=100, 
##'                                      VarImport=0, 
##'                                      models.eval.meth = c('ROC'),
##'                                      do.full.models=FALSE,
##'                                      modeling.id="test")
##' 
##' 
##' #  eval.plot(Biomod.tuning$tune.MAXENT.Phillips at results)
##' par(mfrow=c(1,3))
##' plot(Biomod.tuning$tune.CTA.rpart)
##' plot(Biomod.tuning$tune.CTA.rpart2)
##' plot(Biomod.tuning$tune.RF)
##' }
BIOMOD_tuning <- function(data,
                          models = c('GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','RF','MAXENT.Phillips'),
                          models.options = BIOMOD_ModelingOptions(),
                          method.ANN = 'avNNet',
                          method.RF = 'rf',
                          method.MARS = 'earth',
                          method.GAM = 'gam',
                          method.GLM = 'glmStepAIC',
                          trControl = NULL,
                          metric = 'ROC',
                          ctrl.CTA = NULL,
                          ctrl.RF = NULL,
                          ctrl.ANN = NULL,
                          ctrl.MARS = NULL,                                  
                          ctrl.FDA = NULL,
                          ctrl.GAM = NULL,
                          ctrl.GBM = NULL,
                          ctrl.GLM = NULL,
                          tuneLength = 30,
                          decay.tune.ANN = c(0.001, 0.01, 0.05, 0.1),
                          size.tune.ANN = c(2,4,6,8),
                          maxit.ANN = 500,
                          MaxNWts.ANN = 10 * (ncol(data@data.env.var) + 1) + 10 + 1,
                          type.GLM = c('simple','quadratic','polynomial','s_smoother'),
                          interaction.GLM = c(0,1),
                          cvmethod.ME = 'randomkfold',
                          overlap.ME = FALSE,
                          kfolds.ME = 10,
                          n.bg.ME = 10000,
                          env.ME = NULL,
                          metric.ME = 'ROC',
                          clamp.ME = TRUE,
                          parallel.ME = FALSE,
                          numCores.ME = NULL,
                          Yweights = NULL){
  
  ## MAXENT: http://cran.r-project.org/web/packages/ENMeval/ENMeval.pdf --> ENMevaluate()
  ## or:    http://cran.r-project.org/web/packages/maxent/maxent.pdf -->  tune.maxent()
  #packages <- NULL
  if(sum(c('GLM','GBM','GAM','CTA','ANN','FDA','MARS','RF','MAXENT.Phillips','SRE') %in% models)>0){if(!isNamespaceLoaded("caret")){requireNamespace("caret", quietly = TRUE)}; 
    if(!isNamespaceLoaded('dplyr')){requireNamespace("dplyr", quietly = TRUE)}; 
    if(is.null(trControl)){trControl <- caret::trainControl(method="cv",summaryFunction=caret::twoClassSummary,classProbs=T, returnData = F)}}
  if("MAXENT.Phillips" %in% models){if(!isNamespaceLoaded('ENMeval')){requireNamespace("ENMeval", quietly = TRUE)}}#;packages<-c(packages,"ENMeval")}  
  # if("MAXENT.Tsuruoka" %in% models){if(!isNamespaceLoaded('maxent')){requireNamespace("maxent", quietly = TRUE)}}#;packages<-c(packages,"maxent")}  
  
  tune.SRE <- tune.GLM <- tune.MAXENT.Phillips <- tune.GAM <- tune.GBM <- tune.CTA.rpart <- tune.CTA.rpart2 <- tune.RF <- tune.ANN <- tune.MARS <- tune.FDA <- NULL
    # tune.MAXENT.Tsuruoka <- NULL
    
  
  resp <- data@data.species
  
  if('SRE' %in% models){
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n",
              "Start tuning SRE\n"))
    #if(trControl$method=="cv"){
    tune.SRE <- NULL
    for (rep in 1:trControl$repeats){
      fold <- dismo::kfold(resp, by = resp, 
                           k = trControl$number)
      for(quant in c(0,0.0125,0.025,0.05,0.1)){
        for (i in 1:trControl$number) {
          DATA <- cbind(1:sum(fold==i),
                        resp[fold==i],
                        sre(Response = resp[fold!=i], 
                            Explanatory = data@data.env.var[fold!=i,], 
                            NewData = data@data.env.var[fold==i,], 
                            Quant=quant, 
                            return_extremcond = FALSE))
          tune.SRE <- rbind(tune.SRE,cbind(
            presence.absence.accuracy(DATA, 
                                      threshold=as.vector(PresenceAbsence::optimal.thresholds(DATA,opt.methods=3)[2],mode="numeric")),quant))
        }
      }}
    t<-aggregate(tune.SRE,by=list(quant = tune.SRE$quant),mean)
    if(metric == 'ROC'){models.options@SRE$quant<-t[which.max(t$AUC),"quant"]} 
    if(metric == 'TSS'){models.options@SRE$quant<-t[which.max(t$sensitivity+t$specificity-1),"quant"]} 
    cat(paste("Finished tuning SRE\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }#}
  
  if(metric == 'ROC' | metric == 'TSS'){resp <- as.factor(ifelse(resp == 1 & !is.na(resp), "Presence", "Absence"))}
  
  if('GBM' %in% models){  
    
    if(is.null(ctrl.GBM)){ctrl.GBM <- trControl}
    cat(paste("\n-=-=-=-=-=-=-=-=-=-=\n",
              "Start tuning GBM. Start coarse tuning\n"))
    
    tune.grid <- expand.grid(.interaction.depth = seq(2, 8, by = 3),
                             .n.trees = c(500, 1000, 2500),
                             .shrinkage = c(0.001, 0.01, 0.1),
                             .n.minobsinnode = 10)
    
    try(tune.GBM <- caret::train(data@data.env.var, resp,
                                 method = "gbm",
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.GBM,
                                 verbose = FALSE,
                                 weights = Yweights))
    cat("Best optimization of coarse tuning:\n")
    cat(paste(tune.GBM$bestTune,"\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.GBM)){
      cat("Start fine tuning\n")
      
      if(tune.GBM$bestTune$n.trees==2500){
        cat("Best optimization with large trees! Tuning GBM will take a while.\n")
        n.trees <- seq(2500, 10000, by = 2500)}
      
      if(tune.GBM$bestTune$n.trees==1000){n.trees <- seq(750, 2000, by = 250)}
      
      if(tune.GBM$bestTune$n.trees==500){n.trees <- seq(100, 1000, by = 50)}
      
      tune.grid <- expand.grid(.interaction.depth = c(tune.GBM$bestTune$interaction.depth-1,tune.GBM$bestTune$interaction.depth,tune.GBM$bestTune$interaction.depth+1),
                               .n.trees = n.trees,
                               .shrinkage = c(tune.GBM$bestTune$shrinkage/2,tune.GBM$bestTune$shrinkage,tune.GBM$bestTune$shrinkage*5),
                               .n.minobsinnode = 10)
      tune.GBM <- NULL
      try(tune.GBM <- caret::train(data@data.env.var, resp,
                                   method = "gbm",
                                   tuneGrid = tune.grid,
                                   trControl = ctrl.GBM,
                                   verbose = FALSE,
                                   weights = Yweights))  
    }
    cat(paste("\n Finished tuning GBM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.GBM)){
      if(metric == 'TSS'){
        models.options@GBM$n.trees <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"n.trees"]
        models.options@GBM$interaction.depth <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"interaction.depth"]    
        models.options@GBM$shrinkage <- tune.GBM$results[which.max(apply(tune.GBM$results[,c("Sens","Spec")],1,sum)-1),"shrinkage"] 
      }else{
        models.options@GBM$n.trees <- tune.GBM$bestTune$n.trees
        models.options@GBM$interaction.depth <- tune.GBM$bestTune$interaction.depth    
        models.options@GBM$shrinkage <- tune.GBM$bestTune$shrinkage 
      }}else{ if('GBM' %in% models){cat("Tuning GBM failed!"); tune.GBM <- "FAILED"}}
    
  }
  
  if('RF' %in% models){
    cat("Start tuning RF\n")
    
    if(is.null(ctrl.RF)){ctrl.RF <- trControl}
    tuneLength.rf <- min(tuneLength,ncol(data@data.env.var))
    
    ## give both mtry as bestTune
    try(tune.RF <- caret::train(data@data.env.var, resp,
                                method = method.RF,
                                tuneLength = tuneLength.rf,
                                trControl = ctrl.RF,
                                metric = metric,
                                weights = Yweights))
    cat(paste("Finished tuning RF\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    if(!is.null(tune.RF)){
      if(metric == 'TSS'){
        models.options@RF$mtry <- tune.RF$results[which.max(apply(tune.RF$results[,c("Sens","Spec")],1,sum)-1),"mtry"]
      }else{
        models.options@RF$mtry <- tune.RF$bestTune$mtry
      }}else{ if('RF' %in% models){cat("Tuning RF failed!"); tune.RF <- "FAILED"}}
  }
  
  if('ANN' %in% models){
    cat("Start tuning ANN\n")
    
    if(is.null(ctrl.ANN)){ctrl.ANN <- trControl}
    ## already tuning: 
    # size: optimised by cross validation based on model AUC (NbCv cross validation; tested size will be the following c(2,4,6, 8))
    # decay: optimised by cross validation on model AUC (NbCv cross validation; tested decay will be the following c(0.001, 0.01, 0.05, 0.1)).
    # could increase maxit from 200 to 500
    # a nice option would be to use model averaging for ann: avNNet in package(caret)
    
    ## Create a specific candidate set of models to evaluate:
    tune.grid <- expand.grid(.decay = decay.tune.ANN,
                             .size = size.tune.ANN,
                             .bag = FALSE)
    
    try(tune.ANN <- caret::train(data@data.env.var, resp, 
                                 method = method.ANN,
                                 tuneGrid = tune.grid,
                                 trControl = ctrl.ANN,
                                 ## Automatically standardize data prior to modeling
                                 ## and prediction
                                 preProc = c("center", "scale"),
                                 linout = TRUE,
                                 trace = FALSE,
                                 MaxNWts.ANN = MaxNWts.ANN,
                                 maxit = maxit.ANN,
                                 metric = metric,
                                 weights = Yweights))
    cat(paste("Finished tuning ANN\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    if(!is.null(tune.ANN)){
      if(metric == 'TSS'){
        models.options@ANN$size <- tune.ANN$results[which.max(apply(tune.ANN$results[,c("Sens","Spec")],1,sum)-1),"size"]    
        models.options@ANN$decay <- tune.ANN$results[which.max(apply(tune.ANN$results[,c("Sens","Spec")],1,sum)-1),"decay"]    
        models.options@ANN$maxit <- maxit.ANN
      }else{
        models.options@ANN$size <- tune.ANN$bestTune$size
        models.options@ANN$decay <- tune.ANN$bestTune$decay 
        models.options@ANN$maxit <- maxit.ANN
      }}else{ if('ANN' %in% models){cat("Tuning ANN failed!"); tune.ANN <- "FAILED"}}
    
  }
  
  if('GAM' %in% models){
    cat("Start tuning GAM\n")
    
    if(is.null(ctrl.GAM)){ctrl.GAM <- trControl}
    
    try(tune.GAM <-   caret::train(data@data.env.var, resp, 
                                   method = method.GAM,
                                   trControl = ctrl.GAM,
                                   weights = Yweights))
    cat(paste("Finished tuning GAM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.GAM)){
      if(metric == 'TSS'){
        models.options@GAM$select <- tune.GAM$results[which.max(apply(tune.GAM$results[,c("Sens","Spec")],1,sum)-1),"select"]    
        models.options@GAM$method <- as.character(tune.GAM$results[which.max(apply(tune.GAM$results[,c("Sens","Spec")],1,sum)-1),"method"])    
      }else {
        models.options@GAM$select <- tune.GAM$bestTune$select
        models.options@GAM$method <- as.character(tune.GAM$bestTune$method)
      }}else{ if('GAM' %in% models){cat("Tuning GAM failed!"); tune.GAM <- "FAILED"}}
  }
  
  if('MARS' %in% models){
    cat("Start tuning MARS\n")
    
    if(is.null(ctrl.MARS)){ctrl.MARS <- trControl}
    
    if(is.null(models.options@MARS$nk)){nprune <- 2:max(21, 2 * ncol(data@data.env.var) + 1)
    }else{
      nprune <- 2:min(models.options@MARS$nk,38)}
    tune.grid <- expand.grid(.degree = 1:2, .nprune = nprune)
    try(tune.MARS <-   caret::train(data@data.env.var, resp, 
                                    method = method.MARS,
                                    tuneGrid = tune.grid,
                                    trControl = ctrl.MARS,
                                    weights = Yweights))
    
    cat(paste("Finished tuning MARS\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.MARS)){
      if(metric == 'TSS'){
        if("degree" %in% names(models.options@MARS)){
          models.options@MARS$degree <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"degree"]    
        }else{
          models.options@MARS$interaction.level <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"degree"]-1    
        }
        models.options@MARS$nprune <- tune.MARS$results[which.max(apply(tune.MARS$results[,c("Sens","Spec")],1,sum)-1),"nprune"]    
      }else{
        if("degree" %in% names(models.options@MARS)){
          models.options@MARS$degree <- tune.MARS$bestTune$degree
        }else{
          models.options@MARS$interaction.level <- tune.MARS$bestTune$degree-1    
        }
        models.options@MARS$nprune <- tune.MARS$bestTune$nprune
      }}else{ if('MARS' %in% models){cat("Tuning MARS failed!"); tune.MARS <- "FAILED"}}
  }
  
  
  if('GLM' %in% models){
    cat("Start tuning GLM\n")
    
    if(is.null(ctrl.GLM)){ctrl.GLM <- trControl}
    if("s_smoother" %in% type.GLM){requireNamespace("gam", quietly = TRUE)}
    
    fm<-list()
    GLM.results<-NULL 
    i<-0
    for(type in type.GLM){
      for(IA in interaction.GLM){
        i<-i+1
        try(tune.GLM <-   caret::train( makeFormula("resp",data@data.env.var, type= type,interaction.level = IA),
                                        data=cbind(data@data.env.var,resp=resp),
                                        method = method.GLM,
                                        trControl = ctrl.GLM,
                                        weights = Yweights))  
        try(GLM.results <-  rbind(GLM.results,cbind(tune.GLM$results,il=IA,type=type)))
        try(fm[[i]] <- formula(tune.GLM$finalModel))
      }
    } 
    
    glm.best<-which.max(GLM.results$ROC)
    models.options@GLM$interaction.level <- GLM.results[glm.best,"il"]     
    models.options@GLM$type <- as.character(GLM.results[glm.best,"type"])
    models.options@GLM$myFormula <- formula(paste(data@sp.name,"~",gsub("`","",as.character(fm[[glm.best]])[3])))
    models.options@GLM$test <- "none" 
    
    cat(paste("Finished tuning GLM\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  }      
  
  
  if('FDA' %in% models){
 
    cat("Start tuning FDA\n")
    
    if(is.null(ctrl.FDA)){ctrl.FDA <- trControl}
    
    tune.grid <- expand.grid(.degree = 1:2, .nprune = 2:38)
    try(tune.FDA <- caret::train(data@data.env.var, factor(resp), 
                                 method = "fda",
                                 tuneGrid = tune.grid,                  
                                 trControl = ctrl.FDA,
                                 weights = Yweights))
    cat(paste("Finished tuning FDA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.FDA)){
      #models.options@FDA$method <- "earth"   
      if(metric == 'TSS'){
        models.options@FDA$add_args <- list(degree=tune.FDA$results[which.max(apply(tune.FDA$results[,c("Sens","Spec")],1,sum)-1),"degree"],   
                                            nprune=tune.FDA$results[which.max(apply(tune.FDA$results[,c("Sens","Spec")],1,sum)-1),"nprune"])    
      }else{
        models.options@FDA$add_args <- list(degree=tune.FDA$bestTune$degree,nprune=tune.FDA$bestTune$nprune)
      }}else{ if('FDA' %in% models){cat("Tuning FDA failed!"); tune.FDA <- "FAILED"}}
  }
  
  if('CTA' %in% models){
    cat("Start tuning CTA\n")
    
    if(is.null(ctrl.CTA)){ctrl.CTA <- trControl}    
    
    cat("Tuning Complexity Parameter")    
    try(tune.CTA.rpart <- caret::train(data@data.env.var, resp, 
                                       method = "rpart",
                                       tuneLength = tuneLength,
                                       trControl = ctrl.CTA,
                                       metric=metric,
                                       weights = Yweights))
    
    cat("Tuning Max Tree Depth")
    try(tune.CTA.rpart2 <-  caret::train(data@data.env.var, resp,
                                         method = "rpart2",
                                         tuneLength = tuneLength,
                                         trControl = ctrl.CTA,
                                         metric=metric,
                                         weights = Yweights))
    cat(paste("Finished tuning CTA\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.CTA.rpart)){
      if(metric == 'TSS'){
        models.options@CTA$control$cp <- tune.CTA.rpart$results[which.max(apply(tune.CTA.rpart$results[,c("Sens","Spec")],1,sum)-1),"cp"]    
      }else{
        models.options@CTA$control$cp <- tune.CTA.rpart$bestTune      
      }}else{ if('CTA' %in% models){cat("Tuning CTA cp failed!"); tune.CTA.rpart <- "FAILED"}}
    
    if(!is.null(tune.CTA.rpart2)){
      if(metric == 'TSS'){
        models.options@CTA$control$maxdepth <- tune.CTA.rpart2$results[which.max(apply(tune.CTA.rpart2$results[,c("Sens","Spec")],1,sum)-1),"maxdepth"]    
      }else{
        models.options@CTA$control$maxdepth <- tune.CTA.rpart2$bestTune      
      }}else{ if('CTA' %in% models){cat("Tuning CTA maxdepth failed!"); tune.CTA.rpart2 <- "FAILED"}}
  }
  
  if('MAXENT.Phillips' %in% models){
    cat("Start tuning MAXENT.Phillips\n")
    if(cvmethod.ME != 'randomkfold'){kfolds.ME <- NA}
    try(tune.MAXENT.Phillips <- tuning.maxent(pres=data@data.env.var[data@data.species==1 & !is.na(data@data.species),],
                                              bg= data@data.env.var[data@data.species==0 | is.na(data@data.species),],
                                              method=cvmethod.ME, kfolds = kfolds.ME,#env.ME,
                                              bin.output=TRUE, clamp=clamp.ME, parallel = parallel.ME, numCores = numCores.ME,
                                              categoricals=NULL))
    cat(paste("Finished tuning MAXENT.Phillips\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
    
    if(!is.null(tune.MAXENT.Phillips)){
      if(metric.ME=="ROC"){metric.ME <- "Mean.AUC"}
      if(!metric.ME %in% c("Mean.AUC", "Mean.AUC.DIFF", "Mean.ORmin", "Mean.OR10", "AICc")){metric.ME <- "Mean.AUC"; cat("Invalid metric.ME argument! metric.ME was set to Mean.AUC")}
      if(metric.ME == 'Mean.AUC'){
        models.options@MAXENT.Phillips$linear <- grepl("L",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$quadratic <- grepl("Q",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$hinge <- grepl("H",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$product <- grepl("P",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$threshold <- grepl("T",tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$betamultiplier <- tune.MAXENT.Phillips@results[which.max(tune.MAXENT.Phillips@results[,metric.ME]),"rm"]  
      }else {       
        models.options@MAXENT.Phillips$linear <- grepl("L",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$quadratic <- grepl("Q",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$hinge <- grepl("H",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$product <- grepl("P",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$threshold <- grepl("T",tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"features"]) 
        models.options@MAXENT.Phillips$betamultiplier <- tune.MAXENT.Phillips@results[which.min(tune.MAXENT.Phillips@results[,metric.ME]),"rm"]   
      }}else{ if('MAXENT.Phillips' %in% models){cat("Tuning MAXENT.Phillips failed!"); tune.MAXENT.Phillips <- "FAILED"}}
  }
  
  
  # if('MAXENT.Tsuruoka' %in% models){
  #   cat("Start tuning MAXENT.Tsuruoka\n")
  #   try(tune.MAXENT.Tsuruoka <- as.data.frame(tune.maxent(data@data.env.var,data@data.species,nfold=kfolds.ME,showall=T)))
  #   cat(paste("Finished tuning MAXENT.Tsuruoka\n","\n-=-=-=-=-=-=-=-=-=-=\n"))
  #   
  #   if(!is.null(tune.MAXENT.Tsuruoka)){
  #     models.options@MAXENT.Tsuruoka$l1_regularizer <- tune.MAXENT.Tsuruoka$l1_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     models.options@MAXENT.Tsuruoka$l2_regularizer <- tune.MAXENT.Tsuruoka$l2_regularizer[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #     models.options@MAXENT.Tsuruoka$use_sgd <- ifelse(tune.MAXENT.Tsuruoka[which.max(tune.MAXENT.Tsuruoka$accuracy),]$use_sgd==0,F,T)
  #     models.options@MAXENT.Tsuruoka$set_heldout <- tune.MAXENT.Tsuruoka$set_heldout[which.max(tune.MAXENT.Tsuruoka$accuracy)]
  #   } else { if('MAXENT.Tsuruoka' %in% models){cat("Tuning MAXENT.Tsuruoka failed!"); tune.MAXENT.Tsuruoka <- "FAILED"}}
  # }  
  
  
  return(list(models.options=models.options, tune.SRE =tune.SRE,  tune.CTA.rpart = tune.CTA.rpart, tune.CTA.rpart2 = tune.CTA.rpart2,
              tune.RF = tune.RF, tune.ANN = tune.ANN,  tune.MARS = tune.MARS, tune.FDA = tune.FDA, tune.GBM=tune.GBM,
              tune.GAM = tune.GAM, tune.MAXENT.Phillips = tune.MAXENT.Phillips, 
              # tune.MAXENT.Tsuruoka = tune.MAXENT.Tsuruoka, 
              tune.GLM=tune.GLM))
}


##
#### Modified tuning function from the ENMeval package to tune MAXENT.Phillips (internal function for BIOMOD_tuning)

tuning.maxent <-
  function (occ, env=NULL,pres=NULL, bg=NULL, bg.coords=NULL, occ.grp=NULL, bg.grp=NULL, method=NULL, maxent.args, 
            args.lab, categoricals=NULL, aggregation.factor=c(2,2), kfolds=NA, bin.output=FALSE, 
            clamp, rasterPreds=FALSE, parallel=FALSE, numCores=NULL,RMvalues=seq(0.5, 4, 0.5), fc=c("L", "LQ", "H", "LQH", "LQHP", "LQHPT")) 
  {
    ###NEW
    requireNamespace("ENMeval", quietly = TRUE)
    maxent.args <- ENMeval::make.args(RMvalues, fc)
    args.lab <- ENMeval::make.args(RMvalues, fc, labels = TRUE)
    if(!is.null(pres)){occ<-pres}
    if(!is.null(bg)){bg.coords<-bg}
    #####
    noccs <- nrow(occ)
    if (method == "checkerboard1") 
      group.data <- ENMeval::get.checkerboard1(occ, env, bg.coords, 
                                               aggregation.factor)
    if (method == "checkerboard2") 
      group.data <- ENMeval::get.checkerboard2(occ, env, bg.coords, 
                                               aggregation.factor)
    if (method == "block") 
      group.data <- ENMeval::get.block(occ, bg.coords)
    if (method == "jackknife") 
      group.data <- ENMeval::get.jackknife(occ, bg.coords)
    if (method == "randomkfold") 
      group.data <- ENMeval::get.randomkfold(occ, bg.coords, kfolds)
    if (method == "user") 
      group.data <- ENMeval::get.user(occ.grp, bg.grp)
    nk <- length(unique(group.data$occ.grp))
    ###NEW
    if(is.null(pres)){
      pres <- as.data.frame(extract(env, occ))}
    if(is.null(bg)){
      bg <- as.data.frame(extract(env, bg.coords))}
    #####
    if (any(is.na(colSums(pres)))) {
      message("Warning: some predictors variables are NA at some occurrence points")
    }
    if (any(is.na(colSums(bg)))) {
      message("Warning: some predictors variables are NA at some background points")
    }
    if (!is.null(categoricals)) {
      for (i in 1:length(categoricals)) {
        pres[, categoricals[i]] <- as.factor(pres[, categoricals[i]])
        bg[, categoricals[i]] <- as.factor(bg[, categoricals[i]])
      }
    }
    tune <- function() {
      if (length(maxent.args) > 1 & !parallel) {
        setTxtProgressBar(pb, i)
      }
      x <- rbind(pres, bg)
      p <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
      tmpfolder <- tempfile()
      full.mod <- dismo::maxent(x, p, args = maxent.args[[i]], factors = categoricals, 
                                path = tmpfolder)
      pred.args <- c("outputformat=raw", ifelse(clamp == TRUE, 
                                                "doclamp=true", "doclamp=false"))
      if (rasterPreds == TRUE) {
        predictive.map <- predict(full.mod, env, args = pred.args)
      }
      else {
        predictive.map <- stack()
      }
      AUC.TEST <- double()
      AUC.DIFF <- double()
      OR10 <- double()
      ORmin <- double()
      for (k in 1:nk) {
        train.val <- pres[group.data$occ.grp != k, ]
        test.val <- pres[group.data$occ.grp == k, ]
        bg.val <- bg[group.data$bg.grp != k, ]
        x <- rbind(train.val, bg.val)
        p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))
        mod <- dismo::maxent(x, p, args = maxent.args[[i]], factors = categoricals, 
                             path = tmpfolder)
        ### Specify dismo!!! Problems with biomod!
        AUC.TEST[k] <- dismo::evaluate(test.val, bg, mod)@auc
        AUC.DIFF[k] <- max(0, dismo::evaluate(train.val, bg, mod)@auc - 
                             AUC.TEST[k])
        ########
        p.train <- predict(mod, train.val, args = pred.args)
        p.test <- predict(mod, test.val, args = pred.args)
        if (nrow(train.val) < 10) {
          n90 <- floor(nrow(train.val) * 0.9)
        }
        else {
          n90 <- ceiling(nrow(train.val) * 0.9)
        }
        train.thr.10 <- rev(sort(p.train))[n90]
        OR10[k] <- mean(p.test < train.thr.10)
        train.thr.min <- min(p.train)
        ORmin[k] <- mean(p.test < train.thr.min)
      }
      unlink(tmpfolder, recursive = TRUE)
      stats <- c(AUC.DIFF, AUC.TEST, OR10, ORmin)
      return(list(full.mod, stats, predictive.map))
    }
    if (parallel == TRUE) {
      requireNamespace("foreach")
      allCores <- detectCores()
      if (is.null(numCores)) {
        numCores <- allCores
      }
      c1 <- makeCluster(numCores)
      doParallel::registerDoParallel(c1)
      numCoresUsed <- foreach::getDoParWorkers()
      message(paste("Of", allCores, "total cores using", numCoresUsed))
      message("Running in parallel...")
      out <- foreach::foreach(i = seq_len(length(maxent.args)), .packages = c("dismo", 
                                                                              "raster", "ENMeval")) %dopar% {
                                                                                tune()
                                                                              }
      stopCluster(c1)
    }
    else {
      pb <- txtProgressBar(0, length(maxent.args), style = 3)
      out <- list()
      for (i in 1:length(maxent.args)) {
        out[[i]] <- tune()
      }
      close(pb)
    }
    full.mods <- sapply(out, function(x) x[[1]])
    statsTbl <- as.data.frame(t(sapply(out, function(x) x[[2]])))
    if (rasterPreds) {
      predictive.maps <- stack(sapply(out, function(x) x[[3]]))
    }
    else {
      predictive.maps <- stack()
    }
    AUC.DIFF <- statsTbl[, 1:nk]
    AUC.TEST <- statsTbl[, (nk + 1):(2 * nk)]
    OR10 <- statsTbl[, ((2 * nk) + 1):(3 * nk)]
    ORmin <- statsTbl[, ((3 * nk) + 1):(4 * nk)]
    names(AUC.DIFF) <- paste("AUC.DIFF_bin", 1:nk, sep = ".")
    Mean.AUC.DIFF <- rowMeans(AUC.DIFF)
    Var.AUC.DIFF <- ENMeval::corrected.var(AUC.DIFF, noccs)
    names(AUC.TEST) <- paste("AUC_bin", 1:nk, sep = ".")
    Mean.AUC <- rowMeans(AUC.TEST)
    Var.AUC <- ENMeval::corrected.var(AUC.TEST, noccs)
    names(OR10) <- paste("OR10_bin", 1:nk, sep = ".")
    Mean.OR10 <- rowMeans(OR10)
    Var.OR10 <- apply(OR10, 1, var)
    names(ORmin) <- paste("ORmin_bin", 1:nk, sep = ".")
    Mean.ORmin <- rowMeans(ORmin)
    Var.ORmin <- apply(ORmin, 1, var)
    full.AUC <- double()
    for (i in 1:length(full.mods)) full.AUC[i] <- full.mods[[i]]@results[5]
    nparm <- numeric()
    for (i in 1:length(full.mods)) nparm[i] <- ENMeval::get.params(full.mods[[i]])
    if (rasterPreds == TRUE) {
      aicc <- ENMeval::calc.aicc(nparm, occ, predictive.maps)
    }
    else {
      aicc <- rep(NaN, length(full.AUC))
    }
    features <- args.lab[[1]]
    rm <- args.lab[[2]]
    settings <- paste(args.lab[[1]], args.lab[[2]], sep = "_")
    res <- data.frame(settings, features, rm, full.AUC, Mean.AUC, 
                      Var.AUC, Mean.AUC.DIFF, Var.AUC.DIFF, Mean.OR10, Var.OR10, 
                      Mean.ORmin, Var.ORmin, nparm, aicc)
    if (bin.output == TRUE) {
      res <- as.data.frame(cbind(res, AUC.TEST, AUC.DIFF, OR10, 
                                 ORmin))
    }
    if (rasterPreds == TRUE) {
      names(predictive.maps) <- settings
    }
    results <- ENMeval::ENMevaluation(results = res, predictions = predictive.maps, 
                                      models = full.mods, partition.method = method, occ.pts = occ, 
                                      occ.grp = group.data[[1]], bg.pts = bg.coords, bg.grp = group.data[[2]])
    return(results)
  }




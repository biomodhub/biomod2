# BIOMOD_RangeSize documentation ----------------------------------------------
##' @name bm_ModelAnalysis
##' @author Hélène Blancheteau
##' 
##' @title Analyze the residuals of the single models 
##' 
##' @description This function return differents graphs to help analyse the single models.
##' 
##' 
##' @param bm.mod a \code{BIOMOD.models.out}
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' @param do.plot a \code{logical}, print the plot or not 
##' @param color.by a \code{character} between "full.name", "species", "PA", "RUN" or "algo" 
##' @param do.Rpredicted a \code{logical} : do Rpredicted have to be computed or no ?  
##' 
##' @return
##' \enumerate{
##'   \item a plot residuals ~ observations to detect outliers
##'   \item a histogramm of the residuals 
##'   \item a qqplot of the residuals
##'   \item plot residuals ~ fitted
##'   \item a plot of the different Rsquared available (Rsquared, Rsquared_aj, Psquared)
##' }
##' 
##' 
##' @details 
##' Only the models coherent with a residuals analysis will be selected for some plots.
##' 
##' 
##' @keywords analyze models residuals
##' 
##' 
##' @importFrom foreach foreach %do%
##' @importFrom ggplot2 ggplot geom_hline geom_point stat_qq stat_qq_line vars geom_histogram
##' @importFrom tidyr separate
##' 
##' @export
##' 
##' 

bm_ModelAnalysis <- function(bm.mod, 
                             models.chosen = "all",
                             do.plot = TRUE, 
                             color.by = "full.name",
                             do.Rpredicted = FALSE){
  
  .bm_cat("Model analysis")
  args <- .bm_ModelAnalysis.check.args(bm.mod = bm.mod, models.chosen = models.chosen)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  models.coherent <- grep("GAM|MARS|GLM|CTA", models.chosen, value = T)

  res <- foreach(mod.name = models.coherent, .combine = rbind) %dopar% {
    cat("\n\t> Analysis", mod.name, "...")
    mod <- get(BIOMOD_LoadModels(bm.out = bm.mod, full.name = mod.name))
    infos <- .bm_return_info_analysys(mod)
    return(data.frame("full.name" = mod.name, infos))
  }
  res <- tidyr::separate(res, col = "full.name", into = c("species", "PA", "RUN", "algo"), remove = FALSE)
  cat("\n")
  
  if (do.Rpredicted){
    Psquared <- foreach(mod.name = models.chosen, .combine = c) %dopar% {
      cat("\n\t> Psquared", mod.name, "...")
      PRESS <- .bm_PRESS(bm.mod, mod.name, perc = 1) ## attention perc ? 
      TSS <- .bm_TotalSumSquares(bm.mod, mod.name)
      pred_Rsquared <- 1 - PRESS/TSS
      return(pred_Rsquared)
    }
  } else {
    Psquared = NA
  }

  ## Outliers
  plot_outliers <- ggplot(res, aes(y = residuals, x = obs, color = .data[[color.by]]))+
    geom_point() +
    geom_hline(yintercept = 0, color = "darkred") +
    xlab("Observations") +
    ylab("Residuals") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))

  ## QQ plot of residuals

  plot_qqplot <- ggplot(res, aes(sample = residuals, color = .data[[color.by]])) +
    stat_qq()+
    stat_qq_line(color = "darkred") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))

  ## Histogramme of residuals

  plot_hist <- ggplot(res, aes(x = residuals, color = .data[[color.by]])) +
    facet_wrap(vars(full.name))+
    geom_histogram(fill = NA, position = "dodge", bins = 30) + #ça va pas
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))

  ## Fitted vs resid

  plot_fitted <- ggplot(res, aes(y = residuals, x = fitted, color = .data[[color.by]])) +
    geom_point()+
    geom_hline(yintercept = 0, color = "darkred") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## Psquared
  Psquared <- data.frame(full.name = models.chosen, metric = "Psquared", value = Psquared)
  eval <- get_evaluations(bm.mod)
  if ("Rsquared" %in% unique(eval$metric.eval)){
    eval_Rsq <- eval[eval$full.name %in% models.chosen & eval$metric.eval == "Rsquared", c("full.name", "calibration", "validation")]
    Psquared <- rbind(Psquared, 
                      data.frame(full.name = eval_Rsq$full.name, metric = "Rsquared_calibration", value = eval_Rsq$calibration),
                      data.frame(full.name = eval_Rsq$full.name, metric = "Rsquared_validation", value = eval_Rsq$validation)) ##C'est tordu mais c'est pour pas appeler melt
  }
  if ("Rsquared_aj" %in% unique(eval$metric.eval)){
    eval_Rsq_aj <- eval[eval$full.name %in% models.chosen & eval$metric.eval == "Rsquared_aj", c("full.name", "calibration", "validation")]
    Psquared <- rbind(Psquared, 
                      data.frame(full.name = eval_Rsq_aj$full.name, metric = "Rsquared_aj_calibration", value = eval_Rsq_aj$calibration),
                      data.frame(full.name = eval_Rsq_aj$full.name, metric = "Rsquared_aj_validation", value = eval_Rsq_aj$validation)) 
  }
  Psquared <- tidyr::separate(Psquared, col = "full.name", into = c("species", "PA", "RUN", "algo"), remove = FALSE)
  
  plot_Psquared <- ggplot(Psquared, aes(y= value, x = full.name, color = .data[[color.by]], shape = metric))+
    geom_point()+
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))

  if (do.plot){
    plot(plot_outliers)
    plot(plot_qqplot)
    plot(plot_hist)
    plot(plot_fitted)
    plot(plot_Psquared)
  }
  
  .bm_cat("Done")
  return(res)
}



# .bm_return_info_analysys
# Functions to get residuals and fitted values of a model
# 
# param mod a \code{biomod2_model}
# importFrom stats fitted residuals
# 
# export

setGeneric(".bm_return_info_analysys", function(mod) { standardGeneric(".bm_return_info_analysys")})

setMethod(".bm_return_info_analysys", signature("CTA_biomod2_model"), function(mod){
  resids <- as.vector(residuals(mod@model)) 
  fit <- as.vector(mod@model$y) #[,"y"]
  return(data.frame("obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
})

setMethod(".bm_return_info_analysys", signature("GAM_biomod2_model"), function(mod){
  resids <- as.vector(residuals(mod@model)) 
  fit <- as.vector(fitted(mod@model))
  return(data.frame("obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
})

setMethod(".bm_return_info_analysys", signature("GLM_biomod2_model"), function(mod){
  resids <- as.vector(residuals(mod@model)) 
  fit <- as.vector(fitted(mod@model))
  return(data.frame("obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
})

setMethod(".bm_return_info_analysys", signature("MARS_biomod2_model"), function(mod){
  resids <- as.vector(residuals(mod@model)) 
  fit <- as.vector(fitted(mod@model))
  return(data.frame( "obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
})



# Argument Check ---------------------------------------------------------------

.bm_ModelAnalysis.check.args <- function(bm.mod, models.chosen) {
  
  .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
  
  ## Stop if ordinal -------------------------------------------------------
  if (bm.mod@data.type %in% c("ordinal", "binary")){
    stop("Model analysis is not ready for binary or ordinal data...")
  }
  
  ## Check models.chosen ---------------------------------------------------
  if (models.chosen[1] == 'all') {
    models.chosen <- bm.mod@models.computed
  } else {
    models.chosen <- intersect(models.chosen, bm.mod@models.computed)
  }
  if (length(models.chosen) < 1) {
    stop('No models selected')
  }
  
  ## check that given models exist
  files.check <- paste0(bm.mod@dir.name, "/", bm.mod@sp.name, "/models/",
                        bm.mod@modeling.id, "/", models.chosen)

  not.checked.files <- grep('MAXENT|SRE', files.check)
  if (length(not.checked.files) > 0) {
    files.check <- files.check[-not.checked.files]
  }
  missing.files <- files.check[!file.exists(files.check)]
  if (length(missing.files) > 0) {
    stop(paste0("Models files missing : ", toString(missing.files)))
    if (length(missing.files) == length(files.check)) {
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }
  return(list(models.chosen = models.chosen))
}



.bm_PRESS <- function(bm.mod, mod.name, perc = 0.5) { #perc --> percentage of points evaluate for overfitting 
  
  ## Data
  dataset <- paste0("_", unlist(strsplit(mod.name, "_"))[2:3], collapse = "")
  algo <- unlist(strsplit(mod.name, "_"))[4]
  
  bm.format <- get_formal_data(bm.mod)
  mySpExpl <- get_species_data(bm.format)
  calib.lines <- as.data.frame(get_calib_lines(bm.mod))
  mySpExpl <- mySpExpl[which(calib.lines[,dataset] == TRUE), ]
  
  bm.options <- get_options(bm.mod)
  bm.options <- bm.options@options[[grep(algo, names(bm.options@options))]] ### Warning grep !! 
  
  ## Boucle 
  sample_points <- sample(1:nrow(mySpExpl), nrow(mySpExpl)*perc)
  model.call <- paste0(bm.options@package, "::", bm.options@func)
  
  residuals_predicted <- foreach(point = sample_points, .combine = c) %do% {
    argstmp <- bm.options@args.values[[dataset]]
    argstmp$data <- mySpExpl[-point, ]
    if (algo %in% c("MARS", "RF")){
      formu <- bm_MakeFormula(resp.name = names(mySpExpl)[1], expl.var = mySpExpl[1, 4:ncol(mySpExpl)], type = "simple", interaction.level = 0)
      argstmp$formula <- formu
    }
    # if(algo == "XGBOOST"){
    #   argstmp$label <- argstmp$data[,1]
    #   argstmp$data <- as.matrix(argstmp$data[, 4:ncol(mySpExpl)])
    #   argstmp <- argstmp[c("data", names(argstmp)[which(!(names(argstmp) %in% c("formula", "data")))])]
    # } 
    argstmp <- argstmp[c("formula","data", names(argstmp)[which(!(names(argstmp) %in% c("formula", "data")))])]
    
    
    new.model <- try(do.call(eval(parse(text = model.call)), argstmp), silent = TRUE)
    new.res <- NA
    if (inherits(new.model, 'try-error') != T){
      new.pred <- predict(new.model, mySpExpl[point, ])
      new.res <- mySpExpl[point, 1] - new.pred
    } 
    
    return(new.res)
  }## Return residuals x(-i)
  
  ## Calcul PRESS + R predicted
  PRESS <- sum(residuals_predicted ^ 2)
  return(PRESS)
}
  
  
## TSS with bm_FindOptimStat MAIS on a déjà TSS dans bm_FOS !!
## Pourquoi tout ces acronymes sont les mêmes ?!?

.bm_TotalSumSquares <- function(bm.mod, mod.name){
  dataset <- paste0("_", unlist(strsplit(mod.name, "_"))[2:3], collapse = "")
  bm.format <- get_formal_data(bm.mod)
  mySpExpl <- get_species_data(bm.format)
  calib.lines <- as.data.frame(get_calib_lines(bm.mod))
  mySpExpl <- mySpExpl[which(calib.lines[,dataset] == TRUE), ]
  
  #fit <- tab_res[tab_res$full.name == mod.name, "fitted"]
  y <- mySpExpl[, 1]
  y_mean <- mean(y, na.rm = T)
  TSS <- sum((y - y_mean) ^ 2)
  return(TSS)
}



## Il y a un problème 
## Peut être pas

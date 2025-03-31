###################################################################################################
##' @name bm_ModelAnalysis
##' @author Hélène Blancheteau
##' 
##' @title Analyze the residuals of the single models 
##' 
##' @description This function return several graphs to help analyse the single models.
##' 
##' 
##' @param bm.mod a \code{BIOMOD.models.out}
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function
##' @param do.plot a \code{logical}, print the plots or not 
##' @param color.by a \code{character} between "full.name", "species", "PA", "RUN" or "algo" 
##' to select the color parameter of the plots.
##' 
##' @return
##' A list with: 
##' \enumerate{
##'   \item a dataframe with the observation data, the fitted values and the residuals for all models
##'   \item a dataframe with the compilation of the R scores
##'   \item A to E plots (see details)
##' }
##' 
##' 
##' @details 
##' All the plots will be made for all the models, independently to the different models assumptions. 
##' It is up to the user to interpret the graphs in the light of the model assumptions.
##' 
##' \enumerate{
##'   \item Plot A: residuals ~ observations number. This plot helps to detect one or several outliers.
##'   The x-axis only helps to find the outlier number. 
##'   \item Plot B and C: These are two representation of the distribution of the residuals. If your residuals must follow a normal distribution, 
##'   the points should follow the black line in the Q-Q plot and present a gaussian distribution on the histogram. 
##'   \item Plot D residuals ~ fitted values. This plot helps to detect an heteroscedasticity of the residuals. 
##'   \item Plot E : a plot of the different Rsquared values available (Rsquared, Rsquared_aj).
##'    An big gap between the calibration and validation values can be the sign of an overfitting.
##' }
##' 
##' 
##' 
##' @keywords analyze models residuals
##' 
##' @examples
##' library(terra)
##' 
##' # Load species occurrences (6 species available)
##' data(DataSpecies)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- 'VulpesVulpes'
##' 
##' # Get corresponding presence/absence data
##' myResp <- as.numeric(DataSpecies[, myRespName])
##' 
##' # Get corresponding XY coordinates
##' myRespXY <- DataSpecies[, c('X_WGS84', 'Y_WGS84')]
##' 
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' data(bioclim_current)
##' myExpl <- terra::rast(bioclim_current)
##' 
##' # Transform binary data as count data
##' poissonDistri <- rpois(sum(myResp), 10)
##' myResp[myResp == 1] <- poissonDistri
##' 
##' # ---------------------------------------------------------------
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
##'                                        resp.var = myResp,
##'                                        resp.xy = myRespXY,
##'                                        expl.var = myExpl)
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('GAM', 'GLM', 'MARS'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('Rsquared', 'Rsquared_aj'),
##'                                       seed.val = 42)
##' }
##' 
##' 
##' # ---------------------------------------------------------------
##' # bm_ModelAnalysis
##' analysis <- bm_ModelAnalysis(myBiomodModelOut, color.by = "RUN")
##' 
##' plot(analysis$plotA)
##' 
##' unlink(myRespName, recursive = TRUE)
##' 
##' 
##' @importFrom foreach foreach %do%
##' @importFrom ggplot2 ggplot geom_hline geom_point stat_qq stat_qq_line vars geom_histogram
##' 
##' @export
##' 
##' 
###################################################################################################


bm_ModelAnalysis <- function(bm.mod, 
                             models.chosen = 'all',
                             do.plot = TRUE, 
                             color.by = 'full.name')
{
  .bm_cat("Model analysis")
  cat("> Warning : All the plots will be made for all the models, independently to the different models assumptions. 
      It is up to you to interpret the graphs in the light of the model assumptions.\n")
  
  args <- .bm_ModelAnalysis.check.args(bm.mod = bm.mod, models.chosen = models.chosen,
                                       color.by = color.by, do.plot = do.plot)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  models.coherent <- grep("GAM|MARS|GLM|CTA", models.chosen, value = TRUE)
  
  cat("\n> Residuals computation")
  predictions <- get_predictions(bm.mod, full.name = models.chosen)
  observations <- data.frame(points = 1:max(predictions$points),
                             obs = get_formal_data(bm.mod, subinfo = "resp.var"))
  
  if(bm.mod@data.type == "binary"){
    predictions$pred <- predictions$pred/1000
  }
  
  if(bm.mod@data.type == "ordinal"){
    names_levels <- levels(observations$obs)
    
    predictions$pred <- as.character(predictions$pred)
    to_change <- predictions$pred[predictions$pred %in% names_levels] 
    to_change <- factor(to_change, levels = levels(observations$obs), ordered = T)
    predictions$pred[predictions$pred %in% names_levels] <- as.character(as.numeric(to_change))
    
    predictions$pred <- as.numeric(predictions$pred)
    observations$obs <- as.numeric(observations$obs)
  }
  
  res <- merge(predictions, observations)
  res$residuals <- res$obs - res$pred
  res <- tidyr::separate(res, col = "full.name", into = c("species", "PA", "RUN", "algo"), remove = FALSE)
  
  

  ## Outliers
  plot_outliers <- ggplot(res, aes(y = residuals, x = points, color = .data[[color.by]]))+
    geom_point(size = 1) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    xlab("Observations number") +
    ylab("Residuals") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(colour = guide_legend(override.aes = list(size = 2)))+
    paletteer::scale_color_paletteer_d(palette)+
    ggtitle("Visualising residuals: are there any outliers?")

  ## QQ plot of residuals

  plot_qqplot <- ggplot(res, aes(sample = residuals, color = .data[[color.by]])) +
    stat_qq(size = 1)+
    stat_qq_line(color = "black", size = 1) +
    xlab("Theoretical quantiles") +
    ylab("Standardized residuals") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(colour = guide_legend(override.aes = list(size = 2)))+
    paletteer::scale_color_paletteer_d(palette)+
    ggtitle("Q-Q plot of residuals")

  ## Histogramme of residuals

  plot_hist <- ggplot(res, aes(x = residuals, color = .data[[color.by]])) +
    facet_wrap(vars(full.name))+
    geom_histogram(fill = NA, position = "dodge", bins = 30) + 
    xlab("Residuals") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1),
          , axis.title.y = element_blank(),
          , axis.text.y = element_blank(),
          , axis.ticks.y = element_blank())+
    paletteer::scale_color_paletteer_d(palette)+
    ggtitle("Distribution of residuals")

  ## Fitted vs resid

  plot_fitted <- ggplot(res, aes(y = residuals, x = pred, color = .data[[color.by]])) +
    geom_point(size = 1)+
    geom_hline(yintercept = 0, color = "black", size = 1) +
    xlab("Fitted") +
    ylab("Residuals") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))+
    guides(colour = guide_legend(override.aes = list(size = 2)))+
    paletteer::scale_color_paletteer_d(palette)+
    ggtitle("Residuals ~ Fitted values (Tukey-Anscombe plot)")
  
  if(!bm.mod@data.type %in% c("binary", "ordinal")){
    ## Rscore
    cat("\n> Rscore compilation")
    Rscore <- data.frame(full.name = models.chosen, metric = NA, value = NA)
    eval <- get_evaluations(bm.mod)
    if ("Rsquared" %in% unique(eval$metric.eval)){
      eval_Rsq <- eval[eval$full.name %in% models.chosen & eval$metric.eval == "Rsquared", c("full.name", "calibration", "validation")]
      Rscore <- rbind(Rscore, 
                        data.frame(full.name = eval_Rsq$full.name, metric = "Rsquared_calibration", value = eval_Rsq$calibration),
                        data.frame(full.name = eval_Rsq$full.name, metric = "Rsquared_validation", value = eval_Rsq$validation)) ##C'est tordu mais c'est pour pas appeler melt
    }
    if ("Rsquared_aj" %in% unique(eval$metric.eval)){
      eval_Rsq_aj <- eval[eval$full.name %in% models.chosen & eval$metric.eval == "Rsquared_aj", c("full.name", "calibration", "validation")]
      Rscore <- rbind(Rscore, 
                        data.frame(full.name = eval_Rsq_aj$full.name, metric = "Rsquared_aj_calibration", value = eval_Rsq_aj$calibration),
                        data.frame(full.name = eval_Rsq_aj$full.name, metric = "Rsquared_aj_validation", value = eval_Rsq_aj$validation)) 
    }
    Rscore <- tidyr::separate(Rscore, col = "full.name", into = c("species", "PA", "RUN", "algo"), remove = FALSE)
    Rscore <- Rscore[!is.na(Rscore$metric),]
    
    plot_Rscore <- ggplot(Rscore, aes(y= value, x = full.name, color = .data[[color.by]], shape = metric))+
      geom_point(size = 2)+
      ylim(0,1)+
      xlab("Models")+
      theme(legend.title = element_blank()
            , legend.key = element_rect(fill = "white")
            , axis.text.x = element_text(angle = 45, hjust = 1),
            , axis.title.y = element_blank()) +
      paletteer::scale_color_paletteer_d(palette)+
      ggtitle("R scores")

  } else {
    Rscore <- NULL
  }
  

  if (do.plot){
    if(!bm.mod@data.type %in% c("binary", "ordinal")){print(plot_Rscore)}
    print(plot_outliers + plot_qqplot + plot_hist + plot_fitted + plot_annotation(tag_levels = 'A',
                                                                                  tag_prefix = "(",
                                                                                  tag_suffix = ")") )
  }
  
  .bm_cat("Done")
  return(list(residuals = res,
              R_scores = Rscore,
              plotA = plot_outliers,
              plotB = plot_qqplot,
              plotC = plot_hist,
              plotD = plot_fitted,
              plotE = plot_Rscore))
}





###################################################################################################

.bm_ModelAnalysis.check.args <- function(bm.mod, models.chosen, color.by, do.plot)
{
  .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
  
  ## Stop if ordinal -------------------------------------------------------
  # if (bm.mod@data.type %in% c("ordinal", "binary")){
  #   stop("Model analysis is not ready for binary or ordinal data...")
  # }
  
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
  
  ##check_color_by
  .fun_testIfIn(TRUE, "color.by", color.by, c("full.name", "species", "PA", "RUN", "algo" ))
  
  if (!isNamespaceLoaded("paletteer")) {
    if(!requireNamespace('paletteer', quietly = TRUE)) stop("Package 'paletteer' not found")
  }
  
  if(color.by == "species"){
    palette <- "ggthemes::Green_Orange_Teal"
  } else if (color.by == "algo"){
    palette <- "tvthemes::parksAndRec"
  } else if (color.by == "RUN"){
    palette <- "ggthemes::Classic_20"
  } else if (color.by == "PA"){
    palette <- "ggthemes::Tableau_20"
  } else {
    palette <- "ggsci::default_igv"
  }
  
  
  if(do.plot){
    if (!isNamespaceLoaded("patchwork")) {
      if(!requireNamespace('patchwork', quietly = TRUE)) stop("Package 'patchwork' not found")
    }
    library(patchwork)
  }
  
  return(list(models.chosen = models.chosen,
              palette = palette))
}


###################################################################################################
## Old codes and associated functions

# res <- foreach(mod.name = models.coherent, .combine = rbind) %dopar% {
#   cat("\n\t> Analysis", mod.name, "...")
#   mod <- get(BIOMOD_LoadModels(bm.out = bm.mod, full.name = mod.name))
#   infos <- .bm_return_info_analysys(mod)
#   return(data.frame("full.name" = mod.name, infos))
# }

# if (do.Rpredicted){
#   Psquared <- foreach(mod.name = models.chosen, .combine = c) %dopar% {
#     cat("\n\t> Psquared", mod.name, "...")
#     PRESS <- .bm_PRESS(bm.mod, mod.name, perc = 1) ## attention perc ? 
#     TSS <- .bm_TotalSumSquares(bm.mod, mod.name)
#     pred_Rsquared <- 1 - PRESS/TSS
#     return(pred_Rsquared)
#   }
# } else {
#   Psquared = NA
# }

# .bm_PRESS <- function(bm.mod, mod.name, perc = 0.5) { #perc --> percentage of points evaluate for overfitting 
#   
#   ## Data
#   dataset <- paste0("_", unlist(strsplit(mod.name, "_"))[2:3], collapse = "")
#   algo <- unlist(strsplit(mod.name, "_"))[4]
#   
#   bm.format <- get_formal_data(bm.mod)
#   mySpExpl <- get_species_data(bm.format)
#   calib.lines <- as.data.frame(get_calib_lines(bm.mod))
#   mySpExpl <- mySpExpl[which(calib.lines[,dataset] == TRUE), ]
#   
#   bm.options <- get_options(bm.mod)
#   bm.options <- bm.options@options[[grep(algo, names(bm.options@options))]] ### Warning grep !! 
#   
#   ## Boucle 
#   sample_points <- sample(1:nrow(mySpExpl), nrow(mySpExpl)*perc)
#   model.call <- paste0(bm.options@package, "::", bm.options@func)
#   
#   residuals_predicted <- foreach(point = sample_points, .combine = c) %do% {
#     argstmp <- bm.options@args.values[[dataset]]
#     argstmp$data <- mySpExpl[-point, ]
#     if (algo %in% c("MARS", "RF")){
#       formu <- bm_MakeFormula(resp.name = names(mySpExpl)[1], expl.var = mySpExpl[1, 4:ncol(mySpExpl)], type = "simple", interaction.level = 0)
#       argstmp$formula <- formu
#     }
#     # if(algo == "XGBOOST"){
#     #   argstmp$label <- argstmp$data[,1]
#     #   argstmp$data <- as.matrix(argstmp$data[, 4:ncol(mySpExpl)])
#     #   argstmp <- argstmp[c("data", names(argstmp)[which(!(names(argstmp) %in% c("formula", "data")))])]
#     # } 
#     argstmp <- argstmp[c("formula","data", names(argstmp)[which(!(names(argstmp) %in% c("formula", "data")))])]
#     
#     
#     new.model <- try(do.call(eval(parse(text = model.call)), argstmp), silent = TRUE)
#     new.res <- NA
#     if (inherits(new.model, 'try-error') != T){
#       new.pred <- predict(new.model, mySpExpl[point, ])
#       new.res <- mySpExpl[point, 1] - new.pred
#     } 
#     
#     return(new.res)
#   }## Return residuals x(-i)
#   
#   ## Calcul PRESS + R predicted
#   PRESS <- sum(residuals_predicted ^ 2)
#   return(PRESS)
# }
  
  
## TSS with bm_FindOptimStat MAIS on a déjà TSS dans bm_FOS !!
## Pourquoi tout ces acronymes sont les mêmes ?!?

# .bm_TotalSumSquares <- function(bm.mod, mod.name){
#   dataset <- paste0("_", unlist(strsplit(mod.name, "_"))[2:3], collapse = "")
#   bm.format <- get_formal_data(bm.mod)
#   mySpExpl <- get_species_data(bm.format)
#   calib.lines <- as.data.frame(get_calib_lines(bm.mod))
#   mySpExpl <- mySpExpl[which(calib.lines[,dataset] == TRUE), ]
#   
#   #fit <- tab_res[tab_res$full.name == mod.name, "fitted"]
#   y <- mySpExpl[, 1]
#   y_mean <- mean(y, na.rm = T)
#   TSS <- sum((y - y_mean) ^ 2)
#   return(TSS)
# }


###################################################################################################

# # .bm_return_info_analysys
# # Functions to get residuals and fitted values of a model
# # 
# # param mod a \code{biomod2_model}
# # importFrom stats fitted residuals
# # 
# # export
# 
# setGeneric(".bm_return_info_analysys", function(mod) { standardGeneric(".bm_return_info_analysys")})
# 
# setMethod(".bm_return_info_analysys", signature("CTA_biomod2_model"), function(mod) {
#   resids <- as.vector(residuals(mod@model)) 
#   fit <- as.vector(mod@model$y) #[,"y"]
#   fit <- fit-1
#   return(data.frame("obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
# })
# 
# setMethod(".bm_return_info_analysys", signature("GAM_biomod2_model"), function(mod) {
#   resids <- as.vector(residuals(mod@model)) 
#   fit <- as.vector(fitted(mod@model))
#   return(data.frame("obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
# })
# 
# setMethod(".bm_return_info_analysys", signature("GLM_biomod2_model"), function(mod) {
#   resids <- as.vector(residuals(mod@model)) 
#   fit <- as.vector(fitted(mod@model))
#   return(data.frame("obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
# })
# 
# setMethod(".bm_return_info_analysys", signature("MARS_biomod2_model"), function(mod) {
#   resids <- as.vector(residuals(mod@model)) 
#   fit <- as.vector(fitted(mod@model))
#   return(data.frame( "obs" = 1:length(resids), "residuals" = resids, "fitted" = fit))
# })




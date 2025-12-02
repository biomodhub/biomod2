###################################################################################################
##' @name bm_ModelAnalysis
##' @author Hélène Blancheteau
##' 
##' @title Analyze the residuals of the single models 
##' 
##' @description This function return several graphs to help analyse the single models.
##' 
##' 
##' @param bm.mod a \code{\link{BIOMOD.models.out}} object returned by the 
##' \code{\link{BIOMOD_Modeling}} function
##' @param models.chosen a \code{vector} containing model names to be kept, must be either 
##' \code{all} or a sub-selection of model names that can be obtained with the 
##' \code{\link{get_built_models}} function applied to \code{bm.mod}
##' 
##' @param color.by a \code{character} corresponding to the way plots will be colored, must be 
##' among \code{full.name}, \code{PA}, \code{run} or \code{algo} 
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' 
##' 
##' @return
##' 
##' A \code{list} containing :
##' \enumerate{
##'   \item a \code{data.frame} with variables, predicted values and model residuals
##'   \item a \code{data.frame} with the R score values
##'   \item 5 \code{ggplot} objects (from A to E) representing the different analyses.
##' }
##' 
##' 
##' @details 
##' 
##' 5 plots can be obtained with this function :
##' \describe{
##'   \item{residuals ~ observations number}{to detect outliers. The x-axis only helps to find 
##'   the outlier number.}
##'   \item{residuals Q-Q plot}{to check if residuals follow a normal distribution. Points 
##'   should follow the black line.}
##'   \item{residuals histogram}{to check if residuals follow a normal distribution. Histogram 
##'   should represent a gaussian distribution.}
##'   \item{residuals ~ fitted values}{to detect an heteroscedasticity of the residuals. It 
##'   corresponds to the Tukey-Anscombe plot.}
##'   \item{R scores}{to detect overfitting. Representing \code{Rsquared} and 
##'   \code{Rsquared_aj} values, there should not be a big gap between calibration and 
##'   validation values.}
##' }
##' 
##' \emph{Please not the all plots are made for all models, independently to the different 
##' models assumptions. It is up to the user to interpret the graphics.}
##' 
##' 
##' @keywords analyze models residuals
##' 
##' @family Secondary functions
##' @family Plot functions
##' 
##' 
##' @examples
##' library(terra)
##' 
##' # Load species occurrences (6 species available)
##' data(DataSpecies)
##' head(DataSpecies)
##' 
##' # Select the name of the studied species
##' myRespName <- 'GuloGulo'
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
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExpl <- terra::crop(myExpl, myExtent)
##' }
##' 
##' # ---------------------------------------------------------------#
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
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS','AUCroc'),
##'                                       var.import = 3,
##'                                       seed.val = 42)
##' }
##' 
##' 
##' # ---------------------------------------------------------------#
##' # Explore single models
##' myBiomodAnalysis <- bm_ModelAnalysis(bm.mod = myBiomodModelOut, color.by = "run")
##' 
##' plot(myBiomodAnalysis$plot.outliers)
##' 
##' 
##' 
##' @importFrom foreach foreach %do%
##' @importFrom ggplot2 ggplot geom_hline geom_point geom_histogram stat_qq stat_qq_line vars 
##' 
##' @export
##' 
##' 
###################################################################################################


bm_ModelAnalysis <- function(bm.mod, 
                             models.chosen = 'all',
                             color.by = 'full.name',
                             do.plot = TRUE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_ModelAnalysis.check.args(bm.mod = bm.mod, models.chosen = models.chosen,
                                       color.by = color.by, do.plot = do.plot)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Compute RESIDUALS -------------------------------------------------------------------------
  tab.pred <- get_predictions(bm.mod, full.name = models.chosen)
  tab.obs <- data.frame(points = 1:max(tab.pred$points),
                        obs = get_formal_data(bm.mod, subinfo = "resp.var"))
  
  if (bm.mod@data.type == "binary") {
    tab.pred$pred <- tab.pred$pred / 1000
  }
  
  if (bm.mod@data.type == "ordinal") {
    names_levels <- levels(tab.obs$obs)
    
    tab.pred$pred <- as.character(tab.pred$pred)
    to_change <- tab.pred$pred[tab.pred$pred %in% names_levels] 
    to_change <- factor(to_change, levels = levels(tab.obs$obs), ordered = TRUE)
    tab.pred$pred[tab.pred$pred %in% names_levels] <- as.character(as.numeric(to_change))
    
    tab.pred$pred <- as.numeric(tab.pred$pred)
    tab.obs$obs <- as.numeric(tab.obs$obs)
  }
  
  ggdat <- merge(tab.pred, tab.obs, by = "points")
  ggdat$residuals <- ggdat$obs - ggdat$pred
  
  
  ## 2. Create PLOTS ------------------------------------------------------------------------------
  
  ## a. residuals ~ observations ----------------------------------------------
  gg.outliers <- ggplot(ggdat, aes(y = residuals, x = points, color = .data[[color.by]])) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    labs(x = "Observations number", y = "Residuals"
         , title = "Visualizing residuals: are there any outliers?") +
    paletteer::scale_color_paletteer_d(palette) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ## b. residuals Q-Q plot ----------------------------------------------------
  gg.qqplot <- ggplot(ggdat, aes(sample = residuals, color = .data[[color.by]])) +
    stat_qq(size = 1) +
    stat_qq_line(color = "black", size = 1) +
    labs(x = "Theoretical quantiles", y = "Standardized residuals"
         , title = "Q-Q plot of residuals") +
    paletteer::scale_color_paletteer_d(palette) +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ## c. residuals histogram ---------------------------------------------------
  gg.hist <- ggplot(ggdat, aes(x = residuals, color = .data[[color.by]])) +
    geom_histogram(fill = NA, position = "dodge", bins = 30) +
    facet_wrap(vars(full.name)) +
    labs(x = "Residuals", ylab = "", title = "Distribution of residuals") +
    paletteer::scale_color_paletteer_d(palette)+
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1),
          , axis.text.y = element_blank(),
          , axis.ticks.y = element_blank())


  ## d. residuals ~ fitted values ---------------------------------------------
  gg.fitted <- ggplot(ggdat, aes(y = residuals, x = pred, color = .data[[color.by]])) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    labs(x = "Fitted values", y = "Residuals"
         , title = "Residuals ~ Fitted values (Tukey-Anscombe plot)") +
    paletteer::scale_color_paletteer_d(palette) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  ## e. R scores---------------------------------------------------------------
  if (!bm.mod@data.type %in% c("binary", "ordinal")) {
    
    tab.eval <- get_evaluations(bm.mod)
    tab.eval <- tab.eval[which(tab.eval$full.name %in% models.chosen), ]
    
    Rscore <- foreach (met = c("Rsquared", "Rsquared_aj"), .combine = "rbind") %do%
      {
        if (met %in% unique(tab.eval$metric.eval)) {
          tab <- tab.eval[which(tab.eval$metric.eval == met), c("full.name", "calibration", "validation")]
          return(rbind(data.frame(full.name = tab$full.name
                                  , metric.eval = met
                                  , dataset = "calibration"
                                  , value = tab$calibration)
                       , data.frame(full.name = tab$full.name
                                    , metric.eval = met
                                    , dataset = "validation"
                                    , value = tab$validation)))
        }
      }
    # Rscore <- Rscore[which(!is.na(Rscore$value)), ]
    
    gg.Rscore <- ggplot(Rscore, aes(y = value, x = full.name, color = .data[[color.by]], shape = metric.eval)) +
      geom_point(size = 2) +
      ylim(0, 1) +
      labs(x = "Models", y = "", title = "R scores") +
      paletteer::scale_color_paletteer_d(palette) +
      theme(legend.title = element_blank()
            , legend.key = element_rect(fill = "white")
            , axis.text.x = element_text(angle = 45, hjust = 1),
            , axis.title.y = element_blank())

  } else {
    gg.Rscore <- Rscore <- NULL
  }
  
  ## RETURN PLOTS
  if (do.plot) {
    print(gg.outliers + gg.qqplot + gg.hist + gg.fitted + 
            patchwork::plot_annotation(tag_levels = 'A', tag_prefix = "(", tag_suffix = ")"))
    if (!bm.mod@data.type %in% c("binary", "ordinal")) {
      print(gg.Rscore)
    }
  }
  
  return(list(tab.residuals = ggdat,
              tab.Rscore = Rscore,
              plot.outliers = invisible(gg.outliers),
              plot.qqplot = invisible(gg.qqplot),
              plot.hist = invisible(gg.hist),
              plot.fitted = invisible(gg.fitted),
              plot.Rscore = invisible(gg.Rscore)))
}




###################################################################################################

.bm_ModelAnalysis.check.args <- function(bm.mod, models.chosen, color.by, do.plot)
{
  ## check namespace ----------------------------------------------------------
  if (!isNamespaceLoaded("paletteer")) {
    if (!requireNamespace('paletteer', quietly = TRUE)) stop("Package 'paletteer' not found")
  }
  if(do.plot) {
    if (!isNamespaceLoaded("patchwork")) {
      if (!requireNamespace('patchwork', quietly = TRUE)) stop("Package 'patchwork' not found")
    }
  }
  
  
  ## 1. Check bm.mod argument -------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.mod", bm.mod, "BIOMOD.models.out")
  
  if (bm.mod@data.type == "multiclass") {
    stop("bm_ModelAnalysis is not available for multiclass data models.")
  }
  
  ## 2. Check models.chosen argument ------------------------------------------
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
    stop(paste0("Projection files missing : ", toString(missing.files)))
    if (length(missing.files) == length(files.check)) {
      stop("Impossible to find any models, might be a problem of working directory")
    }
  }
  
  ## 3. Check color.by argument -----------------------------------------------
  .fun_testIfIn(TRUE, "color.by", color.by, c("full.name", "PA", "run", "algo" ))
  palette <- switch(color.by
                    , 'algo' = "tvthemes::parksAndRec"
                    , 'run' = "ggthemes::Classic_20"
                    , 'PA' = "ggthemes::Tableau_20"
                    , 'full.name' = "ggsci::default_igv")
  
  return(list(models.chosen = models.chosen
              , palette = palette))
}


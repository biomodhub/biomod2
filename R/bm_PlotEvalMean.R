###################################################################################################
##' @name bm_PlotEvalMean
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot mean evaluation scores
##' 
##' @description This function represents mean evaluation scores (and their standard deviation) 
##' of species distribution models, from \code{\link{BIOMOD.models.out}} or 
##' \code{\link{BIOMOD.ensemble.models.out}} objects that can be obtained from 
##' \code{\link{BIOMOD_Modeling}} or \code{\link{BIOMOD_EnsembleModeling}} functions. Scores are 
##' represented according to 2 different evaluation methods, and models can be grouped 
##' (see Details).
##' 
##' 
##' @param bm.out a \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' object that can be obtained with the \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param metric.eval a \code{vector} containing evaluation metric names to be used, must 
##' be among \code{ROC}, \code{TSS}, \code{KAPPA}, \code{ACCURACY}, \code{BIAS}, \code{POD}, 
##' \code{FAR}, \code{POFD}, \code{SR}, \code{CSI}, \code{ETS}, \code{HK}, \code{HSS}, \code{OR}, 
##' \code{ORSS}
##' @param group.by a \code{character} corresponding to the way kept models will be combined to 
##' compute mean and sd evaluation scores, must be among \code{model}, \code{algo}, \code{run}, 
##' \code{dataset}
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see Details)
##' 
##' 
##' @return  
##' 
##' A \code{list} containing a \code{data.frame} with mean and standard deviation of evaluation 
##' scores and the corresponding \code{ggplot} object representing them according to 2 different 
##' evaluation methods.
##' 
##' 
##' @details
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item{\code{xlim}}{ : an \code{integer} corresponding to the x maximum limit to represent}
##'   \item{\code{ylim}}{ : an \code{integer} corresponding to the y maximum limit to represent}
##'   \item{\code{main}}{ : a \code{character} corresponding to the graphic title}
##'   \item{\code{col}}{ : a \code{vector} containing new color values}
##' }
##' 
##' 
##' @keywords evaluation ggplot
##' 
##' 
##' @seealso \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.ensemble.models.out}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{get_evaluations}}
##' @family Secundary functions
##' @family Plot functions
##' 
##' 
##' @examples
##' 
##' # Load species occurrences (6 species available)
##' myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
##' DataSpecies <- read.csv(myFile, row.names = 1)
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
##' myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl <- raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' \dontshow{
##' myExtent <- raster::extent(0,30,45,70)
##' myExpl <- raster::stack(raster::crop(myExpl, myExtent))
##' }
##' 
##' # ---------------------------------------------------------------
##' file.out <- paste0(myRespName, "/", myRespName, ".AllModels.models.out")
##' if (file.exists(file.out)) {
##'   myBiomodModelOut <- get(load(file.out))
##' } else {
##' 
##'   # Format Data with true absences
##'   myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                        expl.var = myExpl,
##'                                        resp.xy = myRespXY,
##'                                        resp.name = myRespName)
##' 
##'   # Create default modeling options
##'   myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       bm.options = myBiomodOptions,
##'                                       nb.rep = 2,
##'                                       data.split.perc = 80,
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
##'                                       do.full.models = FALSE,
##'                                       seed.val = 42)
##' }
##' 
##' 
##' # ---------------------------------------------------------------
##' # Get evaluation scores
##' get_evaluations(myBiomodModelOut)
##' 
##' # Represent mean evaluation scores
##' bm_PlotEvalMean(bm.out = myBiomodModelOut)
##' 
##' 
##' @importFrom ggplot2 ggplot aes_string geom_point geom_errorbarh geom_errorbar xlab ylab
##' theme element_blank element_rect coord_cartesian labs
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotEvalMean <- function(bm.out, metric.eval = NULL, group.by = 'algo', do.plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotEvalMean.check.args(bm.out, metric.eval, group.by, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  tmp = strsplit(group.by, '')[[1]]
  group.by <- paste0(toupper(tmp[1]), paste0(tmp[2:length(tmp)], collapse = ''))
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get evaluation values
  scores <- get_evaluations(bm.out, as.data.frame = TRUE)
  scores$Eval.metric <- as.character(scores$Eval.metric)
  
  ## Choose which dataset (calibration or validation) should be used
  eval.data <- ifelse(all(is.na(scores$Evaluating.data)), "Testing.data", "Evaluating.data")
  
  ## Compute mean and sd evaluation scores
  models_mean = tapply(X = scores[, eval.data]
                       , INDEX = list(scores$Eval.metric, scores[, group.by])
                       , FUN = mean, na.rm = TRUE)
  models_sd = tapply(X = scores[, eval.data]
                     , INDEX = list(scores$Eval.metric, scores[, group.by])
                     , FUN = sd, na.rm = TRUE)
  
  ## Prepare data table for graphic
  ggdat <- merge(data.frame(name = colnames(models_mean), t(models_mean)),
                 data.frame(name = colnames(models_sd), t(models_sd)), 
                 by = "name" )
  colnames(ggdat) <- c("name", "mean1", "mean2", "sd1", "sd2")
  
  limits1 <- aes_string(xmax = "mean1 + sd1", xmin = "mean1 - sd1", fill = NULL)
  limits2 <- aes_string(ymax = "mean2 + sd2", ymin = "mean2 - sd2", fill = NULL)
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  gg <- ggplot(ggdat, aes_string(x = "mean1", y = "mean2", colour = "name", fill = NULL)) +
    geom_point() + ## add mean points
    geom_errorbarh(limits1, height = 0) + ## add horizontal error bars
    geom_errorbar(limits2, width = 0) + ## add vertical error bars
    xlab(metric.eval[1]) +
    ylab(metric.eval[2]) +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white"))
  
  if (length(ylim) | length(xlim)) { ## fix scale
    gg <- gg + coord_cartesian(ylim = ylim, xlim = xlim)
  }
  
  if (length(main)) { ## add title
    gg <- gg + labs(title = main)
  }
  
  if (do.plot){ print(gg) }
  return(list(tab = ggdat, plot = invisible(gg)))
}


###################################################################################################

.bm_PlotEvalMean.check.args <- function(bm.out, metric.eval = NULL, group.by = 'Algo', ...)
{
  args <- list(...)
  
  ## 1. Check bm.out argument -------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check metric.eval argument --------------------------------------------
  if (is.null(metric.eval)) {
    scores <- get_evaluations(bm.out, as.data.frame = TRUE)
    metric.eval <- sort(unique(as.character(scores$Eval.metric)))[1:2]
    warnings(toString(metric.eval), " evaluation metric.eval automatically selected")
  } else {
    metric.eval = sort(unique(as.character(metric.eval)))
    if (length(metric.eval) < 2) {
      stop("2 different evaluations metric.eval needed")
    } else if (length(metric.eval) > 2) {
      metric.eval = metric.eval[1:2]
      warning("2 different evaluations metric.eval needed, only the first 2 will be kept")
    }
  }
  
  ## 3. Check group.by argument -----------------------------------------------
  .fun_testIfIn(TRUE, "group.by", group.by, c('model', 'algo', 'run', 'dataset'))
  
  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('xlim', 'ylim', 'main'))
  
  
  return(list(metric.eval = metric.eval,
              xlim = args$xlim,
              ylim = args$ylim,
              main = args$main))
} 


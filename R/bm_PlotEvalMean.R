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
##' be among \code{POD}, \code{FAR}, \code{POFD}, \code{SR}, \code{ACCURACY}, \code{BIAS}, 
##' \code{ROC}, \code{TSS}, \code{KAPPA}, \code{OR}, \code{ORSS}, \code{CSI}, \code{ETS}, 
##' \code{BOYCE}, \code{MPA}
##' @param dataset a \code{character} corresponding to the dataset upon which evaluation metrics 
##' have been calculated and that is to be represented, must be among \code{calibration}, 
##' \code{validation}, \code{evaluation}
##' @param group.by a \code{character} corresponding to the way kept models will be combined to 
##' compute mean and sd evaluation scores, must be among \code{full.name}, \code{PA}, \code{run}, 
##' \code{algo} (if \code{bm.out} is a \code{\link{BIOMOD.models.out}} object), or 
##' \code{full.name}, \code{merged.by.PA}, \code{merged.by.run}, \code{merged.by.algo} 
##' (if \code{bm.out} is a \code{\link{BIOMOD.ensemble.models.out}} object)
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
##'   \item \code{xlim} : an \code{integer} corresponding to the x maximum limit to represent
##'   \item \code{ylim} : an \code{integer} corresponding to the y maximum limit to represent
##'   \item \code{main} : a \code{character} corresponding to the graphic title
##'   \item \code{col} : a \code{vector} containing new color values
##' }
##' 
##' 
##' @keywords evaluation ggplot
##' 
##' 
##' @seealso \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.ensemble.models.out}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{get_evaluations}}
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
##'   # Model single models
##'   myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                       modeling.id = 'AllModels',
##'                                       models = c('RF', 'GLM'),
##'                                       CV.strategy = 'random',
##'                                       CV.nb.rep = 2,
##'                                       CV.perc = 0.8,
##'                                       OPT.strategy = 'bigboss',
##'                                       metric.eval = c('TSS','ROC'),
##'                                       var.import = 3,
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


bm_PlotEvalMean <- function(bm.out, metric.eval = NULL, dataset = 'calibration', group.by = 'algo', do.plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotEvalMean.check.args(bm.out, metric.eval, dataset, group.by, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  if (!is.null.out) {
    ## 1. Get data for graphic ----------------------------------------------------------------------
    ## Get evaluation values
    scores <- get_evaluations(bm.out)
    scores <- scores[scores$metric.eval %in% metric.eval, ]
    
    ## Compute mean and sd evaluation scores
    models_mean = tapply(X = scores[, dataset]
                         , INDEX = list(scores$metric.eval, scores[, group.by])
                         , FUN = mean, na.rm = TRUE)
    models_sd = tapply(X = scores[, dataset]
                       , INDEX = list(scores$metric.eval, scores[, group.by])
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
    
    if (length(ylim) > 0 | length(xlim) > 0) { ## fix scale
      gg <- gg + coord_cartesian(ylim = ylim, xlim = xlim)
    }
    
    if (length(main) > 0) { ## add title
      gg <- gg + labs(title = main)
    }
    
    if (do.plot){ print(gg) }
    return(list(tab = ggdat, plot = invisible(gg)))
  }
}


###################################################################################################

.bm_PlotEvalMean.check.args <- function(bm.out, metric.eval = NULL, dataset, group.by, ...)
{
  args <- list(...)
  
  ## 1. Check bm.out argument -------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check metric.eval argument --------------------------------------------
  scores <- get_evaluations(bm.out)
  
  if (!is.null(scores))
  {
    avail.metrics <- sort(unique(as.character(scores$metric.eval)))
    if (is.null(metric.eval) && length(avail.metrics) > 1) {
      metric.eval <- sort(unique(as.character(scores$metric.eval)))[1:2]
      warning(toString(metric.eval), " evaluation metric.eval automatically selected")
    } else {
      metric.eval = sort(unique(as.character(metric.eval)))
      if (length(metric.eval) < 2) {
        stop("2 different evaluations metric.eval needed")
      } else if (length(metric.eval) > 2) {
        metric.eval = metric.eval[1:2]
        warning("2 different evaluations metric.eval needed, only the first 2 will be kept")
      }
    }
    
    ## 2. Check dataset argument ------------------------------------------------
    .fun_testIfIn(TRUE, "dataset", dataset, c("calibration", "validation", "evaluation"))
    
    ## 3. Check group.by argument -----------------------------------------------
    if (inherits(bm.out, "BIOMOD.models.out")) {
      .fun_testIfIn(TRUE, "group.by", group.by, c("full.name", "PA", "run", "algo"))
    } else if (inherits(bm.out, "BIOMOD.ensemble.models.out")) {
      .fun_testIfIn(TRUE, "group.by", group.by, c("full.name", "merged.by.PA", "merged.by.run", "algo"))
    } 
    if (length(group.by) > 1) {
      group.by = group.by[1]
      warning("`group.by` must contain only one value, only the first one will be kept")
    }
    
    ## 4. Check extra args argument ---------------------------------------------
    .fun_testIfIn(TRUE, "names(args)", names(args), c('xlim', 'ylim', 'main'))
  }
  
  return(list(is.null.out = is.null(scores),
              metric.eval = metric.eval,
              group.by = group.by,
              xlim = args$xlim,
              ylim = args$ylim,
              main = args$main))
} 


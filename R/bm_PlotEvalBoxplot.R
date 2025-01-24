###################################################################################################
##' @name bm_PlotEvalBoxplot
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot boxplot of evaluation scores
##' 
##' @description This function represents boxplot of evaluation scores of species distribution 
##' models, from \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' objects that can be obtained from \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions. Scores are represented according to 2 
##' grouping methods (see Details).
##' 
##' 
##' @param bm.out a \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' object that can be obtained with the \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param dataset a \code{character} corresponding to the dataset upon which evaluation metrics 
##' have been calculated and that is to be represented, must be among \code{calibration}, 
##' \code{validation}, \code{evaluation}
##' @param group.by a 2-length \code{vector} containing the way kept models will be represented,
##' must be among \code{full.name}, \code{PA}, \code{run}, \code{algo} (if \code{bm.out} is a 
##' \code{\link{BIOMOD.models.out}} object), or \code{full.name}, \code{merged.by.PA}, 
##' \code{merged.by.run}, \code{merged.by.algo} (if \code{bm.out} is a 
##' \code{\link{BIOMOD.ensemble.models.out}} object)
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see Details)
##' 
##' 
##' @return  
##' 
##' A \code{list} containing a \code{data.frame} with evaluation scores and the corresponding 
##' \code{ggplot} object representing them in boxplot.
##' 
##' 
##' @details
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item \code{main} : a \code{character} corresponding to the graphic title
##'   \item \code{scales} : a \code{character} corresponding to the \code{scales} argument of 
##'   the \code{\link[ggplot2]{facet_wrap}} function, must be either \code{fixed}, \code{free_x}, 
##'   \code{free_y} or \code{free}
##' }
##' 
##' 
##' @keywords evaluation ggplot boxplot
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
##' # Represent evaluation scores
##' bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'run'))
##' 
##' 
##' @importFrom ggplot2 ggplot aes_string geom_boxplot facet_wrap xlab 
##' theme element_blank element_rect element_text labs
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotEvalBoxplot <- function(bm.out, dataset = 'calibration', group.by = c('algo', 'run'), do.plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotEvalBoxplot.check.args(bm.out, dataset, group.by, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get evaluation values
  scores <- get_evaluations(bm.out)
  
  if (!is.null(scores))
  {
    ## Prepare data table for graphic
    ggdat = scores
    
    ## 2. PLOT graphic ------------------------------------------------------------------------------
    gg <- ggplot(ggdat, aes_string(x = group.by[1], y = dataset, fill = group.by[2])) +
      geom_boxplot() + ## add boxplot
      facet_wrap("metric.eval", scales = scales) +
      xlab("") +
      theme(legend.title = element_blank()
            , legend.key = element_rect(fill = "white")
            , axis.text.x = element_text(angle = 45, hjust = 1))
    
    if (length(main) > 0) { ## add title
      gg <- gg + labs(title = main)
    }
    
    if (do.plot){ print(gg) }
    return(list(tab = ggdat, plot = invisible(gg)))
  }
}


###################################################################################################

.bm_PlotEvalBoxplot.check.args <- function(bm.out, dataset, group.by, ...)
{
  args <- list(...)
  
  ## 1. Check bm.out argument -------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 2. Check dataset argument ------------------------------------------------
  .fun_testIfIn(TRUE, "dataset", dataset, c("calibration", "validation", "evaluation"))
  
  ## 3. Check group.by argument -----------------------------------------------
  if (length(group.by) != 2) { stop("2 group values needed") }
  if (inherits(bm.out, "BIOMOD.models.out")) {
    for (i in 1:length(group.by)) {
      .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c("full.name", "PA", "run", "algo"))
    }
  } else if (inherits(bm.out, "BIOMOD.ensemble.models.out")) {
    for (i in 1:length(group.by)) {
      .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c("full.name", "merged.by.PA", "merged.by.run", "algo"))
    }
  } 
  
  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('main', 'scales'))
  if ("scales" %in% names(args)) {
    .fun_testIfIn(TRUE, "args$scales", args$scales, c('fixed', 'free_x', 'free_y', 'free'))
  } else {
    args$scales = "fixed"
  }
  
  
  return(list(main = args$main,
              scales = args$scales))
} 


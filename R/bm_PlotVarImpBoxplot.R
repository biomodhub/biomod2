###################################################################################################
##' @name bm_PlotVarImpBoxplot
##' @author Damien Georges, Maya Gueguen
##' 
##' @title Plot boxplot of variables importance
##' 
##' @description This function represents boxplot of variables importance of species distribution 
##' models, from \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' objects that can be obtained from \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions. Scores are represented according to 3 
##' grouping methods (see Details).
##' 
##' 
##' @param bm.out a \code{\link{BIOMOD.models.out}} or \code{\link{BIOMOD.ensemble.models.out}} 
##' object that can be obtained with the \code{\link{BIOMOD_Modeling}} or 
##' \code{\link{BIOMOD_EnsembleModeling}} functions
##' @param group.by a 3-length \code{vector} containing the way kept models will be represented,
##' must be among \code{full.name}, \code{PA}, \code{run}, \code{algo}, \code{expl.var} (if 
##' \code{bm.out} is a \code{\link{BIOMOD.models.out}} object), or \code{full.name}, 
##' \code{merged.by.PA}, \code{merged.by.run}, \code{merged.by.algo}, \code{expl.var} 
##' (if \code{bm.out} is a \code{\link{BIOMOD.ensemble.models.out}} object)
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plot is to be rendered or not
##' @param \ldots some additional arguments (see Details)
##' 
##' 
##' @return  
##' 
##' 
##' A \code{list} containing a \code{data.frame} with variables importance and the corresponding 
##' \code{ggplot} object representing them in boxplot.
##' 
##' 
##' @details
##' 
##' \code{...} can take the following values :
##' 
##' \itemize{
##'   \item \code{main} : a \code{character} corresponding to the graphic title
##' }
##' 
##' 
##' @keywords evaluation ggplot boxplot
##' 
##' 
##' @seealso \code{\link{BIOMOD.models.out}}, \code{\link{BIOMOD.ensemble.models.out}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{BIOMOD_EnsembleModeling}}, 
##' \code{\link{get_variables_importance}}
##' @family Secundary functions
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
##' # Get variables importance
##' get_variables_importance(myBiomodModelOut)
##' 
##' # Represent variables importance
##' bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
##' bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'PA'))
##' bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'PA'))
##' 
##' 
##' @importFrom ggplot2 ggplot aes_string geom_boxplot facet_wrap xlab ylab labs 
##' theme element_blank element_rect element_text scale_y_continuous
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotVarImpBoxplot <- function(bm.out, group.by = c('run', 'expl.var', 'algo'), do.plot = TRUE, ...)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_PlotVarImpBoxplot.check.args(bm.out, group.by, ...)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get variables importance values & Prepare data table for graphic
  ggdat = get_variables_importance(bm.out)
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  gg <- ggplot(ggdat, aes_string(x = group.by[1], y = "var.imp", fill = group.by[2])) +
    geom_boxplot() + ## add boxplot
    facet_wrap(group.by[3], scales = "free_x") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%")) + 
    xlab("") +
    ylab("") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (length(main) > 0) { ## add title
    gg <- gg + labs(title = main)
  }
  
  if (do.plot){ print(gg) }
  return(list(tab = ggdat, plot = invisible(gg)))
}


###################################################################################################

.bm_PlotVarImpBoxplot.check.args <- function(bm.out, group.by, ...)
{
  args <- list(...)
  
  ## 1. Check bm.out argument -------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 3. Check group.by argument -----------------------------------------------
  if (length(group.by) != 3) { stop("3 group values needed") }
  if (inherits(bm.out, "BIOMOD.models.out")) {
    for (i in 1:length(group.by)) {
      .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c("full.name", "PA", "run", "algo", "expl.var"))
    }
  } else if (inherits(bm.out, "BIOMOD.ensemble.models.out")) {
    for (i in 1:length(group.by)) {
      .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c("full.name", "merged.by.PA", "merged.by.run", "algo", "expl.var"))
    }
  } 
  
  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('main'))
  
  return(list(main = args$main))
} 


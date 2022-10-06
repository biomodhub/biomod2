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
##' must be among \code{model}, \code{algo}, \code{run}, \code{dataset}, \code{expl.var}
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
##'   \item{\code{main}}{ : a \code{character} corresponding to the graphic title}
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
##' # Get variables importance
##' get_variables_importance(myBiomodModelOut, as.data.frame = TRUE)
##' 
##' # Represent variables importance
##' bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
##' bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('expl.var', 'algo', 'dataset'))
##' bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c('algo', 'expl.var', 'dataset'))
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
  
  
  for (i in 1:length(group.by)) {
    tmp = strsplit(group.by[i], '')[[1]]
    group.by[i] <- paste0(toupper(tmp[1]), paste0(tmp[2:length(tmp)], collapse = ''))
  }
  
  ## 1. Get data for graphic ----------------------------------------------------------------------
  ## Get variables importance values
  scores <- get_variables_importance(bm.out, as.data.frame = TRUE)
  
  ## Prepare data table for graphic
  ggdat = scores
  
  ## 2. PLOT graphic ------------------------------------------------------------------------------
  gg <- ggplot(ggdat, aes_string(x = group.by[1], y = "Var.imp", fill = group.by[2])) +
    geom_boxplot() + ## add boxplot
    facet_wrap(group.by[3], scales = "free_x") +
    scale_y_continuous(breaks = seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%")) + 
    xlab("") +
    ylab("") +
    theme(legend.title = element_blank()
          , legend.key = element_rect(fill = "white")
          , axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (length(main)) { ## add title
    gg <- gg + labs(title = main)
  }
  
  if (do.plot){ print(gg) }
  return(list(tab = ggdat, plot = invisible(gg)))
}


###################################################################################################

.bm_PlotVarImpBoxplot.check.args <- function(bm.out, group.by = 'Algo', ...)
{
  args <- list(...)
  
  ## 1. Check bm.out argument -------------------------------------------------
  .fun_testIfInherits(TRUE, "bm.out", bm.out, c("BIOMOD.models.out", "BIOMOD.ensemble.models.out"))
  
  ## 3. Check group.by argument -----------------------------------------------
  if (length(group.by) != 3) { stop("3 group values needed") }
  for (i in 1:length(group.by)) {
    .fun_testIfIn(TRUE, paste0("group.by[", i, "]"), group.by[i], c('model', 'algo', 'run', 'dataset', 'expl.var'))
  }
  
  ## 4. Check extra args argument ---------------------------------------------
  .fun_testIfIn(TRUE, "names(args)", names(args), c('main'))
  
  
  return(list(main = args$main))
} 


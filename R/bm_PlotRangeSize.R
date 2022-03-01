###################################################################################################
##' @name bm_PlotRangeSize
##' @author Maya Gueguen
##' 
##' @title Plot species range change
##' 
##' @description This function represents species range change from object that can be obtained 
##' from \code{\link{BIOMOD_RangeSize}} function. Several graphics can be obtained, representing 
##' global counts or proportions of gains / losses, as well as spatial representations (see 
##' \href{bm_PlotRangeSize.html#details}{Details}).
##' 
##' 
##' @param bm.range an object returned by the \code{\link{BIOMOD_RangeSize}} function
##' @param do.count (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the count plot is to be computed or not
##' @param do.perc (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the percentage plot is to be computed or not
##' @param do.maps (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the maps plot is to be computed or not
##' @param do.mean (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the mean maps plot is to be computed or not
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether the plots are to be rendered or not
##' 
##' 
##' @return  
##' 
##' A \code{list} containing \code{ggplot} objects representing species range change.
##' 
##' 
##' @details
##' 
##' 4 plots can be obtained with this function :
##' \describe{
##'   \item{Count barplot}{representing absolute number of locations (pixels) lost, stable and 
##'   gained}
##'   \item{Percentage barplot}{representing percentage of locations (pixels) lost, stable, and 
##'   the corresponding Species Range Change (\code{PercGain - PercLoss})}
##'   \item{SRC models maps}{representing spatially locations (pixels) lost, stable and 
##'   gained for each single distribution model}
##'   \item{SRC community averaging maps}{representing spatially locations (pixels) lost, stable 
##'   and gained, taking the majoritary value across single distribution models (and representing 
##'   the percentage of models' agreement)}
##' }
##' \emph{Please see \code{\link{BIOMOD_RangeSize}} function for more details about the values.}
##' 
##' 
##' @keywords ggplot, species range change, projections, gain, loss
##' 
##' 
##' @seealso \code{\link{BIOMOD_RangeSize}}
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
##' 
##' # ---------------------------------------------------------------
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' 
##' # Create default modeling options
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##' 
##' # Model single models
##' myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
##'                                     modeling.id = 'AllModels',
##'                                     bm.options = myBiomodOptions,
##'                                     nb.rep = 2,
##'                                     data.split.perc = 80,
##'                                     metric.eval = c('TSS','ROC'),
##'                                     var.import = 3,
##'                                     do.full.models = FALSE)
##' 
##' # Project single models
##' myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                   proj.name = 'Current',
##'                                   new.env = myExpl,
##'                                   models.chosen = 'all',
##'                                   metric.binary = 'all',
##'                                   metric.filter = 'all',
##'                                   build.clamping.mask = TRUE)
##' 
##' 
##' # ---------------------------------------------------------------
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' myFiles = paste0('external/bioclim/future/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExplFuture = raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' # Project onto future conditions
##' myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                               proj.name = 'Future',
##'                                               new.env = myExplFuture,
##'                                               models.chosen = 'all',
##'                                               metric.binary = 'TSS',
##'                                               build.clamping.mask = TRUE)
##' 
##' # Load current and future binary projections
##' CurrentProj <- stack("GuloGulo/proj_Current/proj_Current_GuloGulo_TSSbin.grd")
##' FutureProj <- stack("GuloGulo/proj_Future/proj_Future_GuloGulo_TSSbin.grd")
##' 
##' # Compute differences
##' myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj)
##' 
##' 
##' # ---------------------------------------------------------------
##' myBiomodRangeSize$Compt.By.Models
##' plot(myBiomodRangeSize$Diff.By.Pixel)
##' 
##' # Represent main results 
##' bm_PlotRangeSize(bm.range = myBiomodRangeSize)
##' 
##' 
##' @importFrom reshape2 melt
##' @importFrom patchwork plot_layout
##' @importFrom ggplot2 ggplot aes_string geom_col geom_raster facet_wrap xlab ylab labs 
##' theme element_blank element_rect scale_fill_manual
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotRangeSize <- function(bm.range, do.count = TRUE, do.perc = TRUE, do.maps = TRUE, do.mean = TRUE, do.plot = TRUE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  if (!is.list(bm.range) ||
      (is.list(bm.range) && length(bm.range) != 2) ||
      (is.list(bm.range) && length(bm.range) == 2 && sum(names(bm.range) == c("Compt.By.Models", "Diff.By.Pixel")) != 2)) {
    stop("'bm.range' must be an object obtained by the BIOMOD_RangeSize function")
  }
  
  ## 1. Create PLOTS for Compt.By.Models ----------------------------------------------------------
  if (do.count || do.perc)
  {
    ggdat = as.data.frame(bm.range$Compt.By.Models)
    
    ## Get models information
    ggdat$full.name = rownames(ggdat)
    ggdat$species = .extract_modelNamesInfo(ggdat$full.name, "species")
    ggdat$data.set = .extract_modelNamesInfo(ggdat$full.name, "data.set")
    ggdat$run.eval = .extract_modelNamesInfo(ggdat$full.name, "run.eval")
    ggdat$models = .extract_modelNamesInfo(ggdat$full.name, "models")
    
    ## Rearrange data
    ggdat = melt(ggdat, measure.vars = c("species", "data.set", "run.eval", "models"), variable.name = "group.level", value.name = "group.value")
    ggdat$group.level = sapply(ggdat$group.level, function(x) 
      switch(x, 'species' = 'Species', 'data.set' = 'Dataset', 'run.eval' = 'Run', 'models' = 'Algo'))
    ggdat = melt(ggdat, measure.vars = c("Loss", "Stable0", "Stable1", "Gain"), variable.name = "count.level", value.name = "count.value")
    ggdat$PercLoss = ggdat$PercLoss * (-1)
    ggdat = melt(ggdat, measure.vars = c("PercLoss", "PercGain", "SpeciesRangeChange"), variable.name = "perc.level", value.name = "perc.value")
    ggdat = melt(ggdat, measure.vars = c("CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
                 , variable.name = "range.level", value.name = "range.value")
    
    ## a. Count plot ----------------------------------------------------------
    if (do.count) {
      gg.count = ggplot(ggdat[which(ggdat$count.level != "Stable0"), ]
                        , aes_string(x = "group.value", y = "count.value", fill = "count.level")) +
        geom_col(position = "stack") +
        facet_wrap("group.level", scales = "free") +
        scale_fill_manual("", values = c("Loss" = "#fc8d62"
                                         , "Gain" = "#66c2a5"
                                         , "Stable1" = "grey")) +
        xlab("") +
        ylab("Number of pixels\n") +
        theme(legend.title = element_blank()
              , legend.key = element_rect(fill = "white")
              , legend.position = "top")
    } else { gg.count = NULL }
    
    ## b. Percentage plot -----------------------------------------------------
    if (do.perc) {
      gg.perc = ggplot(ggdat, aes_string(x = "group.value", y = "perc.value", fill = "perc.level")) +
        geom_col(position = "dodge") +
        facet_wrap("group.level", scales = "free_x") +
        scale_fill_manual("", values = c("PercLoss" = "#fc8d62"
                                         , "PercGain" = "#66c2a5"
                                         , "SpeciesRangeChange" = "#8da0cb")) +
        xlab("") +
        ylab("Percentage (%)\n") +
        theme(legend.title = element_blank()
              , legend.key = element_rect(fill = "white")
              , legend.position = "top")
    } else { gg.perc = NULL }
  }
  
  
  ## 2. Create PLOTS for Diff.By.Pixel ------------------------------------------------------------
  if (do.maps || do.mean)
  {
    ggdat = bm.range$Diff.By.Pixel
    ggdat = as.data.frame(rasterToPoints(ggdat))
    
    ## Get models information
    ggdat = melt(ggdat, id.vars = c("x", "y"), variable.name = "full.name", value.name = "SRC")
    ggdat$full.name = as.character(ggdat$full.name)
    ggdat$species = .extract_modelNamesInfo(ggdat$full.name, "species")
    ggdat$data.set = .extract_modelNamesInfo(ggdat$full.name, "data.set")
    ggdat$run.eval = .extract_modelNamesInfo(ggdat$full.name, "run.eval")
    ggdat$models = .extract_modelNamesInfo(ggdat$full.name, "models")
    
    ## Rearrange data
    ggdat = melt(ggdat, measure.vars = c("species", "data.set", "run.eval", "models"), variable.name = "group.level", value.name = "group.value")
    ggdat$group.level = sapply(ggdat$group.level, function(x) 
      switch(x, 'species' = 'Species', 'data.set' = 'Dataset', 'run.eval' = 'Run', 'models' = 'Algo'))
    
    ## c. SRC maps per model --------------------------------------------------
    if (do.maps) {
      gg.maps = ggplot(ggdat, aes_string(x = "x", y = "y", fill = "as.factor(SRC)")) +
        geom_raster() +
        facet_wrap("full.name") +
        scale_fill_manual("", values = c("-2" = "#fc8d62"
                                         , "-1" = "grey"
                                         , "0" = "white"
                                         , "1" = "#66c2a5")
                          , labels = c("-2" = "Loss"
                                       , "-1" = "Stable1"
                                       , "0" = ""
                                       , "1" = "Gain")) +
        xlab("") +
        ylab("") +
        theme(legend.title = element_blank()
              , legend.key = element_rect(fill = "white")
              , legend.position = "top")
    } else { gg.maps = NULL }
    
    ## d. SRC mean maps per group.level ---------------------------------------
    if (do.mean) {
      
      tab1 = tapply(X = ggdat$SRC
                    , INDEX = list(paste0(ggdat$x, "_", ggdat$y)
                                   , ggdat$group.level
                                   , ggdat$group.value)
                    , FUN = function(x) names(table(x))[which.max(table(x))])
      tab1 = na.exclude(melt(tab1))
      tab1$x = sapply(tab1$Var1, function(x) as.numeric(strsplit(as.character(x), "_")[[1]][1]))
      tab1$y = sapply(tab1$Var1, function(x) as.numeric(strsplit(as.character(x), "_")[[1]][2]))
      colnames(tab1) = c("x_y", "group.level", "group.value", "SRC", "x", "y")
      
      tab2 = tapply(X = ggdat$SRC
                     , INDEX = list(paste0(ggdat$x, "_", ggdat$y)
                                    , ggdat$group.level
                                    , ggdat$group.value)
                     , FUN = function(x) max(table(x)) / sum(table(x)))
      tab2 = na.exclude(melt(tab2))
      tab2$x = sapply(tab2$Var1, function(x) as.numeric(strsplit(as.character(x), "_")[[1]][1]))
      tab2$y = sapply(tab2$Var1, function(x) as.numeric(strsplit(as.character(x), "_")[[1]][2]))
      colnames(tab2) = c("x_y", "group.level", "group.value", "PERC", "x", "y")
      
      tab2$PERC[which(tab2$PERC == 1 & tab1$SRC == 0)] = NA
      

      
      gg.ca1 = ggplot(tab1, aes_string(x = "x", y = "y", fill = "as.factor(SRC)")) +
        geom_raster() +
        facet_wrap("group.level") +
        scale_fill_manual("", values = c("-2" = "#fc8d62"
                                         , "-1" = "grey"
                                         , "0" = "white"
                                         , "1" = "#66c2a5")
                          , labels = c("-2" = "Loss"
                                       , "-1" = "Stable1"
                                       , "0" = ""
                                       , "1" = "Gain")) +
        xlab("") +
        ylab("") +
        labs(title = "Community averaging value across models") +
        theme(legend.title = element_blank()
              , legend.key = element_rect(fill = "white")
              , legend.position = "top")
      
      gg.ca2 = ggplot(tab2, aes_string(x = "x", y = "y", fill = "PERC")) +
        geom_raster() +
        facet_wrap("group.level") +
        scale_fill_viridis_c(""
                             , direction = -1
                             , limits = c(0, 1)
                             , na.value = "white"
                             , breaks = seq(0, 1, 0.5)
                             , labels = paste0(seq(0, 100, 50), "%")) +
        xlab("") +
        ylab("") +
        labs(title = "Percentage of models' agreement") +
        theme(legend.key = element_rect(fill = "white")
              , legend.position = "top")
      
      gg.ca = gg.ca1 | gg.ca2
    } else { gg.ca = NULL }
  }
  
  
  ## RETURN PLOTS
  if (do.plot) { 
    print(gg.count)
    print(gg.perc)
    print(gg.maps)
    print(gg.ca)
  }
  return(list(count = invisible(gg.count)
              , perc = invisible(gg.perc)
              , maps = invisible(gg.maps)
              , ca = invisible(gg.ca)))
}


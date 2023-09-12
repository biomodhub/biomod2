# bm_PlotRangeSize documentation ------------------------------------------------
##' @name bm_PlotRangeSize
##' @author Maya Gueguen
##' 
##' @title Plot species range change
##' 
##' @description This function represents species range change from object that can be obtained 
##' from \code{\link{BIOMOD_RangeSize}} function. Several graphics can be obtained, representing 
##' global counts or proportions of gains / losses, as well as spatial representations (see Details).
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
##' @param row.names (\emph{optional, default} \code{c('Species', 'Dataset', 'Run', 'Algo')}) \cr 
##' A \code{vector} containing tags matching \code{bm.range$Compt.By.Models} rownames splitted by 
##' '_' character
##' 
##' 
##' @return  
##' 
##' A \code{list} containing one or several \code{data.frame} and the corresponding 
##' \code{ggplot} object representing species range change.
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
##' @keywords ggplot "species range change" projections gain loss
##' 
##' 
##' @seealso \code{\link{BIOMOD_RangeSize}}
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
##' # ---------------------------------------------------------------#
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
##' models.proj <- get_built_models(myBiomodModelOut, algo = "RF")
##'   # Project single models
##'   myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                     proj.name = 'CurrentRangeSize',
##'                                     new.env = myExpl,
##'                                     models.chosen = models.proj,
##'                                     metric.binary = 'all')
##' 
##' 
##' 
##' # ---------------------------------------------------------------#
##' # Load environmental variables extracted from BIOCLIM (bio_3, bio_4, bio_7, bio_11 & bio_12)
##' data(bioclim_future)
##' myExplFuture <- terra::rast(bioclim_future)
##' \dontshow{
##' myExtent <- terra::ext(0,30,45,70)
##' myExplFuture <- terra::crop(myExplFuture, myExtent)
##' }
##' 
##' # Project onto future conditions
##' myBiomodProjectionFuture <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
##'                                               proj.name = 'FutureRangeSize',
##'                                               new.env = myExplFuture,
##'                                               models.chosen = models.proj,
##'                                               metric.binary = 'TSS')
##' 
##' # Load current and future binary projections
##' CurrentProj <- get_predictions(myBiomodProj,
##'                                metric.binary = "TSS",
##'                                model.as.col = TRUE)
##' FutureProj <- get_predictions(myBiomodProjectionFuture,
##'                                metric.binary = "TSS",
##'                                model.as.col = TRUE)
##' 
##' # Compute differences
##' myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = CurrentProj, proj.future = FutureProj)
##' 
##' 
##' # ---------------------------------------------------------------#
##' myBiomodRangeSize$Compt.By.Models
##' plot(myBiomodRangeSize$Diff.By.Pixel)
##' 
##' # Represent main results 
##' bm_PlotRangeSize(bm.range = myBiomodRangeSize)
##' 
##' 
##' @importFrom graphics plot.new
##' @importFrom reshape2 melt
##' @importFrom foreach foreach %do%
##' @importFrom terra rast which.max nlyr  classify plot
##' @importFrom ggplot2 ggplot aes_string geom_col geom_tile facet_wrap xlab ylab labs scale_fill_viridis_c
##' theme element_blank element_rect scale_fill_manual
##' 
##' @export
##' 
##' 
#------------------------------------------------------------------------------#


bm_PlotRangeSize <- function(bm.range, do.count = TRUE, do.perc = TRUE
                             , do.maps = TRUE, do.mean = TRUE, do.plot = TRUE
                             , row.names = c('Species', 'Dataset', 'Run', 'Algo'))
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  if (!is.list(bm.range) ||
      (is.list(bm.range) && length(bm.range) != 2) ||
      (is.list(bm.range) && length(bm.range) == 2 && !all(c("Compt.By.Models", "Diff.By.Pixel") %in% names(bm.range)))) {
    stop("'bm.range' must be an object obtained by the BIOMOD_RangeSize function")
  }
  out = list()
  
  ## 1. Create PLOTS for Compt.By.Models ----------------------------------------------------------
  if (do.count || do.perc)
  {
    ggdat = as.data.frame(bm.range$Compt.By.Models)
    
    ## Get models information
    ggdat$full.name = rownames(ggdat)
    {
      tmp = strsplit(ggdat$full.name[1], "_")[[1]][1:length(row.names)]
      names(tmp) = row.names
      warning(paste0("Please check that rownames(bm.range$Compt.By.Models) match 'row.names' argument :\n"
                     , paste0("\t", paste0(tmp, " : ", names(tmp)), collapse = "\n")))
    }
    
    for (ii in 1:length(row.names)) {
      ggdat[[row.names[ii]]] = sapply(ggdat$full.name, function(x) strsplit(x, "_")[[1]][ii])
    }
    
    ## Rearrange data
    ggdat = melt(ggdat, measure.vars = row.names, variable.name = "group.level", value.name = "group.value")
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
        # geom_boxplot() +
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
    out$tab.count = ggdat
    out$plot.count = invisible(gg.count)
    out$plot.perc = invisible(gg.perc)
  }
  
  
  ## 2. Create PLOTS for Diff.By.Pixel ------------------------------------------------------------
  if (do.maps || do.mean)
  {
    ggdat = bm.range$Diff.By.Pixel
    
    ## c. SRC maps per model --------------------------------------------------
    if (do.maps) {
      
      # gg.maps = 'terra::plot(ggdat, col = c("-2" = "#fc8d62", "-1" = "grey", "0" = "white", "1" = "#66c2a5")
      #      , legend.width = 2, legend.shrink = 0.7
      #      , axis.args = list(at = c(-2, -1, 0, 1), labels = c("Loss", "Stable1", "Stable0", "Gain"), cex.axis = 1))'
      gg.maps <- 'plot(ggdat,
           col = data.frame(
             value = c(-2, -1, 0, 1),
             color = c("#fc8d62", "lightgoldenrod2", "grey", "#66c2a5")),
           colNA = "white")'
      # ggdat = as.data.frame(rasterToPoints(ggdat))
      # 
      # ## Get models information
      # ggdat = melt(ggdat, id.vars = c("x", "y"), variable.name = "full.name", value.name = "SRC")
      # for (ii in 1:length(row.names)) { ggdat[[row.names[ii]]] = NA }
      # for (jj in 1:nrow(corres)) {
      #   ind = which(ggdat$full.name == corres$full.name[jj])
      #   for (ii in 1:length(row.names)) {
      #     ggdat[ind, row.names[ii]] = corres[jj, row.names[ii]]
      #   }
      # }
      # 
      # ## Rearrange data
      # ggdat = melt(ggdat, measure.vars = row.names, variable.name = "group.level", value.name = "group.value")
      # 
      # ## Do plot
      # gg.maps = ggplot(ggdat, aes_string(x = "x", y = "y", fill = "as.factor(SRC)")) +
      #   geom_raster() +
      #   facet_wrap("full.name") +
      #   scale_fill_manual("", values = c("-2" = "#fc8d62"
      #                                    , "-1" = "grey"
      #                                    , "0" = "white"
      #                                    , "1" = "#66c2a5")
      #                     , labels = c("-2" = "Loss"
      #                                  , "-1" = "Stable1"
      #                                  , "0" = ""
      #                                  , "1" = "Gain")) +
      #   xlab("") +
      #   ylab("") +
      #   theme(legend.title = element_blank()
      #         , legend.key = element_rect(fill = "white")
      #         , legend.position = "top")
      
      out$tab.maps = ggdat
      out$plot.maps = invisible(gg.maps)
      
    } else { gg.maps = NULL }
    
    ## d. SRC mean maps per group.level ---------------------------------------
    if (do.mean) {
      if(!requireNamespace('ggpubr', quietly = TRUE)) stop("Package 'ggpubr' not found")
      
      corres = data.frame(full.name = names(ggdat))
      for (ii in 1:length(row.names)) {
        corres[[row.names[ii]]] = sapply(corres$full.name, function(x) strsplit(x, "_")[[1]][ii])
      }
      reclass_table = data.frame(is = c(1, 2, 3), becomes = c(1, -1, -2))
      fun_mode = function(x) {
        tmp = table(x)
        return(names(tmp)[which.max(tmp)])
      }
      
      list.cons = list.perc = list()
      for (ii in row.names) {
        for (jj in unique(corres[, ii])) {
          ras = ggdat[[corres$full.name[which(corres[, ii] == jj)]]]
          if (nlyr(ras) > 1) {
            stk = foreach (vali = c(1, -1, -2), .combine = "c") %do% {
              res = ras
              res = classify(res, rcl = matrix(c(vali,1), ncol = 2), others = 0)
              res = sum(res, na.rm = TRUE)
              names(res) = paste0("VAL_", vali)
              res = classify(res, rcl = matrix(c(0,NA), ncol = 2))
              return(res)
            }
            ras1 = which.max(stk)
            ras1 = classify(ras1, reclass_table)
            ras2 = max(stk, na.rm = TRUE) / sum(stk, na.rm = TRUE)
            list.cons[[paste0(ii, "_", jj)]] = ras1
            list.perc[[paste0(ii, "_", jj)]] = ras2
          }
        }
      }
      if (length(list.cons) > 0 && length(list.perc) > 0) {
        stk.cons = rast(list.cons)
        stk.perc = rast(list.perc)
        tab1 = as.data.frame(stk.cons, xy = TRUE)
        tab1 = melt(tab1, id.vars = c("x", "y"))
        tab1$group.level = tab1$group.value = ""
        for (ii in row.names) {
          tab1$group.level[grep(ii, tab1$variable)] = ii
        }
        for (jj in unique(unlist(corres[, 2:ncol(corres)]))) { 
          tab1$group.value[grep(jj, tab1$variable)] = jj 
        }
        tab1$value[which(is.na(tab1$value))] = 0
        
        tab2 = as.data.frame(stk.perc, xy = TRUE)
        tab2 = melt(tab2, id.vars = c("x", "y"))
        tab2$group.level = tab2$group.value = ""
        for (ii in row.names) { tab2$group.level[grep(ii, tab2$variable)] = ii }
        for (jj in unique(unlist(corres[, 2:ncol(corres)]))) { tab2$group.value[grep(jj, tab2$variable)] = jj }
        tab2$value[which(tab2$value == 1 & tab1$value == 0)] = NA
        
        gg.ca1 = ggplot(tab1, aes_string(x = "x", y = "y", fill = "as.factor(value)")) +
          geom_tile() +
          facet_wrap("group.level ~ group.value") +
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
        
        gg.ca2 = ggplot(tab2, aes_string(x = "x", y = "y", fill = "value")) +
          geom_tile() +
          facet_wrap("group.level ~ group.value") +
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
        
        gg.ca = ggpubr::ggarrange(gg.ca1, gg.ca2, ncol = 2)
        out$tab.ca1 = tab1
        out$tab.ca2 = tab2
        out$plot.ca = invisible(gg.ca)
      } else {
        gg.ca = NULL
        warning("'do.mean' is only available if several maps are provided")
      }
    } else { gg.ca = NULL }
  } else {
    gg.maps = NULL
    gg.ca = NULL
  }
  ## RETURN PLOTS
  if (do.plot) { 
    print(gg.count)
    plot.new()
    print(gg.perc, newpage = FALSE)
    eval(parse(text = gg.maps))
    print(gg.ca)
  }
  return(out)
}


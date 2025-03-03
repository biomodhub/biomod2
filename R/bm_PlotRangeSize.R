###################################################################################################
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
##' @param bm.range an \code{BIOMOD.rangesize.out} object returned by the \code{\link{BIOMOD_RangeSize}} function
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
##' @seealso \code{\link{BIOMOD_RangeSize}} \code{\link{bm_RangeSize}}
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
##'                                       metric.eval = c('TSS', 'ROC'),
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
##' # Compute differences
##' myBiomodRangeSize <- BIOMOD_RangeSize(proj.current = myBiomodProj,
##'                                       proj.future = myBiomodProjectionFuture,
##'                                       metric.binary = "TSS")
##' 
##' 
##' # Represent main results
##' bm_PlotRangeSize(bm.range = myBiomodRangeSize)
##' 
##' 
##' 
##' @importFrom graphics plot.new
##' @importFrom reshape2 melt
##' @importFrom foreach foreach %do%
##' @importFrom terra rast which.max nlyr  classify plot
##' @importFrom ggplot2 ggplot geom_col geom_tile geom_label facet_wrap xlab ylab labs scale_fill_viridis_c
##' theme theme_bw element_blank element_rect scale_fill_manual scale_x_discrete guide_legend scale_color_gradientn scale_fill_gradientn
##' @importFrom rlang .data
##' 
##' @export
##' 
##' 
###################################################################################################


bm_PlotRangeSize <- function(bm.range, do.count = TRUE, do.perc = TRUE
                             , do.maps = TRUE, do.mean = TRUE, do.plot = TRUE)
{
  
  args <- .bm_PlotRangeSize.check.args(bm.range, do.count, do.perc, do.maps, do.mean, do.plot)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  out = list()
  
  row.names <- bm.range@row.names
  
  ## 1. Create PLOTS for Compt.By.Models ----------------------------------------------------------
  if (do.count || do.perc)
  {
    ggdat = as.data.frame(bm.range@Compt.By.Models)
    
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
    ggdat = melt(ggdat, measure.vars = c("Loss", "Stable_Abs", "Stable_Pres", "Gain"), variable.name = "count.level", value.name = "count.value")
    ggdat$PercLoss = ggdat$PercLoss * (-1)
    ggdat = melt(ggdat, measure.vars = c("PercLoss", "PercGain", "SpeciesRangeChange"), variable.name = "perc.level", value.name = "perc.value")
    ggdat = melt(ggdat, measure.vars = c("CurrentRangeSize", "FutureRangeSize.NoDisp", "FutureRangeSize.FullDisp")
                 , variable.name = "range.level", value.name = "range.value")
    
    
    
    ## a. Count plot ----------------------------------------------------------
    if (do.count) {
      gg.count = ggplot(ggdat[which(ggdat$count.level != "Stable0"), ]
                        , aes(x = group.value, y = count.value, fill = count.level)) +
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
      gg.perc = ggplot(ggdat, aes(x = group.value, y = perc.value, fill = perc.level)) +
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
    ggdat = bm.range@Diff.By.Pixel
    
    if(type.df){
      ggdf <- cbind(bm.range@coord, bm.range@Diff.By.Pixel)
    }
    
    ## c. SRC maps per model --------------------------------------------------
    if (do.maps) {
      
      if (bm.range@data.type == "binary") {
        if (type.df){
          ggdf <- melt(ggdf, c("x", "y"), variable.name = "models", value.name = "change")
          
          gg.maps <- 'ggplot(ggdf) + 
            facet_wrap(vars(models))+
            geom_point(aes(x = x, y = y, color = as.factor(change))) +
            scale_color_manual(values = c("-2" = "#fc8d62","-1" = "lightgoldenrod2", "0" = "grey", "1" = "#66c2a5"),
                               labels = c("Loss", "Stable_presence", "Stable_Absence", "Gain"),
                               na.translate = F,
                               guide = guide_legend(title = element_blank()))'
        } else {
          gg.maps <- 'plot(ggdat,
                        col = data.frame(
                          value = c(-2, -1, 0, 1),
                          color = c("#fc8d62", "lightgoldenrod2", "grey", "#66c2a5")),
                        colNA = "white")'
        }
        
        
      } else {
        loss.gain <- bm.range@loss.gain
        
        if(type.df){
          min_change <- min(ggdf[,3:ncol(ggdf)], na.rm = T)
          max_change <- max(ggdf[,3:ncol(ggdf)], na.rm = T)
          lim <- max(abs(min_change), abs(max_change))
          if (bm.range@data.type != "ordinal") {lim <- 1}
          
          stable <- loss.gain == "0"
          ggdf[,3:ncol(ggdf)][stable] <- -lim -0.1
          
          ggdf <- melt(ggdf, c("x", "y"), variable.name = "models", value.name = "change")
          
          gg.maps <- 'ggplot(ggdf) + 
  facet_wrap(vars(models))+
  geom_point(aes(x = x, y = y, color = change)) +
  scale_color_gradientn(colours = c("grey", "#d13d04","#fed3c2", "white","#acdece", "#337f67"),
                       values = c(0,0.01,0.45,0.5,0.55,1), limits = c(-lim - 0.1, lim),
                       na.value = "transparent", 
                       breaks = c(-lim/2,0,lim/2),
                       labels = c("Loss", "Stable_presence", "Gain")) + 
  labs("Change intensity")'
          
          
        } else {
          min_change <- min(minmax(ggdat), na.rm = T)
          max_change <- max(minmax(ggdat), na.rm = T)
          lim <- max(abs(min_change), abs(max_change))
          if (bm.range@data.type != "ordinal") {lim <- 1}
          
          for(i in 1:nlyr(ggdat)){
            stable <- loss.gain[[i]][] == "0"
            ggdat[[i]][stable] <- -lim -0.1
          }
          
          gg.maps <- 'ggplot() +
  facet_wrap(~lyr)+
  tidyterra::geom_spatraster(data = ggdat) +
  scale_fill_gradientn(colours = c("grey", "#d13d04","#fed3c2", "white","#acdece", "#337f67"),
                       values = c(0,0.01,0.45,0.5,0.55,1), limits = c(-lim - 0.1, lim),
                       na.value = "transparent", 
                       breaks = c(-lim/2,0,lim/2),
                       labels = c("Loss", "Stable_presence", "Gain")) +
  labs(fill = "Change intensity")'
        }
        
        
      }
      
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
          
          gg.ca1 = ggplot(tab1, aes(x = .data$x, y = .data$y, fill = as.factor(.data$value))) +
            geom_tile() +
            facet_wrap("group.level ~ group.value") +
            scale_fill_manual("", values = c("-2" = "#fc8d62"
                                             , "-1" = "grey"
                                             , "0" = "white"
                                             , "1" = "#66c2a5")
                              , labels = c("-2" = "Loss"
                                           , "-1" = "Stable_Pres"
                                           , "0" = ""
                                           , "1" = "Gain")) +
            xlab("") +
            ylab("") +
            labs(title = "Community averaging value across models") +
            theme(legend.title = element_blank()
                  , legend.key = element_rect(fill = "white")
                  , legend.position = "top")
          
          gg.ca2 = ggplot(tab2, aes(x = .data$x, y = .data$y, fill = .data$value)) +
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
      print(eval(parse(text = gg.maps)))
      print(gg.ca)
  }
  return(gg.maps)
}


###################################################################################################

.bm_PlotRangeSize.check.args <- function(bm.range,
                                         do.count, do.perc, do.maps, do.mean,
                                         do.plot){
  # if (!is.list(bm.range) ||
  #     (is.list(bm.range) && length(bm.range) != 2) ||
  #     (is.list(bm.range) && length(bm.range) == 2 && !all(c("Compt.By.Models", "Diff.By.Pixel") %in% names(bm.range)))) {
  #   stop("'bm.range' must be an object obtained by the BIOMOD_RangeSize function")
  # }
  
  .fun_testIfInherits(TRUE, "bm.range", bm.range, "BIOMOD.rangesize.out")
  
  stopifnot(is.logical(do.count)) ## Useful ? 
  stopifnot(is.logical(do.perc))
  stopifnot(is.logical(do.maps))
  stopifnot(is.logical(do.mean))
  stopifnot(is.logical(do.plot))
  
  type.df <- is.data.frame(bm.range@Diff.By.Pixel)
  
  if((bm.range@data.type != "binary" && do.mean == TRUE) | (bm.range@data.type == "binary" && type.df)){
    do.mean <- FALSE
    warning("'do.mean' is only available for binary data as SpatRaster.")
  } 
  
  
  return(list(do.mean = do.mean,
              type.df = type.df))
}
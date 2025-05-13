

# 1.3 Other Functions -----------------------------------------------------------------------------

### plot.BIOMOD.formated.data.abundance (doc) --------------------------------------------------
##' 
##' @rdname plot
##' @docType methods
##' @author Hélène Blancheteau
##' @title \code{plot} method for \code{BIOMOD.formated.data.abundance} object class
##' 
##' @description Plot the spatial distribution of presences, absences and 
##' pseudo-absences among the different potential dataset (calibration, 
##' validation and evaluation). Available only if coordinates were given to 
##' \code{\link{BIOMOD_FormatingData}}.
##' 
##' 
##' @param x a \code{BIOMOD.formated.data.abundance} 
##' object. Coordinates must be available to be able to use \code{plot}.
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' an \code{data.frame} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{bm_CrossValidation}} functions, to explore the distribution of calibration 
##' and validation datasets
##' @param plot.type a \code{character}, either \code{'points'} (\emph{default}) 
##' or \code{'raster'} (\emph{if environmental variables were given as a raster}). 
##' With \code{plot.type = 'points'} occurrences will be represented as points
##' (better when using fine-grained data). With \code{plot.type = 'raster'}
##' occurrences will be represented as a raster (better when using coarse-grained
##' data)
##' @param plot.output a \code{character}, either \code{'facet'} (\emph{default}) 
##' or \code{'list'}. \code{plot.output} determines whether plots are returned
##' as a single facet with all plots or a \code{list} of individual plots
##' (better when there are numerous graphics)
##' @param PA (\emph{optional, default} \code{'all'}) \cr 
##' If \code{x} is a \code{BIOMOD.formated.data.PA} object, a \code{vector} 
##' containing pseudo-absence set to be represented 
##' @param run (\emph{optional, default} \code{'all'}) \cr 
##' If \code{calib.lines} provided, a \code{vector} containing repetition set to 
##' be represented 
##' @param plot.eval (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether evaluation data should be added to the plot or not
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether the plot is to be rendered or not
##' @param point.size a \code{numeric} to adjust the size of points when
##'  \code{plot.type = 'points'}
##' @param has.mask a \code{logical}. Is a proper mask available ?
##' @param has.mask.eval a \code{logical}. Is a proper mask available for evaluation data ?  
##' 
##' @return a \code{list} with the data used to generate the plot and a
##' \code{ggplot2} object 
##' 
##' @importFrom terra rast minmax crds ext
##' @importFrom ggplot2 ggplot aes scale_color_manual scale_shape_manual scale_fill_manual guides xlim ylim ggtitle facet_wrap theme guide_legend after_stat scale_size scale_alpha_continuous scale_alpha waiver
##' 
##' @export
##' 
##' 
##' @examples
##' 
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
##' ## ----------------------------------------------------------------------- #
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' myBiomodData
##' plot(myBiomodData)
##' 
##' 


.plot.BIOMOD.formated.data.abundance <- function(x,
                                                 calib.lines = NULL,
                                                 plot.type,
                                                 plot.output, 
                                                 run,
                                                 plot.eval,
                                                 point.size = 1.5,
                                                 do.plot = TRUE,
                                                 has.mask = FALSE,
                                                 has.mask.eval = FALSE)
{
  # 1 - extract SpatVector for all required data ----------------------
  
  ## 1.1 - Full Dataset -----------------------------------------------
  filter <- !is.na(x@data.species)
  full.resp <- x@data.species[filter]
  full.xy <- x@coord[filter,]
  full.df <- data.frame(resp = as.numeric(full.resp),
                        x = full.xy[,1],
                        y = full.xy[,2])
  
  full.df.vect <- vect(full.df, geom = c("x","y"))
  names(full.df.vect) <- "resp"
  full.df.vect$dataset <- "Initial dataset"
  full.df.vect$part <- "Initial"
  
  
  ## 1.2 - Eval Dataset -----------------------------------------------
  
  if (plot.eval) {
    eval.resp <- x@eval.data.species
    eval.xy <- x@eval.coord
    eval.df <- data.frame(resp = as.numeric(eval.resp),
                          x = eval.xy[,1],
                          y = eval.xy[,2])
    
    eval.df.vect <- vect(eval.df, geom = c("x","y"))
    names(eval.df.vect) <- "resp"
    eval.df.vect$dataset <- "Initial dataset"
    eval.df.vect$part <- "Evaluation"
    full.df.vect <- rbind(full.df.vect, eval.df.vect)
  }
  
  ## 1.3 - CV dataset -----------------------------
  this_PA <- 'allData'
  if (!is.null(calib.lines)) {
    PA_run.vect <- foreach(this_run = run, .combine = 'rbind') %do%
      {
        if (is.na(this_PA) || this_PA == 'allData') { # run only
          this_name <- paste0("_", this_PA, "_", this_run)
          this_calib <- calib.lines[ , this_name]
          this_valid <- ! calib.lines[ , this_name] 
        }
        else if (is.na(this_run) || this_run == 'allRun') { # PA only
          this_name <- this_PA
          this_calib <- x@PA.table[ , this_PA]
        } else { # PA+run
          this_name <- paste0("_", this_PA, "_", this_run)
          this_calib <- calib.lines[ , this_name] & x@PA.table[ , this_PA]
          this_valid <- ! calib.lines[ , this_name] & x@PA.table[ , this_PA]
        }
        calib.resp <- x@data.species[which(this_calib)]
        # calib.resp <- ifelse(is.na(calib.resp), 30, 
        #                      ifelse(calib.resp == 1, 10, 20))
        calib.xy <- x@coord[which(this_calib),]
        calib.df <- data.frame(resp = calib.resp,
                               x = calib.xy[, 1],
                               y = calib.xy[, 2],
                               part = "Calibration")
        
        if (!is.na(this_run) & this_run != "allRun") { 
          valid.resp <- x@data.species[which(this_valid)]
          # valid.resp <- ifelse(is.na(valid.resp), 31, 
          #                      ifelse(valid.resp == 1, 11, 21))
          valid.xy <- x@coord[which(this_valid),]
          valid.df <- data.frame(resp = valid.resp,
                                 x = valid.xy[, 1],
                                 y = valid.xy[, 2],
                                 part = "Validation")
          calib.df <- rbind(calib.df, valid.df)
        }
        thisdf.vect <- vect(calib.df, geom = c("x","y"))
        names(thisdf.vect) <- c("resp","part")
        thisdf.vect$dataset <- this_name
        thisdf.vect
      }
    full.df.vect <- rbind(full.df.vect, PA_run.vect)
  }
  
  
  # 2- define colors and breaks ------------------------------------
  data_breaks <- c("Initial", "Evaluation", "Calibration", "Validation","Background")              
  data_labels <- data_breaks
  data_labels_facet <- c("Initial", "Evaluation", "Calibration", "Validation", NA) # background
  
  data_colors <- c("Initial" = "#004488",
                   "Evaluation" = "#994455",
                   "Calibration" = "#997700",
                   "Validation" = "#EECC66",
                   "1" = "#D4D4D4")
  
  shape_fit <- 16
  shape_eval <- 17
  data_shape <- c(shape_fit,
                  shape_eval,
                  rep(shape_fit,length(colnames(calib.lines)) ))
  data_alpha <- c()
  data_background <- "#FFFFFF00"
  
  
  
  # 3 - prepare plots -------------------------------------------------------
  this_mask_eval <- rast()
  if(has.mask){
    this_mask <- rast(x@data.mask[["calibration"]])
    this_mask_eval <- this_mask
  } else {
    this_mask <- rast()
  }
  if(has.mask.eval){
    this_mask_eval <- rast(x@data.mask[["evaluation"]])
  }
  if(has.mask | has.mask.eval){
    plot_mask <- foreach(this_dataset = unique(full.df.vect$dataset), 
                         .combine = 'c') %do% {
                           if(this_dataset == "Evaluation dataset"){
                             return(this_mask_eval)
                           } else {
                             return(this_mask)
                           }
                         }
    names(plot_mask) <- unique(full.df.vect$dataset)
  }

  ## 3.1 Raster plot --------------------------------------------------------
  if(plot.type == "raster"){

    rast.plot <- foreach(this_dataset = unique(full.df.vect$dataset), .combine = 'c') %do% {
      this_rast  <-
        rasterize(subset(full.df.vect,
                         full.df.vect$dataset == this_dataset), 
                  plot_mask[[this_dataset]],
                  field = "resp", by = "part", fun = mean, background = 0)
      if (this_dataset == "Initial dataset"){
        names(this_rast) <- this_dataset
      } else {
        names(this_rast) <- paste(this_dataset, c("calibration","validation"), sep = "_")
      }
      this_rast*this_mask
    }
    
    data_colors <- c("#004488")
    if(plot.eval){
      data_colors <- c(data.color,"#994455")
    }
    if(!is.null(calib.lines)){
      nb_run <- ncol(calib.lines)
      data_colors <- c(data_colors, rep(c( "#997700","#EECC66"), nb_run))
    }
    names(data_colors) <- names(rast.plot)

    if(plot.output == "facet"){
      g <- ggplot()+
        tidyterra::geom_spatraster(data = rast.plot, aes(alpha = after_stat(value), fill = lyr))+
        facet_wrap(~lyr)+
        scale_fill_manual(
          NULL,
          breaks = names(rast.plot),
          values = data_colors,
          labels = names(rast.plot),
          drop = FALSE)+
        scale_alpha_continuous(
          NULL,
          range = c(0.1,1),
          na.value = 0)+
        guides(fill = guide_legend(nrow = 2))+
        theme(legend.position = "top",
              legend.key = element_blank(),
              legend.background = element_rect(fill = "grey90"),
              legend.text = ggtext::element_markdown(),
              legend.box = "vertical")
      
    } else {
      g <- lapply(names(rast.plot), function(thisname){
        ggplot()+
          tidyterra::geom_spatraster(data = rast.plot[[thisname]], aes(alpha = after_stat(value)), fill = data_colors[thisname])+
          scale_alpha_continuous(
            NULL,
            range = c(0.1,1),
            na.value = 0)+
          ggtitle(thisname)+
          guides(fill = guide_legend(nrow = 2))+
          theme(legend.position = "top",
                legend.key = element_blank(),
                legend.background = element_rect(fill = "grey90"),
                legend.text = ggtext::element_markdown(),
                legend.box = "vertical")
      })
    }
    if(do.plot){
      print(g)
    }
    return(list("data.vect"  = full.df.vect,
                "data.rast"  = rast.plot,
                "data.label" = data_labels,
                "data.plot"  = g))
  } else {
    
    ## 3.2 Points plot --------------------------------------------------------
    
    data.df <- as.data.frame(full.df.vect, geom = "XY")
    datasets <- unique(data.df$dataset)
    datasets <- sort(datasets)
    datasets <- c("Initial dataset", datasets[-length(datasets)])
    data.df$dataset <- factor(data.df$dataset, datasets)
    
    if (x@data.type %in% c("ordinal", "multiclass")){
      labels_factor <- levels(x@data.species)
      breaks_factor <- 1:length(labels_factor)
    } else {
      labels_factor <- waiver()
      breaks_factor <- waiver()
    }
    
    if(plot.output == "facet"){
      base_g <-  ggplot(data.df)
      if(has.mask){
        base_g <- base_g +
          tidyterra::geom_spatraster(data = this_mask, aes(fill = factor(after_stat(value))))
      }
      
      g <- base_g +      
        geom_point(aes(x = x, y = y, 
                       alpha= resp,
                       color = part,
                       size = resp),
                   shape = 18)+ #size = point.size
        facet_wrap(~dataset)+
        scale_size(
          range =c(0.5,3), 
          labels = labels_factor,
          breaks = breaks_factor
        )+
        scale_color_manual(
          NULL,
          breaks = data_breaks,
          values = data_colors,
          labels = data_labels_facet,
          drop = FALSE)+
        scale_fill_manual(
          guide = "none",
          breaks = data_breaks,
          values = data_colors,
          labels = data_labels,
          na.value = data_background)+
        scale_alpha(
          labels = labels_factor,
          breaks = breaks_factor
        )+
        xlab(NULL)+ ylab(NULL)+
        guides(color = guide_legend(override.aes = list(size = 3)))+
        theme(legend.position = "top",
              legend.key = element_blank(),
              legend.background = element_rect(fill = "grey90"),
              legend.text = ggtext::element_markdown(), 
              legend.box = "vertical")
      
    } else {
      g <- lapply(unique(data.df$dataset), function(thisname){
        base_g <-  ggplot(subset(data.df,
                                 data.df$dataset == thisname))
        if(has.mask){
          base_g <- base_g + 
            tidyterra::geom_spatraster(data = this_mask,
                                       aes(fill = factor(after_stat(value))))
        }
        base_g +      
          geom_point(aes(x = x, y = y, 
                         alpha= resp,
                         color = part,
                         size = resp),
                     shape = 18)+
          scale_size(range =c(0.5,3))+
          scale_color_manual(
            NULL,
            breaks = data_breaks,
            values = data_colors,
            labels = data_labels)+
          scale_fill_manual(
            guide = "none",
            breaks = data_breaks,
            values = data_colors,
            labels = data_labels,
            na.value = data_background)+
          xlab(NULL)+ ylab(NULL)+
          guides(color = guide_legend(override.aes = list(size = 3)))+
          theme(legend.position = "top",
                legend.key = element_blank(),
                legend.background = element_rect(fill = "grey90"),
                legend.box = "vertical")+
          ggtitle(thisname)
      })
      
    }
  }
  if(do.plot){
    print(g)
  }
  return(list("data.vect"  = full.df.vect,
              "data.label" = data_labels,
              "data.plot"  = g))
  
}






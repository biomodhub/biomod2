
## --------------------------------------------------------------------------- #
## 1. BIOMOD.formated.data ---------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.formated.data
##' @aliases BIOMOD.formated.data-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_FormatingData()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
##' \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_CrossValidation}} and 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @param dir.name a \code{character} corresponding to the modeling folder
##' @param sp.name a \code{character} corresponding to the species name
##' 
##' @param sp A \code{vector}, a \code{\link[terra:vect]{SpatVector}} without associated 
##' data (\emph{if presence-only}), or a \code{\link[terra:vect]{SpatVector}}
##' object containing binary data (\code{0} : absence,  \code{1} : presence,
##' \code{NA} : indeterminate) for a single species that will be used to 
##' build the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##'  \code{SpatialPoints}  (if presence-only) or \code{SpatialPointsDataFrame}
##'  object containing binary data.}
##' @param env a \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s).
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' 
##' @param xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to build the 
##' species distribution model(s)
##' @param eval.sp (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector}, a \code{\link[terra:vect]{SpatVector}} without associated 
##' data (\emph{if presence-only}), or a \code{\link[terra:vect]{SpatVector}}
##' object containing binary data (\code{0} : absence, \code{1} : presence, 
##' \code{NA} : indeterminate) for a single species that will be used to
##' evaluate the species distribution model(s) with independent data
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##' \code{SpatialPoints}  (if presence-only) or \code{SpatialPointsDataFrame}
##' object containing binary data.}
##' @param eval.env (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}} or
##' \code{\link[terra:rast]{SpatRaster}} object containing the explanatory
##' variables (in columns or layers) that will be used to evaluate the species
##' distribution model(s) with independent data
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still
##' supported such as \code{RasterStack} and \code{SpatialPointsDataFrame}
##' objects. }
##' @param eval.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or
##' \code{data.frame} containing the corresponding \code{X} and \code{Y}
##' coordinates that will be used to evaluate the species distribution model(s)
##' with independent data
##' 
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether points having one or several missing
##' values for explanatory variables should be removed from the analysis or not
##' @param data.mask (\emph{optional, default} \code{NULL}) \cr 
##' A \code{\link[terra:rast]{SpatRaster}} object containing the mask of the studied area
##' @param shared.eval.env (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether the explanatory variables used for the 
##' evaluation dataset are the same than the ones for calibration (if \code{eval.env} not 
##' provided for example) or not
##' @param filter.raster (\emph{optional, default} \code{FALSE}) \cr 
##' If \code{env} is of raster type, a \code{logical} value defining whether \code{sp} 
##' is to be filtered when several points occur in the same raster cell
##' 
##' @param object a \code{\link{BIOMOD.formated.data}} object
##' 
##' 
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates
##' @slot data.species a \code{vector} containing the species observations (\code{0}, \code{1} or 
##' \code{NA})
##' @slot data.env.var a \code{data.frame} containing explanatory variables
##' @slot data.mask a \code{\link[terra:rast]{SpatRaster}} object containing the mask of the 
##' studied area
##' @slot has.data.eval a \code{logical} value defining whether evaluation data is given
##' @slot eval.coord (\emph{optional, default} \code{NULL}) \cr 
##' A 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates for evaluation data
##' @slot eval.data.species (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the species observations (\code{0}, \code{1} or \code{NA}) for 
##' evaluation data
##' @slot eval.data.env.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{data.frame} containing explanatory variables for evaluation data
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{BIOMOD_Tuning}}, 
##' \code{\link{BIOMOD_CrossValidation}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{bm_RunModelsLoop}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.formated.data")
##' 
##' ## ----------------------------------------------------------------------- #
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
##' summary(myBiomodData)
##' 
##' 
NULL

##' @name BIOMOD.formated.data-class
##' @rdname BIOMOD.formated.data
##' @importFrom terra rast app is.factor subset extract cellFromXY `add<-` 
##' classify rasterize values
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data",
         representation(dir.name = 'character',
                        sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        data.env.var = "data.frame",
                        data.mask = "list",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object){ 
           check.data.mask <- suppressWarnings(
             all(sapply(object@data.mask, function(x) inherits(x, "PackedSpatRaster")))
           )
           if(check.data.mask){
             return(TRUE)
           } else {
             return(FALSE)
           }
         })


# 1.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.formated.data", def = function(sp, env, ...) { standardGeneric("BIOMOD.formated.data") })

## BIOMOD.formated.data(sp = numeric, env = data.frame) -----------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE, data.mask = NULL, shared.eval.env = FALSE
                   , filter.raster = FALSE) 
          {
            
            if (is.null(data.mask)) { 
              suppressWarnings(data.mask <- list("calibration" = wrap(rast())))
            }
            
            if (is.null(eval.sp)) { ## NO EVALUATION DATA
              BFD <- new(
                'BIOMOD.formated.data',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                dir.name = dir.name,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = FALSE
              )
            } else { ## EVALUATION DATA
              BFDeval <- BIOMOD.formated.data(
                sp = eval.sp,
                env = eval.env,
                xy = eval.xy,
                dir.name = dir.name,
                sp.name = sp.name,
                filter.raster = filter.raster
              )
              
              if (rast.has.values(rast(BFDeval@data.mask[[1]]))) {
                data.mask[["evaluation"]] <- BFDeval@data.mask[[1]]
              } else if (shared.eval.env) {
                data.mask[["evaluation"]] <- data.mask[["calibration"]]
              }
              
              BFD <- new(
                'BIOMOD.formated.data',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                dir.name = dir.name,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = TRUE,
                eval.coord = BFDeval@coord,
                eval.data.species = BFDeval@data.species,
                eval.data.env.var = BFDeval@data.env.var
              )
              
              rm('BFDeval')
            }
            
            ## REMOVE NA IF ANY
            if (na.rm) {
              rowToRm <- unique(unlist(lapply(BFD@data.env.var, function(x) { return(which(is.na(x))) })))
              if (length(rowToRm) > 0) {
                cat("\n ! Some NAs have been automatically removed from your data")
                BFD@coord <- BFD@coord[-rowToRm, , drop = FALSE]
                BFD@data.species <- BFD@data.species[-rowToRm]
                BFD@data.env.var <- BFD@data.env.var[-rowToRm, , drop = FALSE]
              }
              if (BFD@has.data.eval) {
                rowToRm <- unique(unlist(lapply(BFD@eval.data.env.var, function(x) { return(which(is.na(x))) })))
                if (length(rowToRm) > 0) {
                  cat("\n ! Some NAs have been automatically removed from your evaluation data")
                  BFD@eval.coord <- BFD@eval.coord[-rowToRm, , drop = FALSE]
                  BFD@eval.data.species <- BFD@eval.data.species[-rowToRm]
                  BFD@eval.data.env.var <- BFD@eval.data.env.var[-rowToRm, , drop = FALSE]
                }
              }
            }
            return(BFD)
          }
)

## BIOMOD.formated.data(sp = data.frame) -----------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL,
                   eval.sp = NULL, eval.env = NULL, eval.xy = NULL,
                   na.rm = TRUE, filter.raster = FALSE)
          {
            if (ncol(sp) > 1) { stop("Invalid response variable") }
            sp <- as.numeric(unlist(sp))
            BFD <- BIOMOD.formated.data(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

## BIOMOD.formated.data(sp = numeric, env = matrix) ----------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'matrix'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL,
                   eval.sp = NULL, eval.env = NULL, eval.xy = NULL,
                   na.rm = TRUE, filter.raster = FALSE)
          {
            env <- as.data.frame(env)
            BFD <- BIOMOD.formated.data(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

## BIOMOD.formated.data(sp = numeric, env = SpatRaster) -----------------------
##' 
##' @rdname BIOMOD.formated.data
##' @export
##' 

setMethod('BIOMOD.formated.data', signature(sp = 'numeric', env = 'SpatRaster'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL,
                   eval.sp = NULL, eval.env = NULL, eval.xy = NULL,
                   na.rm = TRUE, shared.eval.env = FALSE,
                   filter.raster = FALSE)
          {
            
            ## Keep same env variable for eval than calib (+ check for factor)
            if (!is.null(eval.sp) && is.null(eval.env)) {
              output <- check_duplicated_cells(env, eval.xy, eval.sp, filter.raster)
              eval.xy <- output$xy
              eval.sp <- output$sp
              rm(output)
              
              eval.env <- as.data.frame(extract(env, eval.xy, ID = FALSE))
              shared.eval.env <- TRUE
            }
            
            ## Prepare mask of studied area
            data.mask <- prod(classify(env, matrix(c(-Inf,Inf,1), nrow = 1)))
            names(data.mask) <- "Environmental Mask"
            data.mask <- list("calibration" = wrap(data.mask))
            ## Keep same env variable for eval than calib (+ check for factor)
            
            output <- check_duplicated_cells(env, xy, sp, filter.raster)
            xy <- output$xy
            sp <- output$sp
            rm(output)
            
            env <- as.data.frame(extract(env, xy, factors = TRUE, ID = FALSE))
            
            BFD <- BIOMOD.formated.data(sp, env, xy, dir.name, sp.name, 
                                        eval.sp, eval.env, eval.xy,
                                        na.rm = na.rm, data.mask = data.mask,
                                        shared.eval.env = shared.eval.env,
                                        filter.raster = filter.raster)
            return(BFD)
          }
)


# 1.3 Other Functions -----------------------------------------------------------------------------

### plot.BIOMOD.formated.data (doc) --------------------------------------------------
##' 
##' @rdname plot
##' @docType methods
##' @author Remi Patin
##' @title \code{plot} method for \code{\link{BIOMOD.formated.data}} object class
##' 
##' @description Plot the spatial distribution of presences, absences and 
##' pseudo-absences among the different potential dataset (calibration, 
##' validation and evaluation). Available only if coordinates were given to 
##' \code{\link{BIOMOD_FormatingData}}.
##' 
##' 
##' @param x a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}}
##' object. Coordinates must be available to be able to use \code{plot}.
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' an \code{array} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{BIOMOD_CrossValidation}} functions, to explore the distribution of calibration 
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
##' If \code{x} is a \code{\link{BIOMOD.formated.data.PA}} object, a \code{vector} 
##' containing pseudo-absence set to be represented 
##' @param run (\emph{optional, default} \code{'all'}) \cr 
##' If \code{calib.lines} provided, a \code{vector} containing repetition set to 
##' be represented 
##' @param plot.eval (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether evaluation data should be added to the plot or not
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether the plot is to be rendered or not
##' 
##' 
##' @return a \code{list} with the data used to generate the plot and a
##' \code{ggplot2} object 
##' 
##' @importFrom terra rast minmax crds ext
##' @importFrom ggplot2 ggplot aes scale_color_manual scale_shape_manual scale_fill_manual guides xlim ylim ggtitle facet_wrap theme guide_legend after_stat
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


setMethod('plot', signature(x = 'BIOMOD.formated.data', y = "missing"),
          function(x,
                   calib.lines = NULL,
                   plot.type,
                   plot.output, 
                   PA,
                   run,
                   plot.eval,
                   do.plot = TRUE){
            args <- .plot.BIOMOD.formated.data.check.args(x = x,
                                                          calib.lines = calib.lines,
                                                          plot.type = plot.type,
                                                          plot.output = plot.output, 
                                                          PA = PA,
                                                          run = run,
                                                          plot.eval = plot.eval,
                                                          do.plot = do.plot)
            
            for (argi in names(args)) { 
              assign(x = argi, value = args[[argi]]) 
            }
            rm(args)
            
            
            # 1 - extract SpatVector for all required data ----------------------
            
            ## 1.1 - Full Dataset -----------------------------------------------
            filter_PA <- !is.na(x@data.species)
            full.resp <- x@data.species[filter_PA]
            full.resp <- ifelse(full.resp == 1,10,20)
            full.xy <- x@coord[filter_PA,]
            full.df <- data.frame(resp = full.resp,
                                  x = full.xy[,1],
                                  y = full.xy[,2])
            
            full.df.vect <- vect(full.df, geom = c("x","y"))
            names(full.df.vect) <- "resp"
            full.df.vect$dataset <- "Initial dataset"
            
            ## 1.2 - Eval Dataset -----------------------------------------------
            
            if (plot.eval) {
              eval.resp <- x@eval.data.species
              eval.resp <- ifelse(eval.resp == 1,12,22)
              eval.xy <- x@eval.coord
              eval.df <- data.frame(resp = eval.resp,
                                    x = eval.xy[,1],
                                    y = eval.xy[,2])
              
              eval.df.vect <- vect(eval.df, geom = c("x","y"))
              names(eval.df.vect) <- "resp"
              eval.df.vect$dataset <- "Evaluation dataset"
              full.df.vect <- rbind(full.df.vect, eval.df.vect)
            }
            
            ## 1.3 - Pseudo Absences and CV dataset -----------------------------
            if (!is.null(calib.lines) | inherits(x, "BIOMOD.formated.data.PA")) {
              PA_run.vect <-
                foreach(this_PA = PA, .combine = 'rbind') %:%
                foreach(this_run = run, .combine = 'rbind') %do% { # to change if format is updated
                  
                  if( is.na(this_PA) ){ # run only
                    this_name <- this_run
                    this_calib <- calib.lines[ , this_run]
                    this_valid <- ! calib.lines[ , this_run]
                  } else if (is.na(this_run)){ # PA only
                    this_name <- this_PA
                    this_calib <- x@PA.table[ , this_PA]
                  } else { # PA+run
                    this_name <- paste0(this_PA,"_",this_run)
                    this_calib <- calib.lines[ , this_run] & x@PA.table[ , this_PA]
                    this_valid <- !calib.lines[ , this_run] & x@PA.table[ , this_PA]
                  }
                  
                  calib.resp <- x@data.species[this_calib]
                  calib.resp <- ifelse(is.na(calib.resp), 30, 
                                       ifelse(calib.resp == 1, 10, 20))
                  calib.xy <- x@coord[this_calib,]
                  calib.df <- data.frame(resp = calib.resp,
                                         x = calib.xy[, 1],
                                         y = calib.xy[, 2])
                  
                  if (!is.na(this_run)) { 
                    valid.resp <- x@data.species[this_valid]
                    valid.resp <- ifelse(is.na(valid.resp), 31, 
                                         ifelse(valid.resp == 1, 11, 21))
                    valid.xy <- x@coord[this_valid,]
                    valid.df <- data.frame(resp = valid.resp,
                                           x = valid.xy[, 1],
                                           y = valid.xy[, 2])
                    calib.df <- rbind(calib.df, valid.df)
                  }
                  
                  thisdf.vect <- vect(calib.df, geom = c("x","y"))
                  names(thisdf.vect) <- "resp"
                  thisdf.vect$dataset <- this_name
                  thisdf.vect
                }
              full.df.vect <- rbind(full.df.vect, PA_run.vect)
            }
            
            
            # 2- define colors and breaks ------------------------------------
            data_breaks <- c(9,10, 11, 12, # presence
                             19,20, 21, 22, # absence
                             29,30, 31,         # pseudo-absences
                             1)              # background
            data_labels <- c("9" = "**Presences**",
                             "10" = "Presences (calibration)",
                             "11" = "Presences (validation)",
                             "12" = "Presences (evaluation)",
                             "19" = "**True Absences**",
                             "20" = "True Absences (calibration)",
                             "21" = "True Absences (validation)",
                             "22" = "True Absences (evaluation)",
                             "29" = "**Pseudo-Absences**",
                             "30" = "Pseudo-Absences (calibration)",
                             "31" = "Pseudo-Absences (validation)",
                             "1" = "Background")
            data_labels_facet <- c("9" = "**Presences**",
                                   "10" = "calibration",
                                   "11" = "validation",
                                   "12" = "evaluation",
                                   "19" = "**True Absences**",
                                   "20" = "calibration",
                                   "21" = "validation",
                                   "22" = "evaluation",
                                   "29" = "**Pseudo-Absences**",
                                   "30" = "calibration",
                                   "31" = "validation",
                                   "1" = NA) # background
            
            data_colors <- c("9" = NA,
                             "10" = "#004488",
                             "11" = "#6699CC",
                             "12" = "#6699CC",
                             "19" = NA,
                             "20" = "#994455",
                             "21" = "#EE99AA",
                             "22" = "#EE99AA",
                             "29" = NA,
                             "30" = "#997700",
                             "31" = "#EECC66",
                             "1" = "grey70")
            
            shape_fit <- 16
            shape_eval <- 17
            data_shape <- c("9" = NA,
                            "10" = shape_fit,
                            "11" = shape_fit,
                            "12" = shape_eval,
                            "19" = NA,
                            "20" = shape_fit,
                            "21" = shape_fit,
                            "22" = shape_eval,
                            "29" = NA,
                            "30" = shape_fit,
                            "31" = shape_fit,
                            "1" = NA)
            data_alpha <- c("9" = 0,
                            "10" = 1,
                            "11" = 1,
                            "12" = 1,
                            "19" = 0,
                            "20" = 1,
                            "21" = 1,
                            "22" = 1,
                            "29" = 0,
                            "30" = 1,
                            "31" = 1,
                            "1"  = 0)
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
                            field = "resp", background = 1)
                names(this_rast) <- this_dataset
                this_rast*this_mask
              }
              
              if(plot.output == "facet"){
                g <- ggplot()+
                  tidyterra::geom_spatraster(data = rast.plot,
                                             aes(fill = factor(after_stat(value), data_breaks)))+
                  facet_wrap(~lyr)+
                  scale_fill_manual(
                    NULL,
                    breaks = data_breaks,
                    values = data_colors,
                    labels = data_labels_facet,
                    na.value = data_background, 
                    drop = FALSE)+
                  guides(fill = guide_legend(
                    override.aes = list(alpha = data_alpha),
                    ncol = 3))+
                  theme(legend.position = "top",
                        legend.key = element_blank(),
                        legend.background = element_rect(fill = "grey90"),
                        legend.text = ggtext::element_markdown())
                
              } else {
                g <- lapply(names(rast.plot), function(thisname){
                  ggplot()+
                    tidyterra::geom_spatraster(data = rast.plot[[thisname]],
                                               aes(fill = factor(after_stat(value))))+
                    scale_fill_manual(
                      NULL,
                      breaks = data_breaks,
                      values = data_colors,
                      labels = data_labels,
                      na.value = data_background)+
                    ggtitle(thisname)+
                    guides(color = guide_legend(nrow = 2))+
                    theme(legend.position = "top",
                          legend.key = element_blank(),
                          legend.background = element_rect(fill = "grey90"))
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
              if(plot.output == "facet"){
                base_g <-  ggplot(data.df)
                if(has.mask){
                  base_g <- base_g +
                    tidyterra::geom_spatraster(data = this_mask, aes(fill = factor(after_stat(value))))
                }
                
                
                g <- base_g +      
                  geom_point(aes(x = x, y = y, 
                                 color = factor(resp, levels = data_breaks[-12]),
                                 shape = factor(resp, levels = data_breaks[-12])), 
                             alpha = 1, size = 1.5)+
                  facet_wrap(~dataset)+
                  scale_color_manual(
                    NULL,
                    breaks = data_breaks,
                    values = data_colors,
                    labels = data_labels_facet,
                    drop = FALSE)+
                  scale_shape_manual(
                    NULL,
                    breaks = data_breaks,
                    values = data_shape,
                    labels = data_labels_facet,
                    drop = FALSE)+
                  scale_fill_manual(
                    guide = "none",
                    breaks = data_breaks,
                    values = data_colors,
                    labels = data_labels,
                    na.value = data_background)+
                  xlab(NULL)+ ylab(NULL)+
                  guides(color = guide_legend(override.aes = list(size = 3),
                                              ncol = 3))+
                  theme(legend.position = "top",
                        legend.key = element_blank(),
                        legend.background = element_rect(fill = "grey90"),
                        legend.text = ggtext::element_markdown())
                
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
                                   color = as.factor(resp),
                                   shape = as.factor(resp)),
                               alpha = 1, size = 1.5)+
                    scale_color_manual(
                      NULL,
                      breaks = data_breaks,
                      values = data_colors,
                      labels = data_labels)+
                    scale_shape_manual(
                      NULL,
                      breaks = data_breaks,
                      values = data_shape,
                      labels = data_labels)+
                    scale_fill_manual(
                      guide = "none",
                      breaks = data_breaks,
                      values = data_colors,
                      labels = data_labels,
                      na.value = data_background)+
                    xlab(NULL)+ ylab(NULL)+
                    guides(color = guide_legend(override.aes = list(size = 3),
                                                nrow = 2))+
                    theme(legend.position = "top",
                          legend.key = element_blank(),
                          legend.background = element_rect(fill = "grey90"))+
                    ggtitle(thisname)
                })
                
              }
              if(do.plot){
                print(g)
              }
              return(list("data.vect"  = full.df.vect,
                          "data.label" = data_labels,
                          "data.plot"  = g))
            }
          }
)


.plot.BIOMOD.formated.data.check.args <- function(x,
                                                  calib.lines,
                                                  plot.type,
                                                  plot.output, 
                                                  PA,
                                                  run,
                                                  plot.eval,
                                                  do.plot){
  
  
  ## 1 - check x -----------------------------------------
  .fun_testIfInherits(TRUE, "x", x, c("BIOMOD.formated.data", "BIOMOD.formated.data.PA"))
  
  ## 2 - check calib.lines & run -----------------------------------------
  
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- paste0("RUN", seq_len(ncol(calib.lines)))
    
    if (any(colnames(calib.lines) != c(expected_CVnames))) {
      stop("column names for `calib.lines` did not match ", deparse(expected_CVnames))  
    }
    
    if (missing(run)) run <- 'all'
    if ( ! ((length(run) == 1 && run == "all")
            || all(run %in% c(expected_CVnames)))  ) {
      stop("`run` must be 'all' or a combination among ", deparse(expected_CVnames))
    }
    if (length(run) == 1 && run == "all") {
      run <- expected_CVnames
    }
  } else {
    run <- NA
  } 
  
  
  ## 3 - check plot.eval ----------------------
  if (missing(plot.eval)) {
    plot.eval <- x@has.data.eval
  } else {
    stopifnot(is.logical(plot.eval))
    if(plot.eval & !x@has.data.eval){
      plot.eval <- FALSE
      cat('\n  ! Evaluation data are missing and its plot was deactivated')
    }
  }
  
  ## 4 are proper mask available ? -----------------------
  has.mask <- rast.has.values(rast(x@data.mask[["calibration"]]))
  if(plot.eval){
    has.mask.eval <- length(x@data.mask) > 1
  } else {
    has.mask.eval <- FALSE
  }
  if (has.mask | has.mask.eval) {  
    if (!requireNamespace("tidyterra")) {
      stop("Package `tidyterra` is missing. Please install it with `install.packages('tidyterra')`.")
    }
  } 
  ## 5 - check plot.type  ----------------------
  if (missing(plot.type)) {
    plot.type <- "points"
  } else {
    .fun_testIfIn(TRUE, "plot.type", plot.type, c("raster","points"))
    if ( !has.mask & plot.type == "raster") {
      plot.type <- "points"
      cat("\n ! no raster available, `plot.type` automatically set to 'points'\n")
    }
  }
  
  ## 6 - plot.output----------------------
  if (missing(plot.output)) {
    plot.output <- "facet"
  } else {
    .fun_testIfIn(TRUE, "plot.output", plot.output, c("facet","list"))
  }
  
  if(plot.output == "facet"){
    if(!requireNamespace("ggtext")){
      stop("Package `ggtext` is missing. Please install it with `install.packages('ggtext')`.")
    }
  }
  
  
  
  ## 7 - do.plot ----------------------
  # do.plot
  stopifnot(is.logical(do.plot))
  
  ## 8 - check PA ----------------------
  
  if (inherits(x, "BIOMOD.formated.data.PA")) {
    
    if (missing(PA)) PA <- 'all'
    
    expected_PAnames <- colnames(x@PA.table)
    
    if ( ! ((length(PA) == 1 && PA == "all")
            || all(PA %in% c(expected_PAnames)))  ) {
      stop("`PA` must be 'all' or a combination among ", deparse(expected_PAnames))
    }
    if (length(PA) == 1 && PA == "all") {
      PA <- expected_PAnames
    }
    
  } else {
    PA <- NA
  }
  
  ##  9 - check that coordinates are available -------------------------------
  
  if(nrow(x@coord) == 0){
    stop("coordinates are required to plot BIOMOD.formated.data objects")
  }
  

  ## End - return arguments ----------------------------------------------------
  
  
  return(list(x = x,
              calib.lines = calib.lines,
              plot.type = plot.type,
              plot.output = plot.output, 
              PA = PA,
              run = run,
              plot.eval = plot.eval,
              do.plot = do.plot,
              has.mask = has.mask,
              has.mask.eval = has.mask.eval))
}

### show.BIOMOD.formated.data  --------------------------------------------------
##' 
##' @rdname BIOMOD.formated.data
##' @importMethodsFrom methods show
##' @export
##' 

setMethod('show', signature('BIOMOD.formated.data'),
          function(object)
          {
            .bm_cat("BIOMOD.formated.data")
            cat("\ndir.name = ", object@dir.name, fill = .Options$width)
            cat("\nsp.name = ", object@sp.name, fill = .Options$width)
            cat("\n\t",
                sum(object@data.species, na.rm = TRUE),
                'presences, ',
                sum(object@data.species == 0, na.rm = TRUE),
                'true absences and ',
                sum(is.na(object@data.species), na.rm = TRUE),
                'undefined points in dataset',
                fill = .Options$width)
            cat("\n\n\t",
                ncol(object@data.env.var),
                'explanatory variables\n',
                fill = .Options$width)
            print(summary(object@data.env.var))
            
            if (object@has.data.eval) {
              cat("\n\nEvaluation data :", fill = .Options$width)
              cat("\n\t",
                  sum(object@eval.data.species, na.rm = TRUE),
                  'presences, ',
                  sum(object@eval.data.species == 0, na.rm = TRUE),
                  'true absences and ',
                  sum(is.na(object@eval.data.species), na.rm = TRUE),
                  'undefined points in dataset',
                  fill = .Options$width)
              cat("\n\n", fill = .Options$width)
              print(summary(object@eval.data.env.var))
            }
            
            if(inherits(object, "biomod.formated.data.PA")){
              cat(
                "\n\n",
                ncol(object@PA.table),
                'Pseudo Absences dataset available (',
                colnames(object@PA.table),
                ") with ",
                sum(object@PA.table[, 1], na.rm = TRUE) - sum(object@data.species, na.rm = TRUE),
                'absences in each (true abs + pseudo abs)',
                fill = .Options$width
              )
            }
            .bm_cat()
          }
)


### summary.BIOMOD.formated.data  --------------------------------------------------
##' 
##' @rdname summary
##' @docType methods
##' @author Remi Patin
##' 
##' @title \code{summary} method for \code{\link{BIOMOD.formated.data}} object class
##' 
##' @description Summarize the number of presences, absences and 
##' pseudo-absences among the different potential dataset (calibration, 
##' validation and evaluation).  
##' 
##' @param object a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} 
##' object returned by the \code{\link{BIOMOD_FormatingData}} function
##' @param calib.lines (\emph{optional, default} \code{NULL}) \cr
##' an \code{array} object returned by \code{\link{get_calib_lines}} or 
##' \code{\link{BIOMOD_CrossValidation}} functions, to explore the distribution of calibration 
##' and validation datasets
##' 
##' 
##' @return a \code{data.frame} 
##' 
##' @export
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
##' summary(myBiomodData)
##' 
##' 

setMethod('summary', signature(object = 'BIOMOD.formated.data'),
          function(object, calib.lines = NULL)
          {
            args <- .summary.BIOMOD.formated.data.check.args(object = object, calib.lines = calib.lines)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            output <- data.frame("dataset" = "initial",
                                 "run" = NA,
                                 "PA" = NA,
                                 "Presences" = sum(object@data.species, na.rm = TRUE),
                                 "True_Absences" = sum(object@data.species == 0, na.rm = TRUE),
                                 "Pseudo_Absences" = 0,
                                 "Undefined" = sum(is.na(object@data.species), na.rm = TRUE))
            
            if (object@has.data.eval) {
              output <- rbind(output,
                              data.frame("dataset" = "evaluation",
                                         "run" = NA,
                                         "PA" = NA,
                                         "Presences" = sum(object@eval.data.species, na.rm = TRUE),
                                         "True_Absences" = sum(object@eval.data.species == 0, na.rm = TRUE),
                                         "Pseudo_Absences" = 0,
                                         "Undefined" = sum(is.na(object@eval.data.species), na.rm = TRUE)))
            }
            
            PA <- NA
            run <- NA
            if (inherits(object, "BIOMOD.formated.data.PA")) {
              PA <- colnames(object@PA.table)
            }
            if (!is.null(calib.lines)) {
              run <- colnames(calib.lines)
            }
            
            if (!is.null(calib.lines) || inherits(object, "BIOMOD.formated.data.PA")) {
              output <- rbind(output,
                              foreach(this_PA = PA, .combine = 'rbind') %:%
                                foreach(this_run = run, .combine = 'rbind') %do% {
                                  if ( is.na(this_PA) ) { # run only
                                    this_name <- this_run
                                    this_calib <- calib.lines[ , this_run]
                                    this_valid <- ! calib.lines[ , this_run]
                                  } else if (is.na(this_run)) { # PA only
                                    this_name <- this_PA
                                    this_calib <-
                                      object@PA.table[ , this_PA]
                                  } else { # PA+run
                                    this_name <- paste0(this_PA,"_",this_run)
                                    this_calib <- calib.lines[ , this_run] & object@PA.table[ , this_PA]
                                    this_valid <- ! calib.lines[ , this_run] & object@PA.table[ , this_PA]
                                  }
                                  
                                  calib.resp <- object@data.species[this_calib]
                                  tmp <- data.frame("dataset" = "calibration",
                                                    "run" = this_run,
                                                    "PA" = this_PA,
                                                    "Presences" = sum(calib.resp, na.rm = TRUE),
                                                    "True_Absences" = sum(calib.resp == 0, na.rm = TRUE),
                                                    "Pseudo_Absences" = length(which(is.na(calib.resp))),
                                                    "Undefined" = NA)
                                  
                                  if (!is.na(this_run)) { 
                                    valid.resp <- object@data.species[this_valid]
                                    tmp <- rbind(tmp,
                                                 data.frame("dataset" = "validation",
                                                            "run" = this_run,
                                                            "PA" = this_PA,
                                                            "Presences" = sum(valid.resp, na.rm = TRUE),
                                                            "True_Absences" = sum(valid.resp == 0, na.rm = TRUE),
                                                            "Pseudo_Absences" = length(which(is.na(valid.resp))),
                                                            "Undefined" = NA))
                                    
                                  }
                                  return(tmp) # end foreach
                                })
            } 
            output$Total_Absences <- output$True_Absences + output$Pseudo_Absences
            output
          }
)

.summary.BIOMOD.formated.data.check.args <- function(object,
                                                     calib.lines){
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- paste0("RUN", seq_len(ncol(calib.lines)))
    
    if (any(colnames(calib.lines) != c(expected_CVnames))) {
      stop("column names for `calib.lines` did not match ", deparse(expected_CVnames))  
    }
  }
  return(list(object = object,
              calib.lines = calib.lines))
}


## --------------------------------------------------------------------------- #
## 2. BIOMOD.formated.data.PA ------------------------------------------------
## this class inherits from BIOMOD.formated.data and have one more slot 'PA', 
## giving PA selected
## --------------------------------------------------------------------------- #

##' @name BIOMOD.formated.data.PA
##' @aliases BIOMOD.formated.data.PA-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_FormatingData()} output object class (with pseudo-absences)
##' 
##' @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
##' \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_CrossValidation}} and 
##' \code{\link{BIOMOD_Modeling}}
##' 
##' @inheritParams BIOMOD.formated.data
##' @param dir.name a \code{character} corresponding to the modeling folder
##' @param sp.name a \code{character} corresponding to the species name
##' 
##' @param PA.nb.rep (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection, an \code{integer} corresponding to the number of sets 
##' (repetitions) of pseudo-absence points that will be drawn
##' @param PA.strategy (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection, a \code{character} defining the strategy that will be used to 
##' select the pseudo-absence points. Must be \code{random}, \code{sre}, \code{disk} or 
##' \code{user.defined} (see Details)
##' @param PA.nb.absences (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection, and \code{PA.strategy = 'random'} or \code{PA.strategy = 'sre'} 
##' or \code{PA.strategy = 'disk'}, an \code{integer} corresponding to the number of pseudo-absence 
##' points that will be selected for each pseudo-absence repetition (true absences included)
##' @param PA.sre.quant (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'sre'}, a \code{numeric} between \code{0} 
##' and \code{0.5} defining the half-quantile used to make the \code{sre} pseudo-absence selection 
##' (see Details)
##' @param PA.dist.min (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' minimal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in meters, see Details)
##' @param PA.dist.max (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' maximal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in meters, see Details)
##' @param PA.user.table (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'user.defined'}, a \code{matrix} or 
##' \code{data.frame} with as many rows as \code{resp.var} values, as many columns as 
##' \code{PA.nb.rep}, and containing \code{TRUE} or \code{FALSE} values defining which points 
##' will be used to build the species distribution model(s) for each repetition (see Details)
##' 
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether points having one or several missing
##' values for explanatory variables should be removed from the analysis or not
##' @param filter.raster (\emph{optional, default} \code{FALSE}) \cr 
##' If \code{env} is of raster type, a \code{logical} value defining whether \code{sp} 
##' is to be filtered when several points occur in the same raster cell
##' 
##' @slot dir.name a \code{character} corresponding to the modeling folder
##' @slot sp.name a \code{character} corresponding to the species name
##' @slot coord a 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates
##' @slot data.species a \code{vector} containing the species observations (\code{0}, \code{1} or 
##' \code{NA})
##' @slot data.env.var a \code{data.frame} containing explanatory variables
##' @slot data.mask a \code{\link[terra:rast]{SpatRaster}} object containing 
##' the mask of the studied area
##' @slot has.data.eval a \code{logical} value defining whether evaluation data is given
##' @slot eval.coord (\emph{optional, default} \code{NULL}) \cr 
##' A 2-columns \code{data.frame} containing the corresponding \code{X} and \code{Y} 
##' coordinates for evaluation data
##' @slot eval.data.species (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing the species observations (\code{0}, \code{1} or \code{NA}) for 
##' evaluation data
##' @slot eval.data.env.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{data.frame} containing explanatory variables for evaluation data
##' @slot PA.strategy a \code{character} corresponding to the pseudo-absence selection strategy
##' @slot PA.table a \code{data.frame} containing the corresponding table of selected 
##' pseudo-absences (indicated by \code{TRUE} or \code{FALSE}) from the \code{pa.tab} list 
##' element returned by the \code{\link{bm_PseudoAbsences}} function
##' 
##' 
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{bm_PseudoAbsences}}, 
##' \code{\link{BIOMOD_Tuning}}, \code{\link{BIOMOD_CrossValidation}}, 
##' \code{\link{BIOMOD_Modeling}}, \code{\link{bm_RunModelsLoop}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.formated.data.PA")
##' 
##' ## ----------------------------------------------------------------------- #
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
##' # Format Data with pseudo-absences : random method
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName,
##'                                      PA.nb.rep = 0,
##'                                      PA.strategy = 'random',
##'                                      PA.nb.absences = 1000)
##' myBiomodData
##' plot(myBiomodData)
##' 
##' 
NULL

##' @name BIOMOD.formated.data.PA-class
##' @rdname BIOMOD.formated.data.PA
##' 
##' @importFrom terra rast app is.factor subset extract 
##' cellFromXY `add<-` crds vect
##' @export
##' 

# 2.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data.PA",
         contains = "BIOMOD.formated.data",
         representation(PA.strategy = 'character', PA.table = 'data.frame'),
         validity = function(object) { return(TRUE) })


# 2.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.formated.data.PA", 
           def = function(sp, env, ...) { 
             standardGeneric("BIOMOD.formated.data.PA") 
           })

### BIOMOD.formated.data.PA(sp = numeric, env = data.frame) --------------------
##' 
##' @rdname BIOMOD.formated.data.PA
##' @export
##' 

setMethod('BIOMOD.formated.data.PA', signature(sp = 'numeric', env = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , PA.nb.rep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                   , PA.dist.min = 0, PA.dist.max = NULL
                   , PA.sre.quant = 0.025, PA.user.table = NULL
                   , na.rm = TRUE, filter.raster = FALSE) {
            .BIOMOD.formated.data.PA(sp, env, xy, dir.name, sp.name
                                     , eval.sp, eval.env, eval.xy
                                     , PA.nb.rep, PA.strategy, PA.nb.absences
                                     , PA.dist.min, PA.dist.max
                                     , PA.sre.quant, PA.user.table
                                     , na.rm
                                     , filter.raster = filter.raster)
          })

### BIOMOD.formated.data.PA(sp = numeric, env = SpatRaster) -------------------
##' 
##' @rdname BIOMOD.formated.data.PA
##' @export
##' 

setMethod('BIOMOD.formated.data.PA', signature(sp = 'numeric', env = 'SpatRaster'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , PA.nb.rep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                   , PA.dist.min = 0, PA.dist.max = NULL
                   , PA.sre.quant = 0.025, PA.user.table = NULL
                   , na.rm = TRUE, filter.raster = FALSE) {
            .BIOMOD.formated.data.PA(sp, env, xy, dir.name, sp.name
                                     , eval.sp, eval.env, eval.xy
                                     , PA.nb.rep, PA.strategy, PA.nb.absences
                                     , PA.dist.min, PA.dist.max
                                     , PA.sre.quant, PA.user.table
                                     , na.rm
                                     , filter.raster = filter.raster)
          })

### .BIOMOD.formated.data.PA ---------------------------------------------------
.BIOMOD.formated.data.PA <-  function(sp, env, xy, dir.name, sp.name
                                      , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                                      , PA.nb.rep = 1, PA.strategy = 'random', PA.nb.absences = NULL
                                      , PA.dist.min = 0, PA.dist.max = NULL
                                      , PA.sre.quant = 0.025, PA.user.table = NULL
                                      , na.rm = TRUE, filter.raster = FALSE)
{
  

  #### check for categorical vars and filter duplicated data points

  categorical_var <- NULL
  if (inherits(env, 'SpatRaster')) {
    categorical_var <- names(env)[is.factor(env)] 
    
    output <- check_duplicated_cells(env, xy, sp, filter.raster)
    xy <- output$xy
    sp <- output$sp
    rm(output)
    
  }
  

  # Convert sp in SpatVector
  if (is.numeric(sp)) {
    if (nrow(xy) == 0) {
      sp.df <- data.frame(x = 0,
                          y = 0,
                          resp = sp)
    } else {
      sp.df <- data.frame(x = xy[,1],
                          y = xy[,2],
                          resp = sp)
    }
    sp <- vect(sp.df, geom = c("x","y"))
  }
  
  pa.data.tmp <- bm_PseudoAbsences(resp.var = sp,
                                   expl.var = env,
                                   nb.rep = PA.nb.rep,
                                   strategy = PA.strategy,
                                   nb.absences = PA.nb.absences,
                                   sre.quant = PA.sre.quant,
                                   dist.min = PA.dist.min,
                                   dist.max = PA.dist.max,
                                   user.table = PA.user.table)
  
  if (!is.null(pa.data.tmp)) {
    ## Keep same env variable for eval than calib (+ check for factor)
    if (length(categorical_var) > 0) {
      for (cat_var in categorical_var) {
        pa.data.tmp$env[, cat_var] <- as.factor(pa.data.tmp$env[, cat_var])
      }
    }
    
    ## REMOVE NA IF ANY
    if (na.rm) {
      rowToRm <- unique(unlist(lapply(pa.data.tmp$env, function(x) { return(which(is.na(x))) })))
      if (length(rowToRm) > 0) {
        cat("\n ! Some NAs have been automatically removed from your data")
        pa.data.tmp$xy <- pa.data.tmp$xy[-rowToRm, , drop = FALSE]
        pa.data.tmp$sp <- pa.data.tmp$sp[-rowToRm, drop = FALSE]
        pa.data.tmp$env <- pa.data.tmp$env[-rowToRm, , drop = FALSE]
        pa.data.tmp$pa.tab <- pa.data.tmp$pa.tab[-rowToRm, , drop = FALSE]
      }
    }

    if(!inherits(env,"SpatRaster")){
      env <- pa.data.tmp$env
    }
    
    BFD <- BIOMOD.formated.data(sp = pa.data.tmp$sp,
                                env = env,
                                xy = as.data.frame(pa.data.tmp$xy),
                                dir.name = dir.name,
                                sp.name = sp.name, 
                                eval.sp = eval.sp,
                                eval.env = eval.env,
                                eval.xy = eval.xy,
                                na.rm = na.rm,
                                filter.raster = filter.raster)
    
    BFDP <- new('BIOMOD.formated.data.PA',
                dir.name = BFD@dir.name,
                sp.name = BFD@sp.name,
                coord = BFD@coord,
                data.env.var = BFD@data.env.var,
                data.species = BFD@data.species,
                data.mask = BFD@data.mask,
                has.data.eval = BFD@has.data.eval,
                eval.coord = BFD@eval.coord,
                eval.data.species = BFD@eval.data.species,
                eval.data.env.var = BFD@eval.data.env.var,
                PA.strategy = PA.strategy,
                PA.table = as.data.frame(pa.data.tmp$pa.tab))
    
    rm(list = 'BFD')
  } else {
    cat("\n   ! PA selection not done", fill = .Options$width)
    BFDP <- BIOMOD.formated.data(sp = as.vector(values(sp)[,1]),
                                 env = env,
                                 xy = xy,
                                 dir.name = dir.name,
                                 sp.name = sp.name, 
                                 eval.sp = eval.sp,
                                 eval.env = eval.env,
                                 eval.xy = eval.xy)
  }
  rm(list = "pa.data.tmp" )
  return(BFDP)
}


## --------------------------------------------------------------------------- #
## 3. BIOMOD.models.options --------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.models.options
##' @aliases BIOMOD.models.options-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_ModelingOptions()} output object class
##' 
##' @description Class returned by \code{\link{BIOMOD_ModelingOptions}}, and used by 
##' \code{\link{BIOMOD_Tuning}} and \code{\link{BIOMOD_Modeling}}
##' 
##' 
##' @slot GLM a \code{list} containing GLM options
##' @slot GBM a \code{list} containing GBM options
##' @slot GAM a \code{list} containing GAM options
##' @slot CTA a \code{list} containing CTA options
##' @slot ANN a \code{list} containing ANN options
##' @slot SRE a \code{list} containing SRE options
##' @slot FDA a \code{list} containing FDA options
##' @slot MARS a \code{list} containing MARS options
##' @slot RF a \code{list} containing RF options
##' @slot MAXENT a \code{list} containing MAXENT options
##' @slot MAXNET a \code{list} containing MAXNET options
##' 
##' @param object a \code{\link{BIOMOD.models.options}} object
##' 
##' 
##' @seealso \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Tuning}}, 
##' \code{\link{BIOMOD_Modeling}}
##' @family Toolbox objects
##' 
##' 
##' @examples
##' 
##' showClass("BIOMOD.models.options")
##' 
##' ## ----------------------------------------------------------------------- #
##' ## default BIOMOD.models.options object
##' myBiomodOptions <- BIOMOD_ModelingOptions()
##'
##' ## print the object
##' myBiomodOptions
##'
##' 
NULL

##' @name BIOMOD.models.options-class
##' @rdname BIOMOD.models.options
##' @export
##' 

setClass("BIOMOD.models.options",
         representation(GLM = "list",
                        GBM = "list",
                        GAM = "list",
                        CTA = "list",
                        ANN = "list",
                        SRE = "list",
                        FDA = "list",
                        MARS = "list",
                        RF = "list",
                        MAXENT = "list",
                        MAXNET = "list"),
         prototype(GLM = list(type = 'quadratic',
                              interaction.level = 0,
                              myFormula = NULL,
                              test = 'AIC',
                              family = binomial(link = 'logit'),
                              mustart = 0.5,
                              control = glm.control(maxit = 50)),
                   GBM = list(distribution = 'bernoulli',
                              n.trees = 2500,
                              interaction.depth = 7,
                              n.minobsinnode = 5,
                              shrinkage = 0.001,
                              bag.fraction = 0.5,
                              train.fraction = 1,
                              cv.folds = 3,
                              keep.data = FALSE,
                              verbose = FALSE,
                              # class.stratify.cv = 'bernoulli',
                              perf.method = 'cv',
                              n.cores = 1),
                   GAM = list(algo = "GAM_mgcv",
                              type = "s_smoother",
                              k = NULL,
                              interaction.level = 0,
                              myFormula = NULL,
                              family = binomial(link = 'logit'),
                              control = list(epsilon = 1e-06, trace = FALSE, maxit = 100),
                              method = "GCV.Cp",
                              optimizer = c("outer", "newton"),
                              select = FALSE,
                              knots = NULL,
                              paraPen = NULL),
                   CTA = list(method = 'class',
                              parms = 'default',
                              # control = rpart.control(xval = 5, minbucket = 5, minsplit = 5, cp = 0.001, maxdepth = 25),
                              control = list(xval = 5, 
                                             minbucket = 5, 
                                             minsplit = 5,
                                             cp = 0.001, 
                                             maxdepth = 25),
                              cost = NULL),
                   ANN = list(NbCV = 5,
                              size = NULL,
                              decay = NULL,
                              rang = 0.1,
                              maxit = 200),
                   SRE = list(quant = 0.025),
                   FDA = list(method = 'mars', add_args = NULL),
                   MARS = list(type = 'simple',
                               interaction.level = 0,
                               myFormula = NULL,
                               # degree = 1,
                               nk = NULL,
                               penalty = 2,
                               thresh = 0.001,
                               nprune = NULL,
                               pmethod = 'backward'),
                   RF = list(do.classif = TRUE,
                             ntree = 500,
                             mtry = 'default',
                             sampsize = NULL,
                             nodesize = 5,
                             maxnodes = NULL),
                   MAXENT = list(path_to_maxent.jar = getwd(),
                                          memory_allocated = 512,
                                          initial_heap_size = NULL,
                                          max_heap_size = NULL,
                                          background_data_dir = 'default',
                                          maximumbackground = 'default',
                                          maximumiterations = 200,
                                          visible = FALSE,
                                          linear = TRUE,
                                          quadratic = TRUE,
                                          product = TRUE,
                                          threshold = TRUE,
                                          hinge = TRUE,
                                          lq2lqptthreshold = 80,
                                          l2lqthreshold = 10,
                                          hingethreshold = 15,
                                          beta_threshold = -1.0,
                                          beta_categorical = -1.0,
                                          beta_lqp = -1.0,
                                          beta_hinge = -1.0,
                                          betamultiplier = 1,
                                          defaultprevalence = 0.5),
                   MAXNET = list(myFormula = NULL,
                                            regmult = 1,
                                            regfun = maxnet::maxnet.default.regularization)
         ),
         validity = function(object) {
           test <- TRUE
           
           ## GLM ------------------------------------------------------------
           test <- .fun_testIfIn(test, "GLM$type", object@GLM$type,
                                 c("simple", "quadratic", "polynomial", "user.defined"))
           test <- .fun_testIfPosInt(test, "GLM$interaction.level", object@GLM$interaction.level)
           
           if (!is.null(object@GLM$myFormula) && 
               !inherits(object@GLM$myFormula, "formula")) {
             cat("\nGLM$myFormula must be NULL or a formula object")
             test <- FALSE
           }
           
           test <- .fun_testIfIn(test, "GLM$test", object@GLM$test, c("AIC", "BIC", "none"))
           fam <- "none"
           if (!inherits(object@GLM$family, "family")) {
             cat("\nGLM$family must be a valid family object")
             test <- FALSE
           }
           if (!is.list(object@GLM$control)) {
             cat("\nGLM$control must be a list like that returned by glm.control")
             test <- FALSE
           }
           
           ## GBM ------------------------------------------------------------
           test <- .fun_testIfIn(test, "GBM$distribution", object@GBM$distribution, 
                                 c("bernoulli", "huberized", "multinomial", "adaboost"))
           # test <- .fun_testIfPosInt(test, "GBM$n.trees", object@GBM$n.trees)
           if (!is.numeric(object@GBM$n.trees)) {
             cat("\nGBM$n.trees must be a integer")
             test <- FALSE
           } else {
             if (object@GBM$n.trees < 0 | floor(object@GBM$n.trees) != object@GBM$n.trees) {
               cat("\nGBM$n.trees must be a positive integer")
               test <- FALSE
             }
           }
           test <- .fun_testIfPosInt(test, "GBM$interaction.depth", object@GBM$interaction.depth)
           test <- .fun_testIfPosInt(test, "GBM$n.minobsinnode", object@GBM$n.minobsinnode)
           test <- .fun_testIfPosNum(test, "GBM$shrinkage", object@GBM$shrinkage)
           test <- .fun_testIf01(test, "GBM$bag.fraction", object@GBM$bag.fraction)
           test <- .fun_testIf01(test, "GBM$train.fraction", object@GBM$train.fraction)
           test <- .fun_testIfPosInt(test, "GBM$cv.folds", object@GBM$cv.folds)
           if (!is.logical(object@GBM$keep.data)) {
             cat("\nGBM$keep.data must be a logical")
             test <- FALSE
           }
           if (!is.logical(object@GBM$verbose)) {
             cat("\nGBM$verbose must be a logical")
             test <- FALSE
           }
           # test <- .fun_testIfIn(test, "GBM$class.stratify.cv", object@GBM$class.stratify.cv, c('bernoulli', 'multinomial'))
           test <- .fun_testIfIn(test, "GBM$perf.method", object@GBM$perf.method, 
                                 c('OOB', 'test', 'cv'))
           
           ## GAM ------------------------------------------------------------
           test <- .fun_testIfIn(test, "GAM$algo", object@GAM$algo,
                                 c('GAM_mgcv', 'GAM_gam', 'BAM_mgcv'))
           test <- .fun_testIfIn(test, "GAM$type", object@GAM$type,
                                 c('s_smoother', 's', 'lo', 'te'))
           if (!is.null(object@GAM$k)) {
             if (!is.numeric(object@GAM$k)) {
               cat("\nGAM$k must be a integer")
               test <- FALSE
             } else {
               if (object@GAM$k < -1 | object@GAM$k %% 1 != 0) {
                 cat("\nGAM$k must be > -1")
                 test <- FALSE
               }
             }
           }
           test <- .fun_testIfPosInt(test, "GAM$interaction.level", object@GAM$interaction.level)
           if (!is.null(object@GAM$myFormula) &&
               (!inherits(object@GAM$myFormula, "formula"))) {
             cat("\nGAM$myFormula must be NULL or a formula object")
             test <- FALSE
           }
           
           if (!inherits(object@GAM$family, "family")) {
             cat("\nGAM$family must be a valid family object")
             test <- FALSE
           }
           if (!is.list(object@GAM$control)) {
             cat("\nGAM$control must be a list like that returned by gam.control")
             test <- FALSE
           }
           test <- .fun_testIfIn(test, "GAM$method", object@GAM$method,
                                 c("GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", "P-ML"))
           if (any(!object@GAM$optimizer %in% 
                   c("perf", "outer", "newton", "bfgs", "optim", "nlm", "nlm.fd"))) {
             cat("\nGAM$optimizer bad definition (see ?mgcv::gam)")
             test <- FALSE
           }
           if (!is.logical(object@GAM$select)) {
             cat("\nGAM$select must be a logical")
             test <- FALSE
           }
           #            knots=NULL,
           #            paraPen=NULL
           
           ## CTA ------------------------------------------------------------
           test <- .fun_testIfIn(test, "CTA$method", object@CTA$method, c("anova", "poisson", "class", "exp"))
           #parms = 'default',
           if (!is.list(object@CTA$control)) {
             cat("\nCTA$control must be a list like that returned by rpart.control")
             test <- FALSE
           }
           if (length(object@CTA$cost) > 0 && (!is.numeric(object@CTA$cost))) {
             cat("\nCTA$cost must be a non negative cost vector")
             test <- FALSE
           }
           
           
           ## ANN ------------------------------------------------------------
           test <- .fun_testIfPosInt(test, "ANN$NbCV", object@ANN$NbCV)
           if ((is.null(object@ANN$size) || length(object@ANN$size) > 1) &&
               object@ANN$NbCV <= 0) {
             cat("\nANN$size has to be defined as a single integer if ANN$NbCV=0")
             test <- FALSE
           } else {
             if (!is.null(object@ANN$size) &&
                 (!is.numeric(object@ANN$size) || !all(object@ANN$size > 0) || !all(object@ANN$size %% 1 == 0))) {
               cat("\nANN$size must be NULL or a positive (vector of) integer")
               test <- FALSE
             }
           }
           
           
           if ((is.null(object@ANN$decay) | length(object@ANN$decay) > 1) &&
               object@ANN$NbCV <= 0) {
             cat("\nANN$decay has to be defined as a single number if ANN$NbCV=0")
             test <- FALSE
           } else if (!is.null(object@ANN$decay) &&
                      (!is.numeric(object@ANN$decay) || !all(object@ANN$decay > 0))) {
             cat("\nANN$decay must be NULL or a positive (vector of) number")
             test <- FALSE
           }
           
           test <- .fun_testIfPosNum(test, "ANN$rang", object@ANN$rang)
           test <- .fun_testIfPosInt(test, "ANN$maxit", object@ANN$maxit)
           
           ## FDA ------------------------------------------------------------
           test <- .fun_testIfIn(test, "FDA$method", object@FDA$method, c('polyreg', 'mars', 'bruto'))
           if (!is.null(object@FDA$add_args) &&
               !is.list(object@FDA$add_args)) {
             cat("\nFDA$add_args must be a list or NULL"); 
             test <- FALSE
           }
           
           
           ## SRE ------------------------------------------------------------
           if (!is.numeric(object@SRE$quant)) {
             cat("\nSRE$quant must be a numeric"); test <- FALSE
           } else if(object@SRE$quant >= 0.5 | object@SRE$quant < 0){ 
             cat("\nSRE$quant must between 0 and 0.5"); test <- FALSE 
           }
           
           ## MARS ------------------------------------------------------------
           test <- .fun_testIfIn(test, "MARS$type", object@MARS$type, c("simple", "quadratic", "polynomial", "user.defined"))
           test <- .fun_testIfPosInt(test, "MARS$interaction.level", object@MARS$interaction.level)
           if (!is.null(object@MARS$myFormula) && !inherits(object@MARS$myFormula, "formula")) {
             cat("\nMARS$myFormula must be NULL or a formula object")
             test <- FALSE
           }
           # test <- .fun_testIfPosInt(test, "MARS$degree", object@MARS$degree)
           if (!is.null(object@MARS$nk) && 
               (object@MARS$nk < 0 | object@MARS$nk%%1!=0)) {
             cat("\nMARS$nk must be a positive integer or NULL if you want to use default parameter")
             test <- FALSE 
           }
           
           test <- .fun_testIfPosInt(test, "MARS$penalty", object@MARS$penalty)
           test <- .fun_testIfPosNum(test, "MARS$thresh", object@MARS$thresh)
           if (!is.null(object@MARS$nprune) &&
               !is.numeric(object@MARS$nprune)){
             cat("\nMARS$nprune must be a numeric or NULL"); test <- FALSE 
           }
           
           supported.pmethod <- c('backward', 'none', 'exhaustive', 'forward', 'seqrep', 'cv')
           if(!is.element(object@MARS$pmethod, supported.pmethod)){
             cat("\nMARS$pmethod must be a one of", supported.pmethod);
             test <- FALSE 
           }
           
           ## RF ------------------------------------------------------------
           if (!is.logical(object@RF$do.classif)) {
             cat("\nRF$do.classif must be a logical")
             test <- FALSE
           }
           test <- .fun_testIfPosInt(test, "RF$ntree", object@RF$ntree)
           if (object@RF$mtry != "default") {
             test <- .fun_testIfPosInt(test, "RF$mtry", object@RF$mtry)
           }
           if (!is.null(object@RF$sampsize)) {
             test <- .fun_testIfPosInt(test, "RF$sampsize", object@RF$sampsize)
           }
           test <- .fun_testIfPosInt(test, "RF$nodesize", object@RF$nodesize)
           if (length(object@RF$maxnodes) > 0) {
             test <- .fun_testIfPosInt(test, "RF$maxnodes", object@RF$maxnodes)
           }
           
           ## MAXENT ------------------------------------------------------------
           if (!is.character(object@MAXENT$path_to_maxent.jar)) {
             cat("\nMAXENT$path_to_maxent.jar must be a character")
             test <- FALSE
           }
           if (!is.null(object@MAXENT$memory_allocated)) {
             if (!is.numeric(object@MAXENT$memory_allocated)) {
               cat("\nMAXENT$memory_allocated must be a positive integer or NULL for unlimited memory allocation")
               test <- FALSE
             }
           }
           if (!is.character(object@MAXENT$background_data_dir)) {
             cat("\nMAXENT$background_data_dir must be 'default' (=> use the same pseudo absences than other models as background) or a path to the directory where your environmental layer are stored")
             test <- FALSE
           }
           tt <- is.character(object@MAXENT$maximumbackground) | is.numeric(object@MAXENT$maximumbackground)
           if (is.character(object@MAXENT$maximumbackground)) if (object@MAXENT$maximumbackground != "default") tt <- FALSE
           if (!tt) {
             cat("\nMAXENT$maximumbackground must be 'default' or numeric")
             test <- FALSE
           }
           test <- .fun_testIfPosInt(test, "MAXENT$maximumiterations", object@MAXENT$maximumiterations)
           if (!is.logical(object@MAXENT$visible)) {
             cat("\nMAXENT$visible must be a logical")
             test <- FALSE
           }
           if (!is.logical(object@MAXENT$linear)) {
             cat("\nMAXENT$linear must be a logical")
             test <- FALSE
           }
           if (!is.logical(object@MAXENT$quadratic)) {
             cat("\nMAXENT$quadratic must be a logical")
             test <- FALSE
           }
           if (!is.logical(object@MAXENT$product)) {
             cat("\nMAXENT$product must be a logical")
             test <- FALSE
           }
           if (!is.logical(object@MAXENT$threshold)) {
             cat("\nMAXENT$threshold must be a logical")
             test <- FALSE
           }
           if (!is.logical(object@MAXENT$hinge)) {
             cat("\nMAXENT$hinge must be a logical")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$lq2lqptthreshold)) {
             cat("\nMAXENT$lq2lqptthreshold must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$l2lqthreshold)) {
             cat("\nMAXENT$l2lqthreshold must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$lq2lqptthreshold)) {
             cat("\nMAXENT$lq2lqptthreshold must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$hingethreshold)) {
             cat("\nMAXENT$hingethreshold must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$beta_threshold)) {
             cat("\nMAXENT$beta_threshold must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$beta_categorical)) {
             cat("\nMAXENT$beta_categorical must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$beta_lqp)) {
             cat("\nMAXENT$beta_lqp must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$beta_hinge)) {
             cat("\nMAXENT$beta_hinge must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$betamultiplier)) {
             cat("\nMAXENT$betamultiplier must be a numeric")
             test <- FALSE
           }
           if (!is.numeric(object@MAXENT$defaultprevalence)) {
             cat("\nMAXENT$defaultprevalence must be a numeric")
             test <- FALSE
           }
           
           if(!is.null(object@MAXENT$initial_heap_size)){
             test <- .check_bytes_format(test,
                                         object@MAXENT$initial_heap_size,
                                         "initial_heap_size")
           }
           if(!is.null(object@MAXENT$max_heap_size)){
             test <- .check_bytes_format(test,
                                         object@MAXENT$max_heap_size,
                                         "max_heap_size")
           }
           ## MAXNET
           ## As of 2022/11/22 options check for MAXNET are missing 
           
           ## MAXNET (MAXENT.Tsuruoka) --> Obsolete
           ### TO BE DONE ===
           # 		       if(!is.numeric(object@MAXENT.Tsuruoka$l1_regularizer)){ cat("\nMAXENT.Tsuruoka$l1_regularizer must be a numeric"); test <- FALSE }
           # 		       if(!is.numeric(object@MAXENT.Tsuruoka$l2_regularizer)){ cat("\nMAXENT.Tsuruoka$l2_regularizer must be a numeric"); test <- FALSE }
           # 		       if(!is.logical(object@MAXENT.Tsuruoka$use_sgd)){ cat("\nMAXENT.Tsuruoka$use_sgd must be a logical"); test <- FALSE }
           # 		       if(!is.numeric(object@MAXENT.Tsuruoka$set_heldout)){ cat("\nMAXENT.Tsuruoka$set_heldout must be a numeric"); test <- FALSE }
           # 		       if(!is.logical(object@MAXENT.Tsuruoka$verbose)){ cat("\nMAXENT.Tsuruoka$verbose must be a logical"); test <- FALSE }
           
           return(test)
         }
)

### show.BIOMOD.models.options -------------------------------------------------
##'
##' @rdname BIOMOD.models.options
##' @importMethodsFrom methods show
##' @export
##'

setMethod('show', signature('BIOMOD.models.options'),
          function(object)
          {
            .bm_cat("BIOMOD.models.options")
            cat("\n")
            
            ## GLM options
            cat("\nGLM = list( type = '", object@GLM$type, "',", sep = "")
            cat("\n            interaction.level = ", object@GLM$interaction.level, ",", sep = "")
            cat("\n            myFormula = ",  ifelse(length(object@GLM$myFormula) < 1, 'NULL', paste(object@GLM$myFormula[2],
                                                                                                      object@GLM$myFormula[1],
                                                                                                      object@GLM$myFormula[3])), ",", sep = "")
            cat("\n            test = '", object@GLM$test, "',", sep = "")
            cat("\n            family = ", object@GLM$family$family, "(link = '", object@GLM$family$link, "'),", sep = "")
            cat("\n            mustart = ", object@GLM$mustart, ",", sep = "")
            cat("\n            control = glm.control(", .print_control(object@GLM$control), ") ),", sep = "", fill = .Options$width)
            
            ## GBM options
            cat("\n")
            cat("\nGBM = list( distribution = '", object@GBM$distribution, "',", sep = "")
            cat("\n            n.trees = ", object@GBM$n.trees, ",", sep = "")
            cat("\n            interaction.depth = ", object@GBM$interaction.depth, ",", sep = "")
            cat("\n            n.minobsinnode = ", object@GBM$n.minobsinnode, ",", sep = "")
            cat("\n            shrinkage = ", object@GBM$shrinkage, ",", sep = "")
            cat("\n            bag.fraction = ", object@GBM$bag.fraction, ",", sep = "")
            cat("\n            train.fraction = ", object@GBM$train.fraction, ",", sep = "")
            cat("\n            cv.folds = ", object@GBM$cv.folds, ",", sep = "")
            cat("\n            keep.data = ", object@GBM$keep.data, ",", sep = "")
            cat("\n            verbose = ", object@GBM$verbose, ",", sep = "")
            #             cat("\n            class.stratify.cv = '", object@GBM$class.stratify.cv, "',", sep="")
            cat("\n            perf.method = '", object@GBM$perf.method, "',", sep = "")
            cat("\n            n.cores = ", ifelse(length(object@GBM$n.cores), object@GBM$n.cores, 'NULL'), "),", sep = "")
            
            ## GAM options
            cat("\n")
            cat("\nGAM = list( algo = '", object@GAM$algo, "',", sep = "")
            cat("\n            type = '", object@GAM$type, "',", sep = "")
            cat("\n            k = ", ifelse(length(object@GAM$k) < 1, 'NULL', object@GAM$k), ",", sep = "")
            cat("\n            interaction.level = ", object@GAM$interaction.level, ",", sep = "")
            cat("\n            myFormula = ", ifelse(length(object@GAM$myFormula) < 1, 'NULL', paste(object@GAM$myFormula[2],
                                                                                                     object@GAM$myFormula[1],
                                                                                                     object@GAM$myFormula[3])), ",", sep = "")
            cat("\n            family = ", object@GAM$family$family, "(link = '", object@GAM$family$link, "'),", sep = "")
            if (object@GAM$algo == 'GAM_mgcv') {
              cat("\n            method = '", object@GAM$method, "', ", sep = "")
              cat("\n            optimizer = c('", paste(object@GAM$optimizer, collapse = "','"), "'),", sep = "")
              cat("\n            select = ", object@GAM$select, ",", sep = "")
              cat("\n            knots = ",  ifelse(length(object@GLM$knots) < 1, 'NULL', "'user.defined'"), ",", sep = "")
              cat("\n            paraPen = ",  ifelse(length(object@GLM$paraPen) < 1, 'NULL', "'user.defined'"), ",", sep = "")
            }
            cat("\n            control = list(", .print_control(object@GAM$control), ") ),", sep = "", fill = .Options$width)
            
            ## CTA options
            cat("\n")
            cat("\nCTA = list( method = '", object@CTA$method, "',", sep = "")
            cat("\n            parms = '", object@CTA$parms, "',", sep = "")
            cat("\n            cost = ", ifelse(length(object@CTA$cost) < 1, 'NULL', object@CTA$cost), ",", sep = "")
            cat("\n            control = list(", .print_control(object@CTA$control), ") ),", sep = "", fill = .Options$width)
            
            ## ANN options
            cat("\n")
            cat("\nANN = list( NbCV = ", object@ANN$NbCV, ",", sep = "")
            cat("\n            size = ", ifelse(length(object@ANN$size) < 1, 'NULL', object@ANN$size), ",", sep = "")
            cat("\n            decay = ", ifelse(length(object@ANN$decay) < 1, 'NULL', object@ANN$decay), ",", sep = "")
            cat("\n            rang = ", object@ANN$rang, ",", sep = "")
            cat("\n            maxit = ", object@ANN$maxit, "),", sep = "")
            
            ## SRE options
            cat("\n")
            cat("\nSRE = list( quant = ", object@SRE$quant, "),", sep = "")
            
            ## FDA options
            cat("\n")
            cat("\nFDA = list( method = '", object@FDA$method, "',", sep = "")
            cat("\n            add_args = ", ifelse(length(object@FDA$add_args) < 1, 'NULL'
                                                    , paste("list(", paste(.print_control(object@FDA$add_args), collapse = "")
                                                            , ")", sep = "")), "),", sep = "")
            
            ## MARS options
            cat("\n")
            cat("\nMARS = list( type = '", object@MARS$type, "',", sep = "")
            cat("\n             interaction.level = ", object@MARS$interaction.level, ",", sep = "")
            cat("\n             myFormula = ",  ifelse(length(object@MARS$myFormula) < 1, 'NULL', paste(object@GLM$myFormula[2],
                                                                                                        object@GLM$myFormula[1],
                                                                                                        object@GLM$myFormula[3])), ",", sep = "")
            #             cat("\n             degree = ", object@MARS$degree, ",", sep="")
            cat("\n             nk = ", ifelse(length(object@MARS$nk) < 1, 'NULL', object@MARS$nk), ",", sep = "")
            cat("\n             penalty = ", object@MARS$penalty, ",", sep = "")
            cat("\n             thresh = ", object@MARS$thresh, ",", sep = "")
            cat("\n             nprune = ", ifelse(length(object@MARS$nprune) < 1, 'NULL', object@MARS$nprune), ",", sep = "")
            cat("\n             pmethod = '", object@MARS$pmethod, "'),", sep = "")
            
            ## RF options
            cat("\n")
            cat("\nRF = list( do.classif = ", object@RF$do.classif, ",", sep = "")
            cat("\n           ntree = ", object@RF$ntree, ",", sep = "")
            cat("\n           mtry = '", object@RF$mtry, "',", sep = "")
            cat("\n           sampsize = ", ifelse(length(object@RF$sampsize) < 1, 'NULL', object@RF$sampsize), ",", sep = "")
            cat("\n           nodesize = ", object@RF$nodesize, ",", sep = "")
            cat("\n           maxnodes = ", ifelse(length(object@RF$maxnodes) < 1, 'NULL', object@RF$maxnodes),  "),", sep = "")
            
            ## MAXENT options
            cat("\n")
            cat("\nMAXENT = list( path_to_maxent.jar = '", object@MAXENT$path_to_maxent.jar, "', ", sep="")
            cat("\n               memory_allocated = ", ifelse(is.null(object@MAXENT$memory_allocated), 'NULL'
                                                               , object@MAXENT$memory_allocated), ",", sep = "")
            cat("\n               initial heap size = ", 
                ifelse(is.null(object@MAXENT$initial_heap_size), 
                       'NULL'
                       , object@MAXENT$initial_heap_size), ",", sep = "")
            cat("\n               maximum heap size = ", 
                ifelse(is.null(object@MAXENT$max_heap_size), 
                       'NULL', object@MAXENT$max_heap_size), ",", sep = "")
            cat("\n               background_data_dir = ", ifelse(is.character(object@MAXENT$background_data_dir), "'", "")
                , object@MAXENT$background_data_dir, ifelse(is.character(object@MAXENT$background_data_dir), "'", ""), ",", sep = "")
            cat("\n               maximumbackground = ", ifelse(is.character(object@MAXENT$maximumbackground), "'", "")
                , object@MAXENT$maximumbackground, ifelse(is.character(object@MAXENT$maximumbackground), "'", ""), ",", sep = "")
            cat("\n               maximumiterations = ", object@MAXENT$maximumiterations, ",", sep = "")
            cat("\n               visible = ", object@MAXENT$visible, ",", sep = "")
            cat("\n               linear = ", object@MAXENT$linear, ",", sep = "")
            cat("\n               quadratic = ", object@MAXENT$quadratic, ",", sep = "")
            cat("\n               product = ", object@MAXENT$product, ",", sep = "")
            cat("\n               threshold = ", object@MAXENT$threshold, ",", sep = "")
            cat("\n               hinge = ", object@MAXENT$hinge, ",", sep = "")
            cat("\n               lq2lqptthreshold = ", object@MAXENT$lq2lqptthreshold, ",", sep = "")
            cat("\n               l2lqthreshold = ", object@MAXENT$l2lqthreshold, ",", sep = "")
            cat("\n               hingethreshold = ", object@MAXENT$hingethreshold, ",", sep = "")
            cat("\n               beta_threshold = ", object@MAXENT$beta_threshold, ",", sep = "")
            cat("\n               beta_categorical = ", object@MAXENT$beta_categorical, ",", sep = "")
            cat("\n               beta_lqp = ", object@MAXENT$beta_lqp, ",", sep = "")
            cat("\n               beta_hinge = ", object@MAXENT$beta_hinge, ",", sep = "")
            cat("\n               betamultiplier = ", object@MAXENT$betamultiplier, ",", sep = "")
            cat("\n               defaultprevalence = ", object@MAXENT$defaultprevalence, "),", sep = "")
            
            ## MAXNET options
            cat("\n")
            cat("\n MAXNET = list( myFormula = ", .print_formula(object@MAXNET$myFormula), ",", sep = "")
            cat("\n     regmult = ", object@MAXNET$regmult, ",", sep = "")
            cat("\n     regfun = <function> )")
            cat("\n)")
            
            # ## MAXENT.Tsuruoka
            # cat("\n")
            # cat("\nMAXENT.Tsuruoka = list( l1_regularizer = ", object@MAXENT.Tsuruoka$l1_regularizer, ",", sep="")
            # cat("\n                        l2_regularizer = ", object@MAXENT.Tsuruoka$l2_regularizer, ",", sep="")
            # cat("\n                        use_sgd = ", object@MAXENT.Tsuruoka$use_sgd, ",", sep="")
            # cat("\n                        set_heldout = ", object@MAXENT.Tsuruoka$set_heldout, ",", sep="")
            # cat("\n                        verbose = ", object@MAXENT.Tsuruoka$verbose, ")", sep="")
            
            .bm_cat()
          }
)


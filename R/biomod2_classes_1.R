
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
##' \code{\link{bm_Tuning}}, \code{\link{bm_CrossValidation}} and 
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
##' @seealso \code{\link{BIOMOD_FormatingData}}, \code{\link{bm_Tuning}}, 
##' \code{\link{bm_CrossValidation}}, \code{\link{BIOMOD_Modeling}}, 
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

.BIOMOD.formated.data.check.args <- function(sp, env, xy = NULL, eval.sp = NULL, eval.env = NULL
                                             , eval.xy = NULL, filter.raster = FALSE)
{
  ## A.1 Check sp argument --------------------------------------------------------------
  if (inherits(sp, c('Raster','SpatRaster'))) {
    ## sp raster object not supported yet
    stop("Raster response variable not supported yet ! \nPlease extract your presences and your absences by yourself")
    #### TO DO #### extract the 0 and 1 in sp format
  }
  available.types.resp <- c('integer', 'numeric', 'data.frame', 'matrix',
                            'SpatialPointsDataFrame', 'SpatialPoints', 'SpatVector')
  .fun_testIfInherits(TRUE, "sp", sp, available.types.resp)
  
  ## SpatialPoints, SpatialPointsDataFrame, SpatVector
  if (inherits(sp, c('SpatialPoints','SpatVector'))) {
    .tmp <- .check_formating_spatial(resp.var = sp,
                                     expl.var = env,
                                     resp.xy = xy,
                                     eval.data = FALSE)
    sp <- .tmp$resp.var
    xy <- .tmp$resp.xy
    rm(.tmp)
  }
  
  ## data.frame, matrix  : transform into numeric
  if (inherits(sp, c("matrix", "data.frame"))) {
    sp <- .check_formating_table(sp)
  }
  
  ## Check presence/absence
  sp <- .check_formating_resp.var(resp.var = sp, eval.data = FALSE)
  
  ## A.2 Check xy argument --------------------------------------------------------------
  
  if (!is.null(xy)) {
    env.as.df <- inherits(env, c('data.frame'))
    xy <- .check_formating_xy(resp.xy = xy, resp.length = length(sp), env.as.df = env.as.df)
  } else if (inherits(env, c('RasterLayer', 'RasterStack', 'SpatRaster'))) {
    stop("`xy` argument is missing. Please provide `xy` when `env` is a raster.")
  } else {
    xy <- data.frame("x" = numeric(), "y" = numeric())
  }
  
  ## A.3 Check env argument -------------------------------------------------------------
  available.types.expl <- c('integer', 'numeric', 'data.frame', 'matrix',
                            'RasterLayer', 'RasterStack', 'SpatRaster',
                            'SpatialPointsDataFrame', 'SpatVector')
  .fun_testIfInherits(TRUE, "env", env, available.types.expl)
  env <- .check_formating_expl.var(expl.var = env, length.resp.var = length(sp))
  
  ## Check filter.raster argument
  if (inherits(env, "SpatRaster")) {
    stopifnot(is.logical(filter.raster))
  }
  
  #### At this point :
  ####  - sp is a numeric
  ####  - xy is NULL or a data.frame
  ####  - env is a data.frame or a SpatRaster
  
  
  ## DO THE SAME FOR EVALUATION DATA ####################################################
  if (is.null(eval.sp)) {
    cat("\n      ! No data has been set aside for modeling evaluation")
    evaL.env <- eval.xy <- NULL
  } else {
    ## B.1 Check eval.sp argument -------------------------------------------------------
    if (inherits(eval.sp, c('Raster', 'SpatRaster'))) {
      ## eval.sp raster object not supported yet
      stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
      #### TO DO #### extract the 0 and 1 in sp format
    }
    
    .fun_testIfInherits(TRUE, "eval.sp", eval.sp, available.types.resp)
    
    ## SpatialPoints, SpatialPointsDataFrame, SpatVector
    if (inherits(eval.sp, c('SpatialPoints','SpatVector'))) { 
      .tmp <- .check_formating_spatial(resp.var = eval.sp,
                                       expl.var = eval.env, 
                                       resp.xy = eval.xy,
                                       eval.data = TRUE)
      eval.sp <- .tmp$resp.var
      eval.xy <- .tmp$resp.xy
      rm(.tmp)
    }
    ## data.frame, matrix  : transform into numeric
    if (inherits(eval.sp, c("matrix","data.frame"))) {
      eval.sp <- .check_formating_table(eval.sp)
    }
    
    ## Check presence/absence
    eval.sp <- .check_formating_resp.var(resp.var = eval.sp, eval.data = TRUE)
    
    ## B.2 Check eval.xy argument -------------------------------------------------------
    if(!is.null(eval.xy)){
      eval.xy <- .check_formating_xy(resp.xy = eval.xy, resp.length = length(eval.sp))
    }
    
    ## B.3 Check eval.env argument ------------------------------------------------------
    .fun_testIfInherits(TRUE, "eval.env", eval.env, c(available.types.expl, "NULL"))
    if (is.null(eval.env)) {
      if (!(inherits(env, 'SpatRaster'))) {
        stop("If explanatory variable is not a raster and you want to consider evaluation response variable, you have to give evaluation explanatory variables")
      }
    }
    eval.env <- .check_formating_expl.var(expl.var = eval.env, length.resp.var = length(eval.sp))
    
    ## B.4 Remove NA from evaluation data
    if (sum(is.na(eval.sp)) > 0) {
      cat("\n      ! NAs have been automatically removed from Evaluation data")
      if (!is.null(eval.xy)) {
        eval.xy <- eval.xy[-which(is.na(eval.sp)), ]
      }
      eval.sp <- na.omit(eval.sp)
    }
  }
  
  return(list(sp = sp,
              env = env,
              xy = xy,
              eval.sp = eval.sp,
              eval.env = eval.env,
              eval.xy = eval.xy))
}


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
            args <- .BIOMOD.formated.data.check.args(sp, env, xy, eval.sp, eval.env, eval.xy, filter.raster)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            if (!any(sp == 0, na.rm = TRUE) && !any(is.na(sp))) {
              stop("No absences were given and no pseudo-absences were given or configured, at least one of those option is required.")
            }
            
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
            args <- .BIOMOD.formated.data.check.args(sp, env, xy, eval.sp, eval.env, eval.xy, filter.raster)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            if (!any(sp == 0, na.rm = TRUE) && !any(is.na(sp))) {
              stop("No absences were given and no pseudo-absences were given or configured, at least one of those option is required.")
            }
            
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
##' If \code{x} is a \code{\link{BIOMOD.formated.data.PA}} object, a \code{vector} 
##' containing pseudo-absence set to be represented 
##' @param run (\emph{optional, default} \code{'all'}) \cr 
##' If \code{calib.lines} provided, a \code{vector} containing repetition set to 
##' be represented 
##' @param plot.eval (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether evaluation data should be added to the plot or not
##' @param do.plot (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} defining whether the plot is to be rendered or not
##' @param point.size a \code{numeric} to adjust the size of points when
##'  \code{plot.type = 'points'}.
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
                   point.size = 1.5,
                   do.plot = TRUE)
          {
            args <- .plot.BIOMOD.formated.data.check.args(x = x,
                                                          calib.lines = calib.lines,
                                                          plot.type = plot.type,
                                                          plot.output = plot.output, 
                                                          PA = PA,
                                                          run = run,
                                                          plot.eval = plot.eval,
                                                          do.plot = do.plot)
            
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
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
              PA_run.vect <- foreach(this_PA = PA, this_run = run, .combine = 'rbind') %do%
                {
                  if (is.na(this_PA) || this_PA == 'allData') { # run only
                    this_name <- paste0("_", this_PA, "_", this_run)
                    this_calib <- calib.lines[ , this_name]
                    this_valid <- ! calib.lines[ , this_name]
                  } else if (is.na(this_run) || this_run == 'allRun') { # PA only
                    this_name <- this_PA
                    this_calib <- x@PA.table[ , this_PA]
                  } else { # PA+run
                    this_name <- paste0("_", this_PA, "_", this_run)
                    this_calib <- calib.lines[ , this_name] & x@PA.table[ , this_PA]
                    this_valid <- ! calib.lines[ , this_name] & x@PA.table[ , this_PA]
                  }
                  calib.resp <- x@data.species[which(this_calib)]
                  calib.resp <- ifelse(is.na(calib.resp), 30, 
                                       ifelse(calib.resp == 1, 10, 20))
                  calib.xy <- x@coord[which(this_calib),]
                  calib.df <- data.frame(resp = calib.resp,
                                         x = calib.xy[, 1],
                                         y = calib.xy[, 2])
                  
                  if (!is.na(this_run) & this_run != "allRun") { 
                    valid.resp <- x@data.species[which(this_valid)]
                    valid.resp <- ifelse(is.na(valid.resp), 31, 
                                         ifelse(valid.resp == 1, 11, 21))
                    valid.xy <- x@coord[which(this_valid),]
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
                             alpha = 1, size = point.size)+
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
                               alpha = 1, size = point.size)+
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
  
  ## 2 - check PA & run -----------------------------------------
  
  # find possible dataset
  allPA <- allrun <- NA
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun")
    if (inherits(x, "BIOMOD.formated.data.PA")) {
      expected_CVnames <- c(expected_CVnames
                            , sapply(1:ncol(x@PA.table)
                                     , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
                                                           , paste0("_PA", this_PA, "_allRun"))))
    } 
    .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
    
    allPA <- sapply(colnames(calib.lines), function(xx) strsplit(xx, "_")[[1]][2])
    allrun <- sapply(colnames(calib.lines), function(xx) strsplit(xx, "_")[[1]][3])
  } else if (inherits(x, "BIOMOD.formated.data.PA")) {
    allPA <- colnames(x@PA.table)
    allrun <- rep(NA, length(allPA))
  }
  
  # default value for PA and run
  if (missing(PA)) {
    PA <- allPA
  }
  if (missing(run)) {
    run <- allrun
  }
  
  # intersect possible and given dataset and check for PA and run values
  keep <- rep(TRUE, length(allPA))
  if (!is.null(calib.lines)) {
    .fun_testIfIn(TRUE, "run", run, allrun)
    keep[which(!allrun %in% run)] <- FALSE
  }
  if (inherits(x, "BIOMOD.formated.data.PA")) { # PA & CV
    .fun_testIfIn(TRUE, "PA", PA, allPA)
    keep[which(!allPA %in% PA)] <- FALSE
  }
  PA <- allPA[keep]
  run <- allrun[keep]
  
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
            
            if(inherits(object, "BIOMOD.formated.data.PA")){
              PA.length <- sapply(object@PA.table, function(this.pa){
                return(length(which(this.pa))
                       - length(which(object@data.species == 1)))
              })
              PA.unique <- unique(PA.length)
              PA.dataset <- sapply(PA.unique, function(this.PA){
                paste0(names(PA.length)[which(PA.length == this.PA)], collapse = ", ")
              })
              cat(
                "\n\n",
                ncol(object@PA.table),
                'Pseudo Absences dataset available (',
                paste0(colnames(object@PA.table), collapse = ", "),
                ") with ",
                paste0(sapply(seq_len(length(PA.unique)), function(i){
                  paste0(PA.unique[i], " (", PA.dataset[i],")")
                }), collapse = ", "),
                'pseudo absences',
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
##' \code{\link{bm_CrossValidation}} functions, to explore the distribution of calibration 
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
                                 "Presences" = length(which(object@data.species == 1)),
                                 "True_Absences" = length(which(object@data.species == 0)),
                                 "Pseudo_Absences" = 0,
                                 "Undefined" = length(which(is.na(object@data.species))))
            
            if (object@has.data.eval) {
              output <- rbind(output,
                              data.frame("dataset" = "evaluation",
                                         "run" = NA,
                                         "PA" = NA,
                                         "Presences" =  length(which(object@eval.data.species == 1)),
                                         "True_Absences" = length(which(object@eval.data.species == 0)),
                                         "Pseudo_Absences" = 0,
                                         "Undefined" = length(which(is.na(object@eval.data.species)))))
            }
            
            PA <- run <- NA
            if (!is.null(calib.lines)) {
              PA <- sapply(colnames(calib.lines), function(x) strsplit(x, "_")[[1]][2])
              run <- sapply(colnames(calib.lines), function(x) strsplit(x, "_")[[1]][3])
            } else if (inherits(object, "BIOMOD.formated.data.PA")) {
              PA <- colnames(object@PA.table)
              run <- rep(NA, length(PA))
            }
            if (!is.null(calib.lines) || inherits(object, "BIOMOD.formated.data.PA")) {
              output <- 
                rbind(
                  output,
                  foreach(this_run = run, this_PA = PA, .combine = 'rbind')  %do% {
                    if (is.na(this_PA) || this_PA == 'allData') { # run only
                        this_name <- paste0("_", this_PA, "_", this_run)
                        this_calib <- calib.lines[ , this_name]
                        this_valid <- ! calib.lines[ , this_name]
                      } else if (is.na(this_run)) { # PA only
                        this_calib <- ifelse(is.na(object@PA.table[ , this_PA]), FALSE, object@PA.table[ , this_PA])
                      } else { # PA+run
                        this_name <- paste0("_", this_PA, "_", this_run)
                        this_calib <- calib.lines[ , this_name] & object@PA.table[ , this_PA]
                        this_valid <- ! calib.lines[ , this_name] & object@PA.table[ , this_PA]
                      }
                      calib.resp <- object@data.species[which(this_calib)]
                      tmp <- data.frame("dataset" = "calibration",
                                        "run" = this_run,
                                        "PA" = this_PA,
                                        "Presences" = length(which(calib.resp == 1)),
                                        "True_Absences" = length(which(calib.resp == 0)),
                                        "Pseudo_Absences" = 
                                          length(which(is.na(calib.resp))),
                                        "Undefined" = NA)
                      
                      if (!is.na(this_run)) { 
                        valid.resp <- object@data.species[this_valid]
                        tmp <- rbind(tmp,
                                     data.frame("dataset" = "validation",
                                                "run" = this_run,
                                                "PA" = this_PA,
                                                "Presences" = length(which(valid.resp == 1)),
                                                "True_Absences" = length(which(valid.resp == 0)),
                                                "Pseudo_Absences" = 
                                                  length(valid.resp) - 
                                                  length(which(valid.resp == 1)) -
                                                  length(which(valid.resp == 0)),
                                                "Undefined" = NA))
                        
                      }
                      return(tmp) # end foreach
                    })
            } 
            output
          }
)

.summary.BIOMOD.formated.data.check.args <- function(object, calib.lines)
{
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun")
    if (inherits(object, "BIOMOD.formated.data.PA")) {
      if (nrow(calib.lines) != nrow(object@PA.table)) {
        stop("calib.lines and PA.table do not have the same number of rows")
      }
      expected_CVnames <- c(expected_CVnames
                            , sapply(1:ncol(object@PA.table)
                                     , function(this_PA) c(paste0("_PA", this_PA, "_RUN", seq_len(ncol(calib.lines)))
                                                           , paste0("_PA", this_PA, "_allRun"))))
    }
    .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
  }
  return(list(object = object, calib.lines = calib.lines))
}


## ---------------------------------------------------------------------#
## 2. BIOMOD.formated.data.PA ------------------------------------------------
## this class inherits from BIOMOD.formated.data and have one more slot 'PA', 
## giving PA selected
## -------------------------------------------------------------------- #

##' @name BIOMOD.formated.data.PA
##' @aliases BIOMOD.formated.data.PA-class
##' @author Damien Georges
##' 
##' @title \code{BIOMOD_FormatingData()} output object class (with pseudo-absences)
##' 
##' @description Class returned by \code{\link{BIOMOD_FormatingData}}, and used by 
##' \code{\link{bm_Tuning}}, \code{\link{bm_CrossValidation}} and 
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
##' or \code{PA.strategy = 'disk'}, an \code{integer} (or a \code{vector} of \code{integer} the 
##' same size as \code{PA.nb.rep}) corresponding to the number of pseudo-absence points that 
##' will be selected for each pseudo-absence repetition (true absences included)
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
##' \code{\link{bm_Tuning}}, \code{\link{bm_CrossValidation}}, 
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
##' # Keep only presence informations
##' DataSpecies <- DataSpecies[which(DataSpecies[, myRespName] == 1), ]
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
##'                                      PA.nb.rep = 4,
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
  args <- .BIOMOD.formated.data.check.args(sp, env, xy, eval.sp, eval.env, eval.xy, filter.raster)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  if (is.null(PA.strategy) || PA.strategy == 'none' || PA.nb.rep < 1) {
    if (!any(sp == 0, na.rm = TRUE) && !any(is.na(sp))) {
      stop("No absences were given and no pseudo-absences were given or configured, at least one of those option is required.")
    }
  }
  
  ### Check for categorical vars and filter duplicated data points
  categorical_var <- NULL
  if (inherits(env, 'SpatRaster')) {
    categorical_var <- names(env)[is.factor(env)] 
    
    output <- check_duplicated_cells(env, xy, sp, filter.raster, 
                                     PA.user.table = PA.user.table)
    xy <- output$xy
    sp <- output$sp
    PA.user.table <- output$PA.user.table
    rm(output)
  }
  
  ## Convert sp in SpatVector
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
    sp <- vect(sp.df, geom = c("x", "y"))
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


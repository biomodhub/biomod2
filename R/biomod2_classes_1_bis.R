
## --------------------------------------------------------------------------- #
## 1. BIOMOD.formated.data.abundance ---------------------------------------------------
## --------------------------------------------------------------------------- #

##' @name BIOMOD.formated.data.abundance
##' @aliases BIOMOD.formated.data-class
##' @author Blancheteau Hélène 
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
##' @rdname BIOMOD.formated.data.abundance
##' @importFrom terra rast app is.factor subset extract cellFromXY `add<-` 
##' classify rasterize values
##' @export
##' 

# 1.1 Class Definition ----------------------------------------------------------------------------
setClass("BIOMOD.formated.data.abundance",
         representation(data.type = 'character',
                        dir.name = 'character',
                        sp.name = 'character',
                        coord = "data.frame",
                        data.species = "numeric",
                        data.env.var = "data.frame",
                        data.mask = "list",
                        has.data.eval = "logical",
                        eval.coord = "data.frame",
                        eval.data.species = "numeric",
                        eval.data.env.var = "data.frame"),
         validity = function(object) { return(TRUE) })


# 1.2 Constructors -------------------------------------------------------------
setGeneric("BIOMOD.formated.data.abundance", def = function(sp, env, ...) { standardGeneric("BIOMOD.formated.data.abundance") })

.BIOMOD.formated.data.abundance.check.args <- function(sp, env, xy = NULL, eval.sp = NULL, eval.env = NULL
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


## BIOMOD.formated.data.abundance(sp = numeric, env = data.frame) -----------------------------------
##' 
##' @rdname BIOMOD.formated.data.abundance
##' @export
##' 

setMethod('BIOMOD.formated.data.abundance', signature(sp = 'numeric', env = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL
                   , eval.sp = NULL, eval.env = NULL, eval.xy = NULL
                   , na.rm = TRUE, data.mask = NULL, shared.eval.env = FALSE
                   , filter.raster = FALSE) 
          {
            args <- .BIOMOD.formated.data.abundance.check.args(sp, env, xy, eval.sp, eval.env, eval.xy, filter.raster)
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
                'BIOMOD.formated.data.abundance',
                coord = xy,
                data.species = sp,
                data.env.var = env,
                dir.name = dir.name,
                sp.name = sp.name,
                data.mask = data.mask,
                has.data.eval = FALSE,
                data.type = "abundance"
              )
            } else { ## EVALUATION DATA
              BFDeval <- BIOMOD.formated.data.abundance(
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
                'BIOMOD.formated.data.abundance',
                data.type = "abundance",
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

## BIOMOD.formated.data.abundance(sp = data.frame) -----------------------------------
##' 
##' @rdname BIOMOD.formated.data.abundance
##' @export
##' 

setMethod('BIOMOD.formated.data.abundance', signature(sp = 'data.frame'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL,
                   eval.sp = NULL, eval.env = NULL, eval.xy = NULL,
                   na.rm = TRUE, filter.raster = FALSE)
          {
            if (ncol(sp) > 1) { stop("Invalid response variable") }
            sp <- as.numeric(unlist(sp))
            BFD <- BIOMOD.formated.data.abundance(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

## BIOMOD.formated.data.abundance(sp = numeric, env = matrix) ----------------------------
##' 
##' @rdname BIOMOD.formated.data.abundance
##' @export
##' 

setMethod('BIOMOD.formated.data.abundance', signature(sp = 'numeric', env = 'matrix'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL,
                   eval.sp = NULL, eval.env = NULL, eval.xy = NULL,
                   na.rm = TRUE, filter.raster = FALSE)
          {
            env <- as.data.frame(env)
            BFD <- BIOMOD.formated.data.abundance(sp, env, xy, dir.name, sp.name, eval.sp, eval.env, eval.xy, na.rm = na.rm)
            return(BFD)
          }
)

## BIOMOD.formated.data.abundance(sp = numeric, env = SpatRaster) -----------------------
##' 
##' @rdname BIOMOD.formated.data.abundance
##' @export
##' 

setMethod('BIOMOD.formated.data.abundance', signature(sp = 'numeric', env = 'SpatRaster'),
          function(sp, env, xy = NULL, dir.name = '.', sp.name = NULL,
                   eval.sp = NULL, eval.env = NULL, eval.xy = NULL,
                   na.rm = TRUE, shared.eval.env = FALSE,
                   filter.raster = FALSE)
          {
            args <- .BIOMOD.formated.data.abundance.check.args(sp, env, xy, eval.sp, eval.env, eval.xy, filter.raster)
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
            
            BFD <- BIOMOD.formated.data.abundance(sp, env, xy, dir.name, sp.name, 
                                        eval.sp, eval.env, eval.xy,
                                        na.rm = na.rm, data.mask = data.mask,
                                        shared.eval.env = shared.eval.env,
                                        filter.raster = filter.raster)
            return(BFD)
          }
)


# 1.3 Other Functions -----------------------------------------------------------------------------

### plot.BIOMOD.formated.data.abundance (doc) --------------------------------------------------
##' 
##' @rdname plot
##' @docType methods
##' @author Helene Blancheteau
##' @title \code{plot} method for \code{\link{BIOMOD.formated.data.abundance}} object class
##' 
##' @description Plot the spatial distribution of presences, absences and 
##' pseudo-absences among the different potential dataset (calibration, 
##' validation and evaluation). Available only if coordinates were given to 
##' \code{\link{BIOMOD_FormatingData}}.
##' 
##' 
##' @param x a \code{\link{BIOMOD.formated.data.abundance}} or \code{\link{BIOMOD.formated.data.abundance.PA}}
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
##' If \code{x} is a \code{\link{BIOMOD.formated.data.abundance.PA}} object, a \code{vector} 
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


setMethod('plot', signature(x = 'BIOMOD.formated.data.abundance', y = "missing"),
          function(x,
                   calib.lines = NULL,
                   plot.type,
                   plot.output, 
                   run,
                   plot.eval,
                   point.size = 1.5,
                   do.plot = TRUE)
          {
            args <- .plot.BIOMOD.formated.data.abundance.check.args(x = x,
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
            filter <- !is.na(x@data.species)
            full.resp <- x@data.species[filter]
            full.xy <- x@coord[filter,]
            full.df <- data.frame(resp = full.resp,
                                  x = full.xy[,1],
                                  y = full.xy[,2])
            
            full.df.vect <- vect(full.df, geom = c("x","y"))
            names(full.df.vect) <- "resp"
            full.df.vect$dataset <- "Initial dataset"
            
            ## 1.2 - Eval Dataset -----------------------------------------------
            
            if (plot.eval) {
              eval.resp <- x@eval.data.species
              eval.xy <- x@eval.coord
              eval.df <- data.frame(resp = eval.resp,
                                    x = eval.xy[,1],
                                    y = eval.xy[,2])
              
              eval.df.vect <- vect(eval.df, geom = c("x","y"))
              names(eval.df.vect) <- "resp"
              eval.df.vect$dataset <- "Evaluation dataset"
              full.df.vect <- rbind(full.df.vect, eval.df.vect)
            }
            
            ## 1.3 - CV dataset -----------------------------
            if (!is.null(calib.lines)) {
              PA_run.vect <- foreach(this_PA = PA, this_run = run, .combine = 'rbind') %do%
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
                                         y = calib.xy[, 2])
                  
                  if (!is.na(this_run) & this_run != "allRun") { 
                    valid.resp <- x@data.species[which(this_valid)]
                    # valid.resp <- ifelse(is.na(valid.resp), 31, 
                    #                      ifelse(valid.resp == 1, 11, 21))
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
            data_breaks <- c("Initial dataset", "Evaluation dataset", colnames(calib.lines))              
            data_labels <- data_breaks
            data_labels_facet <- c("all data","all data", paste0("RUN", 1:length(colnames(calib.lines)))) # background
            
            data_colors <- c("Initial dataset" = "#004488",
                             "Evaluation dataset" = "#994455",
                             "_allData_RUN1" = "#997700",
                             "_allData_RUN2" = "#997700",
                             "_allData_RUN3" = "#997700")
            
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
                               color = resp,
                               shape = factor(dataset, levels = data_breaks)), 
                           alpha = 1, size = point.size)+
                facet_wrap(~dataset)+
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
                                 color = resp,
                                 shape = as.factor(dataset)),
                             alpha = 1, size = point.size)+
                  # scale_color_manual(
                  #   NULL,
                  #   breaks = data_breaks,
                  #   values = data_colors,
                  #   labels = data_labels)+
                  # scale_shape_manual(
                  #   NULL,
                  #   breaks = data_breaks,
                  #   values = data_shape,
                  #   labels = data_labels)+
                  # scale_fill_manual(
                  #   guide = "none",
                  #   breaks = data_breaks,
                  #   values = data_colors,
                  #   labels = data_labels,
                  #   na.value = data_background)+
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
)


.plot.BIOMOD.formated.data.abundance.check.args <- function(x,
                                                            calib.lines,
                                                            plot.type,
                                                            plot.output, 
                                                            PA,
                                                            run,
                                                            plot.eval,
                                                            do.plot){
  
  
  ## 1 - check x -----------------------------------------
  .fun_testIfInherits(TRUE, "x", x, c("BIOMOD.formated.data.abundance"))
  
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

setMethod('show', signature('BIOMOD.formated.data.abundance'),
          function(object)
          {
            .bm_cat("BIOMOD.formated.data")
            cat("\ndir.name = ", object@dir.name, fill = .Options$width)
            cat("\nsp.name = ", object@sp.name, fill = .Options$width)
            cat("\n\t",
                length(object@data.species),
                'points, with a range from',
                min(object@data.species, na.rm = TRUE),
                'to',
                max(object@data.species, na.rm = TRUE),
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
            .bm_cat()
          }
)


### summary.BIOMOD.formated.data.abundance  --------------------------------------------------
##' 
##' @rdname summary
##' @docType methods
##' @author Helene Blancheteau
##' 
##' @title \code{summary} method for \code{\link{BIOMOD.formated.data}} object class
##' 
##' @description Summarize the number of presences, absences and 
##' pseudo-absences among the different potential dataset (calibration, 
##' validation and evaluation).  
##' 
##' @param object a \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.abundance}} 
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

setMethod('summary', signature(object = 'BIOMOD.formated.data.abundance'),
          function(object, calib.lines = NULL)
          {
            args <- .summary.BIOMOD.formated.data.abundance.check.args(object = object, calib.lines = calib.lines)
            for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
            rm(args)
            
            output <- data.frame("dataset" = "initial",
                                 "run" = NA,
                                 "Points" = length(object@data.species),
                                 "Min" = min(object@data.species, na.rm = T),
                                 "Max" = max(object@data.species, na.rm = T),
                                 "Transformation"  = "Untouched")
            
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
            
            run <- NA
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

.summary.BIOMOD.formated.data.abundance.check.args <- function(object, calib.lines)
{
  if (!is.null(calib.lines)) {
    .fun_testIfInherits(TRUE, "calib.lines", calib.lines, c("matrix"))
    
    expected_CVnames <- c(paste0("_allData_RUN", seq_len(ncol(calib.lines))), "_allData_allRun")
    .fun_testIfIn(TRUE, "colnames(calib.lines)", colnames(calib.lines), expected_CVnames)
  }
  return(list(object = object, calib.lines = calib.lines))
}



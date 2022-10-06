###################################################################################################
##' @name BIOMOD_FormatingData
##' @author Damien Georges, Wilfried Thuiller
##' 
##' @title Format input data, and select pseudo-absences if wanted, for usage in \pkg{biomod2}
##' 
##' @description This function gathers together all input data needed (\emph{xy, 
##' presences/absences, explanatory variables, and the same for evaluation data if available}) to 
##' run \pkg{biomod2} models. It allows to select pseudo-absences if no absence data is available, 
##' with different strategies (see Details).
##' 
##' 
##' @param dir.name (\emph{optional, default} \code{.}) \cr
##' A \code{character} corresponding to the modeling folder
##' @param resp.name a \code{character} corresponding to the species name
##' 
##' @param resp.var a \code{vector}, \code{\link[sp]{SpatialPoints}} (\emph{if presence-only}) or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' build the species distribution model(s)
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s)
##' @param resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to build the 
##' species distribution model(s)
##' 
##' @param eval.resp.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector}, \code{\link[sp]{SpatialPoints}} (\emph{if presence-only}) or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species that will be used to 
##' evaluate the species distribution model(s) with independent data
##' @param eval.expl.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} or 
##' \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables (in 
##' columns or layers) that will be used to evaluate the species distribution model(s) with 
##' independent data
##' @param eval.resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be used to evaluate 
##' the species distribution model(s) with independent data
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
##' A \code{logical} value defining whether points having one or several missing values for 
##' explanatory variables should be removed from the analysis or not
##' 
##' 
##' @return 
##' 
##' A \code{BIOMOD.formated.data} object that can be used to build species distribution model(s) 
##' with the \code{\link{BIOMOD_Modeling}} function. \cr \code{print} and \code{plot} 
##' functions are available to have a summary of the created object.
##' 
##' 
##' @details  
##' 
##' This function gathers and formats all input data needed to run \pkg{biomod2} models. It 
##' supports different kind of inputs (e.g. \code{matrix}, \code{\link[sp]{SpatialPoints}}, 
##' \code{\link[raster:stack]{RasterStack}}) and provides different methods to select 
##' pseudo-absences if needed. \cr \cr
##' 
##' \bold{Concerning explanatory variables and XY coordinates :} 
##' \itemize{
##'   \item if \code{\link[raster:stack]{RasterLayer}} or \code{\link[raster:stack]{RasterStack}} 
##'   provided (for \code{expl.var} or \code{eval.expl.var}), \pkg{biomod2} will extract the 
##'   corresponding values from XY coordinates provided (\code{resp.xy} or \code{eval.resp.xy} 
##'   respectively). \cr \emph{Be sure to give the XY coordinates in the same projection system 
##'   than the raster objects !}
##'   \item if \code{\link[sp]{SpatialPointsDataFrame}} provided (for \code{resp.var} or 
##'   \code{eval.resp.var}), the same rule applies (\emph{same projection system !}).
##'   \item if \code{data.frame} provided (for \code{expl.var} or \code{eval.expl.var}), 
##'   \pkg{biomod2} will simply merge it (\code{cbind}) with \code{resp.var} without 
##'   considering XY coordinates. \cr
##'   \emph{Be sure to give explanatory and response values in the same row order !}
##' }
##' 
##' \bold{Concerning pseudo-absence selection (see \code{\link{bm_PseudoAbsences}}) :}
##' \itemize{
##'   \item{if both presence and absence data are available, and there is enough absences :
##'   set \code{PA.nb.rep = 0} and no pseudo-absence will be selected.
##'   }
##'   \item{if no absence data (or not enough) is available, several pseudo-absence repetitions  
##'   are recommended (to estimate the effect of pseudo-absence selection), as well as high 
##'   number of pseudo-absence points. \cr
##'   \emph{Be sure not to select more pseudo-absence points than maximum number of pixels in 
##'   the studied area ! \cr \cr \cr \cr}
##'   }
##' }
##' 
##' \describe{
##'   \item{Response variable}{
##'   \pkg{biomod2} models single species at a time (no multi-species). Hence, \code{resp.var} 
##'   must be a uni-dimensional object (either a \code{vector}, a one-column \code{matrix}, 
##'   \code{data.frame} or \code{\link[sp]{SpatialPointsDataFrame}}, or a 
##'   \code{\link[sp]{SpatialPoints}} object if presence-only), containing values among :
##'   \itemize{
##'     \item \code{1} : presences
##'     \item \code{0} : true absences (if any)
##'     \item \code{NA} : no information point (might be used to select pseudo-absences if any)
##'   }
##'   If no true absences are available, pseudo-absence selection must be done. \cr
##'   If \code{resp.var} is a non-spatial object (\code{vector}, \code{matrix} or 
##'   \code{data.frame}), XY coordinates must be provided through \code{resp.xy}. \cr 
##'   If pseudo-absence points are to be selected, \code{NA} points must be provided in order to 
##'   select pseudo-absences among them.
##'   }
##'   \item{Explanatory variables}{
##'   Factorial variables are allowed, but might lead to some pseudo-absence strategy or models 
##'   omissions (e.g. \code{sre}).
##'   }
##'   \item{Evaluation data}{
##'   Although \pkg{biomod2} provides tools to automatically divide dataset into calibration and 
##'   validation parts through the modeling process (see \code{NbRunEval} and \code{DataSplit} 
##'   parameters in \code{\link{BIOMOD_Modeling}} function ; or \code{\link{BIOMOD_CrossValidation} 
##'   function}), it is also possible (and strongly advised) to directly provide two independent 
##'   datasets, one for calibration and one for validation.
##'   }
##'   \item{Pseudo-absence selection}{
##'   If no true absences are available, pseudo-absences must be selected from the 
##'   \emph{background data}, meaning data there is no information whether the species of 
##'   interest occurs or not. It corresponds either to the remaining pixels of the \code{expl.var} 
##'   (if provided as a \code{\link[raster:stack]{RasterStack}}) or to the points identified as 
##'   \code{NA} in \code{resp.var} (if \code{expl.var} provided as a \code{matrix} or 
##'   \code{data.frame}). \cr 
##'   Several methods are available to do this selection :
##'   \describe{
##'     \item{random}{all points of initial background are pseudo-absence candidates. 
##'     \code{PA.nb.absences} are drawn randomly, for each \code{PA.nb.rep} requested.
##'     }
##'     \item{sre}{pseudo-absences have to be selected in conditions (combination of explanatory 
##'     variables) that differ in a defined proportion (\code{PA.sre.quant}) from those of 
##'     presence points. A \emph{Surface Range Envelop} model is first run over the species of 
##'     interest, and pseudo-absences are selected outside this envelop. \cr
##'     \emph{This case is appropriate when all the species climatic niche has been sampled, 
##'     otherwise it may lead to over-optimistic model evaluations and predictions !}
##'     }
##'     \item{disk}{pseudo-absences are selected within circles around presence points defined by 
##'     \code{PA.dist.min} and \code{PA.dist.max} distance values (in meters). It allows to select 
##'     pseudo-absence points that are not too close to (avoid same niche and pseudo-replication) 
##'     or too far (localized sampling strategy) from presences.
##'     }
##'     \item{user.defined}{pseudo-absences are defined in advance and given as \code{data.frame} 
##'     through the \code{PA.user.table} parameter.
##'     }
##'   }
##'   }
##' }
##' 
##' 
##' @keywords dataset format evaluation pseudo-absence
##' 
##' 
##' @seealso \code{\link{bm_PseudoAbsences}}, \code{\link{BIOMOD_Modeling}}
##' @family Main functions
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
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
##'                                      expl.var = myExpl,
##'                                      resp.xy = myRespXY,
##'                                      resp.name = myRespName)
##' myBiomodData
##' plot(myBiomodData)
##' 
##' 
##' # ---------------------------------------------------------------
##' # # Transform true absences into potential pseudo-absences
##' # myResp.PA <- ifelse(myResp == 1, 1, NA)
##' # 
##' # # Format Data with pseudo-absences : random method
##' # myBiomodData.r <- BIOMOD_FormatingData(resp.var = myResp.PA,
##' #                                        expl.var = myExpl,
##' #                                        resp.xy = myRespXY,
##' #                                        resp.name = myRespName,
##' #                                        PA.nb.rep = 4,
##' #                                        PA.nb.absences = 1000,
##' #                                        PA.strategy = 'random')
##' # 
##' # # Format Data with pseudo-absences : disk method
##' # myBiomodData.d <- BIOMOD_FormatingData(resp.var = myResp.PA,
##' #                                        expl.var = myExpl,
##' #                                        resp.xy = myRespXY,
##' #                                        resp.name = myRespName,
##' #                                        PA.nb.rep = 4,
##' #                                        PA.nb.absences = 500,
##' #                                        PA.strategy = 'disk',
##' #                                        PA.dist.min = 5,
##' #                                        PA.dist.max = 35)
##' # 
##' # # Format Data with pseudo-absences : SRE method
##' # myBiomodData.s <- BIOMOD_FormatingData(resp.var = myResp.PA,
##' #                                        expl.var = myExpl,
##' #                                        resp.xy = myRespXY,
##' #                                        resp.name = myRespName,
##' #                                        PA.nb.rep = 4,
##' #                                        PA.nb.absences = 1000,
##' #                                        PA.strategy = 'sre',
##' #                                        PA.sre.quant = 0.025)
##' # 
##' # # Format Data with pseudo-absences : user.defined method
##' # myPAtable <- data.frame(PA1 = ifelse(myResp == 1, TRUE, FALSE),
##' #                         PA2 = ifelse(myResp == 1, TRUE, FALSE))
##' # for (i in 1:ncol(myPAtable)) myPAtable[sample(which(myPAtable[, i] == FALSE), 500), i] = TRUE
##' # myBiomodData.u <- BIOMOD_FormatingData(resp.var = myResp.PA,
##' #                                        expl.var = myExpl,
##' #                                        resp.xy = myRespXY,
##' #                                        resp.name = myRespName,
##' #                                        PA.strategy = 'user.defined',
##' #                                        PA.user.table = myPAtable)
##' # 
##' # myBiomodData.r
##' # myBiomodData.d
##' # myBiomodData.s
##' # myBiomodData.u
##' # plot(myBiomodData.r)
##' # plot(myBiomodData.d)
##' # plot(myBiomodData.s)
##' # plot(myBiomodData.u)
##' 
##' 
##' @importFrom sp coordinates
##' @importFrom raster stack
##' 
##' @export
##' 
###################################################################################################

BIOMOD_FormatingData <- function(resp.name,
                                 resp.var,
                                 expl.var,
                                 dir.name = '.',
                                 resp.xy = NULL,
                                 eval.resp.var = NULL,
                                 eval.expl.var = NULL,
                                 eval.resp.xy = NULL,
                                 PA.nb.rep = 0,
                                 PA.nb.absences = 1000,
                                 PA.strategy = 'random',
                                 PA.dist.min = 0,
                                 PA.dist.max = NULL,
                                 PA.sre.quant = 0.025,
                                 PA.user.table = NULL,
                                 na.rm = TRUE)
{
  .bm_cat(paste0(resp.name, " Data Formating"))
  
  ## 1. check args ------------------------------------------------------------
  args <- .BIOMOD_FormatingData.check.args(resp.name,
                                           resp.var,
                                           expl.var,
                                           dir.name,
                                           resp.xy,
                                           eval.resp.var,
                                           eval.expl.var,
                                           eval.resp.xy,
                                           PA.nb.rep,
                                           PA.nb.absences,
                                           PA.strategy,
                                           PA.dist.min,
                                           PA.dist.max,
                                           PA.sre.quant,
                                           PA.user.table)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 2. build BIOMOD.formated.data object -------------------------------------
  out <- NULL
  if( PA.strategy == 'none') { # no Pseudo Absences
    out <- BIOMOD.formated.data(sp = resp.var,
                                xy = resp.xy,
                                env = expl.var,
                                dir.name = dir.name,
                                sp.name = resp.name,
                                eval.sp = eval.resp.var,
                                eval.env = eval.expl.var,
                                eval.xy = eval.resp.xy,
                                na.rm = na.rm)
  } else { # automatic Pseudo Absences selection
    out <- BIOMOD.formated.data.PA(sp = resp.var,
                                   xy = resp.xy,
                                   env = expl.var,
                                   dir.name = dir.name,
                                   sp.name = resp.name,
                                   eval.sp = eval.resp.var,
                                   eval.env = eval.expl.var,
                                   eval.xy = eval.resp.xy,
                                   PA.nb.rep = PA.nb.rep,
                                   PA.strategy = PA.strategy,
                                   PA.nb.absences = PA.nb.absences,
                                   PA.dist.min = PA.dist.min,
                                   PA.dist.max = PA.dist.max,
                                   PA.sre.quant = PA.sre.quant,
                                   PA.user.table = PA.user.table,
                                   na.rm = na.rm)
  }
  .bm_cat("Done")
  return(out)
}


###################################################################################################

.BIOMOD_FormatingData.check.args <- function(resp.name,
                                             resp.var,
                                             expl.var,
                                             dir.name,
                                             resp.xy,
                                             eval.resp.var,
                                             eval.expl.var,
                                             eval.resp.xy,
                                             PA.nb.rep,
                                             PA.nb.absences,
                                             PA.strategy,
                                             PA.dist.min,
                                             PA.dist.max,
                                             PA.sre.quant,
                                             PA.user.table)
{
  ## 0. Checking names (resp.name available ?) --------------------------------
  if (grepl('/', resp.name)) {
    stop(paste0("Response variable name must be a character, and not a file path."
                , "\n Please refer to dir.name parameter to set a modeling folder."))
  }
  if (grepl('_', resp.name) | grepl(' ', resp.name)) {
    resp.name <- paste(unlist(strsplit(resp.name, '_')), collapse = '.')
    resp.name <- paste(unlist(strsplit(resp.name, ' ')), collapse = '.')
    cat('\n      ! Response variable name was converted into', resp.name)
  }
  if (!dir.exists(dir.name)) {
    stop(paste0("Modeling folder '", dir.name, "' does not exist"))
  }
  
  available.types <- c('integer', 'numeric', 'data.frame', 'matrix',
                       'RasterLayer', 'RasterStack',
                       'SpatialPointsDataFrame', 'SpatialPoints')
  
  ## 1. Checking resp.var -----------------------------------------------------
  
  if (!inherits(resp.var, available.types)) { ## resp.var
    stop(paste0("Response variable must be one of ", toString(available.types)))
  }
  
  if (inherits(resp.var, 'Raster')) { ## resp.var raster object not supported yet
    stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
    #### TO DO #### extract the 0 and 1 in sp format
  }
  
  if (!inherits(expl.var, setdiff(available.types, 'SpatialPoints'))) { ## expl.var
    stop(paste0("Explanatory variable must be one of ", toString(available.types)))
  }
  
  if (inherits(resp.var, 'SpatialPoints')) { ## resp.xy
    if (!is.null(resp.xy)) {
      cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
    }
    resp.xy <- data.matrix(coordinates(resp.var))
    if (inherits(resp.var, 'SpatialPointsDataFrame')) {
      resp.var <- resp.var@data
    } else {
      cat("\n      ! Response variable is considered as only presences... Is it really what you want?")
      resp.var <- rep(1, nrow(resp.xy))
    }
  }
  
  ## transforming into numeric if data.frame or matrix
  if (is.matrix(resp.var) | is.data.frame(resp.var)) {
    resp.var = as.data.frame(resp.var)
    if (ncol(resp.var) > 1) {
      stop("You must give a monospecific response variable (1D object)")
    } else {
      resp.var <- as.numeric(resp.var[, 1])
    }
  }
  
  if(is.matrix(expl.var) | is.numeric(expl.var)) {
    expl.var <- as.data.frame(expl.var)
  }
  
  if (inherits(expl.var, 'Raster')) {
    expl.var <- stack(expl.var, RAT = FALSE)
  }
  
  if (inherits(expl.var, 'SpatialPoints')) {
    expl.var <- as.data.frame(expl.var@data)
  }
  
  ## checking xy coordinates validity
  if(!is.null(resp.xy)){
    if(ncol(resp.xy)!=2){
      stop("If given, resp.xy must be a 2 column matrix or data.frame")
    }
    if(nrow(resp.xy) != length(resp.var)){
      stop("Response variable and its coordinates don't match")
    }
    resp.xy <- as.data.frame(resp.xy)
  }
  
  ## converting response var into binary
  resp.var[which(resp.var > 0)] <- 1
  resp.var[which(resp.var <= 0)] <- 0
  
  #### At this point :
  ####  - resp.var is a numeric
  ####  - resp.xy is NULL or a data.frame
  ####  - expl.var is a data.frame or a RasterStack
  ####  - sp.name is a character
  
  ## checking resp and expl var compatibility
  if (is.data.frame(expl.var) && nrow(expl.var) != length(resp.var)) {
    stop("If explanatory variable is not a raster then dimensions of response variable and explanatory variable must match!")
  }
  
  ## PA strategy
  if (is.null(PA.user.table) & PA.nb.rep < 1) {
    cat("\n> No pseudo absences selection !")
    PA.strategy <- "none"
    PA.nb.rep <- 0
  }
  
  if (is.null(PA.strategy) &  PA.nb.rep > 0) {
    cat("\n> Pseudo absences will be selected randomly !")
    PA.strategy <- "random"
  }
  
  if (!is.null(PA.user.table)) {
    cat("\n> Pseudo absences used will be user defined ones !")
    PA.strategy <- "user.defined"
    PA.nb.rep <- 0
  }
  
  if (PA.strategy == "user.defined") {
    if (!(is.matrix(PA.user.table) | is.data.frame(PA.user.table)))
      stop("\n PA.user.table must be a matrix or a data.frame")
    
    if (nrow(PA.user.table) != length(resp.var))
      stop("\n PA.user.table must have as many row than the number of observation of your response variable")
    
    colnames(PA.user.table) <- paste0("PA", 1:ncol(PA.user.table))
  }
  
  ## 2. Checking eval.resp.var ------------------------------------------------
  
  if (!is.null(eval.resp.var)) {
    if (!(class(eval.resp.var)[1] %in% available.types)) { ## eval.respvar
      stop(paste0("Response variable must be one of ", toString(available.types)))
    }
    
    if (inherits(eval.resp.var, 'Raster')) { ## eval.resp.var raster object not supported yet
      stop("Raster response variable not supported yet ! \nPlease extract your Presences and your absences by yourself")
      #### TO DO #### extract the 0 and 1 in sp format
    }
    
    if (!is.null(eval.expl.var)) { ## eval.expl.var
      if (!(class(eval.expl.var)[1] %in% available.types[-which(available.types == 'SpatialPoints')])) {
        stop(paste0("Explanatory variable must be one of ", toString(available.types)))
      }
    } else {
      if (!(inherits(expl.var, 'Raster'))) {
        stop("If explanatory variable is not a raster and you want to consider evaluation response variable, you have to give evaluation explanatory variables")
      }
    }
    
    if (inherits(eval.resp.var, 'SpatialPoints')) { ## eval.resp.xy
      if (!is.null(eval.resp.xy)) {
        cat("\n      ! XY coordinates of response variable will be ignored because spatial response object is given.")
      }
      eval.resp.xy <- data.matrix(coordinates(eval.resp.var))
      if (class(eval.resp.var)[1] == 'SpatialPointsDataFrame') {
        eval.resp.var <- eval.resp.var@data
      } else {
        cat("\n      ! Response variable is considered as only presences... Is it really what you want?")
        eval.resp.var <- rep(1, nrow(eval.resp.xy))
      }
    }
    
    ### transforming into numeric if data.frame or matrix
    if (is.matrix(eval.resp.var) | is.data.frame(eval.resp.var)) {
      if (ncol(eval.resp.var) > 1) {
        stop("You must give a monospecific response variable (1D object)")
      } else {
        eval.resp.var <- as.numeric(eval.resp.var[, 1])
      }
    }
    
    if (is.matrix(eval.expl.var) | is.numeric(eval.expl.var)) {
      eval.expl.var <- as.data.frame(eval.expl.var)
    }
    
    if (inherits(eval.expl.var, 'Raster')) {
      eval.expl.var <- stack(eval.expl.var)
    }
    
    if (inherits(eval.expl.var, 'SpatialPoints')) {
      eval.expl.var <- as.data.frame(eval.expl.var@data)
    }
    
    ## checking xy coordinates validity
    if (!is.null(eval.resp.xy)) {
      if (ncol(eval.resp.xy) != 2) {
        stop("if given, resp.xy must be a 2 column matrix or data.frame")
      }
      if (nrow(eval.resp.xy) != length(eval.resp.var)) {
        stop("Response variable and its coordinates don't match")
      }
      eval.resp.xy <- as.data.frame(eval.resp.xy)
    }
    
    if (is.data.frame(eval.expl.var) && nrow(eval.expl.var) != length(eval.resp.var)) {
      stop("If explanatory variable is not a raster then dimensions of response variable and explanatory variable must match!")
    }
    
    ## remove NA from evaluation data
    if (sum(is.na(eval.resp.var)) > 0) {
      cat("\n      ! NAs have been automatically removed from Evaluation data")
      if (!is.null(eval.resp.xy)) {
        eval.resp.xy <- eval.resp.xy[-which(is.na(eval.resp.var)), ]
      }
      eval.resp.var <- na.omit(eval.resp.var)
    }
    
    ## converting response var into binary
    eval.resp.var[which(eval.resp.var > 0)] <- 1
    eval.resp.var[which(eval.resp.var <= 0)] <- 0
    
    ## check there are both presences and absences in evaluation dataset
    if (sum(eval.resp.var == 1) < 1 | sum(eval.resp.var == 0) < 1) {
      stop("Evaluation response data must have both presences and absences")
    }
  } else {
    cat("\n      ! No data has been set aside for modeling evaluation")
    eval.expl.var <- eval.resp.xy <- NULL
  }
  
  ### PA arguments are not checked here because it will be done later... (may be will do it here later)
  
  return(list(resp.var = resp.var,
              expl.var = expl.var,
              resp.xy = resp.xy,
              resp.name = resp.name,
              dir.name = dir.name,
              eval.resp.var = eval.resp.var,
              eval.expl.var = eval.expl.var,
              eval.resp.xy = eval.resp.xy,
              PA.nb.rep = PA.nb.rep,
              PA.nb.absences = PA.nb.absences,
              PA.strategy = PA.strategy,
              PA.dist.min = PA.dist.min,
              PA.dist.max = PA.dist.max,
              PA.sre.quant = PA.sre.quant,
              PA.user.table = PA.user.table))
}

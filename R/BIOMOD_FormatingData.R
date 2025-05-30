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
##' @param resp.name a \code{character} corresponding to the species name
##' @param resp.var a \code{vector}, a \code{\link[terra:vect]{SpatVector}} without associated 
##' data (\emph{if presence-only}), or a \code{\link[terra:vect]{SpatVector}} object containing 
##' binary data (\code{1} : presence, \code{0} : absence, \code{NA} : indeterminate) or other 
##' data (see \code{data.type} and Details) for a single species that will be used to build the 
##' species distribution model(s)
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##'  \code{SpatialPoints} (if presence-only) or \code{SpatialPointsDataFrame}
##'  object containing binary data or other data.}
##' @param resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be
##' used to build the species distribution model(s)
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables 
##' (in columns or layers) that will be used to build the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }

##' @param dir.name (\emph{optional, default} \code{.}) \cr
##' A \code{character} corresponding to the modeling folder
##' @param data.type a \code{character}, corresponding to the response data type to be used, 
##' must be either \code{binary}, \code{count}, \code{multiclass}, \code{ordinal}, \code{relative}, or 
##' \code{abundance}, and match the data contained in \code{resp.var} 
##' \cr \emph{If not provided, \code{biomod2} will try to guess.}
##' 
##' @param eval.resp.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} or a \code{\link[terra:vect]{SpatVector}} object containing binary data 
##' (\code{1} : presence, \code{0} : absence) for a single species that will be used to evaluate 
##' the species distribution model(s) with independent data
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##'  \code{SpatialPoints} (if presence-only) or \code{SpatialPointsDataFrame}
##'  object containing binary data.}
##' @param eval.resp.xy (\emph{optional, default} \code{NULL}) \cr 
##' If \code{resp.var} is a \code{vector}, a 2-columns \code{matrix} or \code{data.frame} 
##' containing the corresponding \code{X} and \code{Y} coordinates that will be
##' used to evaluate the species distribution model(s) with independent data
##' @param eval.expl.var (\emph{optional, default} \code{NULL}) \cr 
##' A \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables 
##' (in columns or layers) that will be used to evalute the species distribution model(s) with 
##' independent data.
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
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
##' points that will be selected for each pseudo-absence repetition (true absences included). \cr
##' It can also be a \code{vector} of the same length as \code{PA.nb.rep} containing \code{integer} 
##' values corresponding to the different numbers of pseudo-absences to be selected (see Details)
##' @param PA.sre.quant (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'sre'}, a \code{numeric} between \code{0} 
##' and \code{0.5} defining the half-quantile used to make the \code{sre} pseudo-absence selection 
##' (see Details)
##' @param PA.dist.min (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' minimal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in the same projection system units as \code{resp.xy} and \code{expl.var}, see Details)
##' @param PA.dist.max (\emph{optional, default} \code{0}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'disk'}, a \code{numeric} defining the 
##' maximal distance to presence points used to make the \code{disk} pseudo-absence selection 
##' (in the same projection system units as \code{resp.xy} and \code{expl.var}, see Details)
##' @param PA.fact.aggr (\emph{optional, default} \code{NULL}) \cr
##' If pseudo-absence selection and \code{PA.strategy = 'random'} or \code{PA.strategy = 'disk'}, 
##' an \code{integer} defining the factor of aggregation to reduce the spatial resolution of the 
##' environmental variables
##' @param PA.user.table (\emph{optional, default} \code{NULL}) \cr 
##' If pseudo-absence selection and \code{PA.strategy = 'user.defined'}, a \code{matrix} or 
##' \code{data.frame} with as many rows as \code{resp.var} values, as many columns as 
##' \code{PA.nb.rep}, and containing \code{TRUE} or \code{FALSE} values defining which points 
##' will be used to build the species distribution model(s) for each repetition (see Details)
##' 
##' @param na.rm (\emph{optional, default} \code{TRUE}) \cr 
##' A \code{logical} value defining whether points having one or several missing values for 
##' explanatory variables should be removed from the analysis or not
##' @param filter.raster (\emph{optional, default} \code{FALSE}) \cr 
##' If \code{expl.var} is of raster type, a \code{logical} value defining whether \code{resp.var} 
##' is to be filtered when several points occur in the same raster cell
##' @param seed.val (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' 
##' 
##' @return 
##' 
##' A \code{\link{BIOMOD.formated.data}} or \code{\link{BIOMOD.formated.data.PA}} object that can 
##' be used to build species distribution model(s) with the \code{\link{BIOMOD_Modeling}} 
##' function. \cr
##' \href{https://biomodhub.github.io/biomod2/reference/BIOMOD.formated.data.html}{\code{print/show}}, 
##' \href{https://biomodhub.github.io/biomod2/reference/plot.html}{\code{plot}} and 
##' \href{https://biomodhub.github.io/biomod2/reference/summary.html}{\code{summary}} functions 
##' are available to have a summary of the created object. 
##' 
##' 
##' @details  
##' 
##' This function gathers and formats all input data needed to run \pkg{biomod2} models. It 
##' supports different kind of inputs (e.g. \code{matrix},
##' \code{\link[terra:vect]{SpatVector}}, \code{\link[terra:rast]{SpatRaster}})
##' and provides different methods to select pseudo-absences if needed. \cr \cr
##' 
##' \bold{Concerning explanatory variables and XY coordinates :} 
##' \itemize{
##'   \item if \code{\link[terra:rast]{SpatRaster}}, \code{RasterLayer} or \code{RasterStack}
##'   provided for \code{expl.var} or \code{eval.expl.var}, \cr \pkg{biomod2} will extract 
##'   the corresponding values from XY coordinates provided :
##'   \itemize{
##'     \item either through \code{resp.xy} or \code{eval.resp.xy} respectively
##'     \item or through \code{resp.var} or \code{eval.resp.var}, if provided as 
##'     \code{\link[terra:vect]{SpatVector}} or \code{SpatialPointsDataFrame}
##'   }
##'   \emph{Be sure to give the objects containing XY coordinates in the same projection 
##'   system than the raster objects !}
##'    
##'   \item if \code{data.frame} or \code{matrix} provided for \code{expl.var} or
##'    \code{eval.expl.var}, \cr \pkg{biomod2} will simply merge it (\code{cbind}) 
##'    with \code{resp.var} without considering XY coordinates. \cr
##'   \emph{Be sure to give explanatory and response values in the same row order !}
##' }
##' 
##' \bold{Concerning pseudo-absence selection (see \code{\link{bm_PseudoAbsences}}) :}
##' \cr \emph{Only in the case of \code{binary} data !}
##' \itemize{
##'   \item if both presence and absence data are available : \code{PA.nb.rep = 0} and no 
##'   pseudo-absence will be selected.
##'   \item if no absence data is available, several pseudo-absence repetitions  
##'   are recommended (to estimate the effect of pseudo-absence selection), as well as high 
##'   number of pseudo-absence points. \cr
##'   \emph{Be sure not to select more pseudo-absence points than maximum number of pixels in 
##'   the studied area !}
##'   \item it is possible to create several pseudo-absence repetitions \emph{with different 
##'   number of points}, BUT with the same sampling strategy. \code{PA.nb.absences} must contain 
##'   as many values as the number of sets of pseudo-absences (\code{PA.nb.rep}). \cr \cr \cr \cr
##' }
##' 
##' \describe{
##'   \item{Response variable}{
##'   \pkg{biomod2} models single species at a time (no multi-species). \cr
##'   Hence, \code{resp.var} must be an uni-dimensional object, either :
##'   \itemize{
##'     \item a \code{vector}, a one-column \code{matrix} or \code{data.frame}, a 
##'     \code{\link[terra:vect]{SpatVector}} (\emph{without associated data - if presence-only})
##'     \item a \code{SpatialPoints} (\emph{if presence-only})
##'     \item a \code{SpatialPointsDataFrame} or \code{\link[terra:vect]{SpatVector}} object
##'   }
##'   If \code{resp.var} is a non-spatial object (\code{vector}, \code{matrix} or 
##'   \code{data.frame}), XY coordinates must be provided through \code{resp.xy}. \cr \cr
##'   Different data types are available, and require different values :
##'   \describe{
##'     \item{binary}{\code{1} : presences, \code{0} : true absences or \code{NA} : no 
##'     information point (can be used to select pseudo-absences) \cr
##'     \emph{If no true absences are available, pseudo-absence selection must be done.}
##'     }
##'     \item{count}{positive \code{integer} values}
##'     \item{multiclass}{\code{factor} values}
##'     \item{ordinal}{ordered \code{factor} values}
##'     \item{relative}{\code{numeric} values between \code{0} and \code{1}}
##'     \item{abundance}{positive \code{numeric} values}
##'   }
##'   }
##'   \item{Explanatory variables}{
##'   Factorial variables are allowed, but might lead to some pseudo-absence strategy or models 
##'   omissions (e.g. \code{sre}).
##'   }
##'   \item{Evaluation data}{
##'   Although \pkg{biomod2} provides tools to automatically divide dataset into calibration and 
##'   validation parts through the modeling process (see \code{CV.[..]} parameters in 
##'   \code{\link{BIOMOD_Modeling}} function ; or \code{\link{bm_CrossValidation} 
##'   function}), it is also possible (and strongly advised) to directly provide two independent 
##'   datasets, one for calibration/validation and one for evaluation
##'   }
##'   \item{Pseudo-absence selection (see \code{\link{bm_PseudoAbsences}})}{
##'   \emph{Only in the case of \code{binary} data !} \cr
##'   If no true absences are available, pseudo-absences must be selected from the 
##'   \emph{background data}, meaning data there is no information whether the species of 
##'   interest occurs or not. It corresponds either to the remaining pixels of the \code{expl.var} 
##'   (if provided as a \code{\link[terra:rast]{SpatRaster}} or \code{RasterStack})
##'   or to the points identified as  \code{NA} in \code{resp.var} (if \code{expl.var}
##'   provided as a \code{matrix} or \code{data.frame}). \cr \cr
##'   Several methods are available to do this selection :
##'   \describe{
##'     \item{random}{all points of initial background are pseudo-absence candidates. 
##'     \code{PA.nb.absences} are drawn randomly, for each \code{PA.nb.rep} requested.
##'     }
##'     \item{sre}{pseudo-absences have to be selected in conditions (combination of explanatory 
##'     variables) that differ in a defined proportion (\code{PA.sre.quant}) from those of 
##'     presence points. A \emph{Surface Range Envelop} model is first run over the species of 
##'     interest (see \code{\link{bm_SRE}}), and pseudo-absences are selected outside this envelop. \cr
##'     \emph{This case is appropriate when all the species climatic niche has been sampled, 
##'     otherwise it may lead to over-optimistic model evaluations and predictions !}
##'     }
##'     \item{disk}{pseudo-absences are selected within circles around presence points defined by 
##'     \code{PA.dist.min} and \code{PA.dist.max} distance values (in the same projection system 
##'     units as \code{coord} and \code{expl.var}). It allows to select pseudo-absence points that 
##'     are not too close to (avoid same niche and pseudo-replication) or too far (localized 
##'     sampling strategy) from presences.
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
##' # Format Data with true absences
##' myBiomodData <- BIOMOD_FormatingData(resp.name = myRespName,
##'                                      resp.var = myResp,
##'                                      resp.xy = myRespXY,
##'                                      expl.var = myExpl)
##' myBiomodData
##' summary(myBiomodData)
##' plot(myBiomodData)
##' 
##' 
##' # ---------------------------------------------------------------#
##' # # Transform true absences into potential pseudo-absences
##' # myResp.PA <- ifelse(myResp == 1, 1, NA)
##' # 
##' # # Format Data with pseudo-absences : random method
##' # myBiomodData.r <- BIOMOD_FormatingData(resp.name = myRespName,
##' #                                        resp.var = myResp.PA,
##' #                                        resp.xy = myRespXY,
##' #                                        expl.var = myExpl,
##' #                                        PA.nb.rep = 4,
##' #                                        PA.nb.absences = 1000,
##' #                                        PA.strategy = 'random')
##' # 
##' # # Format Data with pseudo-absences : disk method
##' # myBiomodData.d <- BIOMOD_FormatingData(resp.name = myRespName,
##' #                                        resp.var = myResp.PA,
##' #                                        resp.xy = myRespXY,
##' #                                        expl.var = myExpl,
##' #                                        PA.nb.rep = 4,
##' #                                        PA.nb.absences = 500,
##' #                                        PA.strategy = 'disk',
##' #                                        PA.dist.min = 5,
##' #                                        PA.dist.max = 35)
##' # 
##' # # Format Data with pseudo-absences : SRE method
##' # myBiomodData.s <- BIOMOD_FormatingData(resp.name = myRespName,
##' #                                        resp.var = myResp.PA,
##' #                                        resp.xy = myRespXY,
##' #                                        expl.var = myExpl,
##' #                                        PA.nb.rep = 4,
##' #                                        PA.nb.absences = 1000,
##' #                                        PA.strategy = 'sre',
##' #                                        PA.sre.quant = 0.025)
##' # 
##' # # Format Data with pseudo-absences : user.defined method
##' # myPAtable <- data.frame(PA1 = ifelse(myResp == 1, TRUE, FALSE),
##' #                         PA2 = ifelse(myResp == 1, TRUE, FALSE))
##' # for (i in 1:ncol(myPAtable)) myPAtable[sample(which(myPAtable[, i] == FALSE), 500), i] = TRUE
##' # myBiomodData.u <- BIOMOD_FormatingData(resp.name = myRespName,
##' #                                        resp.var = myResp.PA,
##' #                                        resp.xy = myRespXY,
##' #                                        expl.var = myExpl,
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
##' # ---------------------------------------------------------------#
##' # # Select multiple sets of pseudo-absences
##' #
##' # # Transform true absences into potential pseudo-absences
##' # myResp.PA <- ifelse(myResp == 1, 1, NA)
##' # 
##' # # Format Data with pseudo-absences : random method
##' # myBiomodData.multi <- BIOMOD_FormatingData(resp.name = myRespName,
##' #                                            resp.var = myResp.PA,
##' #                                            resp.xy = myRespXY,
##' #                                            expl.var = myExpl,
##' #                                            PA.nb.rep = 4,
##' #                                            PA.nb.absences = c(1000, 500, 500, 200),
##' #                                            PA.strategy = 'random')
##' # myBiomodData.multi
##' # summary(myBiomodData.multi)
##' # plot(myBiomodData.multi)
##' 
##' 
##' @importFrom terra rast crds categories
##' 
##' @export
##' 
##' 
###################################################################################################


BIOMOD_FormatingData <- function(resp.name,
                                 resp.var,
                                 resp.xy = NULL,
                                 expl.var,
                                 dir.name = '.',
                                 data.type = 'binary',
                                 eval.resp.var = NULL,
                                 eval.resp.xy = NULL,
                                 eval.expl.var = NULL,
                                 PA.nb.rep = 0,
                                 PA.nb.absences = 1000,
                                 PA.strategy = NULL,
                                 PA.dist.min = 0,
                                 PA.dist.max = NULL,
                                 PA.sre.quant = 0.025,
                                 PA.fact.aggr = NULL,
                                 PA.user.table = NULL,
                                 na.rm = TRUE,
                                 filter.raster = FALSE,
                                 seed.val = NULL)
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
                                           filter.raster)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 2. build BIOMOD.formated.data object -------------------------------------
  out <- NULL
  if (is.null(PA.strategy) || PA.strategy == 'none') { # no Pseudo Absences
    out <- BIOMOD.formated.data(sp = resp.var,
                                xy = resp.xy,
                                env = expl.var,
                                dir.name = dir.name,
                                sp.name = resp.name,
                                data.type = data.type,
                                eval.sp = eval.resp.var,
                                eval.env = eval.expl.var,
                                eval.xy = eval.resp.xy,
                                na.rm = na.rm,
                                filter.raster = filter.raster)
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
                                   PA.fact.aggr = PA.fact.aggr,
                                   PA.user.table = PA.user.table,
                                   na.rm = na.rm,
                                   filter.raster = filter.raster,
                                   seed.val)
  }
  out@call <- match.call()
  
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
                                             filter.raster)
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
  args <- .BIOMOD.formated.data.check.args(sp = resp.var, env = expl.var, xy = resp.xy
                                           , eval.sp = eval.resp.var, eval.env = eval.expl.var
                                           , eval.xy = eval.resp.xy, filter.raster = filter.raster)
  
  return(list(resp.var = args$sp,
              expl.var = args$env,
              resp.xy = args$xy,
              resp.name = resp.name,
              dir.name = dir.name,
              eval.resp.var = args$eval.sp,
              eval.expl.var = args$eval.env,
              eval.resp.xy = args$eval.xy))
}


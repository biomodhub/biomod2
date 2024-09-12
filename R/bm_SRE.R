###################################################################################################
##' @name bm_SRE
##' @author Wilfried Thuiller, Bruno Lafourcade, Damien Georges 
##' 
##' @title Surface Range Envelope
##' 
##' @description This internal \pkg{biomod2} function allows the user to run a rectilinear surface 
##' range envelop (SRE) (equivalent to 
##' \href{https://caws.org.nz/PPQ567/PPQ\%2006-1\%20pp008-9\%20Busby.pdf}{BIOCLIM}) 
##' using the extreme percentiles (as recommended by Nix or Busby, see References and Details).
##' 
##' @param resp.var a \code{vector}, a \code{\link[terra:vect]{SpatVector}}
##' without associated data (\emph{if presence-only}), 
##' or a \code{\link[terra:vect]{SpatVector}} object containing binary data  
##' (\code{0} : absence, \code{1} : presence, \code{NA} : indeterminate) 
##' for a single species that will be used to build the species distribution model(s)
##' \cr \emph{Note that old format from \pkg{sp} are still supported such as
##'  \code{SpatialPoints}  (if presence-only) or \code{SpatialPointsDataFrame}
##'  object containing binary data.}
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables (in 
##' columns or layers) that will be used to build the SRE model
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' 
##' @param new.env a \code{matrix}, \code{data.frame}, \code{\link[terra:vect]{SpatVector}}
##' or \code{\link[terra:rast]{SpatRaster}} object containing the explanatory variables (in 
##' columns or layers) that will be used to predict the SRE model
##' \cr \emph{Note that old format from \pkg{raster} and \pkg{sp} are still supported such as 
##' \code{RasterStack} and \code{SpatialPointsDataFrame} objects. }
##' 
##' @param quant a \code{numeric} between \code{0} and \code{0.5} defining the half-quantile 
##' corresponding to the most extreme value for each variable not to be taken into account for 
##' determining the tolerance boundaries of the considered species (see Details)
##' @param do.extrem (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether a \code{matrix} containing extreme conditions 
##' supported should be returned or not
##' 
##' 
##' @return 
##' 
##' A \code{vector} or a \code{\link[terra:rast]{SpatRaster}} object, containing binary 
##' (\code{0} or \code{1}) values.
##' 
##' 
##' @details 
##' 
##' \emph{Please refer to References to get more information about surface range envelop models.}
##' 
##' This method is highly influenced by the extremes of the data input. Whereas a linear model 
##' can discriminate the extreme values from the main tendency, the SRE considers them as 
##' important as any other data point leading to changes in predictions. \cr \cr
##' 
##' \emph{The more (non-collinear) variables, the more restrictive the model will be.} \cr \cr
##' 
##' Predictions are returned as binary (\code{0} or \code{1}) values, a site being either 
##' potentially suitable for all the variables, or out of bounds for at least one variable and 
##' therefore considered unsuitable. \cr \cr
##' 
##' \code{quant} determines the threshold from which the data will be taken into account for 
##' calibration. The default value of \code{0.05} induces that the \code{5\%} most extreme values 
##' will be avoided for each variable on each side of its distribution along the gradient, meaning 
##' that a total of \code{10\%} of the data will not be considered.
##' 
##' 
##' @references
##' 
##' \itemize{
##'   \item Nix, H.A., 1986. A biogeographic analysis of Australian elapid snakes. In: 
##'   \emph{Atlas of Elapid Snakes of Australia.} (Ed.) R. Longmore, pp. 4-15. 
##'   \bold{Australian Flora and Fauna Series Number 7.} 
##'   Australian Government Publishing Service: Canberra.
##'   \item Busby, Jeremy. BIOCLIM - a bioclimate analysis and prediction system. 
##'   \emph{Plant protection quarterly} \bold{6} (1991): 8-9.
##' }
##' 
##' @keywords models "surface range envelop" sre quantile
##' 
##' 
##' @seealso \code{\link{bm_PseudoAbsences}}, \code{\link{BIOMOD_FormatingData}}, 
##' \code{\link{bm_ModelingOptions}}, \code{\link{bm_Tuning}}, 
##' \code{\link{bm_RunModelsLoop}}, \code{\link{BIOMOD_Modeling}},
##' @family Secondary functions
##' 
##' @examples
##' 
##' library(terra)
##' ## Load real data
##' data(DataSpecies)
##' myResp.r <- as.numeric(DataSpecies[, 'GuloGulo'])
##' 
##' data(bioclim_current)
##' myExpl.r <- rast(bioclim_current)
##' 
##' myRespXY <- DataSpecies[which(myResp.r == 1), c('X_WGS84', 'Y_WGS84')]
##' myResp.v <- classify(subset(myExpl.r, 1), 
##'                      matrix(c(-Inf, Inf, 0), ncol = 3, byrow = TRUE))
##' myResp.v[cellFromXY(myResp.v, myRespXY)] <- 1
##' 
##' ## Compute SRE for several quantile values
##' sre.100 <- bm_SRE(resp.var = myResp.v,
##'                   expl.var = myExpl.r,
##'                   new.env = myExpl.r,
##'                   quant = 0)
##' sre.095 <- bm_SRE(resp.var = myResp.v,
##'                   expl.var = myExpl.r,
##'                   new.env = myExpl.r,
##'                   quant = 0.025)
##' sre.090 <- bm_SRE(resp.var = myResp.v,
##'                   expl.var = myExpl.r,
##'                   new.env = myExpl.r,
##'                   quant = 0.05)
##' 
##' ## Visualize results
##' res <- c(myResp.v, sre.100, sre.095, sre.090)
##' names(res) <- c("Original distribution", "Full data calibration"
##'                 , "Over 95 percent", "Over 90 percent")
##' plot(res)
##' 
##' 
##' @importFrom terra rast values vect quantile cellFromXY
## quantile classify crds global is.factor mask nlyr subset
##' 
##' @export
##' 
##' 
###################################################################################################

## Remi 20/10/2022
## This function seems to have support for multispecies resp.var although
## the rest of the package do not.

bm_SRE <- function(resp.var = NULL, 
                   expl.var = NULL, 
                   new.env = NULL, 
                   quant = 0.025, 
                   do.extrem = FALSE)
{
  ## 0. Check arguments ---------------------------------------------------------
  args <- .bm_SRE.check.args(resp.var, expl.var, new.env, quant)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)

  ## 1. Determine suitable conditions and make the projection --------------
  lout <- list()
  if (is.data.frame(resp.var) | is.matrix(resp.var)) {
    ### matrix or data.frame ------------
    nb.resp <- ncol(resp.var)
    resp.names <- colnames(resp.var)
    for (j in 1:nb.resp) {
      occ.pts <- which(resp.var[, j] == 1)
      extrem.cond <- t(apply(as.data.frame(expl.var[occ.pts, , drop = FALSE]), 2, quantile
                             , probs = c(0 + quant, 1 - quant), na.rm = TRUE))
      if (!do.extrem) {
        lout[[j]] <- .sre_projection(new.env, extrem.cond)
      }
    }
  } else if (inherits(resp.var, 'SpatRaster')) {
    ## Remi 20/10/2022
    ## this section seems obsolete as resp.var support for Raster/SpatRaster
    ## but also SpatVector/SpatialPoints is currently deactivated.
    ### raster ------------------------------------
    nb.resp <- nlyr(resp.var)
    resp.names <- names(resp.var)
    for (j in 1:nb.resp) {
      extrem.cond <- global(mask(expl.var, subset(resp.var, j), maskvalues = c(0,NA)),
                            fun = quantile,
                            probs = c(0 + quant, 1 - quant),
                            na.rm = TRUE)
      if (!do.extrem) {
        lout[[j]] <- .sre_projection(new.env, extrem.cond)
      }
    }
  } else if (inherits(resp.var, 'SpatVector')) {
    ### SpatVector ----------------------
    nb.resp <- ncol(values(resp.var))
    resp.names <- colnames(values(resp.var))
    for (j in 1:nb.resp) {
      occ.pts <- which(values(resp.var)[, j] == 1)
      if (is.data.frame(expl.var) || is.matrix(expl.var)) {
        extrem.cond <- t(apply(as.data.frame(expl.var[occ.pts, ]), 2, quantile
                               , probs = c(0 + quant, 1 - quant), na.rm = TRUE))
      } else {
        if (inherits(expl.var, 'SpatRaster')) {
          maskTmp <- subset(expl.var, 1)
          maskTmp[] <- NA
          maskTmp[cellFromXY(maskTmp, crds(resp.var)[occ.pts, ])] <- 1
          extrem.cond <- global(mask(expl.var, maskTmp),
                                fun = quantile,
                                probs = c(0 + quant, 1 - quant),
                                na.rm = TRUE)
          
        } else {
          if (inherits(expl.var, 'SpatialPoints')) {
            ## May be good to check corespondances of resp.var and expl.var variables
            extrem.cond <- t(apply(as.data.frame(expl.var[occ.pts, ]), 2, quantile
                                   , probs = c(0 + quant, 1 - quant), na.rm = TRUE))
          } else {
            stop("Unsuported case!")
          }
        }
      }
      if (!do.extrem) {
        lout[[j]] <- .sre_projection(new.env, extrem.cond)
      }
    }
  }
  
  
  ## RETURN results -------------------------------------------------------------------------------
  if (do.extrem) {
    return(as.data.frame(extrem.cond))
  } else {
    # 3. Rearranging the lout object
    if (is.data.frame(new.env)) {
      lout <- simplify2array(lout)
      colnames(lout) <- resp.names
    } else if (inherits(new.env, 'SpatRaster')) {
      lout <- rast(lout)
      if (nlyr(lout) == 1) {
        lout <- subset(lout, 1)
      }
      names(lout) <- resp.names
    }
    return(lout)
  }
}


## SRE Argument check ---------------------------------------------------------

.bm_SRE.check.args <- function(resp.var = NULL, expl.var = NULL, new.env = NULL, quant = 0.025)
{
  ## 0. Check compatibility between resp.var and expl.var arguments -----------
  if (is.vector(resp.var) || inherits(resp.var, c("matrix","data.frame"))) {
    resp.var <- as.data.frame(resp.var)
    if (!is.vector(expl.var) && !inherits(expl.var, c("matrix","data.frame",'SpatVector')))  {
      stop("\n resp.var and expl.var arguments must be of same type (both vector, both matrix, etc)")
    } else {
      if (inherits(expl.var, 'SpatVector')) {
        expl.var <- values(expl.var)
      }
      expl.var <- as.data.frame(expl.var)
      nb.expl.vars <- ncol(expl.var)
      names.expl.vars <- colnames(expl.var)
      if (nrow(resp.var) != nrow(expl.var)) {
        stop("resp.var and expl.var arguments must have the same number of rows")
      }
    }
  }
  
  if (inherits(expl.var, 'SpatVector')) {
    expl.var <- values(expl.var)
    nb.expl.vars <- ncol(expl.var)
    names.expl.vars <- colnames(expl.var)
  }
  
  # back-compatibility with raster package
  if (inherits(resp.var, 'Raster')) {
    if (!inherits(expl.var, 'Raster')) {
      stop("\n resp.var and expl.var arguments must be of same type (both vector, both raster, etc)")
    }
    resp.var <- rast(resp.var)
    expl.var <- rast(expl.var)
  }
  
  if (inherits(resp.var, 'SpatRaster')) {
    if (!inherits(expl.var, 'SpatRaster')) {
      stop("\n resp.var and expl.var arguments must be of same type (both vector, both raster, etc)")
    }
    nb.expl.vars <- nlyr(expl.var)
    names.expl.vars <- names(expl.var)
  }

  ## 1. Check expl.var argument -----------------------------------------------
  if ((inherits(expl.var, 'data.frame') && any(sapply(expl.var, is.factor))) ||
      (inherits(expl.var, c('SpatRaster')) && any(is.factor(expl.var)))) {
    stop("SRE algorithm does not handle factorial variables")
  }

  
  ## 2. Check new.env argument ------------------------------------------------
  if (is.null(new.env)) { ## if no new.env, projection done on expl.var variables
    new.env <- expl.var
  } else { ## check of compatible number of explanatory variables
    if (is.vector(new.env) || is.data.frame(new.env) || is.matrix(new.env))
    {
      new.env <- as.data.frame(new.env)
      if (!all(names.expl.vars %in% colnames(new.env))) {
        stop("expl.var variables names differs in the 2 dataset given")
      }
      new.env <- new.env[,names.expl.vars]
      if (ncol(new.env) != nb.expl.vars) {
        stop("Incompatible number of variables in new.env objects")
      }
    } else if (!inherits(new.env, 'SpatRaster')) {
      new.env <- rast(new.env)
      if (sum(!(names.expl.vars %in% names(new.env))) > 0) {
        stop("expl.var variables names differs in the 2 dataset given")
      }
      new.env <- subset(new.env, names.expl.vars)
      if (nlyr(new.env) != nb.expl.vars) {
        stop("Incompatible number of variables in new.env objects")
      }
    }
  }
  
  ## 3. Check quant argument --------------------------------------------------
  if (quant < 0 || quant >= 0.5) {
    stop("\n quantmust be a 0 to 0.5 numeric")
  }
  
  return(list(resp.var = resp.var,
              expl.var = expl.var,
              new.env = new.env,
              quant = quant))
}


# SRE Projection  ----------------------------------------------------------

.sre_projection <- function(new.env, extrem.cond, mod.name = NULL) {
  if (is.data.frame(new.env) || is.matrix(new.env)) {
    out <- rep(1, nrow(new.env))
    for (this.var in rownames(extrem.cond)) {
      out <- out * as.numeric(new.env[, this.var] >= extrem.cond[this.var, 1] &
                                new.env[, this.var] <= extrem.cond[this.var, 2])
    }
  } else if (inherits(new.env, "SpatRaster")) {
    out <- classify(subset(new.env, 1), 
                    matrix(c(-Inf, Inf, 1), ncol = 3, byrow = TRUE),
                    wopt = list(names = mod.name))
    for (this.var in rownames(extrem.cond)) {
      out <- out * (subset(new.env, this.var) >= extrem.cond[this.var, 1]) * 
        (subset(new.env, this.var) <= extrem.cond[this.var, 2])
    }
    out <- subset(out, 1)
  }
  return(out)
}


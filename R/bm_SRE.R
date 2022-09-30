###################################################################################################
##' @name bm_SRE
##' @author Wilfried Thuiller, Bruno Lafourcade, Damien Georges 
##' 
##' @title Surface Range Envelope
##' 
##' @description This internal \pkg{biomod2} function allows the user to run a rectilinear surface 
##' range envelop (SRE) (equivalent to 
##' \href{https://caws.org.nz/PPQ567/PPQ\%2006-1\%20pp008-9\%20Busby.pdf}{BIOCLIM}) 
##' using the extreme percentiles (as recommended by Nix or Busby, see 
##' \href{https://biomodhub.github.io/biomod2/reference/bm_SRE.html#references}{References} and Details).
##' 
##' 
##' @param resp.var a \code{vector}, \code{\link[sp]{SpatialPoints}} or 
##' \code{\link[sp]{SpatialPointsDataFrame}} object containing binary data (\code{0} : absence, 
##' \code{1} : presence, \code{NA} : indeterminate) for a single species
##' @param expl.var a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables (in 
##' columns or layers) that will be used to build the SRE model
##' @param new.env a \code{matrix}, \code{data.frame}, \code{\link[sp]{SpatialPointsDataFrame}} 
##' or \code{\link[raster:stack]{RasterStack}} object containing the explanatory variables (in 
##' columns or layers) that will be used to predict the SRE model
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
##' A \code{vector} or a \code{\link[raster:raster]{raster}} object, containing binary (\code{0} 
##' or \code{1}) values.
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
##' \emph{The more (non-colinear) variables, the more restrictive the model will be.} \cr \cr
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
##' \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Tuning}}, 
##' \code{\link{bm_RunModelsLoop}}, \code{\link{BIOMOD_Modeling}},
##' @family Secundary functions
##' 
##' 
##' @examples
##' 
##' ## Load real data
##' myFile <- system.file('external/species/mammals_table.csv', package = 'biomod2')
##' DataSpecies <- read.csv(myFile, row.names = 1)
##' myResp.r <- as.numeric(DataSpecies[, 'GuloGulo'])
##' 
##' myFiles <- paste0('external/bioclim/current/bio', c(3, 4, 7, 11, 12), '.grd')
##' myExpl.r <- raster::stack(system.file(myFiles, package = 'biomod2'))
##' 
##' myRespXY <- DataSpecies[which(myResp.r == 1), c('X_WGS84', 'Y_WGS84')]
##' myResp.v <- raster::reclassify(raster::subset(myExpl.r, 1, drop = TRUE), c(-Inf, Inf, 0))
##' myResp.v[raster::cellFromXY(myResp.v, myRespXY)] <- 1
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
##' res <- raster::stack(myResp.v, sre.100, sre.095, sre.090)
##' names(res) <- c("Original distribution", "Full data calibration"
##'                 , "Over 95 percent", "Over 90 percent")
##' plot(res, zlim = c(0, 1))
##' 
##' 
##' @importFrom raster stack subset nlayers mask reclassify coordinates cellFromXY Which 
## quantile
##' 
##' @export
##' 
##' 
###################################################################################################


bm_SRE <- function(resp.var = NULL, 
                   expl.var = NULL, 
                   new.env = NULL, 
                   quant = 0.025, 
                   do.extrem = FALSE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_SRE.check.args(resp.var, expl.var, new.env, quant)
  for (argi in names(args)) { assign(x = argi, value = args[[argi]]) }
  rm(args)
  
  ## 1. Determine suitable conditions and make the projection -------------------------------------
  lout <- list()
  if (is.data.frame(resp.var) | is.matrix(resp.var)) ## matrix or data.frame ------------
  {
    nb.resp <- ncol(resp.var)
    resp.names <- colnames(resp.var)
    for (j in 1:nb.resp) {
      occ.pts <- which(resp.var[, j] == 1)
      extrem.cond <- t(apply(as.data.frame(expl.var[occ.pts, ]), 2, quantile
                             , probs = c(0 + quant, 1 - quant), na.rm = TRUE))
      if (!do.extrem) {
        lout[[j]] <- .sre_projection(new.env, extrem.cond)
      }
    }
  } else if (inherits(resp.var, 'Raster')) ## raster ------------------------------------
  {
    nb.resp <- nlayers(resp.var)
    resp.names <- names(resp.var)
    for (j in 1:nb.resp) {
      occ.pts <- subset(resp.var, j, drop = TRUE)
      x.ooc.pts <- Which(occ.pts != 1, cells = TRUE, na.rm = TRUE)
      occ.pts[x.ooc.pts] <- rep(NA, length(x.ooc.pts))
      extrem.cond <- raster::quantile(mask(expl.var, occ.pts),
                                      probs = c(0 + quant, 1 - quant),
                                      na.rm = TRUE)
      if (!do.extrem) {
        lout[[j]] <- .sre_projection(new.env, extrem.cond)
      }
    }
  } else if (inherits(resp.var, 'SpatialPoints')) ## SpatialPoints ----------------------
  {
    nb.resp <- ncol(resp.var@data)
    resp.names <- colnames(resp.var@data)
    for (j in 1:nb.resp) {
      occ.pts <- which(resp.var@data[, j] == 1)
      if (is.data.frame(expl.var) || is.matrix(expl.var)) {
        extrem.cond <- t(apply(as.data.frame(expl.var[occ.pts, ]), 2, quantile
                               , probs = c(0 + quant, 1 - quant), na.rm = TRUE))
      } else {
        if (inherits(expl.var, 'Raster')) {
          maskTmp <- subset(expl.var, 1, drop = TRUE)
          maskTmp[] <- NA
          maskTmp[cellFromXY(maskTmp, coordinates(resp.var)[occ.pts, ])] <- 1
          extrem.cond <- raster::quantile(mask(expl.var, maskTmp),
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
    } else if (inherits(new.env, 'Raster')) {
      lout <- stack(lout)
      if (nlayers(lout) == 1) {
        lout <- subset(lout, 1, drop = TRUE)
      }
      names(lout) <- resp.names
    }
    return(lout)
  }
}


###################################################################################################

.bm_SRE.check.args <- function(resp.var = NULL, expl.var = NULL, new.env = NULL, quant = 0.025)
{
  ## 0. Check compatibility between resp.var and expl.var arguments -----------
  if (is.vector(resp.var) || is.matrix(resp.var) || is.data.frame(resp.var))
  {
    resp.var <- as.data.frame(resp.var)
    
    if (!is.vector(expl.var) && !is.matrix(expl.var) && !is.data.frame(expl.var) && !inherits(expl.var, 'SpatialPoints'))
    {
      stop("\n resp.var and expl.var arguments must be of same type (both vector, both matrix, etc)")
    } else {
      if (inherits(expl.var, 'SpatialPoints')) {
        expl.var <- as.data.frame(expl.var@data)
      }
      expl.var <- as.data.frame(expl.var)
      nb.expl.vars <- ncol(expl.var)
      names.expl.vars <- colnames(expl.var)
      if (nrow(resp.var) != nrow(expl.var)) {
        stop("resp.var and expl.var arguments must have the same number of rows")
      }
    }
  }
  
  if (inherits(expl.var, 'SpatialPoints')) {
    expl.var <- as.data.frame(expl.var@data)
    nb.expl.vars <- ncol(expl.var)
    names.expl.vars <- colnames(expl.var)
  }
  
  if (inherits(resp.var, 'Raster')) {
    if (!inherits(expl.var, 'Raster')) {
      stop("\n resp.var and expl.var arguments must be of same type (both vector, both raster, etc)")
    }
    nb.expl.vars <- nlayers(expl.var)
    names.expl.vars <- names(expl.var)
  }
  
  
  ## 1. Check expl.var argument -----------------------------------------------
  test_no_factorial_var <- TRUE
  if ((is.data.frame(expl.var) && any(unlist(lapply(expl.var, is.factor)))) ||
      (inherits(expl.var, 'Raster') && any(is.factor(expl.var)))) {
    test_no_factorial_var <- FALSE
  }
  if (!test_no_factorial_var) stop("SRE algorithm does not handle factorial variables")
  
  
  ## 2. Check new.env argument ------------------------------------------------
  if (is.null(new.env)) { ## if no new.env, projection done on expl.var variables
    new.env <- expl.var
  } else { ## check of compatible number of explanatory variables
    if (is.vector(new.env) || is.data.frame(new.env) || is.matrix(new.env))
    {
      new.env <- as.data.frame(new.env)
      if (sum(!(names.expl.vars %in% colnames(new.env))) > 0) {
        stop("expl.var variables names differs in the 2 dataset given")
      }
      new.env <- new.env[,names.expl.vars]
      if (ncol(new.env) != nb.expl.vars) {
        stop("Incompatible number of variables in new.env objects")
      }
    } else if (!inherits(new.env, 'Raster')) {
      new.env <- stack(new.env)
      if (sum(!(names.expl.vars %in% names(new.env))) > 0) {
        stop("expl.var variables names differs in the 2 dataset given")
      }
      new.env <- subset(new.env, names.expl.vars)
      if (nlayers(new.env) != nb.expl.vars) {
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


###################################################################################################

.sre_projection <- function(new.env, extrem.cond)
{
  if (is.data.frame(new.env) || is.matrix(new.env)) {
    out <- rep(1, nrow(new.env))
    for (j in 1:ncol(new.env)) {
      out <- out * as.numeric(new.env[, j] >= extrem.cond[j, 1] &
                                new.env[, j] <= extrem.cond[j, 2])
    }
  } else if (inherits(new.env, "Raster")) {
    out <- reclassify(subset(new.env, 1, drop = TRUE), c(-Inf, Inf, 1))
    for (j in 1:nlayers(new.env)) {
      out <- out * (subset(new.env, j, drop = TRUE) >= extrem.cond[j, 1]) * 
        (subset(new.env, j, drop = TRUE) <= extrem.cond[j, 2])
    }
    out <- subset(out, 1, drop = TRUE)
  }
  return(out)
}


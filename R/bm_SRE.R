##' @name sre
##' @title Surface Range Envelope
##' @description 
##' Run a rectilinear surface range envelop (equivalent to BIOCLIM) using
##' the extreme percentiles as recommended by Nix or Busby.
##' The SRE performs a simple analysis of within which range of each
##' variable the data is recorded and renders predictions.
##' 
##' @param Response  a vector (resp. a matrix/data.frame, a 
##' SpatialPointsDataFrame or a RasterLayer ) giving your species'
##' presences and absences data
##' @param Explanatory a matrix (resp. a data.frame, a 
##' SpatialPointsDataFrame or a RasterStack ) containing the environmental
##' variables for the sites given in Response. It must have as many rows
##' as there are elements in Response.
##' @param NewData The data for which you want to render predictions with
##' the sre. It must be a matrix (resp. a data.frame, a 
##' SpatialPointsDataFrame or a RasterStack ) of the same type as the one
##' given in Explanatory and with precisely the same variable names.
##' @param Quant the value defines the most extreme values for each
##' variable not to be taken into account for determining the tolerance
##' boundaries for the considered species.
##' @param return_extremcond boolean, if TRUE a matrix containing extreme
##' conditions supported is returned
##' 
##' @details 
##' The more variables you put in, the more restrictive your model will 
##' be (if non-colinear variables).
##' 
##' This method is very much influenced by the data input, and more
##' specifically by the extremes.
##' 
##' Where a linear model can discriminate the extreme values from the main
##' tendency, the SRE considers it equal as any other data point which
##' leads to notable changes in predictions.
##' 
##' Note that, as a consequence of its functioning, the predictions are
##' directly given in binary, a site being either potentially suitable for
##' all the variables, either out of bounds for at least one variable and
##' therefore considered unsuitable. 
##' 
##' The quants argument determines the threshold at which the data will be
##' taken into account for calibration : the default of 0.05 induces that
##' the 5\% most extreme values will be avoided for each variable on each
##' side of its distribution along the gradient. So it in fact takes 5\%
##' away at each end of the variables distribution, giving a total of 10\%
##' of data not considered.
##' 
##' @return 
##' A vector (resp. a RasterLayer ) of the same length as there are rows
##' in NewData giving the prediction in binary (1=presence, 0=absence)
##' 
##' @author Wilfried Thuiller, Bruno Lafourcade, Damien Georges 
##' @seealso \code{\link[biomod2]{BIOMOD_Modeling}}, 
##' \code{\link[biomod2]{BIOMOD_ModelingOptions}}, 
##' \code{\link[biomod2]{BIOMOD_Projection}}
##' 
##' @keywords models
##' @keywords multivariate
##' 
##' @examples
##' require(raster)
##' ##' species occurrences
##' DataSpecies <- 
##'   read.csv(
##'     system.file("external/species/mammals_table.csv", package = "biomod2"), 
##'     row.names = 1
##'   )
##' head(DataSpecies)
##' 
##' ##' the name of studied species
##' myRespName <- 'GuloGulo'
##' 
##' ##' the presence/absences data for our species 
##' myResp <- as.numeric(DataSpecies[,myRespName])
##' 
##' ##' the XY coordinates of species data
##' myRespXY <- DataSpecies[which(myResp==1),c("X_WGS84","Y_WGS84")]
##' 
##' ##' Environmental variables extracted from BIOCLIM (bio_3, 
##' ##' bio_4, bio_7, bio_11 & bio_12)
##' myExpl <- 
##'   raster::stack(
##'     system.file("external/bioclim/current/bio3.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio4.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio7.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio11.grd", package = "biomod2"),
##'     system.file("external/bioclim/current/bio12.grd", package = "biomod2")
##'   )

##' myResp <- 
##'   raster::reclassify(
##'     subset(myExpl, 1, drop = TRUE), c(-Inf, Inf, 0)
##'   )
##' myResp[cellFromXY(myResp,myRespXY)] <- 1
##' 
##' ##' Compute some SRE for several quantile values
##' sre.100 <- 
##'   sre(
##'     Response = myResp, 
##'     Explanatory = myExpl, 
##'     NewData=myExpl, 
##'     Quant = 0
##'   )
##'   
##' sre.095 <- 
##'   sre(
##'     Response = myResp, 
##'     Explanatory = myExpl, 
##'     NewData=myExpl, 
##'     Quant = 0.025
##'   )
##' 
##' sre.090 <- 
##'   sre(
##'     Response = myResp, 
##'     Explanatory = myExpl, 
##'     NewData=myExpl, 
##'     Quant = 0.05
##'   )
##'   
##' ##' visualise results
##' par(mfrow=c(2,2),mar=c(6, 5, 5, 3))
##' plot(myResp, main = paste(myRespName, "original distrib."))
##' plot(sre.100, main="full data calibration")
##' plot(sre.095, main="95 %")
##' plot(sre.090, main="90 %")
###################################################################################################


bm_SRE <- function(Response = NULL, 
                Explanatory = NULL, 
                NewData = NULL, 
                Quant = 0.025, 
                return_extremcond = FALSE)
{
  ## 0. Check arguments ---------------------------------------------------------------------------
  args <- .bm_SRE.check.args(Response, Explanatory, NewData, Quant)
  Response <- args$Response
  Explanatory <- args$Explanatory
  NewData <- args$NewData
  Quant <- args$Quant
  rm("args")
  
  ## 1. Determine suitable conditions and make the projection -------------------------------------
  lout <- list()
  if (is.data.frame(Response) | is.matrix(Response)) ## matrix or data.frame ------------
  {
    nb.resp <- ncol(Response)
    resp.names <- colnames(Response)
    for (j in 1:nb.resp) {
      occ.pts <- which(Response[, j] == 1)
      extrem.cond <- t(apply(as.data.frame(Explanatory[occ.pts, ]), 2, quantile
                             , probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))
      if (!return_extremcond) {
        lout[[j]] <- .sre.projection(NewData, extrem.cond)
      }
    }
  } else if (inherits(Response, 'Raster')) ## raster ------------------------------------
  {
    nb.resp <- nlayers(Response)
    resp.names <- names(Response)
    for (j in 1:nb.resp) {
      occ.pts <- subset(Response, j, drop = TRUE)
      x.ooc.pts <- Which(occ.pts != 1, cells = TRUE, na.rm = TRUE)
      occ.pts[x.ooc.pts] <- rep(NA, length(x.ooc.pts))
      extrem.cond <- quantile(mask(Explanatory, occ.pts),
                              probs = c(0 + Quant, 1 - Quant),
                              na.rm = TRUE))
      if (!return_extremcond) {
        lout[[j]] <- .sre.projection(NewData, extrem.cond)
      }
    }
  } else if (inherits(Response, 'SpatialPoints')) ## SpatialPoints ----------------------
  {
    nb.resp <- ncol(Response@data)
    resp.names <- colnames(Response@data)
    for (j in 1:nb.resp) {
      occ.pts <- which(Response@data[, j] == 1)
      if (is.data.frame(Explanatory) || is.matrix(Explanatory)) {
        extrem.cond <- t(apply(as.data.frame(Explanatory[occ.pts, ]), 2, quantile
                               , probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))
      } else {
        if (inherits(Explanatory, 'Raster')) {
          maskTmp <- subset(Explanatory, 1, drop = TRUE)
          maskTmp[] <- NA
          maskTmp[cellFromXY(maskTmp, coordinates(Response)[occ.pts, ])] <- 1
          extrem.cond <- quantile(mask(Explanatory, maskTmp),
                                  probs = c(0 + Quant, 1 - Quant),
                                  na.rm = TRUE)
        } else {
          if (inherits(Explanatory, 'SpatialPoints')) {
            ## May be good to check corespondances of Response and Explanatory variables
            extrem.cond <- t(apply(as.data.frame(Explanatory[occ.pts, ]), 2, quantile
                                   , probs = c(0 + Quant, 1 - Quant), na.rm = TRUE))
          } else {
            stop("Unsuported case!")
          }
        }
      }
      if (!return_extremcond) {
        lout[[j]] <- .sre.projection(NewData, extrem.cond)
      }
    }
  }
  
  
  ## RETURN results -------------------------------------------------------------------------------
  
  if (return_extremcond) {
    return(as.data.frame(extrem.cond))
  } else {
    # 3. Rearranging the lout object
    if (is.data.frame(NewData)) {
      lout <- simplify2array(lout)
      colnames(lout) <- resp.names
    } else if (inherits(NewData, 'Raster')) {
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

.bm_SRE.check.args <- function(Response = NULL, Explanatory = NULL, NewData = NULL, Quant = 0.025)
{
  ## 0. Check compatibility between Response and Explanatory arguments --------
  if (is.vector(Response) || is.matrix(Response) || is.data.frame(Response))
  {
    Response <- as.data.frame(Response)
    
    if (!is.vector(Explanatory) || !is.matrix(Explanatory) || is.data.frame(Explanatory) || !inherits(Explanatory, 'SpatialPoints'))
    {
      stop("\n Response and Explanatory arguments must be of same type (both vector, both matrix, etc)")
    } else {
      if (inherits(Explanatory, 'SpatialPoints')) {
        Explanatory <- as.data.frame(Explanatory@data)
      }
      Explanatory <- as.data.frame(Explanatory)
      nb.expl.vars <- ncol(Explanatory)
      names.expl.vars <- colnames(Explanatory)
      if (nrow(Response) != nrow(Explanatory)) {
        stop("Response and Explanatory arguments must have the same number of rows")
      }
    }
  }
  
  if (inherits(Explanatory, 'SpatialPoints')) {
    Explanatory <- as.data.frame(Explanatory@data)
    nb.expl.vars <- ncol(Explanatory)
    names.expl.vars <- colnames(Explanatory)
  }
  
  if (inherits(Response, 'Raster')) {
    if (!inherits(Explanatory, 'Raster')) {
      stop("\n Response and Explanatory arguments must be of same type (both vector, both raster, etc)")
    }
    nb.expl.vars <- nlayers(Explanatory)
    names.expl.vars <- names(Explanatory)
  }
  
  
  ## 1. Check Explanatory argument --------------------------------------------
  test_no_factorial_var <- TRUE
  if ((is.data.frame(Explanatory) && any(unlist(lapply(Explanatory, is.factor)))) ||
      (inherits(Explanatory, 'Raster') && any(is.factor(Explanatory)))) {
    test_no_factorial_var <- FALSE
  }
  if (!test_no_factorial_var) stop("SRE algorithm does not handle factorial variables")
  
  
  ## 2. Check NewData argument ------------------------------------------------
  if (is.null(NewData)) { ## if no NewData, projection done on Explanatory variables
    NewData <- Explanatory
  } else { ## check of compatible number of explanatory variables
    if (is.vector(NewData) || is.data.frame(NewData) || is.matrix(NewData))
    {
      NewData <- as.data.frame(NewData)
      if (sum(!(names.expl.vars %in% colnames(NewData))) > 0) {
        stop("Explanatory variables names differs in the 2 dataset given")
      }
      NewData <- NewData[,names.expl.vars]
      if (ncol(NewData) != nb.expl.vars) {
        stop("Incompatible number of variables in NewData objects")
      }
    } else if (!inherits(NewData, 'Raster')) {
      NewData <- stack(NewData)
      if (sum(!(names.expl.vars %in% names(NewData))) > 0) {
        stop("Explanatory variables names differs in the 2 dataset given")
      }
      NewData <- subset(NewData, names.expl.vars)
      if (nlayers(NewData) != nb.expl.vars) {
        stop("Incompatible number of variables in NewData objects")
      }
    }
  }
  
  ## 3. Check Quant argument --------------------------------------------------
  if (Quant < 0 || Quant >= 0.5) {
    stop("\n Quantmust be a 0 to 0.5 numeric"))
  }
  
  return(list(Response = Response,
              Explanatory = Explanatory,
              NewData = NewData,
              Quant = Quant))
}


###################################################################################################

.sre.projection <- function(NewData, ExtremCond)
{
  if (is.data.frame(NewData) || is.matrix(NewData)) {
    out <- rep(1, nrow(NewData))
    for (j in 1:ncol(NewData)) {
      out <- out * as.numeric(NewData[, j] >= ExtremCond[j, 1] &
                                NewData[, j] <= ExtremCond[j, 2])
    }
  } else if (inherits(NewData, "Raster")) {
    out <- reclassify(subset(NewData, 1, drop = TRUE), c(-Inf, Inf, 1))
    for (j in 1:nlayers(NewData)) {
      out <- out * (subset(NewData, j, drop = TRUE) >= ExtremCond[j, 1]) * 
        (subset(NewData, j, drop = TRUE) <= ExtremCond[j, 2])
    }
    out <- subset(out, 1, drop = TRUE)
  }
  return(out)
}


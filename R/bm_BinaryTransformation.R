##' @title Convert species' probability of occurrence into binary presence-absence
##'   data using a predefined threshold
##'
##' @description Function that converts an object containing probability values
##'   into a binary presence-absence object according to a pre-defined threshold(s).
##'
##' @param data numeric vector, a \code{matrix}, a \code{data.frame}, a
##'   \code{RasterLayer} or a \code{RasterStack} containing the data to be
##'   converted
##' @param threshold numeric value or a vector containing the threshold to be
##'   used for converting data.
##'
##' @details
##'   If data is a vector or a raster object, then the threshold should be a
##'   numeric value. If data is matrix,dataframe or rasterStack, then the threshold
##'   should have, in theory, as many values as the number of columns or layers
##'   to transform.
##'   In the particular case that the data to convert is a \code{matrix}/\code{data.frame}
##'   with several columns or a \code{RasterStack} with several layers and the
##'   threshold is a single numeric value, the same threshold will be applied
##'   to all columns (resp. layers).
##'
##'
##' @return An object of the same class than \code{data} with binary (0 or 1) values,
##'   usually presence-absence.
##'
##' @author Wilfried Thuiller, Damien Georges
##'
##' @examples
##'   xx <- rnorm(50,10)
##'   yy <- BinaryTransformation(xx, 10)
##'
##'   cbind(xx,yy)
##'
##' @keywords models
##'
##' @export
##' @docType methods
##' @rdname BinaryTransformation-methods


# .convertBin.matrix = function(x, y) {
#   moa <- apply((x > y), 2, as.integer)
#   if (ncol(moa) == 1) { return(moa[, 1]) } else { return(moa) }
# }
# FUN1 = function(data, threshold) {
#   return(sweep(data, 2, threshold, .convertBin.matrix))
# }

.convertBin.matrix = function(data, threshold, doFiltering = FALSE) {
  ind.0 = t(t(data)<threshold)
  data[ind.0] <- 0
  if (!doFiltering) { data[!ind.0] <- 1 }
  if (ncol(data) == 1) { return(data[, 1]) } else { return(data) }
}

.convertBin.array = function(x, y) {
  if (!any(is.na(x))) {
    return(x > y)
  } else {
    return(rep(NA, length(x)))
  }
}

.convertBin.array.filt = function(x, y) {
  if (!any(is.na(x))) {
    return(ifelse(x > y, x, 0))
  } else {
    return(rep(NA, length(x)))
  }
}


setGeneric("bm_BinaryTransformation",
           function(data, threshold, doFiltering = FALSE)
           {
             standardGeneric("bm_BinaryTransformation")
           })

setMethod('bm_BinaryTransformation', signature('data.frame'),
          function(data, threshold, doFiltering = FALSE)
          {
            if (is.numeric(threshold) && !is.na(threshold)) {
              data <- data.matrix(data)
              return(.convertBin.matrix(data, threshold, doFiltering))
            } else { ## return NAs
              return(matrix(NA, ncol = ncol(data), nrow = nrow(data)
                            , dimnames = dimnames(data)))
            }
          })

setMethod('bm_BinaryTransformation', signature('matrix'),
          function(data, threshold, doFiltering = FALSE)
          {
            data <- as.data.frame(data)
            return(bm_BinaryTransformation(data, threshold, doFiltering))
          })

setMethod('bm_BinaryTransformation', signature('numeric'),
          function(data, threshold, doFiltering = FALSE)
          {
            data <- as.data.frame(data)
            return(bm_BinaryTransformation(data, threshold, doFiltering))
          })

setMethod('bm_BinaryTransformation', signature('array'),
          function(data, threshold, doFiltering = FALSE)
          {
            if (length(dim(data)) == length(dim(threshold))) {
              if (sum(dim(data)[-1] != dim(threshold)[-1]) > 0) {
                stop("data and threshold dimensions mismatch")
              }
            } else {
              if (sum(dim(data)[-1] != dim(threshold)) > 0) {
                stop("data and threshold dimensions mismatch")
              }
            }
            
            if (doFiltering) {
              return(sweep(data, 2:length(dim(data)), threshold, .convertBin.array.filt))
            } else {
              return(sweep(data, 2:length(dim(data)), threshold, .convertBin.array))
            }
          })

setMethod('bm_BinaryTransformation', signature('RasterLayer'),
          function(data, threshold, doFiltering = FALSE)
          {
            if(!is.na(threshold)){
              if (doFiltering) {
                return(reclassify(data, c(-Inf, threshold, 0)))
              } else {
                return(reclassify(data,c(-Inf, threshold, 0, threshold, +Inf, 1)))
              }
            } else{ ## return a empty map (NA everywhere)
              return(reclassify(data, c(-Inf, Inf, NA)))
            }
          })

setMethod('bm_BinaryTransformation', signature('RasterStack'),
          function(data, threshold, doFiltering = FALSE)
          {
            if(length(threshold) == 1) {
              threshold <- rep(threshold, raster::nlayers(data))
            }
            StkTmp <- raster::stack()
            for (i in 1:raster::nlayers(data)) {
              ras = bm_BinaryTransformation(raster::subset(data, i, drop = TRUE), threshold[i], doFiltering)
              StkTmp <- raster::addLayer(StkTmp, ras)
            }
            names(StkTmp) <- names(data)
            return(StkTmp)
          })

setMethod('bm_BinaryTransformation', signature('RasterBrick'),
          function(data, threshold, doFiltering = FALSE)
          {
            data <- raster::stack(data, RAT = FALSE)
            return(bm_BinaryTransformation(data, threshold, doFiltering))
          })

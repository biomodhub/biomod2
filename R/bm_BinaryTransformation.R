###################################################################################################
##' @name bm_BinaryTransformation
##' @author Wilfried Thuiller, Damien Georges
##' 
##' @title Convert probability values into binary values using a predefined threshold
##' 
##' @description This function allows to convert probability (not necessary between \code{0} and 
##' \code{1}) values into binary presence-absence (\code{0} or \code{1}) values according to a 
##' predefined threshold (see \href{bm_BinaryTransformation.html#details}{Details}).
##' 
##' @param data a \code{vector}, a \code{matrix} or \code{data.frame}, a 
##' \code{\link[raster:raster]{raster}} or \code{\link[raster:stack]{RasterStack}} containing the 
##' data to be converted
##' @param threshold a \code{numeric} corresponding to the threshold used to convert the given 
##' data
##' @param doFiltering (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether filtered data should be returned, or binary one 
##' (see \href{bm_BinaryTransformation.html#details}{Details})
##' 
##' 
##' @return 
##' 
##' An object of the same class than \code{data} and containing either binary (\code{0} or  
##' \code{1}) values, or filtered values.
##'
##'
##' @details
##' 
##' If \code{data} is a \code{vector} or \code{\link[raster:raster]{raster}}, \code{threshold} 
##' should be a single \code{numeric} value. \cr
##' If \code{data} is a \code{matrix}, \code{data.frame} or 
##' \code{\link[raster:stack]{RasterStack}}, \code{threshold} should be a \code{vector} containing 
##' as many values as the number of columns or layers contained in \code{data}. If only one 
##' \code{numeric} value is given, the same threshold will be applied to all columns or layers. 
##' \cr \cr
##' 
##' If \code{doFiltering = FALSE}, binary (\code{0} or \code{1}) values are returned. \cr 
##' If \code{doFiltering = TRUE}, values will be \emph{filtered} according to \code{threshold}, 
##' meaning that :
##' \itemize{
##'   \item \code{data < threshold} will return 0
##'   \item \code{data >= threshold } will return the actual values of \code{data} (not 
##'   transformed in \code{1})
##' }
##' 
##'
##' @keywords convert, threshold, binary, filter
##' 
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}
##'
##'
##' @examples
##' 
##' ## generate 'probability' vector
##' xx <- rnorm(50, 10)
##' 
##' yy <- bm_BinaryTransformation(data = xx, threshold = 10)
##' yy_filt <- bm_BinaryTransformation(data = xx, threshold = 10, doFiltering = TRUE)
##' 
##' cbind(xx, yy, yy_filt)
##'
##' 
##' @importFrom raster stack subset nlayers addLayer reclassify
##' 
##' @export
##' 
##' 
###################################################################################################

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
              threshold <- rep(threshold, nlayers(data))
            }
            StkTmp <- stack()
            for (i in 1:nlayers(data)) {
              ras = bm_BinaryTransformation(subset(data, i, drop = TRUE), threshold[i], doFiltering)
              StkTmp <- addLayer(StkTmp, ras)
            }
            names(StkTmp) <- names(data)
            return(StkTmp)
          })

setMethod('bm_BinaryTransformation', signature('RasterBrick'),
          function(data, threshold, doFiltering = FALSE)
          {
            data <- stack(data, RAT = FALSE)
            return(bm_BinaryTransformation(data, threshold, doFiltering))
          })


###################################################################################################

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
    return(x >= y)
  } else {
    return(rep(NA, length(x)))
  }
}

.convertBin.array.filt = function(x, y) {
  if (!any(is.na(x))) {
    return(ifelse(x >= y, x, 0))
  } else {
    return(rep(NA, length(x)))
  }
}

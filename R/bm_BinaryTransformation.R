###################################################################################################
##' @name bm_BinaryTransformation
##' @author Wilfried Thuiller, Damien Georges
##' 
##' @title Convert probability values into binary values using a predefined threshold
##' 
##' @description This internal \pkg{biomod2} function allows to convert probability (not necessary 
##' between \code{0} and \code{1}) values into binary presence-absence (\code{0} or \code{1}) values 
##' according to a predefined threshold (see Details).
##' 
##' 
##' @param data a \code{vector}, a \code{matrix}, \code{data.frame}, or a 
##' \code{\link[terra:rast]{SpatRaster}} containing the data to be converted
##' @param threshold a \code{numeric} or a \code{vector} of \code{numeric} corresponding to 
##' the threshold used to convert the given data
##' @param do.filtering (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether filtered data should be returned, or binary one 
##' (see Details)
##' 
##' @return 
##' 
##' An object of the same class than \code{data} and containing either
##'  binary (\code{0} or \code{1}) values, or filtered values.
##'
##' @details
##' 
##' If \code{data} is a \code{vector}, \code{threshold} should be a single 
##' \code{numeric} value. \cr
##' If \code{data} is a \code{matrix}, \code{data.frame} or 
##' \code{\link[terra:rast]{SpatRaster}}, \code{threshold} should be a
##' \code{vector} containing as many values as the number of columns or 
##' layers contained in \code{data}. If only one \code{numeric} value is given,
##'  the same threshold will be applied to all columns or layers. 
##' \cr \cr
##' 
##' If \code{do.filtering = FALSE}, binary (\code{0} or \code{1}) values are returned. \cr 
##' If \code{do.filtering = TRUE}, values will be \emph{filtered} according to \code{threshold}, 
##' meaning that :
##' \itemize{
##'   \item \code{data < threshold} will return 0
##'   \item \code{data >= threshold } will return the actual values of \code{data} (not 
##'   transformed in \code{1})
##' }
##' 
##' @keywords convert threshold binary filter
##' 
##' @seealso \code{\link{BIOMOD_Projection}}, \code{\link{BIOMOD_EnsembleForecasting}}
##' @family Secundary functions
##'
##' 
##' @examples
##' 
##' ## Generate a 0-1000 vector (normal distribution)
##' vec.d <- rnorm(100, 500, 100)
##' 
##' ## From continuous to binary / filtered vector
##' vec.d_bin <- bm_BinaryTransformation(data = vec.d, threshold = 500)
##' vec.d_filt <- bm_BinaryTransformation(data = vec.d, threshold = 500, do.filtering = TRUE)
##' cbind(vec.d, vec.d_bin, vec.d_filt)
##'
##' 
##' @importFrom terra rast nlyr classify subset
##' @export
##' 
##' 
###################################################################################################


## generic methods ----------------------------------------------------------

setGeneric("bm_BinaryTransformation",
           function(data, threshold, do.filtering = FALSE)
           {
             standardGeneric("bm_BinaryTransformation")
           })

## data.frame methods ----------------------------------------------------------
##' @rdname bm_BinaryTransformation
##' @export
##'

setMethod('bm_BinaryTransformation', signature('data.frame'),
          function(data, threshold, do.filtering = FALSE)
          {
            if (length(threshold) > 1) {
              if (length(threshold) != ncol(data)) {
                stop("data and threshold dimensions mismatch")
              } else {
                if (do.filtering) {
                  return(
                    sweep(data, 2, threshold, .convert_bin.array.filt)
                    )
                } else {
                  return(as.data.frame(
                    sweep(data, 2, threshold, .convert_bin.array)
                    ))
                }
              }
            } else {
              if (is.numeric(threshold) && !is.na(threshold)) { #second condition unnecessary?
                data <- data.matrix(data)
                out <- as.data.frame(
                  .convert_bin.matrix(data, threshold, do.filtering)
                )
                colnames(out) <- colnames(data)
                return(out)
              } else { ## return NAs
                return(as.data.frame(
                  matrix(NA, ncol = ncol(data), nrow = nrow(data)
                              , dimnames = dimnames(data))
                ))
              }
            }
          })

## matrix methods ----------------------------------------------------------
##' @rdname bm_BinaryTransformation
##' @export
##'

setMethod('bm_BinaryTransformation', signature('matrix'),
          function(data, threshold, do.filtering = FALSE)
          {
            data <- as.data.frame(data)
            return(data.matrix(
              bm_BinaryTransformation(data, threshold, do.filtering)
              ))
          })

## numeric methods ----------------------------------------------------------
##' @rdname bm_BinaryTransformation
##' @export
##'

setMethod('bm_BinaryTransformation', signature('numeric'),
          function(data, threshold, do.filtering = FALSE)
          {
            data <- as.data.frame(data)
            return(unlist(
              bm_BinaryTransformation(data, threshold, do.filtering)
              ))
          })


## SpatRaster methods ----------------------------------------------------------
##' @rdname bm_BinaryTransformation
##' @importFrom terra `add<-`
##' @export
##'

setMethod('bm_BinaryTransformation', signature('SpatRaster'),
          function(data, threshold, do.filtering = FALSE)
          {
            if(length(threshold) == 1) {
              threshold <- rep(threshold, nlyr(data))
            }
            StkTmp <- rast()
            for (i in 1:nlyr(data)) {
              if(!is.na(threshold[i])){
                if (do.filtering) {
                  ras <- classify(subset(data, i),
                                  matrix(c(-Inf, threshold[i], 0),
                                         ncol = 3, byrow = TRUE),
                                  right = FALSE)
                } else {
                  ras <- classify(subset(data, i), 
                                  matrix(c(-Inf, threshold[i], 0, 
                                           threshold[i], +Inf, 1),
                                         ncol = 3, byrow = TRUE),
                                  right = FALSE)
                }
              } else { ## return a empty map (NA everywhere)
                ras <- classify(subset(data, i),
                                matrix(c(-Inf, Inf, NA),
                                       ncol = 3, byrow = TRUE))
              }
              add(StkTmp) <- ras
            }
            names(StkTmp) <- names(data)
            return(StkTmp)
          })

### .convert_bin.matrix --------------------------------------------------------

.convert_bin.matrix = function(data, threshold, do.filtering = FALSE) {
  ind.0 = t(t(data) < threshold)
  data[ind.0] <- 0
  if (!do.filtering) { 
    data[!ind.0] <- 1 
  }
  if (ncol(data) == 1) { 
    return(data[, 1]) 
  } else { 
    return(data) 
  }
}

.convert_bin.array = function(x, y) {
  return(ifelse(x >= y, 1, 0))
}

.convert_bin.array.filt = function(x, y) {
  x[x <= y] <- 0
  return(x)
}

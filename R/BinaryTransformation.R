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
setGeneric("BinaryTransformation",
           function(data, threshold){
             standardGeneric("BinaryTransformation")
           })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, data.frame-method
setMethod('BinaryTransformation', signature(data='data.frame'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='data.frame')")
  	FUN2 <- function(x,y){
  		moa <- apply((x>y),2,as.integer)
  		if(ncol(moa)==1) return(moa[,1])
  		else return(moa)
  	}
    if(is.numeric(threshold)){
      return(sweep(data.matrix(data), 2, threshold, FUN2))
    } else { ## return NAs
      return( matrix(NA, ncol=ncol(data), nrow=nrow(data), dimnames=dimnames(data)) )
    }

  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, matrix-method
setMethod('BinaryTransformation', signature(data='matrix'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='matrix')")
    data <- as.data.frame(data)
    return(BinaryTransformation(data, threshold))
  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, numeric-method
setMethod('BinaryTransformation', signature(data='numeric'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='numeric')")
    data <- as.data.frame(data)
    return(BinaryTransformation(data, threshold))
  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, array-method
setMethod('BinaryTransformation', signature(data='array'),
          function(data, threshold)
          {
            # cat("\n*** in setMethod('BinaryTransformation', signature(data='array')")
            if(length(dim(data)) == length(dim(threshold))){
              if(sum( dim(data)[-1] != dim(threshold)[-1] ) > 0 ){
                stop("data and threshold dimensions mismatch")
              }
            } else{
              if(sum( dim(data)[-1] != dim(threshold) ) > 0 ){
                stop("data and threshold dimensions mismatch")
              }
            }

            return(sweep(data,2:length(dim(data)),threshold,
                         function(x,y) {
                           if(!any(is.na(x))){
                             return(x>y)
                            } else {
                             return(rep(NA,length(x)) )}
                           } ))
          })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, RasterLayer-method
setMethod('BinaryTransformation', signature(data='RasterLayer'),
  function(data, threshold)
  {
    # cat("\n*** in setMethod('BinaryTransformation', signature(data='RasterLayer')")
    if(!is.na(threshold)){
      return(reclassify(data,c(-Inf,threshold,0, threshold,+Inf,1)))
    } else{ ## return a empty map (NA everywhere)
      return(reclassify(data,c(-Inf,Inf,NA)))
    }

  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, RasterStack-method
setMethod('BinaryTransformation', signature(data='RasterStack'),
  function(data, threshold)
  {
    if(length(threshold) == 1){
      threshold <- rep(threshold, raster::nlayers(data))
    }
    return(calc(data, function(x){x >= threshold}))
#     ## old version
#     StkTmp <- raster::stack()
#     for(i in 1:raster::nlayers(data)){
#       StkTmp <- raster::addLayer(StkTmp, BinaryTransformation(raster::subset(data,i,drop=TRUE), threshold[i]))
#     }
#     names(StkTmp) <- names(data)
#     return(StkTmp)
  })

##' @rdname BinaryTransformation-methods
##' @aliases BinaryTransformation, RasterBrick-method
setMethod('BinaryTransformation', signature(data='RasterBrick'),
  function(data, threshold)
  {
    data <- raster::stack(data)
    return(BinaryTransformation(data, threshold))
  })

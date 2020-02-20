##' @include BiomodClass.R
##' @name FilteringTransformation
##' @aliases FilteringTransformation
##' @aliases FilteringTransformation-methods
##' @aliases FilteringTransformation,data.frame-method
##' @aliases FilteringTransformation,matrix-method
##' @aliases FilteringTransformation,numeric-method
##' @aliases FilteringTransformation,array-method
##' @aliases FilteringTransformation,RasterBrick-method
##' @aliases FilteringTransformation,RasterLayer-method
##' @aliases FilteringTransformation,RasterStack-method
##' 
##' @title Convert species' probability of occurrence into binary 
##' presence-absence data using a predefined threshold
##' 
##' @description
##' Function that converts an object containing probability values into 
##' a filtered object according to a pre-defined threshold(s).
##' 
##' 
##' @param data a numeric vector, a \code{matrix}, a \code{data.frame}, 
##' a \code{RasterLayer} or a \code{RasterStack} containing the data to 
##' be converted
##' @param threshold a numeric value or a vector containing the threshold
##' to be used for converting data.
##' 
##' @details
##' If data is a vector or a raster object, then the threshold should be a
##' numeric value. If data is matrix,dataframe or rasterStack, then the
##' threshold should have, in theory, as many values as the number of
##' columns or layers to transform.
##' In the particular case that the data to convert is a 
##' \code{matrix}/\code{data.frame} with several columns or a 
##' \code{RasterStack} with several layers and the threshold is a single
##' numeric value, the same threshold will be applied to all columns 
##' (resp. layers).  
##' 
##' @return 
##' An object of the same class than \code{data} with the values of data
##' if superior to \code{threshold} and 0 if not.
##' 
##' @author Wilfried Thuiller, Damien Georges
##'
##' @examples
##' xx <- rnorm(50,10)
##' yy <- FilteringTransformation(xx, 10)
##' 
##' cbind(xx,yy)
setGeneric("FilteringTransformation",
           function(data, threshold){
             standardGeneric("FilteringTransformation")
           })

setMethod('FilteringTransformation', signature(data='data.frame'),
  function(data, threshold)
  {
    data <- data.matrix(data)
    data[t(t(data)<threshold)] <-0

    ## check if some thresolds are NAs
    if(any(is.na(threshold))){
      data[,is.na(threshold)] <- NA
    }
    if(ncol(data)==1) data <- data[,1]
  	return(data)

  })

setMethod('FilteringTransformation', signature(data='matrix'),
  function(data, threshold)
  {
    data <- as.data.frame(data)
    return(FilteringTransformation(data, threshold))
  })

setMethod('FilteringTransformation', signature(data='numeric'),
  function(data, threshold)
  {
    data <- as.data.frame(data)
    return(FilteringTransformation(data, threshold))
  })

setMethod('FilteringTransformation', signature(data='array'),
          function(data, threshold)
          {
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
                           if(!is.na(x)){
                             return(ifelse(x>y,x,0))
                           } else {
                             return(rep(NA,length(x)) )}
                         }))
          })


setMethod('FilteringTransformation', signature(data='RasterLayer'),
  function(data, threshold)
  {
    if(!is.na(threshold)){
      return(reclassify(data,c(-Inf,threshold,0)))
    } else{ ## return a empty map (NA everywhere)
      return(reclassify(data,c(-Inf,Inf,NA)))
    }
  })

setMethod('FilteringTransformation', signature(data='RasterStack'),
  function(data, threshold)
  {
    if(length(threshold) == 1){
      threshold <- rep(threshold, raster::nlayers(data))
    }
    StkTmp <- raster::stack()
    for(i in 1:raster::nlayers(data)){
      StkTmp <- raster::addLayer(StkTmp, FilteringTransformation(raster::subset(data,i,drop=TRUE), threshold[i]))
    }
    names(StkTmp) <- names(data)
    return(StkTmp)
  })

setMethod('FilteringTransformation', signature(data='RasterBrick'),
  function(data, threshold)
  {
    data <- raster::stack(data, RAT=FALSE)
    return(FilteringTransformation(data, threshold))
  })

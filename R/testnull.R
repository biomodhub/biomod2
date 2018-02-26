.testnull <-
function(object, Prev = 0.5 , dat){
  if( is.finite(object$deviance) & is.finite(object$null.deviance)){
    if(object$deviance != object$null.deviance){
      if(inherits(dat,'Raster')){
        pred <- predict(object = dat, model=object, type="response")
      } else{
        pred <- predict(object, dat, type="response")
      }
    }
  }
  
  if(!exists('pred')){
    if(inherits(dat,'Raster')){
      pred <- raster::subset(dat,1,drop=TRUE)
      if(Prev < 0.5) pred <- reclassify(x=pred, rcl=c(-Inf,Inf,0))
      if(Prev >= 0.5) pred <- reclassify(x=pred, rcl=c(-Inf,Inf,1))
    } else{
      if(Prev < 0.5) pred <- rep(0, nrow(dat))
      if(Prev >= 0.5) pred <- rep(1, nrow(dat))      
    }
    
  }
    return(pred)
}


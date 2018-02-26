# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Variables Importance tools
# Damien G. - april 2013
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# variables_importance (main function) =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
#
# AIM : give information a


variables_importance <- function(model, data, method="full_rand", nb_rand=1, ...){
  out <- list()
  args <- .variables_importance.check.args(model=model, data=data, method=method, ...)
  
  # test prediction is computable
  ref <- try(predict(args$model, args$data))
  if(inherits(ref,"try-error")) stop("Unable to make model prediction")
  
  # make randomisation
  out$mat <- matrix(0,nrow=length(args$variables), ncol=nb_rand, dimnames=list(args$variables, paste('rand',1:nb_rand,sep="")))
  
  for(r in 1:nb_rand){
    for(v in args$variables){
      data_rand <- randomise_data(args$data,v,method)
      shuffled.pred <- predict(args$model, data_rand)
      out$mat[v,r] <- 1 - max(round(cor(x=ref, y=shuffled.pred, use="pairwise.complete.obs", method="pearson"),digits=6),0,na.rm=T)
    }
  }
  
  class(out) <- "BIOMOD_variables_importances"
  return(out)
}


# variables_importance argument checking =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
.variables_importance.check.args <- function(...){
  args <- list(...)
  
  # test that input data is supported
  supported_models <- c("biomod2_model", "nnet", "rpart", "fda", "gam", "glm", "lm", "gbm", "mars", "randomForest")
  if(!inherits(args$model, supported_models)) stop("Model class unsuported")
  
  # test method is supported
  supported_methods <- c("full_rand")
  if(! args$method %in% supported_methods ) stop("Unknown method")
  
  # get variables names
  if(is.null(args$variables)) args$variables <- colnames(args$data)
  
  return(args)
}

# data_set shuffling =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
randomise_data <- function(data,variable,method){
  if(method=='full_rand'){
    return(full_suffling(data,variable))
  }
  
}


# full shuffling =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
full_suffling <- function(x,id=NULL){
  if(! (is.vector(x) | is.matrix(x) | is.data.frame(x)) ) stop("x must be a 1 or 2 dimention odject")
  
  ## trick to ensure that the randomisation is correctly done was not the case before ##
  set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6"))*1000000)
  out <- NULL
  if(is.null(id)){
    out <- x[sample.int(length(x))]
  } else{
    out <- x
    for(idd in id){
      out[,idd] <- out[sample.int(nrow(x)),idd]
    }
  }
  return(out)
}


##### TEST #####
# setwd("~/__BIOMOD__/DevComputing/")
# x <- rbinom(n=100,size=1,prob=0.3)
# y <- rnorm(100)
# z <- rnorm(100)
# data <- as.data.frame(cbind(x,y,z))
# 
# myGLM <- glm(x~y*z,family='binomial')
# 
# VI <- variables_importance(myGLM, data, nb_rand=10)
# 
# ##
# library(biomod2)
# 
# load("GuloGulo/GuloGulo.test.models.out")
# xx <- BIOMOD_LoadModels(GuloGulo.test.models.out,models='MAXENT.Phillips')
# 
# VI <- variables_importance(get(xx[1]), getModelsInputData(GuloGulo.test.models.out,'expl.var'), nb_rand=10)




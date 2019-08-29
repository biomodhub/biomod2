##' @name variables_importance
##' @aliases variables_importance
##' @title Variables importance calculation
##' @description
##' This function will return a variable importance value for 
##' each variable involved within your model.
##' 
##' @param model the model you want to study variables importance
##'   (one of the models supported within biomod2, ensemble models
##'    are also supported
##' @param data the \code{data.set} on which you want to perform
##'   analyses
##' @param method the randomisation method (only 'full_rand'
##'   available so far)
##' @param nb_rand the number of permutation done for each
##'   variable
##' @param ... additional args (not implemented yet)
##' 
##' @details
##' It's more or less base on the same principle than 
##' \code{\link[randomForest]{randomForest}} variables importance
##' algorithm. The principle is to shuffle a single variable of
##' the given data. Make model prediction with this 'shuffled'
##' data.set. Then we compute a simple correlation (Pearson's by
##' default) between references predictions and the 'shuffled' 
##' one. The return score is 1-cor(pred_ref,pred_shuffled). The
##' highest the value, the more influence the variable has on the
##' model. A value of this 0 assumes no influence of that variable
##' on the model. Note that this technique does not account for
##' interactions between the variables.
##' 
##' @return a \code{list} of class "BIOMOD_variables_importances"
##' which contains:
##' 
##'  - mat: a \code{data.frame} containing variables importance
##'    scores for each permutation run.
##' 
##' @author Damien Georges
##' @seealso \code{\link[biomod2]{randomise_data}}, 
##'   \code{\link[biomod2]{full_suffling}}
##'   
##' @keywords suffle
##' @keywords random
##' @keywords importance
##' 
##' @examples
##' xx <- 
##'   data.frame( 
##'     a = sample(c(0, 1), 100, replace = TRUE),
##'     b = rnorm(100),
##'     c = 1:100
##'   )
##'   
##' mod <- glm(a ~ b + c, data = xx)
##' 
##' variables_importance(
##'   model = mod, 
##'   data = xx[, c('b', 'c')], 
##'   method = "full_rand", 
##'   nb_rand = 3
##' )
##' 
variables_importance <- 
  function(
    model, 
    data, 
    method = "full_rand", 
    nb_rand = 1, 
    ...
  ){
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
  if(!inherits(args$model, supported_models)) stop("Model class unsupported")

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
  if(! (is.vector(x) | is.matrix(x) | is.data.frame(x)) ) stop("x must be a 1 or 2 dimension odject")

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




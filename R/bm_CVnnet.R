###################################################################################################
##' @name bm_CVnnet
##' @author Damien Georges
##' 
##' @title Cross-validation for Neural Networks
##' 
##' @description 
##' 
##' This internal \pkg{biomod2} function allows the user to compute cross-validation for neural 
##' networks in ANN model (see \code{\link[nnet]{nnet}} and \code{\link{BIOMOD_Modeling}}).
##' 
##' @param Input complete dataset with explanatory variables
##' @param Target calibration dataset with observed presence / absence
##' @param size in modeling options : ANN$size
##' @param decay in modeling options ANN$decay
##' @param maxit in modeling options ANN$maxit
##' @param nbCV in modeling options ANN$nbCV
##' @param W weights over calibration lines
##' 
##' 
##' @return  
##' 
##' A \code{data.frame} containing the following elements :
##' \itemize{
##'   \item{\code{Size} : }{the size}
##'   \item{\code{Decay} : }{the decay value}
##'   \item{\code{AUC} : }{the corresponding Area Under Curve}
##' }
##' 
##' 
##' 
##' @seealso \code{\link{BIOMOD_ModelingOptions}}, \code{\link{BIOMOD_Modeling}}, 
##' \code{\link{bm_SampleBinaryVector}}
##' 
##' 
##' @keywords neural networks, cross-validation
##' 
##' 
##' @importFrom pROC auc roc
##' @importFrom nnet nnet
##' 
##' @export
##' 
##'
###################################################################################################


bm_CVnnet = function(Input,
                     Target,
                     size = c(2, 4, 6, 8),
                     decay = c(0.001, 0.01, 0.05, 0.1),
                     maxit = 200,
                     nbCV = 5,
                     W = NULL)
{
  ## Prepare output table
  Eval = data.frame(matrix(0, ncol = 3, nrow = 16, dimnames = list(NULL, c("Size", "Decay", "AUC"))))
  Eval[, 1] = rep(size, 4)
  Eval[, 2] = rep(decay, each = 4)
  
  
  for (i in 1:nbCV) {
    set.seed(555)
    Samp = bm_SampleBinaryVector(ref = Target, ratio = 0.5)
    
    W.tmp = ifelse(is.null(W), rep(1, length(Target)), W)
    
    
    Eval[, 3] = Eval[, 3] + apply(Eval[, 1:2],
                                  1,
                                  Samp,
                                  Target,
                                  Input,
                                  W.tmp,
                                  FUN = function(x, Samp, Target, Input, W.tmp) {
                                    nn = nnet(eval(parse(text = paste("Target[Samp$calibration]",
                                                                      paste(.scope_expSyst(Input[1:10, , drop = FALSE], "GBM"), collapse = "")))),
                                              data = Input[Samp$calibration, , drop = FALSE],
                                              weights = W.tmp[Samp$calibration],
                                              size = x[1],
                                              decay = x[2],
                                              maxit = maxit,
                                              trace = FALSE)
                                    AUC = roc(Target[Samp$evaluation],
                                              as.numeric(predict(nn, Input[Samp$evaluation, , drop = FALSE])),
                                              levels = c(0, 1),
                                              direction = '<')
                                    AUC <- as.numeric(auc(AUC))
                                    return(AUC)
                                  })
  }
  
  Eval[, 3] = Eval[, 3] / nbCV
  z = which.max(Eval[, 3])
  
  return(Eval[z, 1:2])
}

.CV.nnet = function(Input, Target, size=c(2,4,6, 8), decay=c(0.001, 0.01, 0.05, 0.1), maxit=200, nbCV=5, W=NULL){
#   require(pROC, quietly=T)

  Eval = data.frame(matrix(0, ncol=3, nrow=16, dimnames=list(NULL, c("Size", "Decay", "AUC"))))
  Eval[,1] = rep(size,4)
  Eval[,2] = rep(decay, each=4)
  for(i in 1:nbCV){
      set.seed(555)
      Samp = SampleMat2(Target, 0.5)
      
      if(is.null(W)){
          Eval[,3] = Eval[,3] + apply(Eval[,1:2], 1, Samp, Target, Input, FUN=function(x, Samp, Target, Input){
            nn = nnet(eval(parse(text = paste("Target[Samp$calibration]",
                  paste(.scopeExpSyst(Input[1:10, ,drop=FALSE], "GBM"), collapse = "")))),data=Input[Samp$calibration, ,drop=FALSE],
                  size = x[1], decay = x[2], maxit = maxit, trace = FALSE)
            AUC = as.numeric(pROC::auc(pROC::roc(Target[Samp$evaluation], predict(nn, Input[Samp$evaluation,,drop=FALSE]))))
            return(AUC)
          })
      } else{
          Eval[,3] = Eval[,3] + apply(Eval[,1:2], 1, Samp, Target, Input, W, FUN=function(x, Samp, Target, Input, W){
          nn = nnet(eval(parse(text = paste("Target[Samp$calibration]",
                  paste(.scopeExpSyst(Input[1:10, ,drop=FALSE], "GBM"), collapse = "")))),data=Input[Samp$calibration, ,drop=FALSE],
                  weights=W[Samp$calibration], size = x[1], decay = x[2], maxit = maxit, trace = FALSE)
            AUC = as.numeric(pROC::auc(pROC::roc(Target[Samp$evaluation], predict(nn, Input[Samp$evaluation,,drop=FALSE]))))
            return(AUC)
          })
      }
  }
  Eval[,3] = Eval[,3]/nbCV
  z =which.max(Eval[,3])
  return(Eval[z, 1:2])
}

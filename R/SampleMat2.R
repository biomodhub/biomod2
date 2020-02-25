##' @name SampleMat2
##' @title Sample binary vector
##' 
##' @description
##' \code{SampleMat2} is an internal \pkg{biomod2} function that can help
##' user to sample a binary vector keeping the same proportion of 0s and 1s
##' than in the initial vector.
##' 
##' @param ref a binary vector
##' @param ratio the proportion of \code{ref} to sample
##' @param as.logi logical, if FALSE (default) id of cell will be return;
##'   if TRUE, logical vector of same length than ref will be return
##'   
##' @details
##' This function can be useful to help users to select a part of initial
##' dataset that will be only kept for all validation procedures. 
##' 
##' @return 
##' A list of 2 elements is returned :
##' 
##' - `calibration` Ids of cells selected for calibration (1st sample)
##' 
##' - `evaluation` Ids of cells selected for evaluation (1st sample
##'   complementary)
##'   
##' @author Damien Georges
##' 
##' @seealso \code{\link[biomod2]{BIOMOD_FormatingData}}
##' @keywords models
##' @keywords formula
##' @keywords options
##' 
##' @export
##' 
##' @examples
##' a <- sample(c(0,1),100, replace=TRUE)
##' SampleMat2(ref=a, ratio=0.7)
##' 
SampleMat2 <- function(
  ref, 
  ratio, 
  as.logi = FALSE
){
  # set a new random seed to ensure that sampling is random (issue when CTA is involved and seed needs to be set to a fix number)
  set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6"))*1000000)
  
  ntot <- length(ref)
  npres<- sum(ref)    
  ncal <- ceiling(ntot*ratio)

  pres <- sample(which(ref==1), ceiling(npres*ratio))
  absc <- sample(which(ref==0), ncal-length(pres))
  
  if(as.logi){
    calib <- rep(FALSE, ntot)
    calib[c(pres,absc)] <- TRUE
    eval <- !calib
  } else{
    calib <- c(pres,absc)
    eval <- (1:ntot)[-c(pres,absc)]
  }
  
  return(list("calibration"=calib, "evaluation"=eval))
}


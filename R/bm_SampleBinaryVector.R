###################################################################################################
##' @name bm_SampleBinaryVector
##' @author Damien Georges
##' 
##' @title Sample binary vector
##' 
##' @description This internal \pkg{biomod2} function allows the user to sample a binary vector 
##' keeping the same proportion of \code{0} and \code{1} as the initial vector.
##' 
##' 
##' @param obs a \code{vector} containing binary values (either \code{0} or \code{1})
##' @param ratio a \code{numeric} between \code{0} and \code{1} corresponding to the proportion 
##' of \code{obs} values to sample
##' @param as.logical (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether a \code{vector} of \code{TRUE/FALSE} values should 
##' be returned (if \code{as.logical = TRUE}) or a \code{vector} containing the indices of 
##' \code{obs} elements to keep
##' @param seedval (\emph{optional, default} \code{NULL}) \cr 
##' An \code{integer} value corresponding to the new seed value to be set
##' 
##' 
##' @return 
##'  
##' A \code{list} containing to elements is returned :
##' \describe{
##'   \item{calibration}{IDs of elements selected for calibration}
##'   \item{evaluation}{IDs of elements selected for evaluation (complementary to the calibration 
##'   set)}
##' }
##'   
##'   
##' @keywords models formula options
##' 
##' 
##' @seealso \code{\link{bm_CVnnet}}
##' @family Secundary functions
##' 
##' 
##' @examples
##' 
##' ## Generate a binary vector
##' vec.a <- sample(c(0, 1), 100, replace = TRUE)
##' 
##' ## Generate calibration / evaluatation datasets
##' bm_SampleBinaryVector(obs = vec.a, ratio = 0.7)
##' 
##' 
##' @export
##' 
##' 
###################################################################################################


bm_SampleBinaryVector <- function(obs, ratio, as.logical = FALSE, seedval = NULL)
{
  ## Set a new random seed to ensure that sampling is random
  ## (issue when CTA is involved and seed needs to be set to a fix number)
  if (is.null(seedval)) {
    set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6")) * 1000000)
  }
  
  ntot <- length(obs)
  npres <- sum(obs)
  ncal <- ceiling(ntot * ratio)
  
  pres <- sample(which(obs == 1), ceiling(npres * ratio))
  absc <- sample(which(obs == 0), ncal - length(pres))
  
  if (as.logical) {
    calib <- rep(FALSE, ntot)
    calib[c(pres, absc)] <- TRUE
    eval <- !calib
  } else {
    calib <- c(pres, absc)
    eval <- (1:ntot)[-c(pres, absc)]
  }
  
  return(list("calibration" = calib, "evaluation" = eval))
}


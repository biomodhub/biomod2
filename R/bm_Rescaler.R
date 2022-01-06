###################################################################################################
##' @name bm_Rescaler
##' @author Damien Georges
##' 
##' @title Rescale values between \code{0} and \code{1} with a binomial GLM
##'
##' @description 
##' 
##' This internal \pkg{biomod2} function allows the user to rescale predicted values from \code{0} 
##' to \code{1} with a binomial GLM (see \href{bm_Rescaler.html#details}{Details}).
##'
##' @param dataToRescale a \code{vector} containing predicted values to be rescaled
##' @param name a \code{character} corresponding to model name to be used, must be among \code{GLM}, 
##' \code{GBM}, \code{GAM}, \code{CTA}, \code{ANN}, \code{SRE}, \code{FDA}, \code{MARS}, 
##' \code{RF}, \code{MAXENT.Phillips}, \code{MAXENT.Phillips.2}
##' @param doModel (\emph{optional, default} \code{FALSE}) \cr 
##' A \code{logical} value defining whether the scaling model should be re-created or not (if 
##' \code{FALSE}, existing scaling model will be loaded)
##' @param ref (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing observed data to create weights if none given
##' @param weights (\emph{optional, default} \code{NULL}) \cr 
##' A \code{vector} containing values to weight \code{dataToRescale} and influence prevalence
##' 
##' 
##' @return 
##' 
##' A \code{vector} containing rescaled values between \code{0} and \code{1}.
##'
##' @details
##' 
##' \emph{Please refer to parameter \code{rescal.all.models} in \code{\link{BIOMOD_Modeling}} 
##' function to get more information about this procedure.}
##' 
##' Note that scaling model are saved within 
##' \file{resp.name/}\code{name}\file{/models/scaling_models/} folder.
##'
##'
##' @keywords models, projection, rescale, glm
##' 
##'
##' @seealso \code{\link{BIOMOD_Modeling}}, \code{\link{bm_Projection}}, 
##' \code{\link{BIOMOD_Projection}}
##'
##'
##' @examples
##' 
##' ## generate a binary vector
##' a <- sample(c(0, 1), 100, replace = TRUE)
##'
##' ## biased drawing
##' BiasedDrawing <- function(x, m1 = 300, sd1 = 200, m2 = 700, sd2 = 200) {
##'   return(ifelse(x < 0.5, rnorm(1, m1, sd1), rnorm(1, m2, sd2)))
##' }
##' b <- sapply(a, BiasedDrawing)
##' b[which(b < 0)] <- 0
##' b[which(b > 1000)] <- 1000
##' 
##' b_rescaled <- bm_Rescaler(dataToRescale = b, name = "TEST", doModel = TRUE, ref = a)
##' plot(b, b_rescaled)
##' abline(a = 0, b = 0.001, lty = 2)
##'
##' 
###################################################################################################


bm_Rescaler <- function(dataToRescale,
                        name,
                        doModel = FALSE,
                        ref = NULL,
                        weights = NULL)
{
  ## get compressing format depending on OS
  compress.arg = ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  
  ## get output directory name
  name_dir = paste0(getwd(), "/", unlist(strsplit(name, '_'))[1], "/models/scaling_models/")
  name_out = paste0(name_dir, name, "_scaled")
  
  if (doModel) ## create the scaling model -----------------------------------
  {
    if (!file.exists(name_dir)) {
      dir.create(name_dir, showWarnings = FALSE, recursive = TRUE)
    }
    
    Rescaling_GLM = .scaling_model(dataToRescale, ref, prevalence = 0.5, weights)
    save(Rescaling_GLM, file = name_out, compress = compress.arg) 
    
  } else { ## load the scaling model ------------------------------------------
    load(name_out)
  }
  
  ## make the scaling prediction ----------------------------------------------
  if (!inherits(dataToRescale, "Raster")) {
    RescaledData <- predict(Rescaling_GLM, data.frame(pred = as.numeric(dataToRescale)), type = "response")
  } else {
    RescaledData <- predict(dataToRescale, model = Rescaling_GLM, type = 'response')
  }
  
  return(RescaledData)
}


###################################################################################################


.scaling_model <- function(dataToRescale, ref = NULL, ...)
{
  args <- list(...)
  prevalence <- args$prevalence
  weights <- args$weights
  
  ## if no weights given, some are created to rise the define prevalence ------
  if (is.null(weights) && ! is.null(prevalence)) {
    nbPres <- sum(ref, na.rm = TRUE)
    nbAbs <- length(ref) - nbPres
    weights <- rep(1, length(ref))
    
    if (nbAbs > nbPres) { # code absences as 1
      weights[which(ref > 0)] <- (prevalence * nbAbs) / (nbPres * (1 - prevalence))
    } else { # code presences as 1
      weights[which(ref == 0 | is.na(ref))] <- (nbPres * (1 - prevalence)) / (prevalence * nbAbs)
    }
    weights = round(weights[]) # to remove glm & gam warnings
    
  } else if (is.null(weights)) { ## only 1 weights vector ---------------------
    weights <- rep(1, length(ref))
  }
  
  ## define a glm to scale predictions from 0 to 1 ----------------------------
  scaling_model <- glm(ref ~ pred, data = data.frame(ref = as.numeric(ref), pred = as.numeric(dataToRescale))
                       , family = binomial(link = probit), x = TRUE, weights = weights)
  
  return(scaling_model)
}

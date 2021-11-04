.Rescaler5 <- function(dataToRescale,
                       ref = NULL,
                       name,
                       original = FALSE,
                       weights = NULL)
{
  ## get compressing format depending on OS
  compress.arg = ifelse(.Platform$OS.type == 'windows', 'gzip', 'xz')
  
  ## get output directory name
  name_dir = paste0(getwd(), "/", unlist(strsplit(name, '_'))[1], "/models/scaling_models/")
  name_out = paste0(name_dir, name, "_scaled")
  
  if (original) ## create the scaling model -----------------------------------
  {
    if (!file.exists(name_dir)) {
      dir.create(name_dir, showWarnings = FALSE)
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
    weights <- rep(1,length(ref))
  }
  
  ## define a glm to scale predictions from 0 to 1 ----------------------------
  scaling_model <- glm(ref ~ pred, data = data.frame(ref = as.numeric(ref), pred = as.numeric(dataToRescale))
                       , family = binomial(link = probit), x = TRUE, weights = weights)
  
  return(scaling_model)
}

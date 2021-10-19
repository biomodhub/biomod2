
## AUTOMATIC WEIGHTS (BIOMOD_FormatingData) -------------------------------------------------------

.automatic_weights_creation <- function(resp, prev = 0.5, subset = NULL)
{
  if (is.null(subset)) { subset <- rep(TRUE, length(resp)) }
  
  nbPres <- sum(resp[subset], na.rm = TRUE)
  # The number of true absences + pseudo absences to maintain true value of prevalence
  nbAbsKept <- sum(subset, na.rm = T) - sum(resp[subset], na.rm = TRUE)
  Yweights <- rep(1, length(resp))
  
  if (nbAbsKept > nbPres) {
    # code absences as 1
    Yweights[which(resp > 0)] <- (prev * nbAbsKept) / (nbPres * (1 - prev))
  } else {
    # code presences as 1
    Yweights[which(resp == 0 | is.na(resp))] <- (nbPres * (1 - prev)) / (prev * nbAbsKept)
  }
  Yweights = round(Yweights[])
  Yweights[!subset] <- 0
  
  return(Yweights)
}


## TEST PARAMETERS (BIOMOD_ModelingOptions) -------------------------------------------------------

.fun_testIfIn <- function(test, objName, objValue, values)
{
  if (!(objValue %in% values)) {
    cat("\n", paste0(objName, " must be '", paste0(values[1:(length(values) -1)], collapse = "', '"), "' or '", values[length(values)], "'"))
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosNum <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    cat("\n", objName, "must be a numeric")
    test <- FALSE
  } else if (objValue < 0) {
    cat("\n", objName, "must be a positive numeric")
    test <- FALSE
  }
  return(test)
}

.fun_testIf01 <- function(test, objName, objValue)
{
  test <- .fun_testIfPosNum(test, objName, objValue)
  if (test && objValue > 1) {
    cat("\n", objName, "must be a 0 to 1 numeric")
    test <- FALSE
  }
  return(test)
}

.fun_testIfPosInt <- function(test, objName, objValue)
{
  if (!is.numeric(objValue)) {
    cat("\n", objName, "must be a integer")
    test <- FALSE
  } else if (objValue < 0 | objValue %% 1 != 0) {
    cat("\n", objName, "must be a positive integer")
    test <- FALSE
  }
  return(test)
}

## PRINTS (BIOMOD_ModelingOptions) ----------------------------------------------------------------
.print.formula <- function(formula = NULL)
{
  ifelse(length(formula) < 1, 'NULL', paste(formula[2], formula[1], formula[3]))
}

.print.control <- function(ctrl)
{
  out <-  paste0(names(ctrl)[1], " = ", ctrl[[1]])
  if (length(ctrl) > 1)
  {
    i = 2
    while (i <= length(ctrl)) {
      if (is.list(ctrl[[i]])) {
        out <- c(out, paste0(", ", names(ctrl)[i], " = list(",
                             paste0(names(ctrl[[i]]), "=", unlist(ctrl[[i]]), collapse = ", "),
                             ")"))
        i <- i + 1
      } else {
        out <- c(out, paste0(", ", names(ctrl)[i], " = ", ctrl[[i]]))
        i <- i + 1
      }
    }
  }
  return(out)
}
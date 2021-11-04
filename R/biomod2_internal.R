
## BIOMOD SPECIFIC CAT ----------------------------------------------------------------------------

.bmCat <- function(x = NULL, ...)
{
  if (is.null(x))
  {
    cat("\n")
    cat(paste(rep("-=", round(.Options$width / 2)), collapse = ""))
    cat("\n")
  } else
  {
    x.length = nchar(x) + 2
    y.length = (.Options$width - x.length) / 2
    cat("\n")
    cat(paste(rep("-=", round(y.length / 2)), collapse = "")
        , x
        , paste(rep("-=", round(y.length / 2)), collapse = ""))
    cat("\n")
  }
}


## CLEVER CUT (plot BIOMOD_FormatingData) ---------------------------------------------------------
## used in biomod2_classes_1.R file

.clever_cut <- function(x)
{
  nb_col = ceiling(sqrt(x))
  nb_row = ceiling(x / nb_col)
  return(c(nb_row, nb_col))
}


## AUTOMATIC WEIGHTS (BIOMOD_FormatingData) -------------------------------------------------------
## used in biomod2_classes_1.R file

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
## used in biomod2_classes_1.R file

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
## used in biomod2_classes_1.R file

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


## TEMPLATES TO PREDICT MODELS (BIOMOD_Modeling) --------------------------------------------------
## used in biomod2_classes_4.R file

# Fuction to get variables ranges
get_var_type <- function(data) { return(sapply(data, class)) }

get_var_range <- function(data)
{
  get_range <- function(x) {
    if (is.numeric(x)) {
      return(c(min = min(x, na.rm = T), max = max(x, na.rm = T)))
    } else if (is.factor(x)) {
      return(levels(x))
    }
  }
  xx <- lapply(data, get_range)
  names(xx) <- names(data)
  return(xx)
}

# Function to check new data range compatibility with calibrating data #
check_data_range <- function(model, new_data)
{
  ## TODO : remettre en marche cette fonction
  
  #   # get calibration data caracteristics
  #   expl_var_names <- model@expl_var_names
  #   expl_var_type <- model@expl_var_type
  #   expl_var_range <- model@expl_var_range
  #
  #   if(inherits(new_data, "Raster")){ ## raster data case =-=-=-=-=-=-=- #
  #     # check var names compatibility
  #     nd_expl_var_names <- names(new_data)
  #     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
  #       stop("calibration and projections variables names mismatch")
  #     }
  #     # reorder the stack
  #     new_data <- raster::subset(new_data,expl_var_names)
  #     # check var types compatibility (factors)
  #     expl_var_fact <- (expl_var_type=='factor')
  #     nd_expl_var_fact <- is.factor(new_data)
  #     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
  #       stop("calibration and projections variables class mismatch")
  #     }
  #     # check var range compatibility
  #     ### remove all new factors
  #     if(sum(expl_var_fact)>0){ ## there are factorial variables
  #       for(fact_var_id in which(expl_var_fact)){
  #         ## check if new factors occurs
  #         nd_levels <- levels(raster::subset(new_data,fact_var_id))[[1]]
  #         nd_levels <- as.character(nd_levels[,ncol(nd_levels)])
  #         names(nd_levels) <- levels(raster::subset(new_data,fact_var_id))[[1]]$ID
  #         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
  #
  #         ## detect new levels
  #         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
  #
  #         if(length(new_levels)){
  #           for(n_l in new_levels){
  #             # remove points where out of range factors have been detected
  #             new_data[subset(new_data,fact_var_id)[]==as.numeric(names(nd_levels)[which(nd_levels==n_l)])] <- NA
  #           }
  #           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
  #         }
  #       }
  #     }
  #     ## convert data to be sure to get RasterStack output
  # #     new_data <- stack(new_data)
  #   } else{ ## table data case -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
  #     # check var names compatibility
  #     nd_expl_var_names <- colnames(new_data)
  #     if(sum(!(expl_var_names %in% nd_expl_var_names) ) > 0 ){
  #       stop("calibration and projections variables names mismatch")
  #     }
  #     # reorder the stack
  #     new_data <- new_data[,expl_var_names, drop=F]
  #     # check var types compatibility (factors)
  #     expl_var_fact <- (expl_var_type=='factor')
  #     nd_expl_var_fact <- sapply(new_data,is.factor)
  #
  #     if(sum(! (expl_var_fact==nd_expl_var_fact))>0){
  #       stop("calibration and projections variables class mismatch")
  #     }
  #     # check var range compatibility
  #     ### remove all new factors
  #     if(sum(expl_var_fact)>0){ ## there are factorial variables
  #       for(fact_var_id in which(expl_var_fact)){
  #         ## check if new factors occurs
  #         nd_levels <- levels(new_data[,fact_var_id])
  #         cd_levels <- as.character(unlist(expl_var_range[[fact_var_id]]))
  #
  #         ## detect new levels
  #         new_levels <- nd_levels[!(nd_levels %in% cd_levels)]
  #
  #         if(length(new_levels)){
  #           # remove points where out of range factors have been detected
  # #           new_data <- new_data[- which(new_data[,fact_var_id] %in% new_levels),]
  #           new_data[which(new_data[,fact_var_id] %in% new_levels),] <- NA
  #           warning(paste(nd_expl_var_names[fact_var_id]," new levels have been removed from dataset (",toString(new_levels),")",sep=""))
  #         }
  #       }
  #     }
  #   }
  return(new_data)
}

## FOR SINGLE MODELS ----------------------------------------------------------
.template_predict = function(mod, object, newdata, ...)
{
  args <- list(...)
  do_check <- args$do_check
  if (is.null(do_check)) { do_check <- TRUE }
  if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }
  
  if (inherits(newdata, 'Raster')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.RasterStack(object, newdata, ...)")))
    return(res)
  } else if (inherits(newdata, 'data.frame') | inherits(newdata, 'matrix')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.data.frame(object, newdata, ...)")))
    return(res)
  } else {
    stop("invalid newdata input")
  }
}

.template_predict.RasterStack = function(seedval = NULL, predcommand, object, newdata, ...)
{
  args <- list(...)
  filename <- args$filename
  overwrite <- args$overwrite
  on_0_1000 <- args$on_0_1000
  
  if (is.null(overwrite)) { overwrite <- TRUE }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  set.seed(seedval)
  eval(parse(text = paste0("proj <- ", predcommand)))
  
  if (length(get_scaling_model(object))) {
    names(proj) <- "pred"
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5 , dat = proj)
  }
  if (on_0_1000) { proj <- round(proj * 1000) }
  
  # save raster on hard drive ?
  if (!is.null(filename)) {
    cat("\n\t\tWriting projection on hard drive...")
    if (on_0_1000) { ## projections are stored as positive integer
      writeRaster(proj, filename = filename, overwrite = overwrite, datatype = "INT2S", NAflag = -9999)
    } else { ## keep default data format for saved raster
      writeRaster(proj, filename = filename, overwrite = overwrite)
    }
    proj <- raster(filename, RAT = FALSE)
  }
  return(proj)
}

.template_predict.data.frame <- function(seedval = NULL, predcommand, object, newdata, ...)
{
  args <- list(...)
  on_0_1000 <- args$on_0_1000
  omit.na <- args$omit.na
  
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  if (is.null(omit.na)) { omit.na <- FALSE }
  
  ## check if na occurs in newdata cause they are not well supported
  if (omit.na) {
    not_na_rows <- apply(newdata, 1, function(x) { sum(is.na(x)) == 0 })
  } else {
    not_na_rows <- rep(T, nrow(newdata))
  }
  
  set.seed(seedval)
  eval(parse(text = paste0("proj <- ", predcommand)))
  
  ## add original NA from formal dataset
  if (sum(!not_na_rows) > 0) {
    tmp <- rep(NA, length(not_na_rows))
    tmp[not_na_rows] <- proj
    proj <- tmp
    rm('tmp')
  }
  
  if (length(get_scaling_model(object))) {
    proj <- data.frame(pred = proj)
    proj <- .testnull(object = get_scaling_model(object), Prev = 0.5, dat = proj)
  }
  if (on_0_1000) { proj <- round(proj * 1000) }
  
  return(proj)
}

## FOR ENSEMBLE MODELS --------------------------------------------------------
.template_predictEM = function(mod, object, newdata, formal_predictions, ...)
{
  args <- list(...)
  do_check <- args$do_check
  if (is.null(do_check)) { do_check <- TRUE }
  if (do_check) { newdata <- check_data_range(model = object, new_data = newdata) }
  
  if (inherits(newdata, 'Raster') | inherits(formal_predictions, 'Raster')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.RasterStack(object, newdata, formal_predictions, ...)")))
    return(res)
  } else if (inherits(newdata, 'data.frame') | inherits(newdata, 'matrix') | 
             inherits(formal_predictions, 'data.frame') | inherits(formal_predictions, 'matrix')) {
    eval(parse(text = paste0("res = .predict.", mod, "_biomod2_model.data.frame(object, newdata, formal_predictions, ...)")))
    return(res)
  } else {
    stop("invalid newdata input")
  }
}

.template_predictEM.formal_predictions = function(object, newdata, formal_predictions, ...)
{
  args <- list(...)
  filename <- args$filename
  on_0_1000 <- args$on_0_1000
  
  if (is.null(filename)) { filename <- "" }
  if (is.null(on_0_1000)) { on_0_1000 <- FALSE }
  
  if (is.null(formal_predictions)) {
    # make prediction of all models required
    formal_predictions <- raster::stack(lapply(object@model,
                                               function(mod, resp_name, modeling.id)
                                               {
                                                 ## check if model is loaded on memory
                                                 if (is.character(mod)) { mod <- get(load(file.path(resp_name, "models", modeling.id, mod))) }
                                                 return(predict(mod, newdata = newdata, on_0_1000 = on_0_1000))
                                               }, resp_name = object@resp_name, modeling.id = object@modeling.id))
  }
  
  return(formal_predictions)
}


## GET RESIDUAL DEVIANCE AND AIC (in GAM, deprecated ?) -------------------------------------------
## used in Biomod.Models_RE.R file

# .fun_keep <- function(object, AIC)
# {
#   list(df.resid = object$df.resid
#        , deviance = object$deviance
#        , term = as.character(object$formula)[3]
#        , AIC = AIC)
# }


## RANDOMISE DATA (only full shuffling available) -------------------------------------------------
## used in bm_VariableImportance.R file

.randomise_data <- function(data, variable, method)
{
  if (method == 'full_rand') {
    return(.full_shuffling(data, variable))
  }
}

.full_shuffling <- function(x, id = NULL)
{
  if (!(is.vector(x) | is.matrix(x) | is.data.frame(x))) {
    stop("x must be a 1 or 2 dimension odject")
  }
  
  ## Set a new random seed to ensure that sampling is random
  ## (issue when CTA is involved and seed needs to be set to a fix number)
  set.seed(as.double(Sys.time()) + as.numeric(format(Sys.time(), "%OS6")) * 1000000)
  
  out <- NULL
  if (is.null(id)) {
    out <- x[sample.int(length(x))]
  } else {
    out <- x
    for (idd in id) { out[, idd] <- out[sample.int(nrow(x)), idd]  }
  }
  
  return(out)
}
